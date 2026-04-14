from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable
from xml.sax.saxutils import escape

from Bio import SeqIO
from Bio.Restriction import AllEnzymes
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")
PROTEIN_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZ*")


@dataclass(slots=True)
class SequenceSummary:
    sequence: str
    cleaned: str
    molecule_type: str
    invalid_characters: list[str]
    length: int
    gc_percent: float | None


@dataclass(slots=True)
class OrfHit:
    frame_label: str
    start_nt: int
    end_nt: int
    strand: str
    aa_length: int
    protein: str


@dataclass(slots=True)
class SequenceFeature:
    start_nt: int
    end_nt: int
    feature_type: str
    label: str
    strand: str = "+"
    color: str | None = None


@dataclass(slots=True)
class CircularCandidate:
    start_nt: int
    end_nt: int
    overlap_nt: int
    overlap_sequence: str
    mismatch_count: int
    trim_5prime_nt: int
    trim_3prime_nt: int
    score: float
    component_sequence: str
    note: str


@dataclass(slots=True)
class MotifHit:
    start_nt: int
    end_nt: int
    strand: str
    matched_sequence: str


@dataclass(slots=True)
class RestrictionSiteHit:
    enzyme_name: str
    recognition_site: str
    matched_sequence: str
    cut_position: int
    site_start: int
    site_end: int
    cut_index: int
    cut_count: int
    all_cut_positions: list[int]


@dataclass(slots=True)
class SecondaryStructurePrediction:
    input_type: str
    folded_sequence: str
    structure: str
    mfe: float
    model_label: str


def clean_sequence(raw_text: str) -> str:
    return "".join(ch for ch in raw_text.upper() if ch.isalpha() or ch == "*")


def format_fasta(sequence: str, width: int = 70) -> str:
    if not sequence:
        return ""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def detect_molecule_type(sequence: str) -> str:
    if not sequence:
        return "Empty"

    charset = set(sequence)
    if charset <= DNA_ALPHABET:
        return "DNA"
    if charset <= RNA_ALPHABET:
        return "RNA"
    if charset <= PROTEIN_ALPHABET:
        return "Protein"
    return "Mixed/Unknown"


def invalid_characters(sequence: str, molecule_type: str | None = None) -> list[str]:
    if not sequence:
        return []

    detected = molecule_type or detect_molecule_type(sequence)
    if detected == "DNA":
        allowed = DNA_ALPHABET
    elif detected == "RNA":
        allowed = RNA_ALPHABET
    elif detected == "Protein":
        allowed = PROTEIN_ALPHABET
    else:
        allowed = DNA_ALPHABET | RNA_ALPHABET | PROTEIN_ALPHABET

    return sorted({char for char in sequence if char not in allowed})


def calculate_gc(sequence: str) -> float | None:
    if not sequence:
        return None

    molecule_type = detect_molecule_type(sequence)
    if molecule_type not in {"DNA", "RNA"}:
        return None

    gc_count = sum(1 for base in sequence if base in {"G", "C"})
    return (gc_count / len(sequence)) * 100


def summarize_sequence(raw_text: str) -> SequenceSummary:
    cleaned = clean_sequence(raw_text)
    molecule_type = detect_molecule_type(cleaned)
    invalid = invalid_characters(cleaned, molecule_type)
    return SequenceSummary(
        sequence=raw_text,
        cleaned=cleaned,
        molecule_type=molecule_type,
        invalid_characters=invalid,
        length=len(cleaned),
        gc_percent=calculate_gc(cleaned),
    )


def _resolve_restriction_enzyme(name: str):
    normalized = "".join(str(name or "").split())
    if not normalized:
        raise ValueError("Enter the name of a restriction enzyme, for example EcoRI.")
    for enzyme in AllEnzymes:
        if enzyme.__name__.casefold() == normalized.casefold():
            return enzyme
    raise ValueError(f"Restriction enzyme not found: {normalized}")


def resolve_restriction_enzyme_name(name: str) -> str:
    return _resolve_restriction_enzyme(name).__name__


def find_restriction_sites(
    sequence: str,
    *,
    filter_mode: str = "all",
    enzyme_name: str | None = None,
) -> list[RestrictionSiteHit]:
    summary = summarize_sequence(sequence)
    if summary.molecule_type != "DNA":
        raise ValueError("Restriction analysis requires a DNA sequence.")

    mode = str(filter_mode or "all").strip().lower()
    if mode not in {"all", "unique", "double"}:
        raise ValueError(f"Unknown restriction filter: {filter_mode}")

    enzymes = (
        [_resolve_restriction_enzyme(enzyme_name)]
        if enzyme_name is not None
        else sorted(AllEnzymes, key=lambda enzyme: enzyme.__name__.casefold())
    )

    sequence_obj = Seq(summary.cleaned)
    hits: list[RestrictionSiteHit] = []
    for enzyme in enzymes:
        cut_positions = sorted(int(position) for position in enzyme.search(sequence_obj, linear=True))
        cut_count = len(cut_positions)
        if cut_count == 0:
            continue
        if mode == "unique" and cut_count != 1:
            continue
        if mode == "double" and cut_count != 2:
            continue

        recognition_site = str(getattr(enzyme, "site", "") or "")
        recognition_length = len(recognition_site)
        cut_offset = getattr(enzyme, "fst5", 0)
        if not isinstance(cut_offset, int):
            cut_offset = 0

        for cut_index, cut_position in enumerate(cut_positions, start=1):
            site_start = max(1, cut_position - cut_offset)
            site_end = min(len(summary.cleaned), site_start + max(0, recognition_length - 1))
            matched_sequence = summary.cleaned[site_start - 1:site_end]
            hits.append(
                RestrictionSiteHit(
                    enzyme_name=enzyme.__name__,
                    recognition_site=recognition_site,
                    matched_sequence=matched_sequence,
                    cut_position=cut_position,
                    site_start=site_start,
                    site_end=site_end,
                    cut_index=cut_index,
                    cut_count=cut_count,
                    all_cut_positions=cut_positions.copy(),
                )
            )
    return hits


def reverse_complement(sequence: str) -> str:
    summary = summarize_sequence(sequence)
    if summary.molecule_type not in {"DNA", "RNA"}:
        raise ValueError("Reverse complement only applies to DNA or RNA sequences.")
    return str(Seq(summary.cleaned).reverse_complement())


def dna_to_rna(sequence: str) -> str:
    summary = summarize_sequence(sequence)
    if summary.molecule_type not in {"DNA", "RNA"}:
        raise ValueError("DNA to RNA conversion only applies to DNA or RNA sequences.")
    return summary.cleaned.replace("T", "U")


def predict_secondary_structure(
    sequence: str,
    molecule_type: str | None = None,
) -> SecondaryStructurePrediction:
    cleaned = clean_sequence(sequence).replace("*", "")
    detected_type = molecule_type or detect_molecule_type(cleaned)
    if detected_type not in {"DNA", "RNA"}:
        raise ValueError("Secondary structure prediction requires a DNA or RNA sequence.")
    if len(cleaned) < 5:
        raise ValueError("Sequence too short for secondary structure prediction (minimum 5 nt).")

    try:
        import RNA
    except Exception as exc:
        raise RuntimeError("ViennaRNA is not available in this environment.") from exc

    if detected_type == "DNA":
        RNA.params_load_DNA_Mathews2004()
        folded_sequence = cleaned
        model_label = "DNA (Mathews 2004)"
    else:
        RNA.params_load_RNA_Turner2004()
        folded_sequence = cleaned.replace("T", "U")
        model_label = "RNA (Turner 2004)"

    try:
        md = RNA.md()
        fc = RNA.fold_compound(folded_sequence, md)
        structure, mfe = fc.mfe()
    finally:
        RNA.params_load_RNA_Turner2004()

    return SecondaryStructurePrediction(
        input_type=detected_type,
        folded_sequence=folded_sequence,
        structure=structure,
        mfe=float(mfe),
        model_label=model_label,
    )


def export_secondary_structure_svg(
    path: str | Path,
    sequence: str,
    structure: str,
) -> None:
    try:
        import RNA
    except Exception as exc:
        raise RuntimeError("ViennaRNA is not available in this environment.") from exc

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ok = RNA.svg_rna_plot(sequence, structure, str(output_path))
    if not ok or not output_path.exists():
        raise RuntimeError("Could not generate the secondary structure SVG.")


def translate_frames(sequence: str, include_reverse: bool = True) -> list[tuple[str, str]]:
    summary = summarize_sequence(sequence)
    if summary.molecule_type not in {"DNA", "RNA"}:
        raise ValueError("Translation only applies to DNA or RNA sequences.")

    cleaned = summary.cleaned
    frames: list[tuple[str, str]] = []
    for frame in range(3):
        frame_sequence = cleaned[frame:]
        usable_length = len(frame_sequence) - (len(frame_sequence) % 3)
        translated = str(Seq(frame_sequence[:usable_length]).translate(to_stop=False))
        frames.append((f"+{frame + 1}", translated))

    if include_reverse:
        reverse = str(Seq(cleaned).reverse_complement())
        for frame in range(3):
            frame_sequence = reverse[frame:]
            usable_length = len(frame_sequence) - (len(frame_sequence) % 3)
            translated = str(Seq(frame_sequence[:usable_length]).translate(to_stop=False))
            frames.append((f"-{frame + 1}", translated))

    return frames


def find_motif(sequence: str, motif: str) -> list[int]:
    cleaned_sequence = clean_sequence(sequence)
    cleaned_motif = clean_sequence(motif)
    if not cleaned_sequence or not cleaned_motif:
        return []

    positions: list[int] = []
    start = 0
    while True:
        index = cleaned_sequence.find(cleaned_motif, start)
        if index == -1:
            break
        positions.append(index)
        start = index + 1
    return positions


def find_motif_hits(sequence: str, motif: str, include_reverse: bool = True) -> list[MotifHit]:
    cleaned_sequence = clean_sequence(sequence)
    cleaned_motif = clean_sequence(motif)
    if not cleaned_sequence or not cleaned_motif:
        return []

    hits = [
        MotifHit(
            start_nt=position + 1,
            end_nt=position + len(cleaned_motif),
            strand="+",
            matched_sequence=cleaned_sequence[position : position + len(cleaned_motif)],
        )
        for position in find_motif(cleaned_sequence, cleaned_motif)
    ]

    if include_reverse and detect_molecule_type(cleaned_sequence) in {"DNA", "RNA"}:
        reverse_motif = reverse_complement(cleaned_motif)
        if reverse_motif and reverse_motif != cleaned_motif:
            hits.extend(
                MotifHit(
                    start_nt=position + 1,
                    end_nt=position + len(reverse_motif),
                    strand="-",
                    matched_sequence=cleaned_sequence[position : position + len(reverse_motif)],
                )
                for position in find_motif(cleaned_sequence, reverse_motif)
            )

    return sorted(hits, key=lambda hit: (hit.start_nt, hit.end_nt, hit.strand))


def find_circular_candidates(
    sequence: str,
    min_overlap: int = 30,
    max_overlap: int | None = None,
    max_mismatches: int = 0,
    scan_window: int = 300,
    max_candidates: int = 12,
) -> list[CircularCandidate]:
    cleaned = clean_sequence(sequence)
    if len(cleaned) < max(10, min_overlap * 2):
        return []

    if max_overlap is None:
        max_overlap = min(len(cleaned) // 2, max(min_overlap + 50, 80))
    max_overlap = min(max_overlap, len(cleaned) // 2)
    if max_overlap < min_overlap:
        return []
    scan_window = max(0, scan_window)
    max_left_trim = min(scan_window, max(0, len(cleaned) - (min_overlap * 2)))
    max_right_trim = min(scan_window, max(0, len(cleaned) - (min_overlap * 2)))

    seen_ranges: set[tuple[int, int]] = set()
    candidates: list[CircularCandidate] = []

    for overlap in range(max_overlap, min_overlap - 1, -1):
        for left_trim in range(0, max_left_trim + 1):
            left_end = left_trim + overlap
            if left_end > len(cleaned):
                break

            for right_trim in range(0, max_right_trim + 1):
                right_start = len(cleaned) - right_trim - overlap
                if right_start < 0:
                    continue
                if left_end > right_start:
                    continue

                prefix = cleaned[left_trim:left_end]
                suffix = cleaned[right_start : right_start + overlap]
                mismatch_count = sum(1 for left, right in zip(prefix, suffix) if left != right)
                if mismatch_count > max_mismatches:
                    continue

                proposed_regions = [
                    (
                        left_trim + 1,
                        len(cleaned) - right_trim - overlap,
                        cleaned[left_trim : len(cleaned) - right_trim - overlap],
                        "trim_suffix_overlap",
                    ),
                    (
                        left_trim + overlap + 1,
                        len(cleaned) - right_trim,
                        cleaned[left_trim + overlap : len(cleaned) - right_trim],
                        "trim_prefix_overlap",
                    ),
                ]

                for start_nt, end_nt, component_sequence, note in proposed_regions:
                    if start_nt >= end_nt or not component_sequence:
                        continue
                    region_key = (start_nt, end_nt)
                    if region_key in seen_ranges:
                        continue
                    seen_ranges.add(region_key)

                    score = float(overlap * 2) - (mismatch_count * 4) - ((left_trim + right_trim) * 0.05)

                    candidates.append(
                        CircularCandidate(
                            start_nt=start_nt,
                            end_nt=end_nt,
                            overlap_nt=overlap,
                            overlap_sequence=prefix,
                            mismatch_count=mismatch_count,
                            trim_5prime_nt=left_trim,
                            trim_3prime_nt=right_trim,
                            score=score,
                            component_sequence=component_sequence,
                            note=note,
                        )
                    )

    return sorted(
        candidates,
        key=lambda item: (-item.score, -item.overlap_nt, item.start_nt, item.end_nt),
    )[:max_candidates]


def _normalize_feature_color_value(color_value: str | None) -> str | None:
    if not color_value:
        return None
    value = str(color_value).strip()
    if len(value) == 7 and value.startswith("#") and all(ch in "0123456789abcdefABCDEF" for ch in value[1:]):
        return value.upper()
    return None


def _extract_feature_color(qualifiers: dict[str, list[str]]) -> str | None:
    for key in ("ApEinfo_fwdcolor", "ApEinfo_revcolor", "color"):
        values = qualifiers.get(key, [])
        if not values:
            continue
        normalized = _normalize_feature_color_value(values[0])
        if normalized is not None:
            return normalized
    return None


def _color_for_feature(feature: SequenceFeature) -> str:
    manual_color = _normalize_feature_color_value(getattr(feature, "color", None))
    if manual_color is not None:
        return manual_color
    if getattr(feature, "strand", "+") == "-":
        palette = [
            "#fca5a5",
            "#fdba74",
            "#fda4af",
            "#fecaca",
            "#f59e0b",
            "#f9a8d4",
        ]
    else:
        palette = [
            "#84c5f4",
            "#8fd3a6",
            "#93c5fd",
            "#8fded8",
            "#b9b2ff",
            "#60a5fa",
        ]
    seed = sum(ord(char) for char in f"{feature.feature_type}:{feature.label}")
    return palette[seed % len(palette)]


def export_linear_map_svg(
    path: str | Path,
    sequence_length: int,
    features: list[SequenceFeature],
    title: str,
    highlight_range: tuple[int, int] | None = None,
) -> None:
    width = 1400
    height = 440
    left = 110
    right = width - 120
    baseline_y = 220
    feature_y = 184
    usable = max(1, right - left)

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f4f8fb"/>',
        f'<text x="{left}" y="54" font-family="Segoe UI, Arial" font-size="28" font-weight="700" fill="#16324a">{escape(title)}</text>',
        f'<text x="{left}" y="88" font-family="Segoe UI, Arial" font-size="14" fill="#5e7488">Length: {sequence_length} nt/aa</text>',
    ]

    if highlight_range is not None:
        start_nt, end_nt = highlight_range
        x1 = left + ((start_nt - 1) / max(1, sequence_length)) * usable
        x2 = left + (end_nt / max(1, sequence_length)) * usable
        elements.append(
            f'<rect x="{x1:.2f}" y="{feature_y - 26}" width="{max(8, x2 - x1):.2f}" height="86" rx="12" fill="#dbeeff" opacity="0.65"/>'
        )

    elements.extend(
        [
            f'<line x1="{left}" y1="{baseline_y}" x2="{right}" y2="{baseline_y}" stroke="#16324a" stroke-width="4" stroke-linecap="round"/>',
            f'<text x="{left - 20}" y="{baseline_y + 30}" font-family="Segoe UI, Arial" font-size="13" fill="#5e7488">1</text>',
            f'<text x="{right - 20}" y="{baseline_y + 30}" font-family="Segoe UI, Arial" font-size="13" fill="#5e7488">{sequence_length}</text>',
        ]
    )

    for index, feature in enumerate(features):
        start = max(1, min(sequence_length, feature.start_nt))
        end = max(start, min(sequence_length, feature.end_nt))
        x = left + ((start - 1) / max(1, sequence_length)) * usable
        width_feature = max(10, ((end - start + 1) / max(1, sequence_length)) * usable)
        y = feature_y - (index % 2) * 36
        label_y = y - 10 if index % 2 == 0 else y + 48
        color = _color_for_feature(feature)
        elements.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{width_feature:.2f}" height="24" rx="10" fill="{color}" opacity="0.95" stroke="#ffffff" stroke-width="1.5"/>'
        )
        elements.append(
            f'<text x="{x:.2f}" y="{label_y:.2f}" font-family="Segoe UI, Arial" font-size="13" font-weight="600" fill="#16324a">{escape(feature.label)}</text>'
        )

    elements.append("</svg>")
    Path(path).write_text("\n".join(elements), encoding="utf-8")


def _polar(center_x: float, center_y: float, radius: float, angle_deg: float) -> tuple[float, float]:
    radians = math.radians(angle_deg)
    return center_x + radius * math.cos(radians), center_y + radius * math.sin(radians)


def _arc_path(center_x: float, center_y: float, inner_radius: float, outer_radius: float, start_angle: float, end_angle: float) -> str:
    start_outer = _polar(center_x, center_y, outer_radius, start_angle)
    end_outer = _polar(center_x, center_y, outer_radius, end_angle)
    start_inner = _polar(center_x, center_y, inner_radius, start_angle)
    end_inner = _polar(center_x, center_y, inner_radius, end_angle)
    large_arc = 1 if (end_angle - start_angle) % 360 > 180 else 0
    return (
        f"M {start_outer[0]:.2f} {start_outer[1]:.2f} "
        f"A {outer_radius:.2f} {outer_radius:.2f} 0 {large_arc} 1 {end_outer[0]:.2f} {end_outer[1]:.2f} "
        f"L {end_inner[0]:.2f} {end_inner[1]:.2f} "
        f"A {inner_radius:.2f} {inner_radius:.2f} 0 {large_arc} 0 {start_inner[0]:.2f} {start_inner[1]:.2f} Z"
    )


def export_circular_map_svg(
    path: str | Path,
    sequence_length: int,
    features: list[SequenceFeature],
    title: str,
    highlight_range: tuple[int, int] | None = None,
) -> None:
    width = 980
    height = 980
    center_x = width / 2
    center_y = height / 2 + 20
    outer_radius = 280
    inner_radius = 214

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f4f8fb"/>',
        f'<text x="84" y="70" font-family="Segoe UI, Arial" font-size="28" font-weight="700" fill="#16324a">{escape(title)}</text>',
        f'<text x="84" y="102" font-family="Segoe UI, Arial" font-size="14" fill="#5e7488">Length: {sequence_length} nt/aa</text>',
        f'<circle cx="{center_x}" cy="{center_y}" r="{outer_radius}" fill="#ffffff" stroke="#d6e2ec" stroke-width="18"/>',
        f'<circle cx="{center_x}" cy="{center_y}" r="{inner_radius}" fill="#f4f8fb" stroke="#f4f8fb" stroke-width="2"/>',
    ]

    if highlight_range is not None:
        start_nt, end_nt = highlight_range
        start_angle = ((start_nt - 1) / max(1, sequence_length)) * 360 - 90
        end_angle = (end_nt / max(1, sequence_length)) * 360 - 90
        elements.append(
            f'<path d="{_arc_path(center_x, center_y, inner_radius - 24, outer_radius + 18, start_angle, end_angle)}" fill="#dbeeff" opacity="0.7"/>'
        )

    for feature in features:
        start = max(1, min(sequence_length, feature.start_nt))
        end = max(start, min(sequence_length, feature.end_nt))
        start_angle = ((start - 1) / max(1, sequence_length)) * 360 - 90
        end_angle = (end / max(1, sequence_length)) * 360 - 90
        color = _color_for_feature(feature)
        elements.append(
            f'<path d="{_arc_path(center_x, center_y, inner_radius, outer_radius, start_angle, end_angle)}" fill="{color}" opacity="0.96" stroke="#ffffff" stroke-width="1.2"/>'
        )
        mid_angle = (start_angle + end_angle) / 2
        label_x, label_y = _polar(center_x, center_y, outer_radius + 56, mid_angle)
        elements.append(
            f'<text x="{label_x:.2f}" y="{label_y:.2f}" text-anchor="middle" font-family="Segoe UI, Arial" font-size="12" font-weight="600" fill="#16324a">{escape(feature.label)}</text>'
        )

    elements.append(
        f'<text x="{center_x}" y="{center_y - 6}" text-anchor="middle" font-family="Segoe UI, Arial" font-size="18" font-weight="700" fill="#16324a">{sequence_length}</text>'
    )
    elements.append(
        f'<text x="{center_x}" y="{center_y + 20}" text-anchor="middle" font-family="Segoe UI, Arial" font-size="13" fill="#5e7488">nt/aa</text>'
    )
    elements.append("</svg>")
    Path(path).write_text("\n".join(elements), encoding="utf-8")


def _orf_candidates(
    sequence: str,
    strand: str,
    min_amino_acids: int,
) -> Iterable[OrfHit]:
    stop_codons = {"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"}
    start_codons = {"ATG", "AUG"}

    for frame in range(3):
        index = frame
        while index + 3 <= len(sequence):
            codon = sequence[index : index + 3]
            if codon in start_codons:
                scan = index + 3
                while scan + 3 <= len(sequence):
                    stop = sequence[scan : scan + 3]
                    if stop in stop_codons:
                        coding_nt = sequence[index : scan + 3]
                        protein = str(Seq(coding_nt).translate(to_stop=True))
                        if len(protein) >= min_amino_acids:
                            if strand == "+":
                                start_nt = index + 1
                                end_nt = scan + 3
                            else:
                                start_nt = len(sequence) - (scan + 3) + 1
                                end_nt = len(sequence) - index
                            yield OrfHit(
                                frame_label=f"{strand}{frame + 1}",
                                start_nt=start_nt,
                                end_nt=end_nt,
                                strand=strand,
                                aa_length=len(protein),
                                protein=protein,
                            )
                        break
                    scan += 3
            index += 3


def find_orfs(sequence: str, min_amino_acids: int = 30, include_reverse: bool = True) -> list[OrfHit]:
    summary = summarize_sequence(sequence)
    if summary.molecule_type not in {"DNA", "RNA"}:
        raise ValueError("ORF search only applies to DNA or RNA sequences.")

    cleaned = summary.cleaned
    hits = list(_orf_candidates(cleaned, "+", min_amino_acids))
    if include_reverse:
        reverse = str(Seq(cleaned).reverse_complement())
        hits.extend(_orf_candidates(reverse, "-", min_amino_acids))
    return sorted(hits, key=lambda item: (-item.aa_length, item.start_nt, item.frame_label))


def load_first_fasta(path: str | Path) -> tuple[str, str]:
    fasta_path = str(path)
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return record.id, str(record.seq).upper()


def load_sequence_file(path: str | Path) -> tuple[str, str, list[SequenceFeature], str]:
    extension = Path(path).suffix.lower()
    if extension in {".gb", ".gbk", ".genbank"}:
        file_format = "genbank"
    else:
        file_format = "fasta"

    input_path = str(path)
    if file_format == "fasta":
        record = next(SeqIO.parse(input_path, "fasta"))
    else:
        record = next(SeqIO.parse(input_path, file_format))
    features: list[SequenceFeature] = []
    if file_format == "genbank":
        for feature in record.features:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            if start > end:
                continue
            label = (
                feature.qualifiers.get("label", [""])[0]
                or feature.qualifiers.get("gene", [""])[0]
                or feature.qualifiers.get("note", [""])[0]
                or feature.type
            )
            features.append(
                SequenceFeature(
                    start_nt=start,
                    end_nt=end,
                    feature_type=feature.type,
                    label=label,
                    strand="-" if feature.location.strand == -1 else "+",
                    color=_extract_feature_color(feature.qualifiers),
                )
            )

    return record.id, str(record.seq).upper(), features, file_format


_KD_SCALE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9,  "F": 2.8,  "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
_TM_WINDOW   = 20
_TM_THRESHOLD = 1.6   # mean KD score above this → putative TM segment


def _scan_tm_helices(seq: str) -> list[tuple[int, int, float]]:
    """Return list of (start, end, mean_kd) for putative TM windows."""
    hits = []
    for i in range(len(seq) - _TM_WINDOW + 1):
        window = seq[i: i + _TM_WINDOW]
        score = sum(_KD_SCALE.get(aa, 0.0) for aa in window) / _TM_WINDOW
        if score >= _TM_THRESHOLD:
            if hits and i <= hits[-1][1]:
                # extend previous hit if overlapping
                hits[-1] = (hits[-1][0], i + _TM_WINDOW, max(hits[-1][2], score))
            else:
                hits.append((i + 1, i + _TM_WINDOW, score))
    return hits


def _fetch_protein_sol(sequence: str, timeout: int = 20) -> float | None:
    """Query Protein-Sol (Manchester) via form submission.
    Returns scaled-sol score in [0, 1] (population mean ~0.45) or None.
    Current flow redirects to a hosted results.html page."""
    try:
        import re
        import urllib.parse
        import urllib.request

        fasta = f">query\n{sequence}"
        form_data = urllib.parse.urlencode({
            "singleprediction": "Submit",
            "sequence-input": fasta,
        }).encode()
        req = urllib.request.Request(
            "https://protein-sol.manchester.ac.uk/cgi-bin/solubility/sequenceprediction.php",
            data=form_data,
            headers={"Content-Type": "application/x-www-form-urlencoded",
                     "User-Agent": "GeneDraft/1.0"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        redirect_match = re.search(r'window\.location\s*=\s*"([^"]+results\.html)"', html)
        if redirect_match:
            result_url = redirect_match.group(1)
        else:
            job_match = re.search(r"Job id = ([A-Za-z0-9]+)", html)
            if not job_match:
                timestamp_match = re.search(r'var\s+timestamp\s*=\s*"([A-Za-z0-9]+)"', html)
                if not timestamp_match:
                    return None
                job_id = timestamp_match.group(1)
            else:
                job_id = job_match.group(1)
            result_url = (
                "https://protein-sol.manchester.ac.uk/"
                f"results/solubility/run-{job_id}/results.html"
            )

        req2 = urllib.request.Request(result_url, headers={"User-Agent": "GeneDraft/1.0"})
        with urllib.request.urlopen(req2, timeout=timeout) as resp2:
            result_html = resp2.read().decode("utf-8", errors="replace")

        score_match = re.search(
            r"Predicted scaled solubility:</h5>\s*<p[^>]*>\s*([01](?:\.\d+)?)\s*</p>",
            result_html,
            re.I,
        )
        if not score_match:
            score_match = re.search(
                r"scaled solubility value[^0-9]{0,80}([01](?:\.\d+)?)",
                result_html,
                re.I,
            )
        if not score_match:
            return None
        score = float(score_match.group(1))
        if 0.0 <= score <= 1.0:
            return score
    except Exception:
        pass
    return None


def _fetch_soluprot(sequence: str, timeout: int = 30) -> float | None:
    """Query SoluProt (Loschmidt Labs, Masaryk University) via form submission.
    Asynchronous: submits job, polls for result, returns score in [0,1] or None.
    Cutoff: >=0.5 soluble. Endpoint confirmed from page source."""
    try:
        import re
        import time
        import urllib.parse
        import urllib.request

        fasta = f">query\n{sequence}"
        form_data = urllib.parse.urlencode({"fasta": fasta}).encode()
        req = urllib.request.Request(
            "https://loschmidt.chemi.muni.cz/soluprot/",
            data=form_data,
            headers={"Content-Type": "application/x-www-form-urlencoded",
                     "User-Agent": "GeneDraft/1.0"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            html = resp.read().decode("utf-8", errors="replace")
            final_url = resp.geturl()

        # The site redirects to /?job=<id>; older versions embedded the job in HTML.
        job_match = re.search(r'[?&]job=([A-Za-z0-9_\-]+)', final_url)
        if not job_match:
            job_match = re.search(r'[?&]job=([A-Za-z0-9_\-]+)', html)
        if not job_match:
            return None
        job_id = job_match.group(1)

        # Poll for result (max ~25 s)
        result_url = f"https://loschmidt.chemi.muni.cz/soluprot/?job={job_id}"
        for _ in range(5):
            time.sleep(5)
            req2 = urllib.request.Request(result_url,
                                          headers={"User-Agent": "GeneDraft/1.0"})
            with urllib.request.urlopen(req2, timeout=timeout) as resp2:
                result_html = resp2.read().decode("utf-8", errors="replace")

            score_match = re.search(
                r'<table[^>]*class="result_table"[^>]*>.*?<tbody>.*?<tr>.*?'
                r'<td[^>]*>[^<]+</td>\s*<td[^>]*>\s*([01](?:\.\d+)?)\s*</td>',
                result_html,
                re.I | re.S,
            )
            if not score_match:
                score_match = re.search(
                    r'score above 0\.5.*?<td[^>]*>\s*([01](?:\.\d+)?)\s*</td>',
                    result_html,
                    re.I | re.S,
                )
            if not score_match:
                score_match = re.search(
                    r'(?:solubility|score)[^\d]{0,80}([01](?:\.\d+)?)',
                    result_html,
                    re.I | re.S,
                )
            if score_match:
                score = float(score_match.group(1))
                if 0.0 <= score <= 1.0:
                    return score
    except Exception:
        pass
    return None


def fetch_hmmer_pfam(sequence: str, timeout: int = 60) -> list[dict] | None:
    """Run hmmscan against Pfam via EBI HMMER REST API (v1).
    Returns list of domain hit dicts or None if unavailable/no internet.
    Each dict has: name, acc, desc, evalue, start, end, cevalue.
    Safe to call with no internet — returns None silently."""
    try:
        import json
        import time
        import urllib.request

        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        payload = json.dumps({
            "input": f">query\n{sequence}",
            "database": "pfam",
        }).encode()

        # Submit job
        req = urllib.request.Request(
            "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan",
            data=payload,
            headers=headers,
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            job_data = json.loads(resp.read())

        job_id = job_data.get("id") or job_data.get("uuid")
        if not job_id:
            return None

        # Poll for result (max ~55 s total)
        result_url = f"https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{job_id}"
        get_headers = {"Accept": "application/json"}
        deadline = time.time() + timeout
        while time.time() < deadline:
            time.sleep(5)
            req2 = urllib.request.Request(result_url, headers=get_headers)
            with urllib.request.urlopen(req2, timeout=20) as resp2:
                result = json.loads(resp2.read())
            if result.get("status") == "SUCCESS":
                break
        else:
            return None

        # Parse hits — EBI currently returns {"result": {"hits": [...]}}.
        # Keep backward compatibility with older payloads that used "results".
        result_node = result.get("result") or result.get("results") or result
        hits_raw = result_node.get("hits") or []
        hits: list[dict] = []
        for hit in hits_raw:
            for dom in hit.get("domains", []):
                if not dom.get("is_included"):
                    continue
                metadata = hit.get("metadata") or {}
                hits.append({
                    "name":    metadata.get("identifier") or hit.get("name", ""),
                    "acc":     metadata.get("accession") or hit.get("acc",  ""),
                    "desc":    metadata.get("description") or hit.get("desc", ""),
                    "evalue":  hit.get("evalue", None),
                    "start":   dom.get("iali",   dom.get("ienv", 0)),
                    "end":     dom.get("jali",   dom.get("jenv", 0)),
                    "cevalue": dom.get("cevalue", dom.get("ievalue", None)),
                })
        return hits
    except Exception:
        return None


# Sentinel: distinguishes "not yet fetched" from "fetched, returned nothing"
class _HmmerPending:
    pass
_HMMER_PENDING = _HmmerPending()


def analyze_protein(sequence: str, **kwargs) -> str:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    cleaned = clean_sequence(sequence).replace("*", "")
    if detect_molecule_type(cleaned) != "Protein":
        raise ValueError("Protein analysis requires a protein sequence. Translate your DNA/RNA sequence first.")

    _STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
    std_seq = "".join(aa for aa in cleaned if aa in _STANDARD)
    if len(std_seq) < 5:
        raise ValueError("Sequence too short for protein analysis (minimum 5 residues).")

    pa = ProteinAnalysis(std_seq)
    mw          = pa.molecular_weight()
    pi          = pa.isoelectric_point()
    instability = pa.instability_index()
    gravy       = pa.gravy()
    aromaticity = pa.aromaticity()
    helix, turn, sheet = pa.secondary_structure_fraction()
    aa_pct_raw = getattr(pa, "amino_acids_percent", None)
    if callable(aa_pct_raw):
        aa_pct = aa_pct_raw()
    elif isinstance(aa_pct_raw, dict):
        aa_pct = aa_pct_raw
    else:
        aa_pct = {aa: std_seq.count(aa) / len(std_seq) for aa in _STANDARD}

    n_cys     = std_seq.count("C")
    charged   = sum(aa_pct.get(aa, 0) for aa in "DEKR")
    aliphatic = sum(aa_pct.get(aa, 0) for aa in "ILVAM")

    # ── Solubility — local estimate used for expression/localization logic ──
    sol_local = max(0.0, min(1.0, 0.17 * charged - 0.25 * aliphatic
                             - 0.22 * max(0.0, gravy) + 0.27))

    # Use provided remote scores when available; fall back to local
    ps_score = kwargs.get("ps_score")   # Protein-Sol scaled-sol (cutoff 0.45)
    sp_score = kwargs.get("sp_score")   # SoluProt score        (cutoff 0.50)
    fetching = kwargs.get("fetching", False)

    # The score driving expression/localization is Protein-Sol if available,
    # otherwise SoluProt, otherwise local estimate
    sol = ps_score if ps_score is not None else (
          sp_score if sp_score is not None else sol_local)

    def _sol_line(score, label, cutoff, source, pending):
        if pending:
            return f"  {source:<12} fetching..."
        if score is None:
            return f"  {source:<12} unavailable"
        tag = "soluble" if score >= cutoff else "low solubility"
        return f"  {source:<12} {score:.3f}  ({tag}, cutoff {cutoff})"

    # ── Crystallizability (XtalPred-inspired, 0–7) ───────────────────────
    xtal = 0
    if 10_000 <= mw <= 50_000:  xtal += 2
    elif mw <= 100_000:         xtal += 1
    if 5.0 <= pi <= 9.0:        xtal += 2
    else:                       xtal += 1
    if instability < 40:        xtal += 2
    elif instability < 60:      xtal += 1
    if -0.5 <= gravy <= 0.5:    xtal += 1
    xtal_label = ("Favorable" if xtal >= 6 else "Moderate" if xtal >= 4 else "Unfavorable")

    # ── Transmembrane helices (Kyte-Doolittle window scan) ───────────────
    tm_hits = _scan_tm_helices(std_seq)
    is_tm   = len(tm_hits) > 0

    # ── E. coli heterologous expression probability ───────────────────────
    expr_score = 0
    if gravy < 0:          expr_score += 2
    if instability < 40:   expr_score += 2
    if mw < 40_000:        expr_score += 2
    elif mw < 100_000:     expr_score += 1
    if not is_tm:          expr_score += 1
    if charged > 0.25:     expr_score += 1
    if sol > 0.50:         expr_score += 1
    expr_pct = round(expr_score / 9 * 100)
    expr_label = ("Favorable" if expr_pct >= 70
                  else "Moderate" if expr_pct >= 45 else "Unfavorable")

    # ── Predicted localization ────────────────────────────────────────────
    if is_tm:
        loc_label  = "Membrane fraction"
        loc_detail = (f"  {len(tm_hits)} putative TM segment(s) detected.\n"
                      f"  Protein will likely integrate into the lipid bilayer.\n"
                      f"  Consider detergent extraction or cell-free expression.")
    elif sol >= 0.50 and instability < 40 and gravy < 0:
        loc_label  = "Soluble cytoplasm"
        loc_detail = ("  Low hydrophobicity and good stability predict\n"
                      "  soluble expression in the E. coli cytoplasm.")
    elif sol < 0.38 or (gravy > 0.2 and instability >= 40):
        loc_label  = "Inclusion bodies (likely)"
        loc_detail = ("  High hydrophobicity and/or instability predict\n"
                      "  aggregation into inclusion bodies.\n"
                      "  Consider: lower temperature (16-18 °C), fusion tags\n"
                      "  (SUMO, MBP, Trx), or refolding from IB pellet.")
    elif sol < 0.50 or gravy > 0:
        loc_label  = "Inclusion bodies / soluble (mixed)"
        loc_detail = ("  Borderline physicochemical profile — outcome depends\n"
                      "  strongly on expression conditions.\n"
                      "  Recommend: screen multiple temperatures and fusion tags.")
    else:
        loc_label  = "Soluble cytoplasm (moderate confidence)"
        loc_detail = ("  Properties are consistent with soluble expression,\n"
                      "  but experimental validation is recommended.")

    top_aa = sorted(aa_pct.items(), key=lambda x: -x[1])

    sol_note = "  (updating with remote scores...)" if fetching else ""

    lines = [
        "Protein Analysis Report",
        "=" * 46,
        f"Length:               {len(std_seq)} aa",
        f"Molecular weight:     {mw:,.1f} Da  ({mw / 1000:.2f} kDa)",
        f"Isoelectric point:    {pi:.2f}",
        "",
        "Stability",
        f"  Instability index:  {instability:.2f}  ({'stable' if instability < 40 else 'unstable'})",
        f"  GRAVY score:        {gravy:.4f}  ({'hydrophilic' if gravy < 0 else 'hydrophobic'})",
        f"  Aromaticity:        {aromaticity:.4f}",
        "",
        "Predicted secondary structure",
        f"  Alpha-helix:        {helix * 100:.1f}%",
        f"  Beta-turn:          {turn * 100:.1f}%",
        f"  Beta-sheet:         {sheet * 100:.1f}%",
        "",
        f"Cysteines (C):        {n_cys}",
        f"  Possible disulfide: {n_cys // 2} pair(s)",
        "",
        "Solubility" + sol_note,
        _sol_line(ps_score, "Protein-Sol", 0.45, "Protein-Sol", fetching and ps_score is None),
        _sol_line(sp_score, "SoluProt",    0.50, "SoluProt",    fetching and sp_score is None),
        f"  Local est.   {sol_local:.3f}  (linear model, reference only)",
        f"Crystallizability:    {xtal_label}  (score {xtal}/7)",
        "",
        "Transmembrane segments (Kyte-Doolittle, window 20)",
    ]
    if tm_hits:
        for start, end, score in tm_hits:
            lines.append(f"  Segment {start}–{end}  (mean KD {score:.2f})")
    else:
        lines.append("  None detected")

    # ── Pfam domains (HMMER) ─────────────────────────────────────────────
    hmmer_hits: list[dict] | None = kwargs.get("hmmer_hits", _HMMER_PENDING)
    pfam_lines = ["", "─" * 46, "Pfam Domains  (HMMER / EBI)", "─" * 46]
    if hmmer_hits is _HMMER_PENDING:
        pfam_lines.append("  fetching...")
    elif hmmer_hits is None:
        pfam_lines.append("  unavailable (no internet or server error)")
    elif len(hmmer_hits) == 0:
        pfam_lines.append("  No significant Pfam domains found")
    else:
        pfam_lines.append(f"  {'Accession':<12} {'Name':<18} {'Pos':>10}  {'E-value':<12}  Description")
        pfam_lines.append(f"  {'-'*12} {'-'*18} {'-'*10}  {'-'*12}  {'-'*30}")
        for h in hmmer_hits:
            acc_short = h['acc'].split('.')[0]
            pos       = f"{h['start']}-{h['end']}"
            ev        = f"{h['cevalue']:.2e}" if h['cevalue'] is not None else "n/a"
            pfam_lines.append(
                f"  {acc_short:<12} {h['name']:<18} {pos:>10}  {ev:<12}  {h['desc']}"
            )

    lines += pfam_lines
    lines += [
        "",
        "─" * 46,
        "E. coli Heterologous Expression",
        "─" * 46,
        f"  Expression probability:  {expr_pct}%  ({expr_label})",
        f"  Predicted localization:  {loc_label}",
        loc_detail,
        "",
        "Amino acid composition (top 10)",
        f"  {'AA':<5} {'%':>6}",
    ]
    for aa, pct in top_aa[:10]:
        lines.append(f"  {aa:<5} {pct * 100:>5.1f}%")

    return "\n".join(lines)


def fetch_by_accession(
    accession: str,
    email: str,
) -> tuple[str, str, list[SequenceFeature], str]:
    import io as _io
    import time as _time
    from Bio import Entrez

    Entrez.email = email
    acc = accession.strip()

    errors: list[str] = []
    for db in ("nucleotide", "protein"):
        try:
            handle = Entrez.efetch(db=db, id=acc, rettype="gb", retmode="text")
            text = handle.read()
            handle.close()
            record = next(SeqIO.parse(_io.StringIO(text), "genbank"))

            features: list[SequenceFeature] = []
            for feat in record.features:
                if feat.type == "source":
                    continue
                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                if start > end:
                    continue
                label = (
                    feat.qualifiers.get("label",   [""])[0]
                    or feat.qualifiers.get("gene",    [""])[0]
                    or feat.qualifiers.get("product", [""])[0]
                    or feat.qualifiers.get("note",    [""])[0]
                    or feat.type
                )
                features.append(SequenceFeature(
                    start_nt=start,
                    end_nt=end,
                    feature_type=feat.type,
                    label=label,
                    strand="-" if feat.location.strand == -1 else "+",
                    color=_extract_feature_color(feat.qualifiers),
                ))
            return record.id, str(record.seq).upper(), features, db
        except Exception as exc:
            errors.append(f"{db}: {exc}")
            _time.sleep(0.5)

    raise ValueError(f"Could not fetch '{acc}' from NCBI.\n" + "\n".join(errors))


def sanitize_fasta_for_blast(source_path: str | Path, output_path: str | Path) -> tuple[int, str]:
    input_path = str(source_path)
    records = list(SeqIO.parse(input_path, "fasta"))

    if not records:
        raise ValueError("No valid sequences were detected in the FASTA file.")

    sanitized_records: list[SeqRecord] = []
    seen_ids: set[str] = set()
    for index, record in enumerate(records, start=1):
        raw_id = (record.id or f"seq_{index}").strip()
        raw_description = " ".join(str(record.description or "").split()).strip()
        source_label = raw_description or raw_id or f"seq_{index}"
        safe_id_base = "".join(char if (char.isalnum() or char in {"_", "-", "."}) else "_" for char in source_label)
        safe_id_base = safe_id_base[:180].strip("._-") or f"seq_{index}"
        safe_id = safe_id_base
        duplicate_index = 2
        while safe_id in seen_ids:
            suffix = f"_{duplicate_index}"
            allowed_base_length = max(1, 180 - len(suffix))
            safe_id = f"{safe_id_base[:allowed_base_length].rstrip('._-') or 'seq'}{suffix}"
            duplicate_index += 1
        seen_ids.add(safe_id)

        safe_description = raw_description[:240]
        sanitized_records.append(
            SeqRecord(
                Seq(str(record.seq).upper()),
                id=safe_id,
                name=safe_id,
                description=safe_description,
            )
        )

    SeqIO.write(sanitized_records, str(output_path), "fasta")
    sequence_type = detect_molecule_type(str(sanitized_records[0].seq))
    return len(sanitized_records), ("prot" if sequence_type == "Protein" else "nucl")


def save_fasta(path: str | Path, header: str, sequence: str) -> None:
    cleaned = clean_sequence(sequence)
    record = SeqRecord(Seq(cleaned), id=(header or "sequence").strip(), description="")
    SeqIO.write(record, str(path), "fasta")


def save_genbank(
    path: str | Path,
    header: str,
    sequence: str,
    features: list[SequenceFeature] | None = None,
) -> None:
    cleaned = clean_sequence(sequence)
    molecule_type = detect_molecule_type(cleaned)
    if molecule_type not in {"DNA", "RNA", "Protein"}:
        molecule_type = "DNA"

    if molecule_type == "DNA":
        genbank_type = "DNA"
    elif molecule_type == "RNA":
        genbank_type = "RNA"
    else:
        genbank_type = "protein"

    record = SeqRecord(
        Seq(cleaned),
        id=(header or "sequence").strip(),
        name=(header or "sequence").strip(),
        description="",
    )
    record.annotations["molecule_type"] = genbank_type

    record.features = []
    for feature in features or []:
        start = max(0, feature.start_nt - 1)
        end = min(len(cleaned), feature.end_nt)
        if start >= end:
            continue
        qualifiers = {"label": [feature.label or "feature"]}
        feature_color = _normalize_feature_color_value(getattr(feature, "color", None))
        if feature_color is not None:
            qualifiers["ApEinfo_fwdcolor"] = [feature_color]
            qualifiers["ApEinfo_revcolor"] = [feature_color]
            qualifiers["color"] = [feature_color]
        record.features.append(
            SeqFeature(
                FeatureLocation(start, end, strand=(-1 if feature.strand == "-" else 1)),
                type=(feature.feature_type or "misc_feature").strip(),
                qualifiers=qualifiers,
            )
        )

    SeqIO.write(record, str(path), "genbank")
