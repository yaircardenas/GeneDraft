from __future__ import annotations

import ctypes
from dataclasses import dataclass
from itertools import cycle
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
from urllib.parse import urlencode
import webbrowser

from sequence_tools import (
    CircularCandidate,
    MotifHit,
    OrfHit,
    SequenceFeature,
    SequenceSummary,
    clean_sequence,
    export_circular_map_svg,
    export_linear_map_svg,
    find_motif_hits,
    find_orfs,
    find_circular_candidates,
    format_fasta,
    load_sequence_file,
    dna_to_rna,
    export_secondary_structure_svg,
    predict_secondary_structure,
    reverse_complement,
    sanitize_fasta_for_blast,
    save_fasta,
    save_genbank,
    summarize_sequence,
    translate_frames,
    analyze_protein,
    fetch_by_accession,
)


APP_VERSION = "1.0"
APP_RELEASE_YEAR = "2026"
APP_REPOSITORY_URL = "https://github.com/yaircardenas/GeneDraft"
APP_AUTHOR = "Yair Cárdenas-Conejo"
APP_CONTACT_EMAIL = "ycardenasco@secihti.mx"
APP_ORCID = "0000-0002-0190-244X"
APP_AFFILIATION = (
    "1 Secretaría de Ciencia, Humanidades, Tecnología e Innovación, Gobierno de México, "
    "Ciudad de México, México\n"
    "2 Centro Universitario de Investigaciones Biomédicas, Universidad de Colima, "
    "Colima, Colima, México"
)
APP_DESCRIPTION = (
    "GeneDraft is a local sequence workbench for fast sequence editing and lightweight "
    "DNA, RNA, and protein analysis."
)


class _RoundBtn(tk.Button):
    """Toolbar button with a smooth hover effect."""

    def __init__(self, parent, text, command, *,
                 normal_bg: str, hover_bg: str, press_bg: str,
                 fg: str = "#C9D1D9",
                 radius: int = 6, padx: int = 12, pady: int = 5,
                 font: tuple = ("Segoe UI Variable Text Semibold", 9)):
        super().__init__(
            parent,
            text=text,
            command=command,
            bg=normal_bg,
            fg=fg,
            activebackground=hover_bg,
            activeforeground="#E6EDF3",
            relief="flat",
            bd=0,
            padx=padx,
            pady=pady,
            font=font,
            cursor="hand2",
            highlightthickness=1,
            highlightbackground=normal_bg,
            highlightcolor=press_bg,
        )
        self._normal = normal_bg
        self._hover  = hover_bg
        self._press  = press_bg
        self._fg     = fg
        # Lambdas reference self.* so they always use the current theme colors.
        self.bind("<Enter>",           lambda e: self.configure(bg=self._hover, highlightbackground=self._hover))
        self.bind("<Leave>",           lambda e: self.configure(bg=self._normal, fg=self._fg, highlightbackground=self._normal))
        self.bind("<ButtonPress-1>",   lambda e: self.configure(bg=self._press, highlightbackground=self._press))
        self.bind("<ButtonRelease-1>", lambda e: self.configure(bg=self._hover, highlightbackground=self._hover))


@dataclass(slots=True)
class FeatureEntry:
    start_nt: int
    end_nt: int
    feature_type: str
    label: str
    color: str
    strand: str = "+"

    def as_sequence_feature(self) -> SequenceFeature:
        return SequenceFeature(
            start_nt=self.start_nt,
            end_nt=self.end_nt,
            feature_type=self.feature_type,
            label=self.label,
            strand=self.strand,
        )


@dataclass(slots=True)
class SearchMatch:
    start_nt: int
    end_nt: int
    title: str
    details: str
    kind: str
    payload: object | None = None


@dataclass(slots=True)
class MarkerEntry:
    position_nt: int
    label: str


@dataclass(slots=True)
class OperationSnapshot:
    description: str
    sequence_name: str
    editor_text: str
    features: list[FeatureEntry]
    markers: list[MarkerEntry]


@dataclass(slots=True)
class SecondaryStructureExportState:
    suggested_path: Path
    folded_sequence: str
    structure: str


class GeneDraftApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("GeneDraft")
        self.root.geometry("1720x1060")
        self.root.minsize(1380, 860)
        self.root.state("zoomed")

        self.current_theme: str = "night"
        self.colors = self._build_night_colors()
        self._toolbar_btns: list[_RoundBtn] = []
        self._all_menus: list[tk.Menu] = []
        self.app_icon = None
        self.brand_mark_image = None
        self.brand_logo_image = None

        self.editor_font_size = 15
        self.results_font_size = 12
        self.sequence_name = tk.StringVar(value="sequence_1")
        self.length_var    = tk.StringVar(value="Length: 0")
        self.type_var      = tk.StringVar(value="Type: Empty")
        self.gc_var        = tk.StringVar(value="GC: -")
        self.warning_var   = tk.StringVar(value="Ready")
        self.cursor_var    = tk.StringVar(value="Pos: 0")
        self.selection_var = tk.StringVar(value="Selection: 0 nt")
        self.features_var  = tk.StringVar(value="Features: 0")
        self.orf_min_var           = tk.IntVar(value=110)
        self.line_width_var        = tk.IntVar(value=100)
        self.circular_overlap_var    = tk.IntVar(value=15)
        self.circular_mismatch_var   = tk.IntVar(value=2)
        self.circular_scan_var       = tk.IntVar(value=300)
        self.local_blast_source_var = tk.StringVar(value="")
        self.local_blast_db_var     = tk.StringVar(value="")
        self.local_blast_type_var   = tk.StringVar(value="auto")
        self.settings_path = Path.cwd() / "genedraft_settings.json"
        self._recent_files: list[str] = []
        self._ncbi_email: str = ""
        self.session_path = Path.cwd() / "genedraft_session.json"
        self.sequence_name.trace_add("write", self._on_sequence_name_changed)
        self.feature_type_options = (
            "misc_feature", "gene", "CDS", "promoter", "terminator",
            "regulatory", "primer_bind", "repeat_region", "misc_binding",
            "rep_origin", "exon", "intron",
        )
        self.feature_strand_label_to_value = {
            "Forward (+)": "+",
            "Reverse (-)": "-",
        }
        self.feature_strand_value_to_label = {
            value: label for label, value in self.feature_strand_label_to_value.items()
        }

        self.features: list[FeatureEntry] = []
        self.markers: list[MarkerEntry] = []
        self._feature_tags: set[str] = set()
        self._analysis_tags = {
            "motif", "motif_focus", "invalid",
            "orf_plus", "orf_minus", "orf_focus",
            "circular_candidate", "circular_focus", "feature_focus",
        }
        self._feature_colors = cycle([
            "#93c5fd",   # sky blue
            "#6ee7b7",   # emerald
            "#fde68a",   # yellow
            "#fda4af",   # pink
            "#a5f3fc",   # cyan
            "#fdba74",   # orange
            "#c4b5fd",   # lavender
        ])
        self._configure_app_icon()
        self.search_matches: list[SearchMatch] = []
        self.search_index = -1
        self.search_label = ""
        self._internal_edit = False
        self._clean_index_map_cache: list[str] | None = None
        self._history_undo: list[OperationSnapshot] = []
        self._history_redo: list[OperationSnapshot] = []
        self._session_after_id: str | None = None
        self._minimap_after_id: str | None = None
        self._busy_done_after_id: str | None = None
        self._restoring_session = False
        self._session_has_unsaved_work = False
        self._session_had_persistence_event = False
        self._busy_status_message: str | None = None
        self._results_secondary_structure_export: SecondaryStructureExportState | None = None
        self._results_secondary_svg_menu_index: int | None = None

        self._load_app_settings()
        self._configure_styles()
        self._build_menu()
        self._rebuild_recent_menu()
        self._build_ui()
        self._build_context_menus()
        self._bind_shortcuts()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        self.root.after(120, self._apply_initial_layout)
        if not self._try_restore_session():
            self._set_results(
                "Welcome to GeneDraft.\n\n"
                "Paste a sequence, open a FASTA or GenBank file, and use the menu or toolbar.\n"
                "You can search ORFs across the whole sequence or only inside the current selection."
            )
        self._refresh_summary()

    # ------------------------------------------------------------------ palettes

    @staticmethod
    def _build_night_colors() -> dict[str, str]:
        return {
            "bg":           "#0E1117",
            "surface":      "#161C26",
            "surface_alt":  "#1D2535",
            "border":       "#2A3347",
            "border_focus": "#22D3EE",
            "text":         "#E2E8F0",
            "muted":        "#7A8BA0",
            "accent":       "#0EA5E9",
            "accent_dark":  "#38BDF8",
            "accent_soft":  "#0C2233",
            "selection_bg": "#B45309",
            "selection_fg": "#FEF3C7",
            "toolbar_bg":   "#161C26",
            "toolbar_fg":   "#CBD5E1",
            "toolbar_btn":  "#1E2A3A",
            "toolbar_hover":"#263348",
            "toolbar_press":"#0EA5E9",
            "toolbar_sep":  "#2A3347",
            "settings_bg":  "#111822",
            "settings_fg":  "#38BDF8",
            "chrome_text":  "#E2E8F0",
            "chrome_muted": "#7A8BA0",
            "success":      "#052E16",
            "success_fg":   "#4ADE80",
            "warning":      "#2D1B00",
            "warning_fg":   "#FBBF24",
            "danger":       "#2D0A0A",
            "danger_fg":    "#F87171",
            "tag_motif":        "#3D2E00",
            "tag_orf_plus":     "#065F46",
            "tag_orf_minus":    "#881337",
            "tag_orf_focus":    "#1E40AF",
            "tag_circ":         "#0A2540",
            "tag_circ_foc":     "#0E4272",
            "tag_invalid":      "#4A0E0E",
            "tag_feature":      "#0D2B50",
            "tag_feature_focus":"#92400E",
        }

    @staticmethod
    def _build_day_colors() -> dict[str, str]:
        return {
            "bg":           "#F0F4F8",
            "surface":      "#FFFFFF",
            "surface_alt":  "#EBF0F5",
            "border":       "#C8D6E2",
            "border_focus": "#0284C7",
            "text":         "#0F1F2E",
            "muted":        "#5A7388",
            "accent":       "#0284C7",
            "accent_dark":  "#0369A1",
            "accent_soft":  "#E0F2FE",
            "selection_bg": "#0284C7",
            "selection_fg": "#FFFFFF",
            "toolbar_bg":   "#E8EFF5",
            "toolbar_fg":   "#1E3448",
            "toolbar_btn":  "#FFFFFF",
            "toolbar_hover":"#D1E8F5",
            "toolbar_press":"#BAD9EE",
            "toolbar_sep":  "#C8D6E2",
            "settings_bg":  "#DDE8F0",
            "settings_fg":  "#0369A1",
            "chrome_text":  "#0F1F2E",
            "chrome_muted": "#5A7388",
            "success":      "#DCFCE7",
            "success_fg":   "#14532D",
            "warning":      "#FEF9C3",
            "warning_fg":   "#713F12",
            "danger":       "#FEE2E2",
            "danger_fg":    "#7F1D1D",
            "tag_motif":        "#FEF08A",
            "tag_orf_plus":     "#BBF7D0",
            "tag_orf_minus":    "#FECDD3",
            "tag_orf_focus":    "#BAE6FD",
            "tag_circ":         "#E0F2FE",
            "tag_circ_foc":     "#7DD3FC",
            "tag_invalid":      "#FECACA",
            "tag_feature":      "#DBEAFE",
            "tag_feature_focus":"#FBBF24",
        }

    # ------------------------------------------------------------------ styles

    def _configure_styles(self) -> None:
        C = self.colors
        self.root.configure(bg=C["bg"])
        style = ttk.Style()
        style.theme_use("clam")

        # Try modern Windows font; fall back gracefully
        ui_font     = ("Segoe UI Variable Text", 11)
        ui_semi     = ("Segoe UI Variable Text Semibold", 11)
        ui_small    = ("Segoe UI Variable Text", 10)
        ui_semi_sm  = ("Segoe UI Variable Text Semibold", 10)

        style.configure(".", font=ui_font, background=C["bg"], foreground=C["text"])

        # ── Frames ──────────────────────────────────────────────────────────────
        for name, bg in (
            ("App.TFrame",     C["bg"]),
            ("Card.TFrame",    C["surface"]),
            ("Toolbar.TFrame", C["toolbar_bg"]),
            ("Status.TFrame",  C["surface"]),
            ("Sep.TFrame",     C["border"]),
        ):
            style.configure(name, background=bg)

        # ── Labels ──────────────────────────────────────────────────────────────
        style.configure("TLabel",             background=C["bg"],          foreground=C["text"],        font=ui_font)
        style.configure("Card.TLabel",        background=C["surface"],     foreground=C["accent_dark"], font=ui_semi)
        style.configure("Muted.TLabel",       background=C["surface"],     foreground=C["muted"],       font=ui_small)
        style.configure("BrandTitle.TLabel",  background=C["bg"],          foreground=C["chrome_text"], font=ui_semi)
        style.configure("BrandSub.TLabel",    background=C["bg"],          foreground=C["chrome_muted"],font=ui_small)
        style.configure("Status.TLabel",      background=C["surface"],     foreground=C["muted"],       font=ui_small)
        style.configure("SectionHead.TLabel", background=C["surface"],     foreground=C["accent_dark"], font=ui_semi_sm)
        style.configure("Settings.TLabel",    background=C["settings_bg"], foreground=C["settings_fg"], font=ui_semi_sm)

        # ── Action buttons (secondary) ───────────────────────────────────────────
        style.configure(
            "Action.TButton",
            font=ui_semi_sm,
            padding=(10, 5),
            background=C["surface_alt"],
            foreground=C["accent_dark"],
            borderwidth=1,
            relief="flat",
            focusthickness=0,
            bordercolor=C["border"],
            lightcolor=C["border"],
            darkcolor=C["border"],
        )
        style.map(
            "Action.TButton",
            background=[("active", C["toolbar_hover"]), ("pressed", C["accent_soft"])],
            foreground=[("active", C["accent_dark"]), ("pressed", C["accent_dark"])],
            bordercolor=[("active", C["accent"]), ("pressed", C["accent"])],
            relief=[("pressed", "flat"), ("active", "flat")],
        )

        style.configure(
            "PanelAction.TButton",
            font=ui_semi_sm,
            padding=(7, 4),
            background=C["surface_alt"],
            foreground=C["accent_dark"],
            borderwidth=1,
            relief="flat",
            focusthickness=0,
            bordercolor=C["border"],
            lightcolor=C["border"],
            darkcolor=C["border"],
        )
        style.map(
            "PanelAction.TButton",
            background=[("active", C["toolbar_hover"]), ("pressed", C["accent_soft"])],
            foreground=[("active", C["accent_dark"]), ("pressed", C["accent_dark"])],
            bordercolor=[("active", C["accent"]), ("pressed", C["accent"])],
            relief=[("pressed", "flat"), ("active", "flat")],
        )

        # ── Primary button (solid accent fill) ───────────────────────────────────
        style.configure(
            "Primary.TButton",
            font=ui_semi_sm,
            padding=(14, 6),
            background=C["accent"],
            foreground="#ffffff",
            borderwidth=0,
            relief="flat",
            focusthickness=0,
        )
        style.map(
            "Primary.TButton",
            background=[("active", C["accent_dark"]), ("pressed", "#1F6FEB")],
            foreground=[("active", "#ffffff"), ("pressed", "#ffffff")],
            relief=[("pressed", "flat"), ("active", "flat")],
        )

        # ── Entry ────────────────────────────────────────────────────────────────
        style.configure(
            "Modern.TEntry",
            fieldbackground=C["surface_alt"],
            foreground=C["text"],
            insertcolor=C["accent"],
            padding=(8, 5),
            bordercolor=C["border"],
            lightcolor=C["border"],
            darkcolor=C["border"],
        )
        style.map(
            "Modern.TEntry",
            bordercolor=[("focus", C["border_focus"])],
            lightcolor=[("focus", C["border_focus"])],
            darkcolor=[("focus", C["border_focus"])],
        )

        # ── Treeview ─────────────────────────────────────────────────────────────
        style.configure(
            "Modern.Treeview",
            background=C["surface"],
            fieldbackground=C["surface"],
            foreground=C["text"],
            rowheight=28,
            borderwidth=0,
            font=ui_small,
        )
        style.configure(
            "Modern.Treeview.Heading",
            background=C["surface_alt"],
            foreground=C["accent_dark"],
            font=ui_semi_sm,
            relief="flat",
            borderwidth=0,
        )
        style.map(
            "Modern.Treeview",
            background=[("selected", C["accent_soft"])],
            foreground=[("selected", C["accent_dark"])],
        )

        # ── LabelFrame ───────────────────────────────────────────────────────────
        style.configure("Card.TLabelframe",       background=C["surface"], borderwidth=1, relief="groove",
                        bordercolor=C["border"], lightcolor=C["border"], darkcolor=C["border"])
        style.configure("Card.TLabelframe.Label", background=C["surface"], foreground=C["accent_dark"],
                        font=ui_semi_sm)

        # ── Spinbox / Combobox ────────────────────────────────────────────────────
        style.configure("TSpinbox",
                        padding=(6, 4),
                        fieldbackground=C["surface_alt"],
                        foreground=C["text"],
                        bordercolor=C["border"],
                        lightcolor=C["border"],
                        darkcolor=C["border"],
                        arrowcolor=C["muted"],
                        insertcolor=C["accent"])
        style.configure("TCombobox",
                        padding=(6, 4),
                        fieldbackground=C["surface_alt"],
                        foreground=C["text"],
                        selectbackground=C["accent"],
                        selectforeground="#FFFFFF",
                        bordercolor=C["border"],
                        lightcolor=C["border"],
                        darkcolor=C["border"],
                        arrowcolor=C["muted"])
        style.map("TCombobox",
                  fieldbackground=[("readonly", C["surface_alt"])],
                  foreground=[("readonly", C["text"])],
                  selectbackground=[("readonly", C["accent"])],
                  selectforeground=[("readonly", "#FFFFFF")])

        # ── Scrollbars ───────────────────────────────────────────────────────────
        for sb in ("Vertical.TScrollbar", "Horizontal.TScrollbar"):
            style.configure(
                sb,
                gripcount=0,
                background=C["surface_alt"],   # thumb
                troughcolor=C["bg"],            # track
                bordercolor=C["bg"],
                darkcolor=C["surface_alt"],
                lightcolor=C["surface_alt"],
                arrowcolor=C["muted"],
                arrowsize=10,
                relief="flat",
                borderwidth=0,
            )
            style.map(
                sb,
                background=[("active", C["accent"]), ("pressed", C["accent"])],
                arrowcolor=[("active", C["accent_dark"])],
            )

    def _configure_app_icon(self) -> None:
        try:
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("GeneDraft.DNA.Icon")
        except Exception:
            pass

        assets_dir = Path(__file__).resolve().parent / "assets"
        icon_path = assets_dir / "genedraft_app.ico"
        if icon_path.exists():
            try:
                self.root.iconbitmap(default=str(icon_path))
            except Exception:
                pass

        try:
            png_icon_path = assets_dir / "genedraft_app.png"
            if png_icon_path.exists():
                self.app_icon = tk.PhotoImage(file=str(png_icon_path))
            else:
                self.app_icon = self._build_dna_icon()
            self.root.iconphoto(True, self.app_icon)
        except Exception:
            self.app_icon = None

    def _load_brand_mark_image(self) -> tk.PhotoImage | None:
        mark_path = Path(__file__).resolve().parent / "assets" / "genedraft_mark.png"
        if not mark_path.exists():
            return None
        try:
            image = tk.PhotoImage(file=str(mark_path))
            return image.subsample(15, 15)
        except Exception:
            return None

    def _build_dna_icon(self) -> tk.PhotoImage:
        image = tk.PhotoImage(width=64, height=64)

        def fill(x1: int, y1: int, x2: int, y2: int, color: str) -> None:
            image.put(color, to=(x1, y1, x2, y2))

        fill(12, 6, 52, 58, "#0d2137")
        fill(6, 12, 58, 52, "#0d2137")
        fill(14, 10, 50, 54, "#1a6fad")
        fill(10, 14, 54, 50, "#1a6fad")
        fill(16, 14, 48, 50, "#d6eaf8")
        fill(14, 16, 50, 48, "#d6eaf8")

        left_track = [19, 17, 16, 18, 21, 24, 26, 24, 21, 18, 16, 17]
        rung_colors = ("#f59e0b", "#34d399")
        for row, left_x in enumerate(left_track):
            top = 10 + row * 4
            right_x = 64 - left_x - 6
            strand_color = "#0f5ea8" if row % 2 == 0 else "#134e7a"
            fill(left_x, top, left_x + 6, top + 4, strand_color)
            fill(right_x, top, right_x + 6, top + 4, strand_color)
            fill(left_x + 5, top + 1, right_x + 1, top + 3, rung_colors[row % 2])

        fill(26, 18, 38, 22, "#ffffff")
        fill(26, 42, 38, 46, "#ffffff")
        return image

    def _clone_feature(self, feature: FeatureEntry) -> FeatureEntry:
        return FeatureEntry(
            start_nt=feature.start_nt,
            end_nt=feature.end_nt,
            feature_type=feature.feature_type,
            label=feature.label,
            color=feature.color,
            strand=self._normalize_feature_strand(feature.strand),
        )

    def _clone_marker(self, marker: MarkerEntry) -> MarkerEntry:
        return MarkerEntry(position_nt=marker.position_nt, label=marker.label)

    def _create_operation_snapshot(self, description: str) -> OperationSnapshot:
        return OperationSnapshot(
            description=description,
            sequence_name=self.sequence_name.get().strip() or "sequence_1",
            editor_text=self._get_editor_text(),
            features=[self._clone_feature(feature) for feature in self.features],
            markers=[self._clone_marker(marker) for marker in self.markers],
        )

    def _snapshot_signature(self, snapshot: OperationSnapshot) -> tuple[object, ...]:
        return (
            snapshot.sequence_name,
            snapshot.editor_text,
            tuple(
                (feature.start_nt, feature.end_nt, feature.feature_type, feature.label, feature.color, feature.strand)
                for feature in snapshot.features
            ),
            tuple((marker.position_nt, marker.label) for marker in snapshot.markers),
        )

    def _push_operation_snapshot(self, description: str) -> None:
        snapshot = self._create_operation_snapshot(description)
        if self._history_undo and self._snapshot_signature(self._history_undo[-1]) == self._snapshot_signature(snapshot):
            return
        self._history_undo.append(snapshot)
        if len(self._history_undo) > 50:
            self._history_undo = self._history_undo[-50:]
        self._history_redo.clear()

    def _restore_operation_snapshot(self, snapshot: OperationSnapshot) -> None:
        self.sequence_name.set(snapshot.sequence_name or "sequence_1")
        self.features = [self._clone_feature(feature) for feature in snapshot.features]
        self.markers = [self._clone_marker(marker) for marker in snapshot.markers]
        self._reset_search_state()
        self._clear_analysis_highlights()
        self._refresh_feature_table()
        self._refresh_marker_table()
        self._set_editor_text(snapshot.editor_text)

    def undo_biological_operation(self) -> None:
        if not self._history_undo:
            messagebox.showinfo("No history", "There are no biological operations to undo.")
            return
        current_state = self._create_operation_snapshot("current state")
        snapshot = self._history_undo.pop()
        self._history_redo.append(current_state)
        self._restore_operation_snapshot(snapshot)
        self._set_results(f"Operation undone.\n\nRestored action: {snapshot.description}")

    def redo_biological_operation(self) -> None:
        if not self._history_redo:
            messagebox.showinfo("No history", "There are no biological operations to redo.")
            return
        current_state = self._create_operation_snapshot("current state")
        snapshot = self._history_redo.pop()
        self._history_undo.append(current_state)
        self._restore_operation_snapshot(snapshot)
        self._set_results(f"Operation redone.\n\nApplied action: {snapshot.description}")

    def _serialize_session(self) -> dict[str, object]:
        return {
            "sequence_name": self.sequence_name.get().strip() or "sequence_1",
            "editor_text": self._get_editor_text(),
            "has_unsaved_work": self._session_has_unsaved_work,
            "had_persistence_event": self._session_had_persistence_event,
            "features": [
                {
                    "start_nt": feature.start_nt,
                    "end_nt": feature.end_nt,
                    "feature_type": feature.feature_type,
                    "label": feature.label,
                    "color": feature.color,
                    "strand": self._normalize_feature_strand(feature.strand),
                }
                for feature in self.features
            ],
            "markers": [
                {
                    "position_nt": marker.position_nt,
                    "label": marker.label,
                }
                for marker in self.markers
            ],
        }

    def _on_sequence_name_changed(self, *_args) -> None:
        if self._restoring_session:
            return
        self._schedule_session_save()

    def _save_session_now(self) -> None:
        if self._restoring_session:
            return
        data = self._serialize_session()
        try:
            self.session_path.write_text(json.dumps(data, indent=2, ensure_ascii=True), encoding="utf-8")
        except Exception:
            pass

    def _cancel_pending_session_save(self) -> None:
        if self._session_after_id is None:
            return
        try:
            self.root.after_cancel(self._session_after_id)
        except tk.TclError:
            pass
        self._session_after_id = None

    def _mark_session_persisted(self) -> None:
        if self._restoring_session:
            return
        self._session_had_persistence_event = True
        self._session_has_unsaved_work = False
        self._cancel_pending_session_save()
        self._discard_saved_session()

    def _schedule_session_save(self, mark_dirty: bool = True) -> None:
        if self._restoring_session:
            return
        if self._session_had_persistence_event:
            return
        if mark_dirty:
            self._session_has_unsaved_work = True
        self._cancel_pending_session_save()
        self._session_after_id = self.root.after(700, self._run_scheduled_session_save)

    def _run_scheduled_session_save(self) -> None:
        self._session_after_id = None
        self._save_session_now()

    def _schedule_minimap_redraw(self) -> None:
        if not hasattr(self, "minimap_canvas"):
            return
        if self._minimap_after_id is not None:
            try:
                self.root.after_cancel(self._minimap_after_id)
            except tk.TclError:
                pass
        self._minimap_after_id = self.root.after(90, self._run_scheduled_minimap_redraw)

    def _run_scheduled_minimap_redraw(self) -> None:
        self._minimap_after_id = None
        self._redraw_minimap()

    def _discard_saved_session(self) -> None:
        try:
            self.session_path.unlink(missing_ok=True)
        except Exception:
            pass

    def _try_restore_session(self) -> bool:
        if not self.session_path.exists():
            return False
        try:
            data = json.loads(self.session_path.read_text(encoding="utf-8"))
        except Exception:
            return False
        session_dirty = bool(data.get("has_unsaved_work", True))
        had_persistence_event = bool(data.get("had_persistence_event", False))
        editor_text = str(data.get("editor_text", "") or "")
        feature_rows = list(data.get("features", []) or [])
        marker_rows = list(data.get("markers", []) or [])
        if not editor_text and not feature_rows and not marker_rows:
            return False
        if had_persistence_event or not session_dirty:
            self._discard_saved_session()
            return False
        if not self._ask_confirmation_dialog(
            "Recover previous session",
            "An unsaved session from the previous run was found.\n\n"
            "Do you want to recover it now?\n\n"
            "If you choose No, the saved session will be discarded and the app will open clean.",
            subtitle="Session recovery",
        ):
            self._discard_saved_session()
            return False

        self._restoring_session = True
        try:
            self._session_had_persistence_event = had_persistence_event
            self._session_has_unsaved_work = session_dirty
            self._history_undo.clear()
            self._history_redo.clear()
            self.sequence_name.set(str(data.get("sequence_name", "sequence_1") or "sequence_1"))
            self.features = [
                FeatureEntry(
                    start_nt=int(row.get("start_nt", 1)),
                    end_nt=int(row.get("end_nt", 1)),
                    feature_type=str(row.get("feature_type", "misc_feature") or "misc_feature"),
                    label=str(row.get("label", "feature") or "feature"),
                    color=str(row.get("color", next(self._feature_colors)) or next(self._feature_colors)),
                    strand=self._normalize_feature_strand(str(row.get("strand", "+") or "+")),
                )
                for row in feature_rows
            ]
            self.markers = [
                MarkerEntry(
                    position_nt=max(1, int(row.get("position_nt", 1))),
                    label=str(row.get("label", "marker") or "marker"),
                )
                for row in marker_rows
            ]
            self._refresh_feature_table()
            self._refresh_marker_table()
            self._set_editor_text(editor_text)
        finally:
            self._restoring_session = False

        self._set_results(
            "Previous unsaved session recovered.\n\n"
            f"Features: {len(self.features)}\n"
            f"Markers: {len(self.markers)}"
        )
        return True

    def _on_close(self) -> None:
        if self._minimap_after_id is not None:
            try:
                self.root.after_cancel(self._minimap_after_id)
            except tk.TclError:
                pass
        self._cancel_pending_session_save()
        if self._session_had_persistence_event:
            self._discard_saved_session()
        elif self._session_has_unsaved_work:
            self._save_session_now()
        else:
            self._discard_saved_session()
        self.root.destroy()

    # ------------------------------------------------------------------ app settings

    def _load_app_settings(self) -> None:
        if not self.settings_path.exists():
            return
        try:
            data = json.loads(self.settings_path.read_text(encoding="utf-8"))
        except Exception:
            return

        source = self._resolve_config_path(str(data.get("local_blast_source", "") or "").strip())
        db_prefix = self._resolve_config_path(str(data.get("local_blast_db", "") or "").strip())
        db_type = str(data.get("local_blast_type", "auto") or "auto").strip()

        if source:
            self.local_blast_source_var.set(source)
        if db_prefix and self._blast_db_prefix_exists(db_prefix):
            self.local_blast_db_var.set(db_prefix)
        if db_type in {"auto", "nucl", "prot"}:
            self.local_blast_type_var.set(db_type)

        self._recent_files = [p for p in data.get("recent_files", []) if Path(p).exists()]
        self._ncbi_email   = str(data.get("ncbi_email", "") or "").strip()

    def _save_app_settings(self) -> None:
        data = {
            "local_blast_source": self._serialize_config_path(self.local_blast_source_var.get().strip()),
            "local_blast_db": self._serialize_config_path(self.local_blast_db_var.get().strip()),
            "local_blast_type": self.local_blast_type_var.get().strip() or "auto",
            "recent_files": self._recent_files,
            "ncbi_email":   self._ncbi_email,
        }
        try:
            self.settings_path.write_text(
                json.dumps(data, indent=2, ensure_ascii=True),
                encoding="utf-8",
            )
        except Exception:
            pass

    def _resolve_config_path(self, path_text: str) -> str:
        normalized = path_text.strip().strip('"')
        if not normalized:
            return ""
        path = Path(normalized)
        if path.is_absolute():
            return str(path)
        return str((Path.cwd() / path).resolve())

    def _serialize_config_path(self, path_text: str) -> str:
        normalized = path_text.strip().strip('"')
        if not normalized:
            return ""
        path = Path(normalized)
        try:
            resolved = path.resolve()
        except Exception:
            resolved = path
        try:
            return str(resolved.relative_to(Path.cwd().resolve()))
        except Exception:
            return str(resolved)

    def _push_recent_file(self, path: str) -> None:
        path = str(Path(path).resolve())
        self._recent_files = [p for p in self._recent_files if p != path]
        self._recent_files.insert(0, path)
        self._recent_files = self._recent_files[:5]
        self._rebuild_recent_menu()
        self._save_app_settings()

    def _rebuild_recent_menu(self) -> None:
        self._recent_menu.delete(0, "end")
        if not self._recent_files:
            self._recent_menu.add_command(label="(empty)", state="disabled")
            return
        for path in self._recent_files:
            label = Path(path).name
            self._recent_menu.add_command(
                label=label,
                command=lambda p=path: self._open_recent(p),
            )

    def _open_recent(self, path: str) -> None:
        if not Path(path).exists():
            messagebox.showerror("File not found", f"Could not find:\n{path}")
            self._recent_files = [p for p in self._recent_files if p != path]
            self._rebuild_recent_menu()
            self._save_app_settings()
            return
        try:
            header, sequence, loaded_features, file_format = load_sequence_file(path)
        except Exception as exc:
            messagebox.showerror("Could not open", str(exc))
            return
        self._push_operation_snapshot("open file")
        self.sequence_name.set(header)
        self.features = [
            FeatureEntry(
                start_nt=feature.start_nt,
                end_nt=feature.end_nt,
                feature_type=feature.feature_type,
                label=feature.label,
                color=next(self._feature_colors),
                strand=self._normalize_feature_strand(feature.strand),
            )
            for feature in loaded_features
        ]
        self.markers = []
        self._reset_search_state()
        self._refresh_feature_table()
        self._refresh_marker_table()
        self._clear_analysis_highlights()
        self._set_editor_text(self._format_for_editor(sequence))
        self._set_results(
            "File loaded.\n\n"
            f"Name: {Path(path).name}\n"
            f"Format: {file_format}\n"
            f"Imported features: {len(self.features)}"
        )
        self._push_recent_file(path)

    def _normalize_blast_db_prefix(self, path_text: str) -> str:
        normalized = path_text.strip().strip('"')
        if not normalized:
            return ""
        path = Path(normalized)
        if path.suffix.lower() in {".pin", ".psq", ".phr", ".pdb", ".nin", ".nsq", ".nhr", ".ndb"}:
            path = path.with_suffix("")
        return str(path)

    def _ensure_blast_paths_supported(self, *paths: str) -> bool:
        checked = [str(Path.cwd().resolve())]
        checked.extend(path.strip() for path in paths if path and path.strip())
        offenders = [path for path in checked if " " in path]
        if not offenders:
            return True

        details = "\n".join(f"- {path}" for path in offenders)
        messagebox.showerror(
            "BLAST path issue",
            "Local BLAST on Windows does not work reliably when a path contains spaces.\n\n"
            "Please rename or move the project folder, FASTA file, or BLAST database path and try again.\n\n"
            "Recommended: C:\\GeneDraft\n"
            "Avoid: C:\\My Projects\\GeneDraft\n\n"
            f"Paths with spaces detected:\n{details}"
        )
        return False

    def _blast_db_prefix_exists(self, db_prefix: str) -> bool:
        prefix = Path(self._normalize_blast_db_prefix(db_prefix))
        return any(
            prefix.with_suffix(ext).exists()
            for ext in (".pin", ".psq", ".phr", ".pdb", ".nin", ".nsq", ".nhr", ".ndb")
        )

    def _discover_local_blast_db_prefixes(self) -> list[tuple[str, str]]:
        blast_root = Path.cwd() / "blast_db"
        if not blast_root.exists():
            return []

        found: dict[str, str] = {}
        for pattern, db_type in (("*.pin", "prot"), ("*.nin", "nucl")):
            for file_path in blast_root.rglob(pattern):
                prefix = str(file_path.with_suffix(""))
                found[prefix] = db_type
        return sorted(found.items(), key=lambda item: item[0].lower())

    def _apply_local_blast_db_config(self, db_prefix: str, db_type: str, source: str = "") -> None:
        db_prefix = self._normalize_blast_db_prefix(db_prefix)
        self.local_blast_db_var.set(db_prefix)
        self.local_blast_type_var.set(db_type if db_type in {"auto", "nucl", "prot"} else "auto")
        self.local_blast_source_var.set(source.strip())
        self._save_app_settings()

    # ------------------------------------------------------------------ menu bar

    def _build_menu(self) -> None:
        C = self.colors
        menubar = tk.Menu(self.root, bg=C["surface"], fg=C["text"],
                          activebackground=C["accent_soft"], activeforeground=C["accent_dark"],
                          relief="flat", borderwidth=0)
        self.root.configure(menu=menubar)
        self._menubar = menubar
        self._all_menus.append(menubar)

        def menu(label):
            m = tk.Menu(menubar, tearoff=0, bg=C["surface"], fg=C["text"],
                        activebackground=C["accent_soft"], activeforeground=C["accent_dark"],
                        font=("Segoe UI Variable Text", 11))
            menubar.add_cascade(label=label, menu=m)
            self._all_menus.append(m)
            return m

        m = menu("File")
        m.add_command(label="New Window             Ctrl+N", command=self.new_window)
        m.add_separator()
        m.add_command(label="Open...                Ctrl+O", command=self.open_sequence_file)
        m.add_command(label="Fetch by Accession...  Ctrl+Shift+O", command=self.fetch_accession_dialog)
        self._recent_menu = tk.Menu(m, tearoff=0, bg=C["surface"], fg=C["text"],
                                    activebackground=C["accent_soft"], activeforeground=C["accent_dark"],
                                    font=("Segoe UI Variable Text", 11))
        self._all_menus.append(self._recent_menu)
        m.add_cascade(label="Open Recent", menu=self._recent_menu)
        m.add_separator()
        m.add_command(label="Save FASTA...          Ctrl+S", command=self.save_fasta_file)
        m.add_command(label="Save GenBank...        Ctrl+Shift+S", command=self.save_genbank_file)
        m.add_command(label="Copy FASTA             Ctrl+Shift+C", command=self.copy_fasta_to_clipboard)
        m.add_separator()
        m.add_command(label="Clear All              Ctrl+Shift+L", command=self.clear_all)

        m = menu("Edit")
        m.add_command(label="Undo Operation         Ctrl+Alt+Z", command=self.undo_biological_operation)
        m.add_command(label="Redo Operation         Ctrl+Alt+Y", command=self.redo_biological_operation)

        m = menu("Sequence")
        m.add_command(label="Normalize              Ctrl+L", command=self.normalize_sequence)
        m.add_command(label="Validate               F5", command=self.validate_sequence)
        m.add_command(label="Go To Position         Ctrl+G", command=self.go_to_position)
        m.add_separator()
        m.add_command(label="Reverse Complement     Ctrl+R", command=self.replace_with_reverse_complement)
        m.add_command(label="Reverse Complement Selection", command=self.reverse_complement_selection)
        m.add_command(label="Convert DNA To RNA", command=self.convert_dna_to_rna)
        m.add_separator()
        m.add_command(label="Translate              Ctrl+T", command=self.translate_selection)
        m.add_separator()
        m.add_command(label="Add Marker Here", command=self.add_marker_at_cursor)
        m.add_command(label="Copy Selection Coordinates", command=self.copy_selection_coordinates)

        m = menu("Analysis")
        m.add_command(label="Find Motif             Ctrl+F", command=self.search_motif)
        m.add_command(label="Find ORFs              Ctrl+Shift+F", command=self.show_orfs)
        m.add_command(label="Secondary Structure    Ctrl+Alt+S", command=self.show_secondary_structure_analysis)
        m.add_command(label="Protein Properties     Ctrl+P", command=self.show_protein_analysis)
        m.add_separator()
        m.add_command(label="Circularize Selection", command=self.circularize_selection)
        m.add_command(label="Mark Component", command=self.mark_current_component)
        m.add_command(label="Crop Component", command=self.crop_to_current_component)
        m.add_separator()
        m.add_command(label="Clear Highlights       Esc", command=self.clear_analysis_marks)

        m = menu("Features")
        m.add_command(label="Add Feature            Ctrl+E", command=self.add_feature)
        m.add_command(label="Edit Feature", command=self.edit_selected_feature)
        m.add_command(label="Delete Feature", command=self.delete_selected_feature)
        m.add_command(label="Focus Feature", command=self.focus_selected_feature)
        m.add_command(label="Feature Summary", command=self.show_selected_feature_summary)
        m.add_command(label="Translate Feature/CDS", command=self.translate_selected_feature)
        m.add_command(label="Copy Feature Coordinates", command=self.copy_selected_feature_coordinates)

        m = menu("Markers")
        m.add_command(label="Add Marker Here", command=self.add_marker_at_cursor)
        m.add_command(label="Add Marker From Selection", command=self.add_marker_from_selection)
        m.add_command(label="Go To Marker", command=self.focus_selected_marker)
        m.add_command(label="Delete Marker", command=self.delete_selected_marker)

        m = menu("Export")
        m.add_command(label="Linear Map SVG...", command=self.export_linear_map)
        m.add_command(label="Circular Map SVG...", command=self.export_circular_map)
        m.add_separator()
        m.add_command(label="Selection To FASTA...", command=self.save_selected_fasta_file)
        m.add_command(label="Selection To GenBank...", command=self.save_selected_genbank_file)

        m = menu("BLAST")
        m.add_command(label="BLAST Web (Editor)", command=self.run_blast_web_from_selection)
        m.add_command(label="BLAST Web (Results)", command=self.run_blast_web_from_results)
        m.add_separator()
        m.add_command(label="BLAST Local (Editor)", command=self.run_blast_local_from_selection)
        m.add_command(label="BLAST Local (Results)", command=self.run_blast_local_from_results)
        m.add_separator()
        m.add_command(label="Configure BLAST DB...", command=self.configure_local_blast_db)
        m.add_command(label="Build BLAST DB...", command=self.build_local_blast_db)
        m.add_command(label="BLAST DB Info", command=self.show_local_blast_db_info)

        m = menu("About")
        m.add_command(label="Information", command=self.show_about_information)
        m.add_command(label="Help", command=self.show_about_help)
        m.add_command(label="Shortcuts", command=self.show_about_shortcuts)
        m.add_separator()
        m.add_command(label="Version", command=self.show_about_version)
        m.add_command(label="Citation", command=self.show_about_citation)
        m.add_command(label="License and Use", command=self.show_about_license)
        m.add_command(label="Acknowledgments", command=self.show_about_acknowledgments)
        m.add_command(label="Contact", command=self.show_about_contact)

    # ------------------------------------------------------------------ toolbar helpers

    def _toolbar_btn(self, parent, text, command, tooltip="", width=None):
        """Rounded toolbar button."""
        C = self.colors
        btn = _RoundBtn(
            parent, text, command,
            normal_bg=C["toolbar_btn"],
            hover_bg=C["toolbar_hover"],
            press_bg=C["toolbar_press"],
            fg=C["toolbar_fg"],
            radius=6, padx=11, pady=5,
            font=("Segoe UI Variable Text Semibold", 9),
        )
        self._toolbar_btns.append(btn)
        return btn

    def _toolbar_sep(self, parent):
        """Vertical separator for the toolbar."""
        C = self.colors
        sep = tk.Frame(parent, width=1, bg=C["toolbar_sep"])
        sep.pack(side="left", fill="y", padx=6, pady=4)

    # ------------------------------------------------------------------ theme switching

    def _toggle_theme(self) -> None:
        new_theme = "day" if self.current_theme == "night" else "night"
        old_colors = self.colors
        new_colors = self._build_day_colors() if new_theme == "day" else self._build_night_colors()

        # Build old-hex → new-hex mapping (exact lowercase keys)
        color_map: dict[str, str] = {}
        for key in old_colors:
            if key in new_colors:
                ov = old_colors[key].lower()
                nv = new_colors[key].lower()
                if ov != nv:
                    color_map[ov] = nv

        self.current_theme = new_theme
        self.colors = new_colors

        # Re-apply ttk styles
        self._configure_styles()

        # Walk the full widget tree
        self._remap_widget_tree(self.root, color_map)

        # Update menus
        C = new_colors
        menu_kw = dict(bg=C["surface"], fg=C["text"],
                       activebackground=C["accent_soft"], activeforeground=C["accent_dark"])
        for m in self._all_menus:
            try:
                m.configure(**menu_kw)
            except Exception:
                pass

        # Update toolbar buttons
        for btn in self._toolbar_btns:
            btn._normal = C["toolbar_btn"]
            btn._hover  = C["toolbar_hover"]
            btn._press  = C["toolbar_press"]
            btn._fg     = C["toolbar_fg"]
            btn.configure(
                bg=C["toolbar_btn"],
                fg=C["toolbar_fg"],
                activebackground=C["toolbar_hover"],
                activeforeground=C["text"],
                highlightbackground=C["toolbar_btn"],
            )

        # Update theme button label
        self._theme_btn.configure(text="Day" if new_theme == "night" else "Night")

        # Update Text widget highlight tags
        self._reconfigure_text_tags()

        # Redraw minimap with new colors
        self._schedule_minimap_redraw()

    def _remap_widget_tree(self, widget: tk.BaseWidget, color_map: dict[str, str]) -> None:
        """Recursively remap colors on all tk (non-ttk) widgets."""
        props = ("bg", "fg", "background", "foreground",
                 "insertbackground", "selectbackground", "selectforeground",
                 "highlightbackground", "highlightcolor",
                 "activebackground", "activeforeground", "troughcolor")
        for prop in props:
            try:
                val = widget.cget(prop)
                if isinstance(val, str):
                    mapped = color_map.get(val.lower())
                    if mapped:
                        widget.configure(**{prop: mapped})
            except Exception:
                pass
        for child in widget.winfo_children():
            self._remap_widget_tree(child, color_map)

    def _reconfigure_text_tags(self) -> None:
        """(Re-)apply editor highlight tag colors from current palette."""
        C = self.colors
        is_night = self.current_theme == "night"
        bright_fg = "#E2E8F0" if is_night else C["text"]
        self.sequence_text.tag_configure("motif",               background=C["tag_motif"])
        self.sequence_text.tag_configure("motif_focus",         background=C["tag_orf_focus"])
        self.sequence_text.tag_configure("invalid",             background=C["tag_invalid"])
        self.sequence_text.tag_configure("orf_plus",            background=C["tag_orf_plus"],  foreground=bright_fg)
        self.sequence_text.tag_configure("orf_minus",           background=C["tag_orf_minus"], foreground=bright_fg)
        orf_focus_fg = "#FFFFFF" if is_night else C["text"]
        self.sequence_text.tag_configure("orf_focus",           background=C["tag_orf_focus"], foreground=orf_focus_fg)
        self.sequence_text.tag_configure("circular_candidate",  background=C["tag_circ"])
        self.sequence_text.tag_configure("circular_focus",      background=C["tag_circ_foc"])
        self.sequence_text.tag_configure("feature_focus",       background=C["tag_feature_focus"], foreground=C["text"])
        # Keep the built-in selection tag on top so it is always visible over other tags
        self.sequence_text.tag_configure("sel",
                                         background=C["selection_bg"],
                                         foreground=C["selection_fg"])
        self.sequence_text.tag_raise("sel")

    # ------------------------------------------------------------------ main UI

    def _build_ui(self) -> None:
        C = self.colors
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(2, weight=1)

        # Toolbar
        toolbar = tk.Frame(self.root, bg=C["toolbar_bg"], height=40)
        toolbar.grid(row=0, column=0, sticky="ew", padx=18, pady=(0, 1))
        toolbar.columnconfigure(100, weight=1)   # spacer

        def tb(text, cmd, tip=""):
            btn = self._toolbar_btn(toolbar, text, cmd, tip)
            btn.pack(side="left", padx=1, pady=3)
            return btn

        # File and history
        tb("Open", self.open_sequence_file)
        tb("Save FASTA", self.save_fasta_file)
        tb("< Undo", self.undo_biological_operation)
        tb("Redo >", self.redo_biological_operation)
        self._toolbar_sep(toolbar)

        # Core editing
        tb("Normalize", self.normalize_sequence)
        tb("Translate", self.translate_selection)
        self._toolbar_sep(toolbar)

        # Analysis
        tb("Find Motif", self.search_motif)
        tb("Find ORFs", self.show_orfs)
        tb("RNA Fold", self.show_secondary_structure_analysis)
        tb("Protein Analysis", self.show_protein_analysis)
        tb("Circularize", self.circularize_selection)
        self._toolbar_sep(toolbar)

        # Features
        tb("Add Feature", self.add_feature)
        self._toolbar_sep(toolbar)

        # BLAST
        tb("BLAST web",     self.run_blast_web_from_selection)
        tb("BLAST local",   self.run_blast_local_from_selection)
        self._toolbar_sep(toolbar)
        tb("Clear All", self.clear_all)

        # Theme toggle — right side
        self._theme_btn = self._toolbar_btn(toolbar, "Day", self._toggle_theme)
        self._theme_btn.pack(side="right", padx=(0, 8), pady=3)
        self._toolbar_frame = toolbar

        # Settings strip
        settings_bar = ttk.Frame(self.root, style="Card.TFrame", padding=(10, 6))
        settings_bar.grid(row=1, column=0, sticky="ew")  # placeholder; placed below toolbar via pack later

        # Use a secondary frame approach: toolbar + settings in a column
        # Reset grid rows so toolbar=0, settings=1, content=2, status=3/4
        toolbar.grid_forget()
        settings_bar.grid_forget()

        self.root.rowconfigure(0, weight=0)
        self.root.rowconfigure(1, weight=0)
        self.root.rowconfigure(2, weight=1)
        self.root.rowconfigure(3, weight=0)
        self.root.rowconfigure(4, weight=0)

        toolbar.grid(row=0, column=0, sticky="ew", padx=18)

        # Settings strip
        settings_strip = tk.Frame(self.root, bg=C["settings_bg"], padx=12, pady=5)
        settings_strip.grid(row=1, column=0, sticky="ew", padx=18, pady=(2, 6))

        col = 0
        def slabel(text):
            nonlocal col
            ttk.Label(settings_strip, text=text, style="Settings.TLabel").grid(row=0, column=col, sticky="w", padx=(0 if col == 0 else 4, 2))
            col += 1

        def swidget(w):
            nonlocal col
            w.grid(row=0, column=col, sticky="w", padx=(0, 10))
            col += 1

        slabel("Min ORF (aa)")
        self.orf_min_spinbox = ttk.Spinbox(settings_strip, from_=1, to=5000, textvariable=self.orf_min_var, width=6)
        swidget(self.orf_min_spinbox)

        sep_frm = ttk.Frame(settings_strip, style="Sep.TFrame", width=1)
        sep_frm.grid(row=0, column=col, sticky="ns", pady=2, padx=4); col += 1

        slabel("Circular overlap")
        swidget(ttk.Spinbox(settings_strip, from_=10, to=500, textvariable=self.circular_overlap_var, width=6))

        slabel("Circular mismatches")
        swidget(ttk.Spinbox(settings_strip, from_=0, to=5, textvariable=self.circular_mismatch_var, width=5))

        slabel("Circular scan")
        swidget(ttk.Spinbox(settings_strip, from_=20, to=2000, textvariable=self.circular_scan_var, width=6))

        sep_frm2 = ttk.Frame(settings_strip, style="Sep.TFrame", width=1)
        sep_frm2.grid(row=0, column=col, sticky="ns", pady=2, padx=4); col += 1

        slabel("FASTA width")
        swidget(ttk.Spinbox(settings_strip, from_=40, to=200, textvariable=self.line_width_var, width=6))

        # spacer
        settings_strip.columnconfigure(col, weight=1)

        # Main content
        self.content = ttk.Panedwindow(self.root, orient=tk.HORIZONTAL)
        self.content.grid(row=2, column=0, padx=18, pady=(0, 0), sticky="nsew")

        # Editor panel
        editor_panel = ttk.Frame(self.content, style="Card.TFrame", padding=(16, 14, 12, 14))
        editor_panel.columnconfigure(0, weight=1)
        editor_panel.rowconfigure(2, weight=1)

        editor_head = ttk.Frame(editor_panel, style="Card.TFrame")
        editor_head.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        editor_head.columnconfigure(2, weight=1)
        ttk.Label(editor_head, text="Sequence editor", style="Card.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(editor_head, text="Name", style="Muted.TLabel").grid(row=0, column=1, sticky="e", padx=(14, 6))
        ttk.Entry(editor_head, textvariable=self.sequence_name, style="Modern.TEntry", width=24).grid(row=0, column=2, sticky="ew")
        editor_nav = ttk.Frame(editor_head, style="Card.TFrame")
        editor_nav.grid(row=0, column=3, sticky="e", padx=(12, 0))
        ttk.Button(editor_nav, text="< Previous", command=self.show_previous_search_hit, style="Action.TButton").grid(row=0, column=0, padx=(0, 6))
        ttk.Button(editor_nav, text="Next >", command=self.show_next_search_hit, style="Action.TButton").grid(row=0, column=1)

        self.minimap_canvas = tk.Canvas(
            editor_panel,
            height=66,
            bg=self.colors["surface_alt"],
            highlightthickness=1,
            highlightbackground=self.colors["border"],
            relief="flat",
            cursor="hand2",
        )
        self.minimap_canvas.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self.minimap_canvas.bind("<Button-1>", self._on_minimap_click)
        self.minimap_canvas.bind("<B1-Motion>", self._on_minimap_click)
        self.minimap_canvas.bind("<Configure>", lambda _event: self._redraw_minimap())

        editor_frame = ttk.Frame(editor_panel)
        editor_frame.grid(row=2, column=0, sticky="nsew")
        editor_frame.columnconfigure(0, weight=1)
        editor_frame.rowconfigure(0, weight=1)

        self.sequence_text = tk.Text(
            editor_frame,
            wrap="none",
            undo=True,
            font=("Cascadia Code", self.editor_font_size),
            bg=self.colors["surface"],
            fg=self.colors["text"],
            insertbackground=self.colors["accent"],
            selectbackground=self.colors["selection_bg"],
            selectforeground=self.colors["selection_fg"],
            relief="flat",
            highlightthickness=1,
            highlightbackground=self.colors["border"],
            highlightcolor=self.colors["accent"],
            padx=12,
            pady=12,
            insertwidth=2,
        )
        ed_scroll_y = ttk.Scrollbar(editor_frame, orient="vertical",   command=self.sequence_text.yview)
        ed_scroll_x = ttk.Scrollbar(editor_frame, orient="horizontal", command=self.sequence_text.xview)
        self.sequence_text.configure(yscrollcommand=ed_scroll_y.set, xscrollcommand=ed_scroll_x.set)
        self.sequence_text.grid(row=0, column=0, sticky="nsew")
        ed_scroll_y.grid(row=0, column=1, sticky="ns")
        ed_scroll_x.grid(row=1, column=0, sticky="ew")

        self._reconfigure_text_tags()

        self.sequence_text.bind("<<Modified>>",      self._on_text_modified)
        self.sequence_text.bind("<ButtonRelease-1>", self._update_position_info)
        self.sequence_text.bind("<KeyRelease>",      self._update_position_info)
        self.sequence_text.bind("<Control-v>", lambda _event: (self.paste_into_sequence(), "break")[1])
        self.sequence_text.bind("<Control-V>", lambda _event: (self.paste_into_sequence(), "break")[1])

        # Side panel
        side_panel = ttk.Frame(self.content, style="App.TFrame", padding=(10, 0, 0, 0))
        side_panel.columnconfigure(0, weight=1)
        side_panel.rowconfigure(0, weight=1)
        self.side_stack = ttk.Panedwindow(side_panel, orient=tk.VERTICAL)
        self.side_stack.grid(row=0, column=0, sticky="nsew")

        # Results
        results_panel = ttk.Frame(self.side_stack, style="Card.TFrame", padding=(14, 12, 12, 12))
        results_panel.columnconfigure(0, weight=1)
        results_panel.rowconfigure(1, weight=1)

        res_head = ttk.Frame(results_panel, style="Card.TFrame")
        res_head.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        res_head.columnconfigure(0, weight=1)
        ttk.Label(res_head, text="Results", style="Card.TLabel").grid(row=0, column=0, sticky="w")
        res_btns = ttk.Frame(res_head, style="Card.TFrame")
        res_btns.grid(row=0, column=1, sticky="e")
        ttk.Button(res_btns, text="BLAST web",   command=self.run_blast_web_from_results,   style="Action.TButton").grid(row=0, column=0, padx=(0, 4))
        ttk.Button(res_btns, text="BLAST local", command=self.run_blast_local_from_results, style="Action.TButton").grid(row=0, column=1)

        results_frame = ttk.Frame(results_panel)
        results_frame.grid(row=1, column=0, sticky="nsew")
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)

        self.results_text = tk.Text(
            results_frame,
            wrap="word",
            font=("Cascadia Code", self.results_font_size),
            bg=self.colors["surface"],
            fg=self.colors["text"],
            selectbackground=self.colors["selection_bg"],
            selectforeground=self.colors["selection_fg"],
            relief="flat",
            highlightthickness=1,
            highlightbackground=self.colors["border"],
            highlightcolor=self.colors["accent"],
            padx=12,
            pady=12,
            state="disabled",
        )
        res_scroll = ttk.Scrollbar(results_frame, orient="vertical", command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=res_scroll.set)
        self.results_text.grid(row=0, column=0, sticky="nsew")
        res_scroll.grid(row=0, column=1, sticky="ns")

        # Features table
        feat_frame = ttk.LabelFrame(self.side_stack, text="Features", style="Card.TLabelframe", padding=(10, 8))
        feat_frame.columnconfigure(0, weight=1)
        feat_frame.rowconfigure(1, weight=1)

        feat_toolbar = ttk.Frame(feat_frame, style="Card.TFrame")
        feat_toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 4))
        for i in range(4):
            feat_toolbar.columnconfigure(i, weight=1)
        ttk.Button(feat_toolbar, text="Add", command=self.add_feature, style="PanelAction.TButton").grid(row=0, column=0, padx=2, sticky="ew")
        ttk.Button(feat_toolbar, text="Edit", command=self.edit_selected_feature, style="PanelAction.TButton").grid(row=0, column=1, padx=2, sticky="ew")
        ttk.Button(feat_toolbar, text="Delete", command=self.delete_selected_feature, style="PanelAction.TButton").grid(row=0, column=2, padx=2, sticky="ew")
        ttk.Button(feat_toolbar, text="Focus", command=self.focus_selected_feature, style="PanelAction.TButton").grid(row=0, column=3, padx=2, sticky="ew")

        columns = ("type", "start", "end", "strand", "label")
        self.feature_table = ttk.Treeview(feat_frame, columns=columns, show="headings", height=11, style="Modern.Treeview")
        self.feature_table.heading("type",  text="Type")
        self.feature_table.heading("start", text="Start")
        self.feature_table.heading("end",   text="End")
        self.feature_table.heading("strand", text="Strand")
        self.feature_table.heading("label", text="Label")
        self.feature_table.column("type",  width=110, anchor="w")
        self.feature_table.column("start", width=72,  anchor="center")
        self.feature_table.column("end",   width=72,  anchor="center")
        self.feature_table.column("strand", width=78, anchor="center")
        self.feature_table.column("label", width=180, anchor="w")

        feat_scroll = ttk.Scrollbar(feat_frame, orient="vertical", command=self.feature_table.yview)
        self.feature_table.configure(yscrollcommand=feat_scroll.set)
        self.feature_table.grid(row=1, column=0, sticky="nsew")
        feat_scroll.grid(row=1, column=1, sticky="ns")
        self.feature_table.bind("<<TreeviewSelect>>", self._on_feature_table_select)
        self.feature_table.bind("<Button-1>", self._on_feature_table_click)
        self.feature_table.bind("<Double-Button-1>", self._on_feature_table_click)

        marker_frame = ttk.LabelFrame(self.side_stack, text="Markers", style="Card.TLabelframe", padding=(10, 8))
        marker_frame.columnconfigure(0, weight=1)
        marker_frame.rowconfigure(1, weight=1)

        marker_toolbar = ttk.Frame(marker_frame, style="Card.TFrame")
        marker_toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 4))
        for i in range(3):
            marker_toolbar.columnconfigure(i, weight=1)
        ttk.Button(marker_toolbar, text="Add", command=self.add_marker_at_cursor, style="PanelAction.TButton").grid(row=0, column=0, padx=2, sticky="ew")
        ttk.Button(marker_toolbar, text="Go", command=self.focus_selected_marker, style="PanelAction.TButton").grid(row=0, column=1, padx=2, sticky="ew")
        ttk.Button(marker_toolbar, text="Delete", command=self.delete_selected_marker, style="PanelAction.TButton").grid(row=0, column=2, padx=2, sticky="ew")

        self.marker_table = ttk.Treeview(marker_frame, columns=("pos", "label"), show="headings", height=6, style="Modern.Treeview")
        self.marker_table.heading("pos", text="Pos")
        self.marker_table.heading("label", text="Label")
        self.marker_table.column("pos", width=80, anchor="center")
        self.marker_table.column("label", width=220, anchor="w")
        marker_scroll = ttk.Scrollbar(marker_frame, orient="vertical", command=self.marker_table.yview)
        self.marker_table.configure(yscrollcommand=marker_scroll.set)
        self.marker_table.grid(row=1, column=0, sticky="nsew")
        marker_scroll.grid(row=1, column=1, sticky="ns")
        self.marker_table.bind("<<TreeviewSelect>>", self._on_marker_table_select)

        self.side_stack.add(results_panel, weight=2)
        self.side_stack.add(feat_frame, weight=5)
        self.side_stack.add(marker_frame, weight=3)

        self.content.add(editor_panel, weight=8)
        self.content.add(side_panel, weight=3)

        # Status bar
        status_bar = tk.Frame(self.root, bg=C["accent"], height=2)
        status_bar.grid(row=3, column=0, sticky="ew", padx=18)

        status = tk.Frame(self.root, bg=C["surface_alt"], padx=14, pady=7)
        status.grid(row=4, column=0, sticky="ew", padx=18, pady=(0, 18))

        status_items = [
            self.length_var,
            self.type_var,
            self.gc_var,
            self.selection_var,
            self.features_var,
            self.cursor_var,
        ]
        for i, var in enumerate(status_items):
            if i > 0:
                tk.Label(status, text="|", bg=C["surface_alt"],
                         fg=C["border"], font=("Segoe UI Variable Text", 12)).grid(row=0, column=i * 2 - 1, padx=7)
            tk.Label(status, textvariable=var, bg=C["surface_alt"],
                     fg=C["muted"], font=("Segoe UI Variable Text Semibold", 10)).grid(row=0, column=i * 2, sticky="w")

        badge_column = len(status_items) * 2
        status.columnconfigure(badge_column, weight=1)
        self.status_warning_label = tk.Label(
            status,
            textvariable=self.warning_var,
            bg=C["surface"],
            fg=C["muted"],
            font=("Segoe UI Variable Text Semibold", 10),
            padx=12,
            pady=4,
            bd=0,
            relief="flat",
            highlightthickness=1,
            highlightbackground=C["border"],
            highlightcolor=C["border"],
        )
        self.status_warning_label.grid(row=0, column=badge_column + 1, sticky="e")
        self._apply_warning_badge_style("idle")

    def _apply_initial_layout(self, _retry: int = 0) -> None:
        if not hasattr(self, "content") or not hasattr(self, "side_stack"):
            return
        self.root.update_idletasks()
        ready = True
        try:
            total_width = self.content.winfo_width()
            if total_width > 10:
                self.content.sashpos(0, int(total_width * 0.75))
            else:
                ready = False
        except tk.TclError:
            ready = False
        try:
            total_height = self.side_stack.winfo_height()
            if total_height > 10:
                self.side_stack.sashpos(0, int(total_height * 0.20))
                self.side_stack.sashpos(1, int(total_height * 0.72))
            else:
                ready = False
        except tk.TclError:
            ready = False
        if not ready and _retry < 8:
            self.root.after(80, lambda: self._apply_initial_layout(_retry + 1))

    # ------------------------------------------------------------------ context menus

    def _build_context_menus(self) -> None:
        C = self.colors
        cm_kw = dict(tearoff=0, bg=C["surface"], fg=C["text"],
                     activebackground=C["accent_soft"], activeforeground=C["accent_dark"],
                     font=("Segoe UI Variable Text", 11))

        self.sequence_context_menu = tk.Menu(self.root, **cm_kw)
        self.sequence_context_menu.add_command(label="Copy",                              command=self.copy_selected_sequence)
        self.sequence_context_menu.add_command(label="Cut",                              command=self.cut_selected_sequence)
        self.sequence_context_menu.add_command(label="Paste",                               command=self.paste_into_sequence)
        self.sequence_context_menu.add_separator()
        self.sequence_context_menu.add_command(label="Find motif", command=self.search_motif)
        self.sequence_context_menu.add_command(label="Find ORFs", command=self.show_orfs)
        self.sequence_context_menu.add_command(label="Secondary structure", command=self.show_secondary_structure_analysis)
        self.sequence_context_menu.add_command(label="Circularize selection", command=self.circularize_selection)
        self.sequence_context_menu.add_command(label="BLAST web",                           command=self.run_blast_web_from_selection)
        self.sequence_context_menu.add_command(label="BLAST local",                         command=self.run_blast_local_from_selection)
        self.sequence_context_menu.add_command(label="Translate", command=self.translate_selection)
        self.sequence_context_menu.add_command(label="Reverse complement selection", command=self.reverse_complement_selection)
        self.sequence_context_menu.add_command(label="Export selection FASTA", command=self.save_selected_fasta_file)
        self.sequence_context_menu.add_command(label="Export selection GenBank", command=self.save_selected_genbank_file)
        self.sequence_context_menu.add_command(label="Copy selection coords",             command=self.copy_selection_coordinates)
        self.sequence_context_menu.add_command(label="Add marker here", command=self.add_marker_at_cursor)
        self.sequence_context_menu.add_command(label="Create feature from selection", command=self.add_feature)
        self.sequence_context_menu.add_command(label="Clear highlights",              command=self.clear_analysis_marks)
        self.sequence_context_menu.add_separator()
        self.sequence_context_menu.add_command(label="Select all",                    command=lambda: self._select_all_in_widget(self.sequence_text))

        self.results_context_menu = tk.Menu(self.root, **cm_kw)
        self.results_context_menu.add_command(label="Copy",          command=lambda: self._copy_from_widget(self.results_text))
        self.results_context_menu.add_command(
            label="Export secondary structure SVG...",
            command=self.export_secondary_structure_svg_from_results,
        )
        self._results_secondary_svg_menu_index = self.results_context_menu.index("end")
        self.results_context_menu.add_separator()
        self.results_context_menu.add_command(label="BLAST web",       command=self.run_blast_web_from_results)
        self.results_context_menu.add_command(label="BLAST local",       command=self.run_blast_local_from_results)
        self.results_context_menu.add_command(label="Protein Analysis", command=self.protein_analysis_from_results)
        self.results_context_menu.add_command(label="Select all",       command=lambda: self._select_all_in_widget(self.results_text))

        self.feature_context_menu = tk.Menu(self.root, **cm_kw)
        self.feature_context_menu.add_command(label="Edit feature",   command=self.edit_selected_feature)
        self.feature_context_menu.add_command(label="Focus feature",  command=self.focus_selected_feature)
        self.feature_context_menu.add_command(label="Feature summary",  command=self.show_selected_feature_summary)
        self.feature_context_menu.add_command(label="Translate feature/CDS", command=self.translate_selected_feature)
        self.feature_context_menu.add_command(label="Copy coordinates", command=self.copy_selected_feature_coordinates)
        self.feature_context_menu.add_command(label="Delete feature", command=self.delete_selected_feature)
        self.feature_context_menu.add_separator()
        self.feature_context_menu.add_command(label="Add feature",  command=self.add_feature)

        self.marker_context_menu = tk.Menu(self.root, **cm_kw)
        self.marker_context_menu.add_command(label="Go to marker",     command=self.focus_selected_marker)
        self.marker_context_menu.add_command(label="Copy coordinates", command=self.copy_selected_marker_coordinates)
        self.marker_context_menu.add_command(label="Delete marker", command=self.delete_selected_marker)
        self.marker_context_menu.add_separator()
        self.marker_context_menu.add_command(label="Add marker here", command=self.add_marker_at_cursor)

        self.sequence_text.bind("<Button-3>", self._show_sequence_context_menu)
        self.results_text.bind("<Button-3>",  self._show_results_context_menu)
        self.feature_table.bind("<Button-3>", self._show_feature_context_menu)
        self.feature_table.bind("<Double-1>", self._handle_feature_table_double_click)
        self.marker_table.bind("<Button-3>", self._show_marker_context_menu)

        self._all_menus.extend([
            self.sequence_context_menu,
            self.results_context_menu,
            self.feature_context_menu,
            self.marker_context_menu,
        ])

    # ------------------------------------------------------------------ shortcuts

    def _bind_shortcuts(self) -> None:
        self._bind_action("<Control-n>",         self.new_window)
        self._bind_action("<Control-o>",         self.open_sequence_file)
        self._bind_action("<Control-Shift-O>",   self.fetch_accession_dialog)
        self._bind_action("<Control-p>",         self.show_protein_analysis)
        self._bind_action("<Control-s>",       self.save_fasta_file)
        self._bind_action("<Control-Shift-S>", self.save_genbank_file)
        self._bind_action("<Control-Shift-s>", self.save_genbank_file)
        self._bind_action("<Control-Shift-C>", self.copy_fasta_to_clipboard)
        self._bind_action("<Control-Shift-c>", self.copy_fasta_to_clipboard)
        self._bind_action("<Control-l>",       self.normalize_sequence)
        self._bind_action("<Control-Shift-L>", self.clear_all)
        self._bind_action("<Control-Shift-l>", self.clear_all)
        self._bind_action("<F5>",              self.validate_sequence)
        self._bind_action("<Control-r>",       self.replace_with_reverse_complement)
        self._bind_action("<Control-t>",       self.translate_selection)
        self._bind_action("<Control-g>",       self.go_to_position)
        self._bind_action("<Control-G>",       self.go_to_position)
        self._bind_action("<Control-Alt-z>",   self.undo_biological_operation)
        self._bind_action("<Control-Alt-Z>",   self.undo_biological_operation)
        self._bind_action("<Control-Alt-y>",   self.redo_biological_operation)
        self._bind_action("<Control-Alt-Y>",   self.redo_biological_operation)
        self._bind_action("<Control-f>",       self.search_motif)
        self._bind_action("<Control-Shift-F>", self.show_orfs)
        self._bind_action("<Control-Shift-f>", self.show_orfs)
        self._bind_action("<Control-Alt-s>",   self.show_secondary_structure_analysis)
        self._bind_action("<Control-Alt-S>",   self.show_secondary_structure_analysis)
        self._bind_action("<Control-e>",       self.add_feature)
        self._bind_action("<F6>",              self.show_next_search_hit)
        self._bind_action("<Shift-F6>",        self.show_previous_search_hit)
        self._bind_action("<Escape>",          self.clear_analysis_marks)

    def _bind_action(self, sequence: str, callback) -> None:
        self.root.bind(sequence, lambda event: self._invoke_action(callback))

    def _invoke_action(self, callback):
        callback()
        return "break"

    def _show_sequence_context_menu(self, event) -> str:
        self.sequence_text.focus_set()
        try:
            self.sequence_context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.sequence_context_menu.grab_release()
        return "break"

    def _show_results_context_menu(self, event) -> str:
        self.results_text.focus_set()
        if self._results_secondary_svg_menu_index is not None:
            state = "normal" if self._results_secondary_structure_export is not None else "disabled"
            self.results_context_menu.entryconfigure(self._results_secondary_svg_menu_index, state=state)
        try:
            self.results_context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.results_context_menu.grab_release()
        return "break"

    def _show_feature_context_menu(self, event) -> str:
        row_id = self.feature_table.identify_row(event.y)
        if row_id:
            self.feature_table.selection_set(row_id)
            self.feature_table.focus(row_id)
        else:
            current_selection = self.feature_table.selection()
            if current_selection:
                self.feature_table.selection_remove(*current_selection)
        try:
            self.feature_context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.feature_context_menu.grab_release()
        return "break"

    def _handle_feature_table_double_click(self, _event) -> str:
        self.edit_selected_feature()
        return "break"

    def _show_marker_context_menu(self, event) -> str:
        row_id = self.marker_table.identify_row(event.y)
        if row_id:
            self.marker_table.selection_set(row_id)
            self.marker_table.focus(row_id)
        try:
            self.marker_context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.marker_context_menu.grab_release()
        return "break"

    # ------------------------------------------------------------------ clipboard / selection helpers

    def _copy_from_widget(self, widget: tk.Text) -> None:
        try:
            selected_text = widget.get("sel.first", "sel.last")
        except tk.TclError:
            return
        self.root.clipboard_clear()
        self.root.clipboard_append(selected_text)

    def _select_all_in_widget(self, widget: tk.Text) -> None:
        widget.focus_set()
        widget.tag_add(tk.SEL, "1.0", "end-1c")
        widget.mark_set(tk.INSERT, "1.0")
        widget.see(tk.INSERT)

    def _replace_sequence_selection(self, new_text: str) -> bool:
        selection = self._get_selection_indices()
        if selection is None:
            return False
        start_index, end_index = selection
        self._clean_index_map_cache = None
        self.sequence_text.delete(start_index, end_index)
        self.sequence_text.insert(start_index, new_text)
        self.sequence_text.tag_remove(tk.SEL, "1.0", tk.END)
        self._refresh_summary()
        self._redraw_features()
        return True

    def _get_clean_selection_from_widget(self, widget: tk.Text) -> str | None:
        try:
            selected_text = widget.get("sel.first", "sel.last")
        except tk.TclError:
            return None
        cleaned = clean_sequence(selected_text)
        return cleaned or None

    def _get_fasta_line_width(self) -> int:
        try:
            value = int(self.line_width_var.get())
        except (tk.TclError, ValueError):
            self.line_width_var.set(100)
            return 100
        if value < 20:
            self.line_width_var.set(100)
            return 100
        return value

    def _format_for_editor(self, sequence: str) -> str:
        return format_fasta(sequence, width=self._get_fasta_line_width())

    def _get_selected_clean_sequence(self) -> str | None:
        return self._get_clean_selection_from_widget(self.sequence_text)

    def _extract_sequence_from_paste(self, raw_text: str) -> str:
        lines = raw_text.splitlines()
        if any(line.lstrip().startswith(">") for line in lines):
            sequence_lines = [
                line
                for line in lines
                if line.strip() and not line.lstrip().startswith(">") and not line.lstrip().startswith(";")
            ]
            return clean_sequence("".join(sequence_lines))
        return clean_sequence(raw_text)

    def copy_selected_sequence(self) -> None:
        self._copy_from_widget(self.sequence_text)

    def cut_selected_sequence(self) -> None:
        selection = self._get_selection_indices()
        if selection is None:
            return
        self.copy_selected_sequence()
        start_index, end_index = selection
        self._clean_index_map_cache = None
        self.sequence_text.delete(start_index, end_index)
        self._refresh_summary()
        self._redraw_features()

    def paste_into_sequence(self) -> None:
        try:
            clipboard_text = self.root.clipboard_get()
        except tk.TclError:
            return
        pasted_sequence = self._extract_sequence_from_paste(clipboard_text)
        if not pasted_sequence:
            return
        self._push_operation_snapshot("paste sequence")
        self._clean_index_map_cache = None
        selection = self._get_selection_indices()
        if selection is not None:
            start_index, end_index = selection
            self.sequence_text.delete(start_index, end_index)
            self.sequence_text.insert(start_index, pasted_sequence)
        else:
            self.sequence_text.insert(tk.INSERT, pasted_sequence)
        normalized_sequence = self._get_summary().cleaned
        self._set_editor_text(self._format_for_editor(normalized_sequence))

    # ------------------------------------------------------------------ BLAST results window

    def _make_dialog_header(self, dialog: tk.Toplevel, title: str, subtitle: str = "") -> None:
        """Dark header with a thin left accent bar."""
        C = self.colors
        # Outer frame — dark surface with bottom border line
        outer = tk.Frame(dialog, bg=C["surface_alt"])
        outer.grid(row=0, column=0, sticky="ew")
        # Thin left accent bar
        tk.Frame(outer, bg=C["accent"], width=4).pack(side="left", fill="y")
        # Content area
        inner = tk.Frame(outer, bg=C["surface_alt"], padx=16, pady=12)
        inner.pack(side="left", fill="both", expand=True)
        tk.Label(inner, text=title, bg=C["surface_alt"], fg=C["text"],
                 font=("Segoe UI Variable Text Semibold", 13)).pack(anchor="w")
        if subtitle:
            tk.Label(inner, text=subtitle, bg=C["surface_alt"], fg=C["accent_dark"],
                     font=("Segoe UI Variable Text", 10)).pack(anchor="w", pady=(2, 0))
        # Bottom separator line
        tk.Frame(dialog, bg=C["border"], height=1).grid(row=0, column=0, sticky="sew")

    def _ask_confirmation_dialog(
        self,
        title: str,
        message: str,
        *,
        subtitle: str = "",
        yes_label: str = "Yes",
        no_label: str = "No",
    ) -> bool:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)

        self._make_dialog_header(dialog, title, subtitle)

        container = ttk.Frame(dialog, style="Card.TFrame", padding=(18, 16, 18, 16))
        container.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        container.columnconfigure(0, weight=1)

        body = tk.Label(
            container,
            text=message,
            justify="left",
            anchor="w",
            wraplength=480,
            bg=C["surface"],
            fg=C["text"],
            font=("Segoe UI", 10),
        )
        body.grid(row=0, column=0, sticky="ew")

        buttons = ttk.Frame(container, style="Card.TFrame")
        buttons.grid(row=1, column=0, sticky="e", pady=(16, 0))

        result = {"value": False}

        def choose(value: bool) -> None:
            result["value"] = value
            dialog.destroy()

        ttk.Button(buttons, text=no_label, command=lambda: choose(False), style="Action.TButton").grid(
            row=0, column=0, padx=(0, 8)
        )
        ttk.Button(buttons, text=yes_label, command=lambda: choose(True), style="Primary.TButton").grid(
            row=0, column=1
        )

        dialog.protocol("WM_DELETE_WINDOW", lambda: choose(False))
        dialog.bind("<Escape>", lambda _event: choose(False))
        dialog.bind("<Return>", lambda _event: choose(True))
        dialog.update_idletasks()
        dialog.geometry(f"+{self.root.winfo_rootx() + 120}+{self.root.winfo_rooty() + 120}")
        dialog.focus_force()
        self.root.wait_window(dialog)
        return bool(result["value"])

    def _show_blast_results_window(self, title: str, text: str) -> None:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.geometry("1220x760")
        dialog.minsize(900, 560)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)
        dialog.rowconfigure(1, weight=1)

        self._make_dialog_header(dialog, title, "BLAST alignment results")

        frame = ttk.Frame(dialog, style="Card.TFrame", padding=(14, 12, 14, 14))
        frame.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        text_frame = ttk.Frame(frame, style="Card.TFrame")
        text_frame.grid(row=0, column=0, sticky="nsew")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)

        viewer = tk.Text(
            text_frame, wrap="none",
            font=("Cascadia Code", 11),
            bg=C["surface"], fg=C["text"],
            insertbackground=C["accent"],
            selectbackground=C["selection_bg"],
            selectforeground=C["selection_fg"],
            relief="flat",
            highlightthickness=1,
            highlightbackground=C["border"],
            highlightcolor=C["accent"],
            padx=12, pady=12,
        )
        viewer.insert("1.0", text)
        viewer.configure(state="disabled")
        scroll_y = ttk.Scrollbar(text_frame, orient="vertical",   command=viewer.yview)
        scroll_x = ttk.Scrollbar(text_frame, orient="horizontal", command=viewer.xview)
        viewer.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        viewer.grid(row=0, column=0, sticky="nsew")
        scroll_y.grid(row=0, column=1, sticky="ns")
        scroll_x.grid(row=1, column=0, sticky="ew")

        # Right-click menu
        ctx = tk.Menu(dialog, tearoff=0, bg=C["surface"], fg=C["text"],
                      activebackground=C["accent_soft"], activeforeground=C["accent_dark"])

        def _copy_selection():
            try:
                sel = viewer.get("sel.first", "sel.last")
            except tk.TclError:
                return
            dialog.clipboard_clear()
            dialog.clipboard_append(sel)

        ctx.add_command(label="Copy", command=_copy_selection)

        def _show_ctx(event):
            try:
                ctx.tk_popup(event.x_root, event.y_root)
            finally:
                ctx.grab_release()

        viewer.bind("<Button-3>", _show_ctx)

    def _show_document_window(self, title: str, subtitle: str, text: str) -> None:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.geometry("1080x720")
        dialog.minsize(860, 520)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)
        dialog.rowconfigure(1, weight=1)

        self._make_dialog_header(dialog, title, subtitle)

        frame = ttk.Frame(dialog, style="Card.TFrame", padding=(14, 12, 14, 14))
        frame.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        text_frame = ttk.Frame(frame, style="Card.TFrame")
        text_frame.grid(row=0, column=0, sticky="nsew")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)

        viewer = tk.Text(
            text_frame,
            wrap="word",
            font=("Segoe UI", 11),
            bg=C["surface"],
            fg=C["text"],
            insertbackground=C["accent"],
            selectbackground=C["selection_bg"],
            selectforeground=C["selection_fg"],
            relief="flat",
            highlightthickness=1,
            highlightbackground=C["border"],
            highlightcolor=C["accent"],
            padx=16,
            pady=14,
        )
        viewer.insert("1.0", text.rstrip())
        viewer.configure(state="disabled")
        scroll_y = ttk.Scrollbar(text_frame, orient="vertical", command=viewer.yview)
        viewer.configure(yscrollcommand=scroll_y.set)
        viewer.grid(row=0, column=0, sticky="nsew")
        scroll_y.grid(row=0, column=1, sticky="ns")

        ctx = tk.Menu(dialog, tearoff=0, bg=C["surface"], fg=C["text"],
                      activebackground=C["accent_soft"], activeforeground=C["accent_dark"])

        def _copy_selection() -> None:
            try:
                sel = viewer.get("sel.first", "sel.last")
            except tk.TclError:
                return
            dialog.clipboard_clear()
            dialog.clipboard_append(sel)

        def _copy_all() -> None:
            dialog.clipboard_clear()
            dialog.clipboard_append(viewer.get("1.0", "end-1c"))

        ctx.add_command(label="Copy", command=_copy_selection)
        ctx.add_command(label="Copy all", command=_copy_all)

        def _show_ctx(event) -> None:
            try:
                ctx.tk_popup(event.x_root, event.y_root)
            finally:
                ctx.grab_release()

        viewer.bind("<Button-3>", _show_ctx)

    def _about_information_text(self) -> str:
        return (
            "GeneDraft\n"
            f"Version {APP_VERSION}\n\n"
            f"{APP_DESCRIPTION}\n\n"
            "Autor\n"
            f"{APP_AUTHOR}\n\n"
            "Afiliación\n"
            f"{APP_AFFILIATION}\n\n"
            "Autor de correspondencia\n"
            f"{APP_AUTHOR}*\n"
            f"ORCID: {APP_ORCID}\n\n"
            "Primary focus\n"
            "- FASTA and GenBank editing\n"
            "- motif search and ORF review\n"
            "- selection-based translation\n"
            "- DNA/RNA secondary structure prediction\n"
            "- circular component delimitation\n"
            "- feature annotation\n"
            "- protein analysis\n"
            "- BLAST web and local BLAST workflows\n"
        )

    def _about_help_text(self) -> str:
        return (
            "Getting started\n"
            "1. Open a FASTA or GenBank file, or paste a sequence into the editor.\n"
            "2. The editor keeps the sequence normalized in FASTA blocks using the current FASTA width.\n"
            "3. Use Find Motif, Find ORFs, Translate, Circularize, and Add Feature as the main workflow.\n\n"
            "Open and save\n"
            "- Open FASTA or GenBank from File > Open.\n"
            "- Save FASTA or GenBank from the File menu.\n"
            "- Edit the visible Name field in the editor header to control exported sequence names.\n\n"
            "Normalize\n"
            "- Normalize rewrites the current content as a cleaned FASTA sequence.\n"
            "- FASTA width is controlled from the settings strip.\n"
            "- Pasted sequence content is normalized automatically.\n\n"
            "Motif search\n"
            "- Use Find Motif to locate a subsequence.\n"
            "- Matches are highlighted in the editor.\n"
            "- F6 and Shift+F6 move between matches and select the active range.\n\n"
            "ORFs\n"
            "- ORF search scans the full sequence on both strands.\n"
            "- Minimum ORF size is controlled from the top settings strip.\n"
            "- Active ORF results can be converted into features quickly with Ctrl+E.\n\n"
            "Translate\n"
            "- Translate works on the current selection only.\n"
            "- Select the region first, then run Translate.\n"
            "- Translation reports forward and reverse strand frames for the selected region.\n\n"
            "Secondary structure\n"
            "- Secondary Structure works on the current DNA or RNA selection.\n"
            "- DNA uses ViennaRNA DNA parameters; RNA uses ViennaRNA RNA parameters.\n"
            "- The report includes dot-bracket structure and delta G (MFE).\n"
            "- Right-click in the Results panel to export the secondary-structure SVG on demand.\n\n"
            "Circularize\n"
            "- Circularize works on the current selection.\n"
            "- It ranks candidates using terminal overlap, mismatches, and scan window.\n"
            "- Review candidates with F6 and Shift+F6, then mark or crop the selected component.\n\n"
            "Features\n"
            "- Add Feature uses the current selection.\n"
            "- If the selection matches an active motif or ORF hit, the dialog pre-fills useful defaults.\n"
            "- Features can be focused, edited, translated, exported, and reused in maps.\n\n"
            "BLAST web\n"
            "- Select a region in the editor or results panel and launch BLAST web.\n"
            "- No local installation is required for BLAST web.\n\n"
            "BLAST local\n"
            "- Local BLAST uses an already built or newly built BLAST+ database.\n"
            "- Database creation is available from the BLAST menu.\n"
            "- Keep BLAST database paths and project paths simple, especially on Windows.\n\n"
            "Python requirements\n"
            "- Python 3\n"
            "- tkinter\n"
            "- Biopython 1.84 or newer\n"
            "- ViennaRNA Python bindings\n\n"
            "Programs expected in PATH for local BLAST\n"
            "- makeblastdb\n"
            "- blastn\n"
            "- blastp\n"
        )

    def _about_shortcuts_text(self) -> str:
        return (
            "Shortcuts\n\n"
            "General\n"
            "- Ctrl+O: Open file\n"
            "- Ctrl+S: Save FASTA\n"
            "- Ctrl+Shift+S: Save GenBank\n"
            "- Ctrl+Shift+C: Copy current FASTA\n"
            "- Ctrl+Shift+L: Clear All\n\n"
            "Sequence editing\n"
            "- Ctrl+L: Normalize\n"
            "- F5: Validate\n"
            "- Ctrl+R: Reverse complement\n"
            "- Ctrl+T: Translate selection\n"
            "- Ctrl+Alt+S: Secondary structure\n"
            "- Ctrl+G: Go to position\n\n"
            "Analysis\n"
            "- Ctrl+F: Find motif\n"
            "- Ctrl+Shift+F: Find ORFs\n"
            "- F6: Next result\n"
            "- Shift+F6: Previous result\n"
            "- Esc: Clear analysis highlights\n\n"
            "Annotation\n"
            "- Ctrl+E: Add Feature from the current selection\n\n"
            "Workflow note\n"
            "- Search hits selected with F6 or Shift+F6 can be turned into features quickly with Ctrl+E.\n"
        )

    def _about_version_text(self) -> str:
        return (
            "GeneDraft version information\n\n"
            f"Version: {APP_VERSION}\n"
            f"Release year: {APP_RELEASE_YEAR}\n"
            f"Repository: {APP_REPOSITORY_URL}\n"
            "License: MIT\n"
        )

    def _about_citation_text(self) -> str:
        return (
            "Suggested citation\n\n"
            f"Cárdenas-Conejo Y. ({APP_RELEASE_YEAR}). GeneDraft (Version {APP_VERSION}) "
            f"[Computer software]. {APP_REPOSITORY_URL}.\n\n"
            "Autor\n"
            f"{APP_AUTHOR}\n"
            f"ORCID: {APP_ORCID}\n\n"
            "Afiliación\n"
            f"{APP_AFFILIATION}\n\n"
            "Citation note\n"
            "- Update the repository field once the GitHub project is public.\n"
            "- If you release new versions, cite the version actually used in the analysis.\n"
        )

    def _about_license_text(self) -> str:
        return (
            "License and use\n\n"
            "GeneDraft is distributed under the MIT License.\n\n"
            "Permission summary\n"
            "- Use, copy, modify, publish, distribute, sublicense, and sell are permitted.\n"
            "- Keep the copyright notice and license text with redistributed copies.\n\n"
            "Warranty\n"
            "The software is provided 'as is', without warranty of any kind.\n\n"
            "Disclaimer\n"
            "GeneDraft is provided as research-support software. Users are responsible for "
            "sequence interpretation, biological conclusions, and downstream decisions.\n"
        )

    def _about_acknowledgments_text(self) -> str:
        return (
            "Acknowledgments\n\n"
            "GeneDraft builds on and benefits from the following tools and resources:\n"
            "- Biopython\n"
            "- NCBI BLAST and BLAST+\n"
            "- Python and tkinter\n\n"
            "These tools support sequence parsing, analysis, local BLAST workflows, and desktop UI behavior.\n"
        )

    def _about_contact_text(self) -> str:
        return (
            "Contact\n\n"
            f"Autor: {APP_AUTHOR}\n"
            f"Autor de correspondencia: {APP_AUTHOR}*\n"
            f"Email: {APP_CONTACT_EMAIL}\n"
            f"ORCID: {APP_ORCID}\n\n"
            "Afiliación\n"
            f"{APP_AFFILIATION}\n\n"
            "Repositorio\n"
            f"{APP_REPOSITORY_URL}\n"
        )

    def show_about_information(self) -> None:
        self._show_document_window("About GeneDraft", "Project information", self._about_information_text())

    def show_about_help(self) -> None:
        self._show_document_window("GeneDraft Help", "How to use the application", self._about_help_text())

    def show_about_shortcuts(self) -> None:
        self._show_document_window("GeneDraft Shortcuts", "Keyboard shortcuts for efficient use", self._about_shortcuts_text())

    def show_about_version(self) -> None:
        self._show_document_window("GeneDraft Version", "Release details", self._about_version_text())

    def show_about_citation(self) -> None:
        self._show_document_window("How to Cite GeneDraft", "Academic citation guidance", self._about_citation_text())

    def show_about_license(self) -> None:
        self._show_document_window("GeneDraft License and Use", "Academic use, commercial use, and disclaimer", self._about_license_text())

    def show_about_acknowledgments(self) -> None:
        self._show_document_window("GeneDraft Acknowledgments", "External tools and resources", self._about_acknowledgments_text())

    def show_about_contact(self) -> None:
        self._show_document_window("GeneDraft Contact", "Author and contact information", self._about_contact_text())

    # ------------------------------------------------------------------ BLAST dispatch

    def run_blast_web_from_selection(self) -> None:
        sequence = self._get_clean_selection_from_widget(self.sequence_text)
        if not sequence:
            messagebox.showinfo("No selection", "Select a sequence to send to BLAST web.")
            return
        self._run_blast_web(sequence)

    def run_blast_web_from_results(self) -> None:
        sequence = self._get_clean_selection_from_widget(self.results_text)
        if not sequence:
            messagebox.showinfo("No selection", "Select a sequence in Results to send to BLAST web.")
            return
        self._run_blast_web(sequence)

    def run_blast_local_from_selection(self) -> None:
        sequence = self._get_clean_selection_from_widget(self.sequence_text)
        if not sequence:
            messagebox.showinfo("No selection", "Select a sequence to run local BLAST.")
            return
        self._run_blast_local(sequence)

    def run_blast_local_from_results(self) -> None:
        sequence = self._get_clean_selection_from_widget(self.results_text)
        if not sequence:
            messagebox.showinfo("No selection", "Select a sequence in Results to run local BLAST.")
            return
        self._run_blast_local(sequence)

    def _run_blast_web(self, sequence: str) -> None:
        try:
            program = self._detect_query_program(sequence)
        except Exception as exc:
            messagebox.showerror("Could not detect sequence type", str(exc))
            return

        params = {"PAGE_TYPE": "BlastSearch", "PROGRAM": program, "QUERY": sequence}
        url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?{urlencode(params)}"
        try:
            webbrowser.open(url)
        except Exception:
            self.root.clipboard_clear()
            self.root.clipboard_append(url)
            self._set_results(
                "Could not open the browser from the app.\n\n"
                "The BLAST web URL was copied to the clipboard."
            )
            return

        self._set_results(
            "BLAST web opened in the browser.\n\n"
            f"Program: {program}\n"
            f"Query length: {len(sequence)}"
        )

    def _run_blast_local(self, sequence: str) -> None:
        db_prefix = self._normalize_blast_db_prefix(self.local_blast_db_var.get())
        if not db_prefix or not self._blast_db_prefix_exists(db_prefix):
            detected = self._discover_local_blast_db_prefixes()
            if detected:
                if self._ask_confirmation_dialog(
                    "No local database",
                    "There is no valid local BLAST database configured.\n\n"
                    "Existing BLAST databases were detected. Do you want to choose one now?",
                    subtitle="BLAST database configuration",
                ):
                    self.configure_local_blast_db()
            else:
                messagebox.showinfo("No local database", "Create or configure a local BLAST database first.")
            return
        if not self._ensure_blast_paths_supported(db_prefix, self.local_blast_source_var.get().strip()):
            return
        self.local_blast_db_var.set(db_prefix)
        self._save_app_settings()

        configured_db_type = self.local_blast_type_var.get().strip() or "auto"
        db_type = configured_db_type
        if db_type == "auto":
            fasta_source = self.local_blast_source_var.get().strip()
            if fasta_source:
                try:
                    _, source_sequence, _, _ = load_sequence_file(fasta_source)
                    db_type = self._detect_blast_db_type(source_sequence)
                except Exception:
                    db_type = "nucl"
            else:
                db_type = "nucl"

        try:
            program = self._detect_query_program(sequence, db_type=db_type)
        except Exception as exc:
            messagebox.showerror("Incompatible query", str(exc))
            return

        self._set_busy_status(f"Working: running local BLAST with {program}...")
        try:
            with tempfile.TemporaryDirectory(prefix="genedraft_blast_") as temp_dir:
                query_path = Path(temp_dir) / "query.fasta"
                query_path.write_text(f">selection\n{sequence}\n", encoding="utf-8")

                command = [
                    program, "-query", str(query_path), "-db", db_prefix,
                    "-outfmt", "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send stitle",
                    "-max_target_seqs", "50",
                ]
                try:
                    completed = subprocess.run(command, capture_output=True, text=True, check=True)
                except FileNotFoundError:
                    messagebox.showerror("BLAST+ not found", f"Could not find {program} in PATH.")
                    return
                except subprocess.CalledProcessError as exc:
                    messagebox.showerror("Local BLAST failed", (exc.stderr or exc.stdout or str(exc)).strip())
                    return

            output = (completed.stdout or "").strip()
            if not output:
                self._set_results(f"Local BLAST completed with no hits.\n\nProgram: {program}\nDatabase: {db_prefix}")
                return

            lines = [line for line in output.splitlines() if line.strip()]
        finally:
            self._clear_busy_status()

        header = "qseqid\tsseqid\tpident\tlength\tevalue\tbitscore\tqstart\tqend\tsstart\tsend\tstitle"
        display_text = "\n".join([header, *lines])
        top_hit = lines[0].split("\t")
        self._set_results(
            f"Local BLAST completed.\n\nProgram: {program}\nDatabase: {db_prefix}\n"
            f"Hits: {len(lines)}\n"
            f"Top hit: {top_hit[1] if len(top_hit) > 1 else '-'}\n"
            f"Identity: {top_hit[2] if len(top_hit) > 2 else '-'}%\n"
            f"E-value: {top_hit[4] if len(top_hit) > 4 else '-'}"
        )
        self._show_blast_results_window(f"BLAST local: {program}", display_text)

    # ------------------------------------------------------------------ text events

    def _on_text_modified(self, _event) -> None:
        if self._internal_edit:
            self.sequence_text.edit_modified(False)
            return
        if self.sequence_text.edit_modified():
            self.sequence_text.edit_modified(False)
            self._clean_index_map_cache = None
            self._history_redo.clear()
            self._reset_search_state()
            self._clear_analysis_highlights()
            self._refresh_summary()
            self._redraw_features()
            self._schedule_session_save()

    def _refresh_summary(self) -> None:
        summary = self._get_summary()
        self.length_var.set(f"Length: {summary.length}")
        self.type_var.set(f"Type: {summary.molecule_type}")
        self.gc_var.set("GC: -" if summary.gc_percent is None else f"GC: {summary.gc_percent:.1f}%")
        if self._busy_status_message:
            status = self._busy_status_message
            badge_mode = "busy"
        elif not summary.cleaned:
            status = "No sequence"
            badge_mode = "idle"
        elif summary.invalid_characters:
            status = f"Review: {', '.join(summary.invalid_characters)}"
            badge_mode = "review"
        else:
            status = "Valid"
            badge_mode = "idle"
        self.warning_var.set(status)
        self._apply_warning_badge_style(badge_mode)
        self.features_var.set(f"Features: {len(self.features)}")
        self._update_position_info()
        self._schedule_minimap_redraw()

    def _apply_warning_badge_style(self, mode: str) -> None:
        if not hasattr(self, "status_warning_label"):
            return
        C = self.colors
        styles = {
            "busy": (C["warning"], C["warning_fg"], C["warning_fg"]),
            "review": (C["danger"], C["danger_fg"], C["danger_fg"]),
            "done": (C["success"], C["success_fg"], C["success_fg"]),
            "idle": (C["surface"], C["muted"], C["border"]),
        }
        bg, fg, border = styles.get(mode, styles["idle"])
        if mode == "idle":
            self.status_warning_label.grid_remove()
            return
        self.status_warning_label.grid()
        self.status_warning_label.configure(
            bg=bg,
            fg=fg,
            highlightbackground=border,
            highlightcolor=border,
        )

    def _set_busy_status(self, message: str) -> None:
        if self._busy_done_after_id is not None:
            try:
                self.root.after_cancel(self._busy_done_after_id)
            except tk.TclError:
                pass
            self._busy_done_after_id = None
        self._busy_status_message = message.strip() or "Working..."
        self.warning_var.set(self._busy_status_message)
        self._apply_warning_badge_style("busy")
        try:
            self.root.configure(cursor="watch")
            self.root.update_idletasks()
        except tk.TclError:
            pass

    def _clear_busy_status(self) -> None:
        self._busy_status_message = None
        if self._busy_done_after_id is not None:
            try:
                self.root.after_cancel(self._busy_done_after_id)
            except tk.TclError:
                pass
            self._busy_done_after_id = None
        try:
            self.root.configure(cursor="")
            self.root.update_idletasks()
        except tk.TclError:
            pass
        summary = self._get_summary()
        if summary.cleaned and not summary.invalid_characters:
            self.warning_var.set("Done")
            self._apply_warning_badge_style("done")
            try:
                self._busy_done_after_id = self.root.after(2000, self._clear_temporary_done_status)
            except tk.TclError:
                self._busy_done_after_id = None
                self._refresh_summary()
            return
        self._refresh_summary()

    def _clear_temporary_done_status(self) -> None:
        self._busy_done_after_id = None
        self._refresh_summary()

    def _update_position_info(self, _event=None) -> None:
        index = self.sequence_text.index(tk.INSERT)
        prefix = self.sequence_text.get("1.0", index)
        summary = self._get_summary()
        clean_position = len(clean_sequence(prefix))
        if summary.length == 0:
            display_position = 0
        elif clean_position == 0:
            display_position = 1
        else:
            display_position = min(clean_position, summary.length)
        self.cursor_var.set(f"Pos: {display_position}")

        selection_range = self._get_selection_clean_range()
        if selection_range is None:
            self.selection_var.set("Selection: 0 nt")
        else:
            start_clean, end_clean = selection_range
            self.selection_var.set(f"Selection: {end_clean - start_clean} nt")
        self._schedule_minimap_redraw()

    def _get_editor_text(self) -> str:
        return self.sequence_text.get("1.0", "end-1c")

    def _set_editor_text(self, text: str) -> None:
        self._internal_edit = True
        self._clean_index_map_cache = None
        self.sequence_text.delete("1.0", tk.END)
        self.sequence_text.insert("1.0", text)
        self.sequence_text.edit_modified(False)
        self._internal_edit = False
        self._refresh_summary()
        self._redraw_features()
        self._schedule_session_save()

    def _get_summary(self) -> SequenceSummary:
        return summarize_sequence(self._get_editor_text())

    def _set_results(
        self,
        text: str,
        *,
        secondary_structure_export: SecondaryStructureExportState | None = None,
    ) -> None:
        self._results_secondary_structure_export = secondary_structure_export
        self.results_text.configure(state="normal")
        self.results_text.delete("1.0", tk.END)
        self.results_text.insert("1.0", text)
        self.results_text.configure(state="disabled")

    def _get_orf_min_length(self) -> int | None:
        try:
            value = int(self.orf_min_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("Invalid value", "The minimum ORF length must be an integer.")
            self.orf_min_var.set(110)
            return None
        if value < 1:
            messagebox.showwarning("Invalid value", "The minimum ORF length must be at least 1.")
            self.orf_min_var.set(110)
            return None
        return value

    def _detect_blast_db_type(self, sequence_text: str) -> str:
        sequence_type = summarize_sequence(sequence_text).molecule_type
        if sequence_type in {"DNA", "RNA"}:
            return "nucl"
        if sequence_type == "Protein":
            return "prot"
        raise ValueError("Could not determine whether the sequence is nucleotide or protein.")

    def _detect_query_program(self, sequence_text: str, db_type: str | None = None) -> str:
        sequence_type = summarize_sequence(sequence_text).molecule_type
        if sequence_type in {"DNA", "RNA"}:
            if db_type in {None, "nucl"}:
                return "blastn"
        if sequence_type == "Protein":
            if db_type in {None, "prot"}:
                return "blastp"
        raise ValueError("The selection and the database are not compatible for automatic blastn/blastp mode.")

    def select_blast_fasta(self) -> None:
        file_path = filedialog.askopenfilename(
            title="Select FASTA for local database",
            filetypes=[("FASTA", "*.fasta *.fa *.fna *.faa *.fas"), ("All files", "*.*")],
        )
        if not file_path:
            return
        self.local_blast_source_var.set(file_path)
        stem = Path(file_path).stem
        db_dir = Path.cwd() / "blast_db" / stem
        self.local_blast_db_var.set(str(db_dir / stem))
        try:
            _, sequence, _, _ = load_sequence_file(file_path)
            self.local_blast_type_var.set(self._detect_blast_db_type(sequence))
        except Exception:
            self.local_blast_type_var.set("auto")
        self._save_app_settings()

    def _suggest_blast_db_prefix(self, file_path: str) -> str:
        stem = Path(file_path).stem
        db_dir = Path.cwd() / "blast_db" / stem
        return str(db_dir / stem)

    def _coerce_blast_db_prefix(self, source: str, db_prefix: str) -> str:
        source = source.strip()
        db_prefix = db_prefix.strip()
        suggested = self._suggest_blast_db_prefix(source) if source else ""
        if not db_prefix:
            return suggested

        prefix_path = Path(db_prefix)
        if db_prefix.endswith(("\\", "/")) or not prefix_path.name:
            if not suggested:
                return db_prefix
            return str(prefix_path / Path(suggested).name)

        if prefix_path.exists() and prefix_path.is_dir():
            if not suggested:
                return str(prefix_path / prefix_path.name)
            return str(prefix_path / Path(suggested).name)

        return db_prefix

    def configure_local_blast_db(self) -> None:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title("Configure BLAST DB")
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)

        self._make_dialog_header(
            dialog,
            "Configure local BLAST database",
            "Use an existing database without running makeblastdb again",
        )

        container = ttk.Frame(dialog, style="Card.TFrame", padding=(18, 16, 18, 16))
        container.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        container.columnconfigure(1, weight=1)

        discovered = self._discover_local_blast_db_prefixes()
        current_prefix = self._normalize_blast_db_prefix(self.local_blast_db_var.get())
        current_type = self.local_blast_type_var.get().strip() or "auto"

        detected_var = tk.StringVar(value=current_prefix)
        prefix_var = tk.StringVar(value=current_prefix)
        type_var = tk.StringVar(value=current_type)
        source_var = tk.StringVar(value=self.local_blast_source_var.get())

        def sync_from_detected(_event=None) -> None:
            selected = detected_var.get().strip()
            if not selected:
                return
            prefix_var.set(selected)
            for prefix, db_type in discovered:
                if prefix == selected:
                    type_var.set(db_type)
                    break

        def browse_existing_db() -> None:
            file_path = filedialog.askopenfilename(
                title="Select an existing BLAST database file",
                filetypes=[
                    ("BLAST indices", "*.pin *.nin *.psq *.nsq *.phr *.nhr"),
                    ("All files", "*.*"),
                ],
                parent=dialog,
            )
            if not file_path:
                return
            normalized = self._normalize_blast_db_prefix(file_path)
            prefix_var.set(normalized)
            detected_var.set(normalized)
            suffix = Path(file_path).suffix.lower()
            if suffix.startswith(".p"):
                type_var.set("prot")
            elif suffix.startswith(".n"):
                type_var.set("nucl")

        ttk.Label(container, text="Detected databases", style="SectionHead.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 10))
        detected_combo = ttk.Combobox(
            container,
            textvariable=detected_var,
            values=[prefix for prefix, _db_type in discovered],
            state="normal",
            width=48,
        )
        detected_combo.grid(row=0, column=1, sticky="ew")
        detected_combo.bind("<<ComboboxSelected>>", sync_from_detected)
        ttk.Button(container, text="Browse DB", command=browse_existing_db, style="Action.TButton").grid(row=0, column=2, sticky="ew", padx=(8, 0))

        ttk.Label(container, text="DB prefix", style="SectionHead.TLabel").grid(row=1, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        ttk.Entry(container, textvariable=prefix_var, style="Modern.TEntry", width=48).grid(row=1, column=1, sticky="ew", pady=(10, 0))
        ttk.Combobox(container, textvariable=type_var, values=("auto", "nucl", "prot"), state="readonly", width=8).grid(row=1, column=2, sticky="ew", padx=(8, 0), pady=(10, 0))

        ttk.Label(container, text="Source FASTA (optional)", style="SectionHead.TLabel").grid(row=2, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        ttk.Entry(container, textvariable=source_var, style="Modern.TEntry", width=48).grid(row=2, column=1, sticky="ew", pady=(10, 0))
        ttk.Button(
            container,
            text="Browse FASTA",
            command=lambda: source_var.set(
                filedialog.askopenfilename(
                    title="Select source FASTA",
                    filetypes=[("FASTA", "*.fasta *.fa *.fna *.faa *.fas"), ("All files", "*.*")],
                    parent=dialog,
                ) or source_var.get()
            ),
            style="Action.TButton",
        ).grid(row=2, column=2, sticky="ew", padx=(8, 0), pady=(10, 0))

        if discovered:
            help_text = "You can choose a detected database or manually point to the prefix of an existing one."
        else:
            help_text = "No databases were detected in blast_db. You can manually point to an existing prefix."
        ttk.Label(container, text=help_text, style="Muted.TLabel").grid(row=3, column=0, columnspan=3, sticky="w", pady=(10, 16))

        buttons = ttk.Frame(container, style="Card.TFrame")
        buttons.grid(row=4, column=0, columnspan=3, sticky="e")

        def submit() -> None:
            db_prefix = self._normalize_blast_db_prefix(prefix_var.get())
            if not db_prefix:
                messagebox.showinfo("No database", "Enter the BLAST database prefix.", parent=dialog)
                return
            if not self._blast_db_prefix_exists(db_prefix):
                messagebox.showerror(
                    "Missing database path",
                    "No BLAST database files were found for that prefix.\n"
                    "Select a .pin/.nin file or enter the correct prefix.",
                    parent=dialog,
                )
                return
            self._apply_local_blast_db_config(db_prefix, type_var.get().strip() or "auto", source_var.get().strip())
            self._set_results(
                "Local BLAST database configured.\n\n"
                f"DB prefix: {db_prefix}\n"
                f"Type: {self.local_blast_type_var.get() or 'auto'}\n"
                f"Source FASTA: {self.local_blast_source_var.get().strip() or '(not set)'}"
            )
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=dialog.destroy, style="Action.TButton").grid(row=0, column=0, padx=(0, 8))
        ttk.Button(buttons, text="Save DB", command=submit, style="Primary.TButton").grid(row=0, column=1)

        dialog.bind("<Escape>", lambda event: dialog.destroy())
        dialog.bind("<Return>", lambda event: submit())
        dialog.update_idletasks()
        dialog.geometry(f"+{self.root.winfo_rootx() + 120}+{self.root.winfo_rooty() + 120}")
        dialog.focus_force()
        self.root.wait_window(dialog)

    def _build_local_blast_db_from_values(self, source: str, db_prefix: str, requested_type: str) -> bool:
        if not source:
            messagebox.showinfo("No FASTA", "Select the FASTA file you want to use as the source database.")
            return False

        source_path = Path(source)
        if not source_path.exists():
            messagebox.showerror("Missing path", "The selected FASTA file does not exist.")
            return False

        db_prefix = self._coerce_blast_db_prefix(source, db_prefix)
        if not db_prefix:
            messagebox.showerror("Missing database prefix", "Could not determine an output name for the BLAST database.")
            return False
        if not self._ensure_blast_paths_supported(source, db_prefix):
            return False

        requested_type = requested_type.strip() or "auto"
        prefix_path = Path(db_prefix)
        prefix_path.parent.mkdir(parents=True, exist_ok=True)
        sanitized_fasta_path = prefix_path.parent / f"{prefix_path.name}__sanitized_input.fasta"

        self._set_busy_status("Working: building local BLAST database...")
        try:
            try:
                record_count, detected_db_type = sanitize_fasta_for_blast(source_path, sanitized_fasta_path)
            except Exception as exc:
                messagebox.showerror("Could not prepare FASTA", str(exc))
                return False

            db_type = detected_db_type if requested_type == "auto" else requested_type

            command = [
                "makeblastdb", "-in", str(sanitized_fasta_path), "-input_type", "fasta",
                "-dbtype", db_type, "-out", str(prefix_path),
                "-title", prefix_path.name, "-blastdb_version", "5",
            ]

            try:
                completed = subprocess.run(command, capture_output=True, text=True, check=True)
            except FileNotFoundError:
                messagebox.showerror("BLAST+ not found", "makeblastdb was not found in PATH.")
                return False
            except subprocess.CalledProcessError as exc:
                messagebox.showerror("Could not build database", (exc.stderr or exc.stdout or str(exc)).strip())
                return False

            self._apply_local_blast_db_config(str(prefix_path), db_type, str(source_path))
            self._set_results(
                "Local BLAST database created.\n\n"
                f"FASTA: {source_path.name}\n"
                f"Loaded sequences: {record_count}\n"
                f"Type: {db_type}\n"
                f"DB prefix: {prefix_path}\n\n"
                "Note: the database is built from a sanitized FASTA and without parse_seqids to avoid problems "
                "with leading comments or long headers."
            )
            output_text = (completed.stdout or "") + ("\n" + completed.stderr if completed.stderr else "")
            if output_text.strip():
                self._show_blast_results_window("BLAST database build", output_text.strip())
            return True
        finally:
            self._clear_busy_status()

    def build_local_blast_db(self) -> None:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title("Build BLAST DB")
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)

        self._make_dialog_header(dialog, "Build local BLAST database",
                                 "Select a FASTA file and define the output prefix")

        container = ttk.Frame(dialog, style="Card.TFrame", padding=(18, 16, 18, 16))
        container.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        container.columnconfigure(1, weight=1)

        source_var = tk.StringVar(value=self.local_blast_source_var.get())
        db_var     = tk.StringVar(value=self._coerce_blast_db_prefix(self.local_blast_source_var.get(), self.local_blast_db_var.get()))
        type_var   = tk.StringVar(value=self.local_blast_type_var.get() or "auto")

        def browse_fasta() -> None:
            file_path = filedialog.askopenfilename(
                title="Select FASTA for local database",
                filetypes=[("FASTA", "*.fasta *.fa *.fna *.faa *.fas"), ("All files", "*.*")],
                parent=dialog,
            )
            if not file_path:
                return
            source_var.set(file_path)
            db_var.set(self._suggest_blast_db_prefix(file_path))
            try:
                _, sequence, _, _ = load_sequence_file(file_path)
                type_var.set(self._detect_blast_db_type(sequence))
            except Exception:
                type_var.set("auto")

        ttk.Label(container, text="Source FASTA", style="SectionHead.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 10))
        ttk.Entry(container, textvariable=source_var, style="Modern.TEntry", width=46).grid(row=0, column=1, sticky="ew")
        ttk.Button(container, text="Browse", command=browse_fasta, style="Action.TButton").grid(row=0, column=2, sticky="ew", padx=(8, 0))

        ttk.Label(container, text="DB prefix", style="SectionHead.TLabel").grid(row=1, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        ttk.Entry(container, textvariable=db_var, style="Modern.TEntry", width=46).grid(row=1, column=1, sticky="ew", pady=(10, 0))
        ttk.Combobox(container, textvariable=type_var, values=("auto", "nucl", "prot"), state="readonly", width=8).grid(row=1, column=2, sticky="ew", padx=(8, 0), pady=(10, 0))

        ttk.Label(container, text="The app sanitizes the FASTA before running makeblastdb.", style="Muted.TLabel").grid(row=2, column=0, columnspan=3, sticky="w", pady=(10, 16))

        buttons = ttk.Frame(container, style="Card.TFrame")
        buttons.grid(row=3, column=0, columnspan=3, sticky="e")

        def submit() -> None:
            db_var.set(self._coerce_blast_db_prefix(source_var.get().strip(), db_var.get().strip()))
            if self._build_local_blast_db_from_values(source_var.get().strip(), db_var.get().strip(), type_var.get().strip()):
                dialog.destroy()

        ttk.Button(buttons, text="Cancel",     command=dialog.destroy, style="Action.TButton").grid(row=0, column=0, padx=(0, 8))
        ttk.Button(
            buttons,
            text="Use existing",
            command=lambda: (dialog.destroy(), self.configure_local_blast_db()),
            style="Action.TButton",
        ).grid(row=0, column=1, padx=(0, 8))
        ttk.Button(buttons, text="Build DB", command=submit,         style="Primary.TButton").grid(row=0, column=2)

        dialog.bind("<Escape>", lambda event: dialog.destroy())
        dialog.bind("<Return>", lambda event: submit())
        dialog.update_idletasks()
        dialog.geometry(f"+{self.root.winfo_rootx() + 120}+{self.root.winfo_rooty() + 120}")
        dialog.focus_force()
        self.root.wait_window(dialog)

    def show_local_blast_db_info(self) -> None:
        db_prefix = self._normalize_blast_db_prefix(self.local_blast_db_var.get())
        if not db_prefix:
            messagebox.showinfo("No database", "There is no local BLAST database configured yet.")
            return
        status = "Available" if self._blast_db_prefix_exists(db_prefix) else "Not found"
        self._set_results(
            "Local BLAST database configured.\n\n"
            f"DB prefix: {db_prefix}\n"
            f"Type: {self.local_blast_type_var.get() or 'auto'}\n"
            f"Source FASTA: {self.local_blast_source_var.get().strip() or '(not set)'}\n"
            f"Status: {status}\n"
            f"Saved config: {self.settings_path.name}"
        )

    # ------------------------------------------------------------------ search state

    def _reset_search_state(self) -> None:
        self.search_matches = []
        self.search_index   = -1
        self.search_label   = ""

    def _clear_analysis_highlights(self) -> None:
        for tag in self._analysis_tags:
            self.sequence_text.tag_remove(tag, "1.0", tk.END)

    def clear_analysis_marks(self) -> None:
        self._reset_search_state()
        self._clear_analysis_highlights()
        for item in self.feature_table.selection():
            self.feature_table.selection_remove(item)
        self._set_results("Analysis highlights cleared. Features remain in place.")

    def _clear_feature_tags(self) -> None:
        for tag in self._feature_tags:
            self.sequence_text.tag_delete(tag)
        self._feature_tags.clear()

    def _clean_positions_to_text_indices(self) -> list[str]:
        if self._clean_index_map_cache is not None:
            return self._clean_index_map_cache
        indices: list[str] = []
        raw_text = self._get_editor_text()
        for offset, char in enumerate(raw_text):
            if char.isalpha() or char == "*":
                indices.append(self.sequence_text.index(f"1.0+{offset}c"))
        self._clean_index_map_cache = indices
        return indices

    def _highlight_clean_span(self, start_clean: int, end_clean: int, tag_name: str) -> bool:
        if start_clean < 0 or end_clean <= start_clean:
            return False
        clean_index_map = self._clean_positions_to_text_indices()
        if start_clean >= len(clean_index_map):
            return False
        end_clean_position = min(end_clean - 1, len(clean_index_map) - 1)
        start_index = clean_index_map[start_clean]
        end_index   = self.sequence_text.index(f"{clean_index_map[end_clean_position]}+1c")
        self.sequence_text.tag_add(tag_name, start_index, end_index)
        return True

    def _highlight_nt_range(self, start_nt: int, end_nt: int, tag_name: str) -> bool:
        return self._highlight_clean_span(start_nt - 1, end_nt, tag_name)

    def _select_nt_range(self, start_nt: int, end_nt: int) -> bool:
        clean_index_map = self._clean_positions_to_text_indices()
        if not clean_index_map:
            return False
        start_clean = max(0, start_nt - 1)
        end_clean = max(start_clean + 1, end_nt)
        if start_clean >= len(clean_index_map):
            return False
        end_clean_position = min(end_clean - 1, len(clean_index_map) - 1)
        start_index = clean_index_map[start_clean]
        end_index = self.sequence_text.index(f"{clean_index_map[end_clean_position]}+1c")
        self.sequence_text.tag_remove(tk.SEL, "1.0", tk.END)
        self.sequence_text.tag_add(tk.SEL, start_index, end_index)
        self.sequence_text.mark_set(tk.INSERT, start_index)
        self.sequence_text.see(start_index)
        self.sequence_text.focus_set()
        self._update_position_info()
        return True

    def _scroll_to_nt(self, start_nt: int) -> None:
        clean_index_map = self._clean_positions_to_text_indices()
        if not clean_index_map:
            return
        position = max(1, min(start_nt, len(clean_index_map)))
        index = clean_index_map[position - 1]
        self.sequence_text.mark_set(tk.INSERT, index)
        self.sequence_text.see(index)
        self._update_position_info()

    def _get_selection_indices(self) -> tuple[str, str] | None:
        try:
            return self.sequence_text.index("sel.first"), self.sequence_text.index("sel.last")
        except tk.TclError:
            return None

    def _get_selection_clean_range(self) -> tuple[int, int] | None:
        selection = self._get_selection_indices()
        if selection is None:
            return None
        start_index, end_index = selection
        prefix        = self.sequence_text.get("1.0", start_index)
        selected_text = self.sequence_text.get(start_index, end_index)
        start_clean   = len(clean_sequence(prefix))
        selected_clean = clean_sequence(selected_text)
        if not selected_clean:
            return None
        return start_clean, start_clean + len(selected_clean)

    def _get_insert_nt_position(self) -> int:
        summary = self._get_summary()
        if summary.length <= 0:
            return 1
        index = self.sequence_text.index(tk.INSERT)
        prefix = self.sequence_text.get("1.0", index)
        clean_position = len(clean_sequence(prefix))
        return max(1, min(summary.length, clean_position if clean_position > 0 else 1))

    def _copy_text_to_clipboard(self, text: str) -> None:
        self.root.clipboard_clear()
        self.root.clipboard_append(text)

    def _get_selection_sequence(self) -> tuple[int, int, str] | None:
        summary = self._get_summary()
        selection_range = self._get_selection_clean_range()
        if selection_range is None:
            return None
        start_clean, end_clean = selection_range
        return start_clean + 1, end_clean, summary.cleaned[start_clean:end_clean]

    def _collect_features_for_range(self, start_nt: int, end_nt: int) -> list[SequenceFeature]:
        collected: list[SequenceFeature] = []
        for feature in self.features:
            overlap_start = max(feature.start_nt, start_nt)
            overlap_end = min(feature.end_nt, end_nt)
            if overlap_start > overlap_end:
                continue
            collected.append(
                SequenceFeature(
                    start_nt=overlap_start - start_nt + 1,
                    end_nt=overlap_end - start_nt + 1,
                    feature_type=feature.feature_type,
                    label=feature.label,
                    strand=self._normalize_feature_strand(feature.strand),
                )
            )
        return collected

    def _feature_display_color(self, feature: FeatureEntry | SequenceFeature) -> str:
        return "#0EA5E9" if self._normalize_feature_strand(feature.strand) == "+" else "#F87171"

    def _feature_border_color(self, feature: FeatureEntry | SequenceFeature) -> str:
        return "#38BDF8" if self._normalize_feature_strand(feature.strand) == "+" else "#FCA5A5"

    def _format_coordinates(self, start_nt: int, end_nt: int, strand: str | None = None) -> str:
        length = max(0, end_nt - start_nt + 1)
        if strand is None:
            return f"{start_nt}-{end_nt} ({length} nt)"
        return f"{start_nt}-{end_nt} {strand} ({length} nt)"

    def _get_selected_feature(self) -> tuple[int, FeatureEntry] | None:
        feature_index = self._get_selected_feature_index()
        if feature_index is None:
            return None
        return feature_index, self.features[feature_index]

    def _feature_sequence(self, feature: FeatureEntry) -> str:
        summary = self._get_summary()
        if not summary.cleaned or feature.start_nt < 1 or feature.end_nt > summary.length:
            return ""
        sequence = summary.cleaned[feature.start_nt - 1 : feature.end_nt]
        if self._normalize_feature_strand(feature.strand) == "-":
            try:
                return reverse_complement(sequence)
            except Exception:
                return sequence
        return sequence

    def _feature_gc_percent(self, feature: FeatureEntry) -> str:
        sequence = self._feature_sequence(feature)
        if not sequence:
            return "-"
        gc = summarize_sequence(sequence).gc_percent
        return "-" if gc is None else f"{gc:.1f}%"

    def _translate_feature_sequence(self, feature: FeatureEntry) -> str | None:
        sequence = self._feature_sequence(feature)
        if not sequence:
            return None
        molecule_type = summarize_sequence(sequence).molecule_type
        if molecule_type not in {"DNA", "RNA"}:
            return None
        try:
            frames = translate_frames(sequence, include_reverse=False)
        except Exception:
            return None
        return frames[0][1] if frames else None

    def go_to_position(self) -> None:
        summary = self._get_summary()
        if summary.length <= 0:
            messagebox.showinfo("No sequence", "No sequence is loaded.")
            return
        position = simpledialog.askinteger(
            "Go to position",
            f"Target position (1-{summary.length}):",
            initialvalue=self._get_insert_nt_position(),
            minvalue=1,
            maxvalue=summary.length,
        )
        if position is None:
            return
        self._scroll_to_nt(position)
        self._set_results(f"Direct navigation.\n\nCurrent position: {position}")

    def save_selected_fasta_file(self) -> None:
        selection = self._get_selection_sequence()
        if selection is None:
            messagebox.showinfo("No selection", "Select a sequence region first.")
            return
        start_nt, end_nt, selected_sequence = selection
        file_path = filedialog.asksaveasfilename(
            title="Save selection as FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA", "*.fasta"), ("All files", "*.*")],
        )
        if not file_path:
            return
        header = f"{self.sequence_name.get().strip() or 'sequence'}_{start_nt}_{end_nt}"
        try:
            save_fasta(file_path, header, selected_sequence)
        except Exception as exc:
            messagebox.showerror("Could not save", str(exc))
            return
        self._mark_session_persisted()
        self._set_results(
            "Selection exported as FASTA.\n\n"
            f"File: {Path(file_path).name}\n"
            f"Range: {start_nt}-{end_nt}"
        )

    def save_selected_genbank_file(self) -> None:
        selection = self._get_selection_sequence()
        if selection is None:
            messagebox.showinfo("No selection", "Select a sequence region first.")
            return
        start_nt, end_nt, selected_sequence = selection
        file_path = filedialog.asksaveasfilename(
            title="Save selection as GenBank",
            defaultextension=".gbk",
            filetypes=[("GenBank", "*.gbk *.gb"), ("All files", "*.*")],
        )
        if not file_path:
            return
        header = f"{self.sequence_name.get().strip() or 'sequence'}_{start_nt}_{end_nt}"
        try:
            save_genbank(file_path, header, selected_sequence, self._collect_features_for_range(start_nt, end_nt))
        except Exception as exc:
            messagebox.showerror("Could not save", str(exc))
            return
        self._mark_session_persisted()
        self._set_results(
            "Selection exported as GenBank.\n\n"
            f"File: {Path(file_path).name}\n"
            f"Range: {start_nt}-{end_nt}"
        )

    def copy_selection_coordinates(self) -> None:
        selection = self._get_selection_sequence()
        if selection is None:
            messagebox.showinfo("No selection", "Select a sequence region first.")
            return
        start_nt, end_nt, _selected_sequence = selection
        coords = self._format_coordinates(start_nt, end_nt)
        self._copy_text_to_clipboard(coords)
        self._set_results(f"Coordinates copied.\n\nSelection: {coords}")

    def _refresh_marker_table(self) -> None:
        if not hasattr(self, "marker_table"):
            return
        for item_id in self.marker_table.get_children():
            self.marker_table.delete(item_id)
        for index, marker in enumerate(self.markers):
            self.marker_table.insert("", tk.END, iid=str(index), values=(marker.position_nt, marker.label))
        self._schedule_minimap_redraw()
        self._schedule_session_save()

    def _get_selected_marker_index(self) -> int | None:
        if not hasattr(self, "marker_table"):
            return None
        selection = self.marker_table.selection()
        if not selection:
            return None
        return int(selection[0])

    def _on_marker_table_select(self, _event) -> None:
        self.focus_selected_marker()

    def add_marker_at_cursor(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "No sequence is loaded.")
            return
        self._push_operation_snapshot("add marker")
        position = self._get_insert_nt_position()
        self.markers.append(MarkerEntry(position_nt=position, label=f"marker_{len(self.markers) + 1}"))
        self.markers.sort(key=lambda marker: marker.position_nt)
        self._refresh_marker_table()
        self._redraw_minimap()
        self._set_results(f"Marker added.\n\nPosition: {position}")

    def add_marker_from_selection(self) -> None:
        selection = self._get_selection_sequence()
        if selection is None:
            messagebox.showinfo("No selection", "Select a sequence region first.")
            return
        start_nt, _end_nt, _selected_sequence = selection
        self._push_operation_snapshot("add marker from selection")
        self.markers.append(MarkerEntry(position_nt=start_nt, label=f"marker_{len(self.markers) + 1}"))
        self.markers.sort(key=lambda marker: marker.position_nt)
        self._refresh_marker_table()
        self._redraw_minimap()
        self._set_results(f"Marker added from selection.\n\nPosition: {start_nt}")

    def focus_selected_marker(self) -> None:
        marker_index = self._get_selected_marker_index()
        if marker_index is None:
            return
        marker = self.markers[marker_index]
        self._scroll_to_nt(marker.position_nt)
        self._set_results(
            "Marker focused.\n\n"
            f"Label: {marker.label}\n"
            f"Position: {marker.position_nt}"
        )

    def delete_selected_marker(self) -> None:
        marker_index = self._get_selected_marker_index()
        if marker_index is None:
            messagebox.showinfo("No marker", "Select a marker to delete.")
            return
        self._push_operation_snapshot("delete marker")
        removed = self.markers.pop(marker_index)
        self._refresh_marker_table()
        self._redraw_minimap()
        self._set_results(
            "Marker deleted.\n\n"
            f"Label: {removed.label}\n"
            f"Position: {removed.position_nt}"
        )

    def copy_selected_marker_coordinates(self) -> None:
        marker_index = self._get_selected_marker_index()
        if marker_index is None:
            messagebox.showinfo("No marker", "Select a marker first.")
            return
        marker = self.markers[marker_index]
        coords = f"{marker.position_nt}"
        self._copy_text_to_clipboard(coords)
        self._set_results(f"Marker coordinates copied.\n\n{marker.label}: {coords}")

    def _on_minimap_click(self, event) -> str:
        summary = self._get_summary()
        if summary.length <= 0:
            return "break"
        width = max(1, self.minimap_canvas.winfo_width())
        left = 18
        right = max(left + 1, width - 18)
        x = max(left, min(event.x, right))
        ratio = (x - left) / max(1, right - left)
        position = max(1, min(summary.length, int(round(ratio * (summary.length - 1))) + 1))
        self._scroll_to_nt(position)
        return "break"

    def _redraw_minimap(self) -> None:
        if not hasattr(self, "minimap_canvas"):
            return
        canvas = self.minimap_canvas
        canvas.delete("all")
        width = max(1, canvas.winfo_width())
        height = max(1, canvas.winfo_height())
        canvas.create_rectangle(0, 0, width, height, fill=self.colors["surface_alt"], outline="")
        summary = self._get_summary()
        if summary.length <= 0:
            canvas.create_text(width / 2, height / 2, text="Mini map", fill=self.colors["muted"], font=("Segoe UI", 10))
            return

        left = 18
        right = width - 18
        top = 18
        bottom = height - 18
        baseline_y = (top + bottom) / 2
        usable = max(1, right - left)

        def x_for_nt(position_nt: int) -> float:
            return left + ((position_nt - 1) / max(1, summary.length - 1)) * usable

        canvas.create_line(left, baseline_y, right, baseline_y, fill=self.colors["accent_dark"], width=4, capstyle=tk.ROUND)
        for feature in self.features:
            if feature.start_nt < 1 or feature.end_nt > summary.length or feature.start_nt > feature.end_nt:
                continue
            x1 = x_for_nt(feature.start_nt)
            x2 = x_for_nt(feature.end_nt)
            y1 = top + (4 if self._normalize_feature_strand(feature.strand) == "+" else 18)
            y2 = y1 + 12
            canvas.create_rectangle(
                x1,
                y1,
                max(x1 + 4, x2),
                y2,
                fill=self._feature_display_color(feature),
                outline=self._feature_border_color(feature),
                width=1,
            )

        selection = self._get_selection_sequence()
        if selection is not None:
            start_nt, end_nt, _selected_sequence = selection
            canvas.create_rectangle(
                x_for_nt(start_nt),
                top,
                x_for_nt(end_nt),
                bottom,
                fill="#D29922",
                outline="",
                stipple="gray25",
            )

        current_match = self._get_current_circular_match()
        if current_match is not None:
            canvas.create_rectangle(
                x_for_nt(current_match.start_nt),
                top - 2,
                x_for_nt(current_match.end_nt),
                bottom + 2,
                outline="#0EA5E9",
                width=2,
            )

        for marker in self.markers:
            if 1 <= marker.position_nt <= summary.length:
                x = x_for_nt(marker.position_nt)
                canvas.create_line(x, top - 2, x, bottom + 2, fill="#D29922", width=2)

        cursor_nt = self._get_insert_nt_position()
        cursor_x = x_for_nt(cursor_nt)
        canvas.create_polygon(
            cursor_x, top - 4,
            cursor_x - 6, top + 6,
            cursor_x + 6, top + 6,
            fill=self.colors["accent_dark"],
            outline="",
        )
        canvas.create_text(left, height - 6, text="1", anchor="sw", fill=self.colors["muted"], font=("Segoe UI", 9))
        canvas.create_text(right, height - 6, text=str(summary.length), anchor="se", fill=self.colors["muted"], font=("Segoe UI", 9))

    # ------------------------------------------------------------------ features

    def _normalize_feature_strand(self, strand: str) -> str:
        return "-" if str(strand).strip() == "-" else "+"

    def _invert_feature_strand(self, strand: str) -> str:
        return "-" if self._normalize_feature_strand(strand) == "+" else "+"

    def _get_selected_feature_index(self) -> int | None:
        selection = self.feature_table.selection()
        if not selection:
            return None
        return int(selection[0])

    def _refresh_feature_table(self) -> None:
        for item_id in self.feature_table.get_children():
            self.feature_table.delete(item_id)
        for index, feature in enumerate(self.features):
            self.feature_table.insert("", tk.END, iid=str(index),
                values=(
                    feature.feature_type,
                    feature.start_nt,
                    feature.end_nt,
                    self._normalize_feature_strand(feature.strand),
                    feature.label,
                ))
        self.features_var.set(f"Features: {len(self.features)}")
        self._schedule_minimap_redraw()
        self._schedule_session_save()

    def _redraw_features(self) -> None:
        self._clear_feature_tags()
        summary = self._get_summary()
        for index, feature in enumerate(self.features):
            if feature.start_nt < 1 or feature.end_nt > summary.length or feature.start_nt > feature.end_nt:
                continue
            tag_name = f"feature_{index}"
            self._feature_tags.add(tag_name)
            self.sequence_text.tag_configure(
                tag_name,
                background=self._feature_display_color(feature),
                relief="raised",
                borderwidth=1,
            )
            self._highlight_nt_range(feature.start_nt, feature.end_nt, tag_name)

        self.sequence_text.tag_raise("motif")
        self.sequence_text.tag_raise("motif_focus")
        self.sequence_text.tag_raise("orf_plus")
        self.sequence_text.tag_raise("orf_minus")
        self.sequence_text.tag_raise("circular_candidate")
        self.sequence_text.tag_raise("circular_focus")
        self.sequence_text.tag_raise("orf_focus")
        self.sequence_text.tag_raise("invalid")
        self.sequence_text.tag_raise("feature_focus")
        self.sequence_text.tag_raise("sel")

    def _set_feature_focus(self, feature_index: int) -> None:
        self.sequence_text.tag_remove("feature_focus", "1.0", tk.END)
        if feature_index < 0 or feature_index >= len(self.features):
            return
        feature = self.features[feature_index]
        if self._highlight_nt_range(feature.start_nt, feature.end_nt, "feature_focus"):
            self._scroll_to_nt(feature.start_nt)

    def _on_feature_table_select(self, _event) -> None:
        self.focus_selected_feature()

    def _on_feature_table_click(self, event) -> None:
        """Clicking on empty space in the feature table deselects everything."""
        if not self.feature_table.identify_row(event.y):
            for item in self.feature_table.selection():
                self.feature_table.selection_remove(item)
            self.sequence_text.tag_remove("feature_focus", "1.0", tk.END)

    def focus_selected_feature(self) -> None:
        feature_index = self._get_selected_feature_index()
        if feature_index is None:
            return
        self._set_feature_focus(feature_index)
        self.show_selected_feature_summary()

    def show_selected_feature_summary(self) -> None:
        selected = self._get_selected_feature()
        if selected is None:
            messagebox.showinfo("No feature", "Select a feature first.")
            return
        _feature_index, feature = selected
        feature_sequence = self._feature_sequence(feature)
        summary_lines = [
            "Feature summary",
            "",
            f"Label: {feature.label}",
            f"Type: {feature.feature_type}",
            f"Range: {feature.start_nt}-{feature.end_nt}",
            f"Strand: {self._normalize_feature_strand(feature.strand)}",
            f"Length: {max(0, feature.end_nt - feature.start_nt + 1)} nt",
            f"Local GC: {self._feature_gc_percent(feature)}",
        ]
        if feature_sequence:
            summary_lines.extend([
                "",
                "Oriented sequence:",
                format_fasta(feature_sequence, width=60),
            ])
        if feature.feature_type.strip().lower() == "cds":
            protein = self._translate_feature_sequence(feature)
            if protein is not None:
                summary_lines.extend([
                    "",
                    f"Estimated translation: {max(0, len(protein))} aa",
                    format_fasta(protein, width=60),
                ])
        self._set_results("\n".join(summary_lines).rstrip())

    def translate_selected_feature(self) -> None:
        selected = self._get_selected_feature()
        if selected is None:
            messagebox.showinfo("No feature", "Select a feature first.")
            return
        _feature_index, feature = selected
        protein = self._translate_feature_sequence(feature)
        if protein is None:
            messagebox.showinfo("Not translatable", "The selected feature could not be translated as a CDS.")
            return
        self._set_results(
            "Feature/CDS translation\n\n"
            f"Label: {feature.label}\n"
            f"Range: {feature.start_nt}-{feature.end_nt}\n"
            f"Strand: {self._normalize_feature_strand(feature.strand)}\n"
            f"Protein length: {len(protein)} aa\n\n"
            f"{format_fasta(protein, width=60)}"
        )

    def copy_selected_feature_coordinates(self) -> None:
        selected = self._get_selected_feature()
        if selected is None:
            messagebox.showinfo("No feature", "Select a feature first.")
            return
        _feature_index, feature = selected
        coords = self._format_coordinates(feature.start_nt, feature.end_nt, self._normalize_feature_strand(feature.strand))
        self._copy_text_to_clipboard(coords)
        self._set_results(f"Feature coordinates copied.\n\n{feature.label}: {coords}")

    # ------------------------------------------------------------------ search navigation

    def _set_search_matches(self, matches: list[SearchMatch], label: str) -> None:
        self.search_matches = matches
        self.search_index   = 0 if matches else -1
        self.search_label   = label

    def _render_current_search_hit(self) -> None:
        if not self.search_matches or self.search_index < 0:
            return
        match = self.search_matches[self.search_index]
        position_line = f"Result {self.search_index + 1} of {len(self.search_matches)}"
        text = (
            f"{self.search_label}\n\n"
            f"{position_line}\n"
            f"{match.title}\n"
            f"Range: {match.start_nt}-{match.end_nt}\n\n"
            f"{match.details}"
        )
        self._set_results(text.rstrip())

    def _get_current_search_match(self) -> SearchMatch | None:
        if not self.search_matches or self.search_index < 0:
            return None
        if self.search_index >= len(self.search_matches):
            return None
        return self.search_matches[self.search_index]

    def _get_feature_defaults_from_active_match(self) -> tuple[str, str, str]:
        default_label = f"feature_{len(self.features) + 1}"
        default_type = "misc_feature"
        default_strand = "+"
        selection_range = self._get_selection_clean_range()
        match = self._get_current_search_match()
        if selection_range is None or match is None:
            return default_label, default_type, default_strand

        start_clean, end_clean = selection_range
        if (start_clean + 1, end_clean) != (match.start_nt, match.end_nt):
            return default_label, default_type, default_strand

        if match.kind == "motif":
            motif_name = match.title.removeprefix("Motif ").split(" (", 1)[0].strip()
            return (
                f"motif_{motif_name}" if motif_name else default_label,
                "misc_feature",
                self._normalize_feature_strand(str(getattr(match.payload, "strand", "+"))),
            )

        if match.kind == "orf" and isinstance(match.payload, OrfHit):
            frame_label = match.payload.frame_label.strip() or "orf"
            return (
                f"ORF_{frame_label}",
                "CDS",
                self._normalize_feature_strand(match.payload.strand),
            )

        return default_label, default_type, default_strand

    def _focus_current_search_hit(self) -> None:
        if not self.search_matches or self.search_index < 0:
            messagebox.showinfo("No results", "There is no active search to navigate.")
            return
        match = self.search_matches[self.search_index]
        self._clear_analysis_highlights()

        if match.kind == "motif":
            for current in self.search_matches:
                self._highlight_nt_range(current.start_nt, current.end_nt, "motif")
            self._highlight_nt_range(match.start_nt, match.end_nt, "motif_focus")
        elif match.kind == "orf":
            for current in self.search_matches:
                hit = current.payload
                if isinstance(hit, OrfHit):
                    tag_name = "orf_plus" if hit.strand == "+" else "orf_minus"
                    self._highlight_nt_range(hit.start_nt, hit.end_nt, tag_name)
            self._highlight_nt_range(match.start_nt, match.end_nt, "orf_focus")
        elif match.kind == "circular":
            for current in self.search_matches:
                self._highlight_nt_range(current.start_nt, current.end_nt, "circular_candidate")
            self._highlight_nt_range(match.start_nt, match.end_nt, "circular_focus")

        self._select_nt_range(match.start_nt, match.end_nt)
        self._scroll_to_nt(match.start_nt)
        self._redraw_features()
        self._render_current_search_hit()

    def show_next_search_hit(self) -> None:
        if not self.search_matches:
            messagebox.showinfo("No results", "Run a motif search, ORF search, or circularization first.")
            return
        self.search_index = (self.search_index + 1) % len(self.search_matches)
        self._focus_current_search_hit()

    def show_previous_search_hit(self) -> None:
        if not self.search_matches:
            messagebox.showinfo("No results", "Run a motif search, ORF search, or circularization first.")
            return
        self.search_index = (self.search_index - 1) % len(self.search_matches)
        self._focus_current_search_hit()

    # ------------------------------------------------------------------ file operations

    def new_window(self) -> None:
        subprocess.Popen([sys.executable, __file__])

    # ------------------------------------------------------------------ fetch accession
    def _get_ncbi_email(self) -> str:
        if getattr(self, "_ncbi_email", ""):
            return self._ncbi_email
        email = simpledialog.askstring(
            "NCBI Email",
            "NCBI requires an email for database access.\nEnter your email (saved for future use):",
            parent=self.root,
        )
        if email and email.strip():
            self._ncbi_email = email.strip()
            self._save_app_settings()
            return self._ncbi_email
        return ""

    def fetch_accession_dialog(self) -> None:
        accession = simpledialog.askstring(
            "Fetch by Accession",
            "Enter NCBI accession number:\n(e.g.  NM_007294  ·  NP_001254  ·  EU490707  ·  P04637)",
            parent=self.root,
        )
        if not accession or not accession.strip():
            return
        email = self._get_ncbi_email()
        if not email:
            return
        self._set_busy_status("Working: fetching accession from NCBI...")
        try:
            try:
                header, sequence, loaded_features, db = fetch_by_accession(accession.strip(), email)
            except Exception as exc:
                messagebox.showerror("Fetch failed", str(exc))
                return
            self._push_operation_snapshot("fetch accession")
            self.sequence_name.set(header)
            self.features = [
                FeatureEntry(
                    start_nt=feature.start_nt,
                    end_nt=feature.end_nt,
                    feature_type=feature.feature_type,
                    label=feature.label,
                    color=next(self._feature_colors),
                    strand=self._normalize_feature_strand(feature.strand),
                )
                for feature in loaded_features
            ]
            self.markers = []
            self._reset_search_state()
            self._refresh_feature_table()
            self._refresh_marker_table()
            self._clear_analysis_highlights()
            self._set_editor_text(self._format_for_editor(sequence))
            self._set_results(
                f"Fetched from NCBI.\n\n"
                f"Accession: {header}\n"
                f"Database:  {db}\n"
                f"Imported features: {len(self.features)}"
            )
        finally:
            self._clear_busy_status()

    # ------------------------------------------------------------------ protein analysis
    def protein_analysis_from_results(self) -> None:
        sequence = self._get_clean_selection_from_widget(self.results_text)
        if not sequence:
            messagebox.showwarning("No selection", "Select a protein sequence in the Results panel first.")
            return
        summary = summarize_sequence(sequence)
        if summary.molecule_type != "Protein":
            messagebox.showwarning(
                "Protein required",
                "The selected text does not look like a protein sequence.",
            )
            return
        self._run_protein_analysis(summary.cleaned)

    def _run_protein_analysis(self, sequence: str) -> None:
        """Show local results immediately, then update with remote scores in background."""
        import threading
        from sequence_tools import _fetch_protein_sol, _fetch_soluprot, fetch_hmmer_pfam

        # 1. Instant local report — all remote fields show "fetching..."
        try:
            report = analyze_protein(sequence, fetching=True)
            self._set_results(report)
        except Exception as exc:
            messagebox.showerror("Analysis error", str(exc))
            return

        # 2. Background thread fetches all three remote sources in parallel
        self._set_busy_status("Working: fetching remote protein annotations...")

        def _fetch():
            results: dict = {}

            def _get_ps():
                results["ps"] = _fetch_protein_sol(sequence)
            def _get_sp():
                results["sp"] = _fetch_soluprot(sequence)
            def _get_hmmer():
                results["hmmer"] = fetch_hmmer_pfam(sequence)

            workers = [
                threading.Thread(target=_get_ps,    daemon=True),
                threading.Thread(target=_get_sp,    daemon=True),
                threading.Thread(target=_get_hmmer, daemon=True),
            ]
            for w in workers:
                w.start()
            for w in workers:
                w.join()

            self.root.after(0, lambda: _update(
                results.get("ps"), results.get("sp"), results.get("hmmer")
            ))

        def _update(ps, sp, hmmer):
            try:
                report = analyze_protein(
                    sequence,
                    ps_score=ps, sp_score=sp,
                    hmmer_hits=hmmer,
                    fetching=False,
                )
                self._set_results(report)
            except Exception:
                pass
            finally:
                self._clear_busy_status()

        threading.Thread(target=_fetch, daemon=True).start()

    def show_protein_analysis(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showwarning("No sequence", "There is no sequence to analyze.")
            return
        if summary.molecule_type != "Protein":
            messagebox.showwarning(
                "Protein required",
                "Protein Analysis requires a protein sequence.\n\n"
                "Paste a protein sequence, or translate a CDS first\n"
                "(select a region and use Ctrl+T).",
            )
            return
        self._run_protein_analysis(summary.cleaned)

    def _safe_name_fragment(self, text: str) -> str:
        cleaned = "".join(ch if (ch.isalnum() or ch in {"_", "-"}) else "_" for ch in text.strip())
        return cleaned.strip("._-")[:80] or "sequence"

    def _secondary_structure_svg_path(self, molecule_type: str, start_nt: int, end_nt: int) -> Path:
        output_dir = Path.cwd() / "generated_structures"
        base_name = self._safe_name_fragment(self.sequence_name.get().strip() or "sequence")
        return output_dir / f"{base_name}_{molecule_type.lower()}_{start_nt}_{end_nt}_secondary_structure.svg"

    def export_secondary_structure_svg_from_results(self) -> None:
        export_state = self._results_secondary_structure_export
        if export_state is None:
            messagebox.showinfo(
                "No secondary structure",
                "The current Results panel does not contain a secondary structure prediction to export.",
            )
            return

        export_state.suggested_path.parent.mkdir(parents=True, exist_ok=True)
        file_path = filedialog.asksaveasfilename(
            title="Export secondary structure SVG",
            defaultextension=".svg",
            initialdir=str(export_state.suggested_path.parent),
            initialfile=export_state.suggested_path.name,
            filetypes=[("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not file_path:
            return

        try:
            export_secondary_structure_svg(
                file_path,
                export_state.folded_sequence,
                export_state.structure,
            )
        except Exception as exc:
            messagebox.showerror("Could not export", str(exc))
            return

        self.warning_var.set(f"Secondary structure SVG exported: {Path(file_path).name}")

    def show_secondary_structure_analysis(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showwarning("No sequence", "There is no sequence to analyze.")
            return
        selection_range = self._get_selection_clean_range()
        if selection_range is None:
            messagebox.showinfo("No selection", "Select the DNA or RNA region you want to fold.")
            return

        start_clean, end_clean = selection_range
        selected_sequence = summary.cleaned[start_clean:end_clean]
        selected_type = summary.molecule_type
        if selected_type not in {"DNA", "RNA"}:
            selected_type = summarize_sequence(selected_sequence).molecule_type
        if selected_type not in {"DNA", "RNA"}:
            messagebox.showwarning(
                "DNA or RNA required",
                "Secondary structure prediction requires a DNA or RNA selection.",
            )
            return

        self._set_busy_status("Working: predicting secondary structure...")
        try:
            try:
                prediction = predict_secondary_structure(selected_sequence, molecule_type=selected_type)
                svg_path = self._secondary_structure_svg_path(prediction.input_type, start_clean + 1, end_clean)
            except Exception as exc:
                messagebox.showerror("Secondary structure error", str(exc))
                return

            lines = [
                "Secondary Structure Prediction",
                "=" * 46,
                f"Selected range:       {start_clean + 1}-{end_clean}",
                f"Input type:           {prediction.input_type}",
                f"Folding model:        {prediction.model_label}",
                f"Selected length:      {len(prediction.folded_sequence)} nt",
                f"Delta G (MFE):        {prediction.mfe:.2f} kcal/mol",
                "",
                "Sequence",
                format_fasta(prediction.folded_sequence, width=60),
                "",
                "Dot-bracket structure",
                format_fasta(prediction.structure, width=60),
                "",
                "SVG export: Right-click this Results panel and choose",
                "'Export secondary structure SVG...'",
            ]
            self._set_results(
                "\n".join(lines).rstrip(),
                secondary_structure_export=SecondaryStructureExportState(
                    suggested_path=svg_path,
                    folded_sequence=prediction.folded_sequence,
                    structure=prediction.structure,
                ),
            )
        finally:
            self._clear_busy_status()

    def convert_dna_to_rna(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to convert.")
            return
        if summary.molecule_type == "RNA":
            self._set_results("The current sequence is already RNA.")
            return
        if summary.molecule_type != "DNA":
            messagebox.showwarning("DNA required", "DNA to RNA conversion only applies to DNA sequences.")
            return
        try:
            converted = dna_to_rna(summary.cleaned)
        except Exception as exc:
            messagebox.showerror("Conversion error", str(exc))
            return

        self._push_operation_snapshot("convert DNA to RNA")
        self._reset_search_state()
        self._clear_analysis_highlights()
        self._set_editor_text(self._format_for_editor(converted))
        self._set_results(
            "DNA converted to RNA.\n\n"
            f"Length: {len(converted)} nt\n"
            "All thymine residues were replaced with uracil."
        )

    def open_sequence_file(self) -> None:
        file_path = filedialog.askopenfilename(
            title="Open sequence file",
            filetypes=[
                ("Sequence files", "*.fasta *.fa *.fna *.faa *.fas *.gb *.gbk *.genbank"),
                ("FASTA", "*.fasta *.fa *.fna *.faa *.fas"),
                ("GenBank", "*.gb *.gbk *.genbank"),
                ("All files", "*.*"),
            ],
        )
        if not file_path:
            return

        self._set_busy_status(f"Working: opening {Path(file_path).name}...")
        try:
            try:
                header, sequence, loaded_features, file_format = load_sequence_file(file_path)
            except Exception as exc:
                messagebox.showerror("Could not open", str(exc))
                return

            self._push_operation_snapshot("open file")
            self.sequence_name.set(header)
            self.features = [
                FeatureEntry(
                    start_nt=feature.start_nt,
                    end_nt=feature.end_nt,
                    feature_type=feature.feature_type,
                    label=feature.label,
                    color=next(self._feature_colors),
                    strand=self._normalize_feature_strand(feature.strand),
                )
                for feature in loaded_features
            ]
            self.markers = []
            self._reset_search_state()
            self._refresh_feature_table()
            self._refresh_marker_table()
            self._clear_analysis_highlights()
            self._set_editor_text(self._format_for_editor(sequence))
            self._set_results(
                "File loaded.\n\n"
                f"Name: {Path(file_path).name}\n"
                f"Format: {file_format}\n"
                f"Imported features: {len(self.features)}"
            )
            self._push_recent_file(file_path)
        finally:
            self._clear_busy_status()

    def save_fasta_file(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showwarning("No sequence", "There is no sequence to save.")
            return
        file_path = filedialog.asksaveasfilename(
            title="Save FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA", "*.fasta"), ("All files", "*.*")],
        )
        if not file_path:
            return
        try:
            save_fasta(file_path, self.sequence_name.get(), summary.cleaned)
        except Exception as exc:
            messagebox.showerror("Could not save", str(exc))
            return
        self._mark_session_persisted()
        self._set_results(
            "FASTA file saved.\n\n"
            f"Name: {Path(file_path).name}\n"
            f"Exported length: {summary.length}"
        )

    def save_genbank_file(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showwarning("No sequence", "There is no sequence to save.")
            return
        file_path = filedialog.asksaveasfilename(
            title="Save GenBank",
            defaultextension=".gbk",
            filetypes=[("GenBank", "*.gbk *.gb"), ("All files", "*.*")],
        )
        if not file_path:
            return
        try:
            save_genbank(
                file_path, self.sequence_name.get(), summary.cleaned,
                [feature.as_sequence_feature() for feature in self.features],
            )
        except Exception as exc:
            messagebox.showerror("Could not save", str(exc))
            return
        self._mark_session_persisted()
        self._set_results(
            "GenBank file saved.\n\n"
            f"Name: {Path(file_path).name}\n"
            f"Exported features: {len(self.features)}"
        )

    def copy_fasta_to_clipboard(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showwarning("No sequence", "There is no sequence to copy.")
            return
        fasta_text = f">{self.sequence_name.get().strip() or 'sequence'}\n{format_fasta(summary.cleaned)}"
        self.root.clipboard_clear()
        self.root.clipboard_append(fasta_text)
        self._set_results("The current FASTA sequence was copied to the clipboard.")

    # ------------------------------------------------------------------ sequence operations

    def normalize_sequence(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence content to normalize.")
            return
        self._push_operation_snapshot("normalize sequence")
        self._reset_search_state()
        width = self._get_fasta_line_width()
        self._set_editor_text(self._format_for_editor(summary.cleaned))
        self._set_results(
            "The sequence was cleaned and formatted in FASTA blocks.\n\n"
            f"Current width: {width} nt/aa per line"
        )

    def clear_all(self) -> None:
        self._push_operation_snapshot("clear all")
        self.sequence_name.set("sequence_1")
        self.features = []
        self.markers = []
        self._reset_search_state()
        self._refresh_feature_table()
        self._refresh_marker_table()
        self._clear_analysis_highlights()
        self._clear_feature_tags()
        self._set_editor_text("")
        self._set_results("Editor content, highlights, and features were cleared.")

    def validate_sequence(self) -> None:
        summary = self._get_summary()
        self._clear_analysis_highlights()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to validate.")
            return
        self._set_busy_status("Working: validating sequence...")
        try:
            raw_text = self._get_editor_text()
            invalid_count = 0
            for offset, char in enumerate(raw_text):
                if not (char.isalpha() or char == "*"):
                    continue
                if char.upper() in summary.invalid_characters:
                    self.sequence_text.tag_add("invalid", f"1.0+{offset}c", f"1.0+{offset + 1}c")
                    invalid_count += 1
            if summary.invalid_characters:
                self._set_results(
                    "Validation completed.\n\n"
                    f"Detected type: {summary.molecule_type}\n"
                    f"Invalid characters: {', '.join(summary.invalid_characters)}\n"
                    f"Positions marked in red: {invalid_count}"
                )
            else:
                self._set_results(
                    "Validation completed.\n\n"
                    f"Detected type: {summary.molecule_type}\n"
                    f"Length: {summary.length}\n"
                    "No invalid characters were detected."
                )
            self._redraw_features()
        finally:
            self._clear_busy_status()

    def replace_with_reverse_complement(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to transform.")
            return
        self._push_operation_snapshot("reverse complement")
        try:
            reversed_sequence = reverse_complement(summary.cleaned)
        except Exception as exc:
            messagebox.showerror("Could not compute", str(exc))
            return
        length = summary.length
        self.features = [
            FeatureEntry(
                start_nt=length - feature.end_nt + 1,
                end_nt=length - feature.start_nt + 1,
                feature_type=feature.feature_type,
                label=feature.label,
                color=feature.color,
                strand=self._invert_feature_strand(feature.strand),
            )
            for feature in self.features
        ]
        self.markers = [
            MarkerEntry(position_nt=length - marker.position_nt + 1, label=marker.label)
            for marker in self.markers
            if 1 <= marker.position_nt <= length
        ]
        self.markers.sort(key=lambda marker: marker.position_nt)
        self._reset_search_state()
        self._refresh_feature_table()
        self._refresh_marker_table()
        self._clear_analysis_highlights()
        self._set_editor_text(self._format_for_editor(reversed_sequence))
        self._set_results(
            "Reverse complement applied.\n\n"
            f"New length: {len(reversed_sequence)}\n"
            f"Adjusted features: {len(self.features)}"
        )

    def _prompt_motif_search(self) -> str | None:
        """Custom motif search dialog."""
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title("Find motif")
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)

        self._make_dialog_header(
            dialog,
            "Find motif or subsequence",
            "Exact search on the forward and reverse strands when applicable.",
        )

        container = ttk.Frame(dialog, style="Card.TFrame", padding=(22, 18, 22, 18))
        container.grid(row=1, column=0, sticky="nsew", padx=14, pady=14)
        container.columnconfigure(1, weight=1)

        ttk.Label(container, text="Motif", style="SectionHead.TLabel").grid(
            row=0, column=0, sticky="w", padx=(0, 12))

        motif_var = tk.StringVar()
        entry = tk.Entry(
            container,
            textvariable=motif_var,
            font=("Cascadia Code", 14),
            bg=C["surface"],
            fg=C["text"],
            insertbackground=C["accent"],
            relief="flat",
            highlightthickness=2,
            highlightbackground=C["border"],
            highlightcolor=C["accent"],
            width=28,
        )
        entry.grid(row=0, column=1, sticky="ew", ipady=6)

        ttk.Label(
            container,
            text="Example: ATGCGT  |  TATAAA  |  GCGC",
            style="Muted.TLabel",
        ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(8, 18))

        buttons = ttk.Frame(container, style="Card.TFrame")
        buttons.grid(row=2, column=0, columnspan=2, sticky="e")
        result: list[str] = []

        def submit() -> None:
            val = motif_var.get().strip()
            if val:
                result.append(val)
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=dialog.destroy,
                   style="Action.TButton").grid(row=0, column=0, padx=(0, 8))
        ttk.Button(buttons, text="Search", command=submit,
                   style="Primary.TButton").grid(row=0, column=1)

        dialog.bind("<Escape>", lambda e: dialog.destroy())
        dialog.bind("<Return>",  lambda e: submit())
        dialog.update_idletasks()
        dialog.geometry(f"+{self.root.winfo_rootx() + 140}+{self.root.winfo_rooty() + 140}")
        entry.focus_set()
        dialog.lift()
        self.root.wait_window(dialog)

        return result[0] if result else None

    def search_motif(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to analyze.")
            return
        motif = self._prompt_motif_search()
        if motif is None:
            return
        cleaned_motif = clean_sequence(motif)
        if not cleaned_motif:
            messagebox.showwarning("Empty motif", "Enter a valid motif.")
            return
        self._set_busy_status(f"Working: searching motif {cleaned_motif}...")
        try:
            hits = find_motif_hits(summary.cleaned, cleaned_motif, include_reverse=True)
            self._clear_analysis_highlights()
            if not hits:
                self._reset_search_state()
                self._redraw_features()
                self._set_results(f"No matches were found for: {cleaned_motif}")
                return
            matches = [
                SearchMatch(
                    start_nt=hit.start_nt,
                    end_nt=hit.end_nt,
                    title=f"Motif {cleaned_motif} ({'+' if hit.strand == '+' else '-'})",
                    details=(
                        f"Strand: {'forward' if hit.strand == '+' else 'reverse'}\n"
                        f"Match position: {hit.start_nt}-{hit.end_nt}\n"
                        f"Detected sequence: {hit.matched_sequence}\n"
                        f"Motif length: {len(cleaned_motif)} nt/aa\n\n"
                        "Use F6 and Shift+F6 to move between matches."
                    ),
                    kind="motif",
                    payload=hit,
                )
                for hit in hits
            ]
            self._set_search_matches(matches, f"Motif search: {cleaned_motif}")
            self._focus_current_search_hit()
        finally:
            self._clear_busy_status()

    def show_orfs(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to analyze.")
            return
        min_aa = self._get_orf_min_length()
        if min_aa is None:
            return
        self._set_busy_status(f"Working: scanning ORFs (minimum {min_aa} aa)...")
        try:
            try:
                hits = find_orfs(summary.cleaned, min_amino_acids=min_aa, include_reverse=True)
            except Exception as exc:
                messagebox.showerror("Could not analyze", str(exc))
                return
            if not hits:
                self._reset_search_state()
                self._clear_analysis_highlights()
                self._redraw_features()
                self._set_results(f"No ORFs were found with a minimum length of {min_aa} aa.")
                return
            matches = [
                SearchMatch(
                    start_nt=hit.start_nt,
                    end_nt=hit.end_nt,
                    title=f"ORF frame {hit.frame_label}",
                    details=(
                        f"Strand: {'forward' if hit.strand == '+' else 'reverse'}\n"
                        f"Length: {hit.aa_length} aa\n"
                        f"Protein:\n{format_fasta(hit.protein, width=60)}\n\n"
                        "Use F6 and Shift+F6 to navigate between ORFs."
                    ),
                    kind="orf",
                    payload=hit,
                )
                for hit in hits
            ]
            self._set_search_matches(matches, f"ORFs detected across the full sequence (minimum {min_aa} aa)")
            self._focus_current_search_hit()
        finally:
            self._clear_busy_status()

    def translate_selection(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to translate.")
            return
        selection_range = self._get_selection_clean_range()
        if selection_range is None:
            messagebox.showinfo("No selection", "Select the sequence region you want to translate.")
            return
        start_clean, end_clean = selection_range
        selected_sequence = summary.cleaned[start_clean:end_clean]
        try:
            frames = translate_frames(selected_sequence, include_reverse=True)
        except Exception as exc:
            messagebox.showerror("Could not translate", str(exc))
            return
        self._reset_search_state()
        lines = [
            "Translation", "",
            f"Selected range: {start_clean + 1}-{end_clean}",
            f"Selected length: {len(selected_sequence)} nt/aa", "",
            "Forward strand", "",
        ]
        for frame_label, protein in frames[:3]:
            lines.append(f"{frame_label}:")
            lines.append(format_fasta(protein, width=60))
            lines.append("")
        if len(frames) > 3:
            lines.extend(["Reverse strand", ""])
            for frame_label, protein in frames[3:]:
                lines.append(f"{frame_label}:")
                lines.append(format_fasta(protein, width=60))
                lines.append("")
        self._set_results("\n".join(lines).rstrip())

    def reverse_complement_selection(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to transform.")
            return
        selection = self._get_selection_indices()
        if selection is None:
            messagebox.showinfo("No selection", "Select a sequence region first.")
            return
        selected_text = self.sequence_text.get(selection[0], selection[1])
        cleaned_selection = clean_sequence(selected_text)
        if not cleaned_selection:
            messagebox.showinfo("No selection", "The selection does not contain valid nucleotides.")
            return
        try:
            transformed = reverse_complement(cleaned_selection)
        except Exception as exc:
            messagebox.showerror("Could not compute", str(exc))
            return
        self._push_operation_snapshot("reverse complement selection")
        if not self._replace_sequence_selection(transformed):
            return
        self._set_results(
            "Reverse complement applied to the selection.\n\n"
            f"Transformed length: {len(cleaned_selection)} nt"
        )

    # ------------------------------------------------------------------ circularize

    def _get_circular_scan_window(self) -> int | None:
        try:
            value = int(self.circular_scan_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("Invalid value", "Circular scan length must be an integer.")
            self.circular_scan_var.set(300)
            return None
        if value < 20:
            messagebox.showwarning("Invalid value", "Circular scan length must be at least 20 nt.")
            self.circular_scan_var.set(300)
            return None
        return value

    def _get_circular_max_mismatches(self) -> int | None:
        try:
            value = int(self.circular_mismatch_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("Invalid value", "Circular mismatches must be an integer.")
            self.circular_mismatch_var.set(2)
            return None
        if not 0 <= value <= 5:
            messagebox.showwarning("Invalid value", "Circular mismatches must be between 0 and 5.")
            self.circular_mismatch_var.set(2)
            return None
        return value

    def _get_circular_overlap_min(self) -> int | None:
        try:
            value = int(self.circular_overlap_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("Invalid value", "Circular overlap must be an integer.")
            self.circular_overlap_var.set(15)
            return None
        if value < 5:
            messagebox.showwarning("Invalid value", "Circular overlap must be at least 5.")
            self.circular_overlap_var.set(15)
            return None
        return value

    def _get_current_circular_match(self) -> SearchMatch | None:
        if not self.search_matches or self.search_index < 0:
            return None
        match = self.search_matches[self.search_index]
        return match if match.kind == "circular" else None

    def circularize_selection(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to analyze.")
            return
        selection_range = self._get_selection_clean_range()
        if selection_range is None:
            messagebox.showinfo("No selection", "Select an approximate component region first.")
            return
        min_overlap = self._get_circular_overlap_min()
        if min_overlap is None:
            return
        max_mismatches = self._get_circular_max_mismatches()
        if max_mismatches is None:
            return
        scan_window = self._get_circular_scan_window()
        if scan_window is None:
            return

        start_clean, end_clean = selection_range
        selected_sequence = summary.cleaned[start_clean:end_clean]
        self._set_busy_status("Working: scanning circularization candidates...")
        try:
            candidates = find_circular_candidates(
                selected_sequence,
                min_overlap=min_overlap,
                max_mismatches=max_mismatches,
                scan_window=scan_window,
                max_candidates=12,
            )

            if not candidates:
                self._reset_search_state()
                self._clear_analysis_highlights()
                self._redraw_features()
                self._set_results(
                    "Circularization: no candidates found.\n\n"
                    f"Analyzed range: {start_clean + 1}-{end_clean}\n"
                    f"Minimum overlap: {min_overlap}\n"
                    f"Maximum mismatches: {max_mismatches}\n"
                    f"Scan window from each end: {scan_window} nt\n"
                    "\n"
                    "No compatible terminal overlap was detected in the selection."
                )
                return

            matches: list[SearchMatch] = []
            for index, candidate in enumerate(candidates, start=1):
                global_start = start_clean + candidate.start_nt
                global_end = start_clean + candidate.end_nt
                details = (
                    f"Candidate {index}\n"
                    f"Terminal overlap: {candidate.overlap_nt} nt\n"
                    f"Overlap mismatches: {candidate.mismatch_count}\n"
                    f"5' trim: {candidate.trim_5prime_nt} nt\n"
                    f"3' trim: {candidate.trim_3prime_nt} nt\n"
                    f"Proposed component length: {len(candidate.component_sequence)} nt\n"
                    f"Trim origin: {candidate.note}\n"
                    f"Score: {candidate.score:.1f}\n\n"
                    "Use F6 and Shift+F6 to review candidates. "
                    "You can mark or crop the component using the buttons."
                )
                matches.append(
                    SearchMatch(
                        start_nt=global_start,
                        end_nt=global_end,
                        title=f"Circular candidate {index}",
                        details=details,
                        kind="circular",
                        payload=candidate,
                    )
                )
            self._set_search_matches(matches, "Circularize selection")
            self._focus_current_search_hit()
        finally:
            self._clear_busy_status()

    def mark_current_component(self) -> None:
        match = self._get_current_circular_match()
        if match is None:
            messagebox.showinfo("No candidate", "Run circularization first and choose a candidate.")
            return
        self._push_operation_snapshot("mark component as feature")
        feature = FeatureEntry(
            start_nt=match.start_nt, end_nt=match.end_nt,
            feature_type="misc_feature",
            label=f"circular_component_{len(self.features) + 1}",
            color=next(self._feature_colors),
            strand="+",
        )
        self.features.append(feature)
        self._refresh_feature_table()
        self._redraw_features()
        self.feature_table.selection_set(str(len(self.features) - 1))
        self._set_feature_focus(len(self.features) - 1)
        self._set_results(
            "Component marked as feature.\n\n"
            f"Label: {feature.label}\n"
            f"Range: {feature.start_nt}-{feature.end_nt}"
        )

    def crop_to_current_component(self) -> None:
        match = self._get_current_circular_match()
        if match is None:
            messagebox.showinfo("No candidate", "Run circularization first and choose a candidate.")
            return
        self._push_operation_snapshot("crop to current component")
        summary = self._get_summary()
        new_sequence = summary.cleaned[match.start_nt - 1 : match.end_nt]
        if not new_sequence:
            return
        new_features: list[FeatureEntry] = []
        for feature in self.features:
            overlap_start = max(feature.start_nt, match.start_nt)
            overlap_end   = min(feature.end_nt,   match.end_nt)
            if overlap_start > overlap_end:
                continue
            new_features.append(FeatureEntry(
                start_nt=overlap_start - match.start_nt + 1,
                end_nt=overlap_end   - match.start_nt + 1,
                feature_type=feature.feature_type,
                label=feature.label,
                color=feature.color,
                strand=self._normalize_feature_strand(feature.strand),
            ))
        self.features = new_features
        self.markers = [
            MarkerEntry(position_nt=marker.position_nt - match.start_nt + 1, label=marker.label)
            for marker in self.markers
            if match.start_nt <= marker.position_nt <= match.end_nt
        ]
        self._refresh_feature_table()
        self._refresh_marker_table()
        self._reset_search_state()
        self._clear_analysis_highlights()
        self.sequence_name.set(f"{self.sequence_name.get().strip() or 'sequence'}_component")
        self._set_editor_text(self._format_for_editor(new_sequence))
        self._set_results(
            "Sequence cropped to the current component.\n\n"
            f"New range: 1-{len(new_sequence)}\n"
            f"Retained features: {len(self.features)}"
        )

    # ------------------------------------------------------------------ map export

    def _collect_export_features(self) -> list[SequenceFeature]:
        return [
            SequenceFeature(
                start_nt=feature.start_nt, end_nt=feature.end_nt,
                feature_type=feature.feature_type, label=feature.label,
                strand=self._normalize_feature_strand(feature.strand),
            )
            for feature in self.features
        ]

    def export_linear_map(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to export.")
            return
        file_path = filedialog.asksaveasfilename(
            title="Export linear map", defaultextension=".svg",
            filetypes=[("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not file_path:
            return
        highlight_range = None
        current_match = self._get_current_circular_match()
        if current_match is not None:
            highlight_range = (current_match.start_nt, current_match.end_nt)
        try:
            export_linear_map_svg(
                file_path, summary.length, self._collect_export_features(),
                self.sequence_name.get().strip() or "GeneDraft",
                highlight_range=highlight_range,
            )
        except Exception as exc:
            messagebox.showerror("Could not export", str(exc))
            return
        self._set_results(f"Linear map exported.\n\nFile: {Path(file_path).name}")

    def export_circular_map(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to export.")
            return
        file_path = filedialog.asksaveasfilename(
            title="Export circular map", defaultextension=".svg",
            filetypes=[("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not file_path:
            return
        highlight_range = None
        current_match = self._get_current_circular_match()
        if current_match is not None:
            highlight_range = (current_match.start_nt, current_match.end_nt)
        try:
            export_circular_map_svg(
                file_path, summary.length, self._collect_export_features(),
                self.sequence_name.get().strip() or "GeneDraft",
                highlight_range=highlight_range,
            )
        except Exception as exc:
            messagebox.showerror("Could not export", str(exc))
            return
        self._set_results(f"Circular map exported.\n\nFile: {Path(file_path).name}")

    # ------------------------------------------------------------------ feature add/delete

    def _get_default_feature_range(self, summary: SequenceSummary) -> tuple[int, int]:
        selection_range = self._get_selection_clean_range()
        if selection_range is not None:
            start_clean, end_clean = selection_range
            return start_clean + 1, end_clean
        if summary.length <= 0:
            return 1, 1
        return 1, 1

    def _prompt_feature_details(self, start_nt: int, end_nt: int) -> tuple[str, str] | None:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title("New feature")
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)

        self._make_dialog_header(dialog, "Create feature", f"Range: {start_nt} - {end_nt}")

        container = ttk.Frame(dialog, style="Card.TFrame", padding=(18, 14, 18, 16))
        container.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        container.columnconfigure(1, weight=1)

        ttk.Label(container, text="Label", style="SectionHead.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 10))
        label_var = tk.StringVar(value=f"feature_{len(self.features) + 1}")
        label_entry = ttk.Entry(container, textvariable=label_var, style="Modern.TEntry", width=28)
        label_entry.grid(row=0, column=1, sticky="ew", pady=(0, 10))

        ttk.Label(container, text="Type", style="SectionHead.TLabel").grid(row=1, column=0, sticky="w", padx=(0, 10))
        type_var = tk.StringVar(value="misc_feature")
        type_combo = ttk.Combobox(container, textvariable=type_var, values=self.feature_type_options, state="normal", width=26)
        type_combo.grid(row=1, column=1, sticky="ew")

        ttk.Label(container, text="Choose a common type or enter a custom one.", style="Muted.TLabel").grid(row=2, column=0, columnspan=2, sticky="w", pady=(8, 16))

        buttons = ttk.Frame(container, style="Card.TFrame")
        buttons.grid(row=3, column=0, columnspan=2, sticky="e")
        result: dict[str, str] = {}

        def submit() -> None:
            result["label"] = label_var.get().strip() or f"feature_{len(self.features) + 1}"
            result["type"]  = type_var.get().strip() or "misc_feature"
            dialog.destroy()

        def cancel() -> None:
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=cancel, style="Action.TButton").grid(row=0, column=0, padx=(0, 8))
        ttk.Button(buttons, text="Save feature", command=submit, style="Primary.TButton").grid(row=0, column=1)

        dialog.protocol("WM_DELETE_WINDOW", cancel)
        dialog.bind("<Escape>", lambda event: cancel())
        dialog.bind("<Return>", lambda event: submit())
        dialog.update_idletasks()
        dialog.geometry(f"+{self.root.winfo_rootx() + 120}+{self.root.winfo_rooty() + 120}")
        dialog.lift()
        dialog.focus_force()
        label_entry.focus_set()
        label_entry.selection_range(0, tk.END)
        self.root.wait_window(dialog)

        if not result:
            return None
        return result["label"], result["type"]

    def _prompt_feature_editor_dialog(
        self,
        *,
        summary: SequenceSummary,
        dialog_title: str,
        dialog_heading: str,
        initial_label: str,
        initial_type: str,
        initial_start: int,
        initial_end: int,
        initial_strand: str = "+",
    ) -> dict[str, object] | None:
        C = self.colors
        dialog = tk.Toplevel(self.root)
        dialog.title(dialog_title)
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.configure(bg=C["bg"])
        dialog.columnconfigure(0, weight=1)

        self._make_dialog_header(
            dialog,
            dialog_heading,
            f"Available range: 1 - {max(1, summary.length)}",
        )

        container = ttk.Frame(dialog, style="Card.TFrame", padding=(18, 14, 18, 16))
        container.grid(row=1, column=0, sticky="nsew", padx=12, pady=12)
        container.columnconfigure(1, weight=1)

        ttk.Label(container, text="Label", style="SectionHead.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 10))
        label_var = tk.StringVar(value=initial_label.strip() or f"feature_{len(self.features) + 1}")
        label_entry = ttk.Entry(container, textvariable=label_var, style="Modern.TEntry", width=28)
        label_entry.grid(row=0, column=1, sticky="ew", pady=(0, 10))

        ttk.Label(container, text="Type", style="SectionHead.TLabel").grid(row=1, column=0, sticky="w", padx=(0, 10))
        type_var = tk.StringVar(value=initial_type.strip() or "misc_feature")
        ttk.Combobox(container, textvariable=type_var, values=self.feature_type_options, state="normal", width=26).grid(row=1, column=1, sticky="ew")

        ttk.Label(container, text="Start", style="SectionHead.TLabel").grid(row=2, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        start_var = tk.StringVar(value=str(initial_start))
        ttk.Entry(container, textvariable=start_var, style="Modern.TEntry", width=12).grid(row=2, column=1, sticky="w", pady=(10, 0))

        ttk.Label(container, text="End", style="SectionHead.TLabel").grid(row=3, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        end_var = tk.StringVar(value=str(initial_end))
        ttk.Entry(container, textvariable=end_var, style="Modern.TEntry", width=12).grid(row=3, column=1, sticky="w", pady=(10, 0))

        ttk.Label(container, text="Strand", style="SectionHead.TLabel").grid(row=4, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        strand_var = tk.StringVar(
            value=self.feature_strand_value_to_label.get(self._normalize_feature_strand(initial_strand), "Forward (+)")
        )
        ttk.Combobox(
            container,
            textvariable=strand_var,
            values=tuple(self.feature_strand_label_to_value.keys()),
            state="readonly",
            width=22,
        ).grid(row=4, column=1, sticky="w", pady=(10, 0))

        ttk.Label(
            container,
            text="Choose a common type or enter a custom one. Strand is stored as + or -.",
            style="Muted.TLabel",
        ).grid(row=5, column=0, columnspan=2, sticky="w", pady=(10, 16))

        buttons = ttk.Frame(container, style="Card.TFrame")
        buttons.grid(row=6, column=0, columnspan=2, sticky="e")
        result: dict[str, object] = {}

        def submit() -> None:
            try:
                start_nt = int(start_var.get().strip())
                end_nt = int(end_var.get().strip())
            except ValueError:
                messagebox.showwarning("Invalid value", "Start and end must be integers.", parent=dialog)
                return
            if start_nt < 1 or end_nt < 1:
                messagebox.showwarning("Invalid value", "Start and end must be at least 1.", parent=dialog)
                return
            if start_nt > end_nt:
                messagebox.showwarning("Invalid value", "Start cannot be greater than end.", parent=dialog)
                return
            if end_nt > summary.length:
                messagebox.showwarning("Out of range", f"End cannot exceed the current length ({summary.length}).", parent=dialog)
                return
            result["label"] = label_var.get().strip() or f"feature_{len(self.features) + 1}"
            result["type"] = type_var.get().strip() or "misc_feature"
            result["start_nt"] = start_nt
            result["end_nt"] = end_nt
            result["strand"] = self.feature_strand_label_to_value.get(strand_var.get(), "+")
            dialog.destroy()

        def cancel() -> None:
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=cancel, style="Action.TButton").grid(row=0, column=0, padx=(0, 8))
        ttk.Button(buttons, text="Save feature", command=submit, style="Primary.TButton").grid(row=0, column=1)

        dialog.protocol("WM_DELETE_WINDOW", cancel)
        dialog.bind("<Escape>", lambda event: cancel())
        dialog.bind("<Return>", lambda event: submit())
        dialog.update_idletasks()
        dialog.geometry(f"+{self.root.winfo_rootx() + 120}+{self.root.winfo_rooty() + 120}")
        dialog.lift()
        dialog.focus_force()
        label_entry.focus_set()
        label_entry.selection_range(0, tk.END)
        self.root.wait_window(dialog)

        if not result:
            return None
        return result

    def add_feature(self) -> None:
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "There is no sequence to annotate.")
            return
        start_nt, end_nt = self._get_default_feature_range(summary)
        initial_label, initial_type, initial_strand = self._get_feature_defaults_from_active_match()
        feature_details = self._prompt_feature_editor_dialog(
            summary=summary,
            dialog_title="New feature",
            dialog_heading="Create feature",
            initial_label=initial_label,
            initial_type=initial_type,
            initial_start=start_nt,
            initial_end=end_nt,
            initial_strand=initial_strand,
        )
        if feature_details is None:
            return
        self._push_operation_snapshot("add feature")
        feature = FeatureEntry(
            start_nt=int(feature_details["start_nt"]),
            end_nt=int(feature_details["end_nt"]),
            feature_type=str(feature_details["type"]).strip() or "misc_feature",
            label=str(feature_details["label"]).strip() or f"feature_{len(self.features) + 1}",
            color=next(self._feature_colors),
            strand=self._normalize_feature_strand(str(feature_details["strand"])),
        )
        self.features.append(feature)
        self._refresh_feature_table()
        self._redraw_features()
        self.feature_table.selection_set(str(len(self.features) - 1))
        self._set_feature_focus(len(self.features) - 1)
        self._set_results(
            "Feature added.\n\n"
            f"Label: {feature.label}\n"
            f"Type: {feature.feature_type}\n"
            f"Range: {feature.start_nt}-{feature.end_nt}\n"
            f"Strand: {feature.strand}"
        )

    def edit_selected_feature(self) -> None:
        feature_index = self._get_selected_feature_index()
        if feature_index is None:
            messagebox.showinfo("No feature", "Select a feature to edit.")
            return
        summary = self._get_summary()
        if not summary.cleaned:
            messagebox.showinfo("No sequence", "No sequence is loaded.")
            return

        feature = self.features[feature_index]
        updated = self._prompt_feature_editor_dialog(
            summary=summary,
            dialog_title="Edit feature",
            dialog_heading="Edit feature",
            initial_label=feature.label,
            initial_type=feature.feature_type,
            initial_start=feature.start_nt,
            initial_end=feature.end_nt,
            initial_strand=feature.strand,
        )
        if updated is None:
            return

        self._push_operation_snapshot("edit feature")
        feature.start_nt = int(updated["start_nt"])
        feature.end_nt = int(updated["end_nt"])
        feature.feature_type = str(updated["type"]).strip() or "misc_feature"
        feature.label = str(updated["label"]).strip() or feature.label
        feature.strand = self._normalize_feature_strand(str(updated["strand"]))
        self._refresh_feature_table()
        self._redraw_features()
        self.feature_table.selection_set(str(feature_index))
        self._set_feature_focus(feature_index)
        self._set_results(
            "Feature updated.\n\n"
            f"Label: {feature.label}\n"
            f"Type: {feature.feature_type}\n"
            f"Range: {feature.start_nt}-{feature.end_nt}\n"
            f"Strand: {feature.strand}"
        )

    def delete_selected_feature(self) -> None:
        feature_index = self._get_selected_feature_index()
        if feature_index is None:
            messagebox.showinfo("No feature", "Select a feature to delete.")
            return
        self._push_operation_snapshot("delete feature")
        removed = self.features.pop(feature_index)
        self._refresh_feature_table()
        self._redraw_features()
        self._set_results(
            "Feature deleted.\n\n"
            f"Label: {removed.label}\n"
            f"Type: {removed.feature_type}\n"
            f"Range: {removed.start_nt}-{removed.end_nt}\n"
            f"Strand: {removed.strand}"
        )


def main() -> None:
    root = tk.Tk()
    try:
        root.call("tk", "scaling", 1.2)
    except tk.TclError:
        pass
    GeneDraftApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
