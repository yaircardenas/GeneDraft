# GeneDraft

GeneDraft: a local desktop workbench for DNA, RNA, and protein sequence editing and exploratory analysis

Yair Cárdenas-Conejo*1,2

Secretaría de Ciencia, Humanidades, Tecnología e Innovación, Gobierno de México, Ciudad de México, México

Centro Universitario de Investigaciones Biomédicas, Universidad de Colima, Colima, Colima, México

Corresponding author: Yair Cárdenas-Conejo*

ORCID: 0000-0002-0190-244X

## Abstract

GeneDraft is a local desktop workbench for small-to-medium molecular sequence tasks that are often split across separate editors, web services, and ad hoc scripts. The software integrates sequence editing, FASTA and GenBank import/export, motif search, open reading frame (ORF) detection, selection-based translation, local DNA/RNA secondary-structure prediction, feature annotation, heuristic circularization review, and BLAST entry points within a single graphical interface. In its current implementation, GeneDraft supports DNA, RNA, and protein sequences and combines local analysis with optional network-enabled functions including accession retrieval, web BLAST launching, and remote protein-annotation enrichment. Protein-oriented reports include physicochemical summaries, transmembrane-segment inspection, and heuristic expression-oriented interpretation, with optional enrichment through Protein-Sol, SoluProt, and Pfam/HMMER. GeneDraft is implemented in Python with a tkinter-based user interface, uses Biopython for sequence handling, and uses ViennaRNA for local thermodynamics-based folding of selected nucleotide regions. The software is intended for rapid exploratory work, construct review, teaching, and bench-adjacent sequence inspection rather than for genome-scale or high-throughput analysis. Its practical value lies not in the novelty of any single function, but in the combination of capabilities it brings into one local workspace. Its main contribution is therefore workflow integration: a local, editor-centered environment that reduces tool switching during routine sequence handling and preliminary interpretation.

## Keywords

sequence analysis; sequence editing; protein analysis; annotation; ORF detection; BLAST; FASTA; GenBank; desktop software; molecular biology

## Introduction

Routine molecular sequence work often consists of many small but essential operations: opening a FASTA or GenBank record, cleaning and validating sequence text, checking sequence type and composition, reviewing candidate motifs or open reading frames, translating a selected region, annotating features, and exporting the result for downstream analysis. Integrated desktop platforms have shown the value of bringing these tasks into a single environment [1]. At the same time, many everyday workflows remain distributed across separate editors, browser-based services, command-line tools, and one-off scripts. This fragmentation slows exploratory work, complicates teaching demonstrations, and increases the likelihood of generating inconsistent intermediate files.

Existing software covers many parts of this landscape. Mature integrated platforms such as Geneious and DNASTAR Lasergene/EditSeq are widely used in molecular biology [2,3], and open-source alternatives such as UGENE are available [4]. However, for users whose main need is rapid editor-centered inspection rather than maximal analytical breadth, even capable general-purpose environments may involve more interface overhead than is desirable for quick bench-adjacent tasks. In that context, the practical gap is often not the absence of algorithms but the absence of a lightweight local workspace that keeps frequent sequence operations close together.

GeneDraft was developed to address that narrower but common use case. The goal of the software is not to replace specialized genome browsers, large-scale annotation systems, comparative genomics environments, or production bioinformatics pipelines. Instead, GeneDraft is intended as a local-first desktop workbench for rapid sequence inspection, lightweight annotation, exploratory translation and folding, quick protein review, and practical export in a single interface. The central premise is that the combination of these routine capabilities within one environment can be useful even when none of the individual operations is itself methodologically novel.

The current version extends this editor-centered workflow beyond nucleotide handling by supporting protein-oriented review in the same environment. This is useful when users move directly from nucleotide inspection to translation-derived protein interpretation, or when they need a compact local report before deciding whether a sequence warrants more specialized downstream analysis.

## Software architecture and implementation

GeneDraft is implemented in Python and uses tkinter for the desktop graphical interface. Core sequence parsing, type-aware editing, translation, sequence summaries, GenBank feature handling, and file export are implemented with support from Biopython [5]. Local secondary-structure prediction for selected DNA or RNA regions is implemented through the ViennaRNA library [6]. The software follows a local-first design in which the main editing, annotation, translation, folding, and map-export functions operate without an internet connection.

From a user perspective, the application is organized as a persistent desktop workspace rather than as a collection of separate utilities. FASTA and GenBank records can be loaded directly, pasted into the editor, curated interactively, and exported back to standard formats. Imported GenBank features remain editable, and local state such as markers, recent files, settings, and unsaved-session recovery is preserved to support iterative exploratory work rather than a formal project-database model.

Optional network-enabled modules are reserved for tasks that benefit from external resources. These include retrieval of records by NCBI accession, web BLAST launching, and remote enrichment of protein analysis through public online services. Local BLAST workflows are also supported when BLAST+ executables are installed in the user environment and a compatible database has been configured or built from FASTA input [7].

## Functional scope

Functionally, GeneDraft combines three user-facing layers within the same interface: sequence handling, exploratory analysis, and annotation/output. The sequence-handling layer is centered on direct interaction with DNA, RNA, or protein text, including import from FASTA or GenBank, manual editing, normalization and validation, reverse-complement or DNA-to-RNA transformations where relevant, and continuous display of contextual metrics such as sequence length, GC content, cursor position, and selection size. This makes the application useful not only for file conversion, but for iterative curation of working sequences.

On top of that workspace, the exploratory-analysis layer addresses common bench-adjacent questions about a sequence rather than large-scale computation. Users can inspect motifs, identify candidate ORFs, translate selected regions across reading frames, and fold selected nucleotide segments locally to obtain dot-bracket structure, minimum free energy, and an SVG secondary-structure rendering through ViennaRNA [6]. A parallel annotation/output layer allows features and markers to be created, reviewed, refined, exported, and reused during GenBank writing, while heuristic circularization review, coordinate-oriented operations, linear and circular SVG map export, BLAST integration, accession retrieval, and operation history provide downstream paths from the same working context. In this sense, the scope of GeneDraft is not merely to open or translate sequences, but to support a coherent small-scale workflow from sequence intake to preliminary interpretation and export.

## Protein analysis workflow

One of the main additions in the current implementation is a dedicated protein-analysis workflow. This module accepts protein sequences directly or sequences obtained by translating nucleotide regions and generates a compact report intended for fast exploratory interpretation. The local analysis includes molecular weight, isoelectric point, instability index, GRAVY score, aromaticity, amino-acid composition, secondary-structure fractions, cysteine counts, and a simple transmembrane-segment scan based on a Kyte-Doolittle window [8].

GeneDraft also derives heuristic expression-oriented summaries from these local properties. These include a local solubility estimate, a coarse crystallizability score, a simplified estimate of heterologous expression favorability in *Escherichia coli*, and an interpretation of likely soluble, membrane-associated, or aggregation-prone behavior. These outputs are intended as lightweight decision-support cues for exploratory use and not as substitutes for dedicated protein-expression or structural-biology pipelines.

When network access is available, the protein-analysis report can be enriched asynchronously with remote predictors and domain annotation services. In the current version, these optional enrichments include Protein-Sol [9] and SoluProt [10] solubility scores and Pfam-domain detection through the EBI HMMER service [11], which builds on accelerated profile hidden Markov model search methods implemented in HMMER3 [12]. The local report is shown immediately, and remote annotations are inserted when available. This design keeps the workflow responsive while allowing best-effort access to external resources.

## Positioning and intended use

GeneDraft is intended for rapid exploratory sequence work, especially in settings where users need to move quickly between inspection, annotation, translation, lightweight protein review, and export. Representative use cases include small-to-medium scale molecular biology projects, construct review, preliminary annotation of cloned or assembled regions, classroom demonstrations, and quick inspection of viral or plasmid-like components before downstream analysis in more specialized platforms.

The software is particularly useful when the main need is to reduce transitions between independent tools rather than to maximize analytical depth in a single domain. Its value lies in integration, responsiveness, and low setup burden for common sequence tasks.

## Limitations

GeneDraft is not intended to replace genome-scale analysis environments, high-throughput annotation systems, formal comparative genomics workflows, or specialized structural-biology software. Circularization support is heuristic and based on terminal overlap rather than full-sequence alignment. Protein-expression and solubility outputs combine simple local rules with optional web-derived annotations, so they should be interpreted as exploratory decision-support cues rather than validated predictive endpoints. Likewise, local secondary-structure predictions provide rapid thermodynamics-based guidance for selected regions rather than comprehensive structural modeling. Local BLAST requires BLAST+ command-line tools to be installed and configured. Some optional functions depend on external services and therefore inherit the availability, response behavior, and interface stability of those services.

Although the software is designed for cross-platform desktop use, testing to date has been concentrated on Windows. As with any exploratory sequence-analysis application, biological interpretation and downstream analytical decisions remain the responsibility of the user.

## Conclusion

GeneDraft contributes a local, editor-centered environment for small-to-medium sequence tasks. Its strength lies less in introducing new algorithms than in bringing routine sequence editing, exploratory analysis, annotation, protein follow-up, and export into a single low-friction workspace. The software is therefore best viewed as a practical desktop workbench for rapid sequence review, construct-oriented work, and teaching contexts in which the combination of capabilities is more useful than maximizing analytical breadth in any single domain.

## Availability and source code

Project name: GeneDraft

Programming language: Python

Main dependencies: tkinter, Biopython, ViennaRNA, optional NCBI BLAST+

License: MIT License

Version described here: 1.0

Source code repository: https://github.com/yaircardenas/GeneDraft

Archived release DOI: https://doi.org/10.5281/zenodo.19445034

Software concept DOI: https://doi.org/10.5281/zenodo.19445033

## References

1. Kearse M, Moir R, Wilson A, Stones-Havas S, Cheung M, Sturrock S, et al. Geneious Basic: an integrated and extendable desktop software platform for the organization and analysis of sequence data. Bioinformatics. 2012;28(12):1647-1649. doi:10.1093/bioinformatics/bts199.
2. Geneious Prime pricing. Geneious. https://www.geneious.com/pricing. Accessed April 6, 2026.
3. Lasergene Molecular Biology sequence analysis software. DNASTAR. https://www.dnastar.com/software/lasergene/molecular-biology/. Accessed April 6, 2026.
4. Okonechnikov K, Golosova O, Fursov M, UGENE team. Unipro UGENE: a unified bioinformatics toolkit. Bioinformatics. 2012;28(8):1166-1167. doi:10.1093/bioinformatics/bts091.
5. Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-1423. doi:10.1093/bioinformatics/btp163.
6. Lorenz R, Bernhart SH, Honer Zu Siederdissen C, Tafer H, Flamm C, Stadler PF, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011;6:26. doi:10.1186/1748-7188-6-26.
7. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: architecture and applications. BMC Bioinformatics. 2009;10:421. doi:10.1186/1471-2105-10-421.
8. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 1982;157(1):105-132. doi:10.1016/0022-2836(82)90515-0.
9. Hebditch M, Carballo-Amador MA, Charonis S, Curtis R, Warwicker J. Protein-Sol: a web tool for predicting protein solubility from sequence. Bioinformatics. 2017;33(19):3098-3100. doi:10.1093/bioinformatics/btx345.
10. Hon J, Marusiak M, Martinek T, Kunka A, Zendulka J, Bednar D, et al. SoluProt: prediction of soluble protein expression in *Escherichia coli*. Bioinformatics. 2021;37(1):23-28. doi:10.1093/bioinformatics/btaa1102.
11. Potter SC, Luciani A, Eddy SR, Park Y, Lopez R, Finn RD. HMMER web server: 2018 update. Nucleic Acids Res. 2018;46(W1):W200-W204. doi:10.1093/nar/gky448.
12. Eddy SR. Accelerated Profile HMM Searches. PLoS Comput Biol. 2011;7(10):e1002195. doi:10.1371/journal.pcbi.1002195.

## Declarations

### Conflicts of interest

The author declares no conflicts of interest.

### Author contributions

Y.C.-C. conceived the project, developed the software, curated the documentation, and wrote the manuscript.

### Operating system support

Designed for Windows, Linux, and macOS; currently tested on Windows.

### Funding

This study was supported by a grant from the Secretaría de Ciencia, Humanidades, Tecnología e Innovación (CBF-2025-G-902).

### Data and software availability

GeneDraft source code is available under the MIT License at https://github.com/yaircardenas/GeneDraft. Archived releases are available through Zenodo at https://doi.org/10.5281/zenodo.19445034.
