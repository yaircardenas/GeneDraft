# GitHub and Zenodo Checklist for GeneDraft

This checklist is for shipping GeneDraft as a small software preprint plus a citable software release.

## 1. Make the repository publication-ready

- Choose the final repository name, ideally `GeneDraft`.
- Add a `LICENSE` file in the repo root.
- Finalize `README.md` with scope, installation, intended use, and limitations.
- Rename `CITATION.cff.draft` to `CITATION.cff` after metadata is final.
- Confirm that `requirements.txt` is current.
- Confirm that the main launch command in the README actually works on a clean environment.
- Add at least one screenshot of the application for the README and release page.

## 2. Freeze the first citable software version

- Decide which version the preprint will describe, for example `v1.0.0`.
- Write release notes for that version.
- Summarize the main capabilities in the release notes, not every UI detail.
- Tag the release only after the repository metadata and citation files are in place.

## 3. Prepare the software metadata

- Finalize author name format.
- Finalize affiliations.
- Finalize corresponding author email.
- Confirm license choice.
- Add the final GitHub repository URL.
- Add keywords consistently across `README`, `CITATION.cff`, and the preprint.
- Decide whether the public-facing project URL will be the GitHub repo or a separate site.

## 4. Connect GitHub and Zenodo

- Create or log into Zenodo.
- Enable GitHub archiving for the GeneDraft repository in Zenodo.
- Publish the GitHub release only after Zenodo integration is active.
- Verify that Zenodo archives the tagged release and mints a version-specific DOI.
- Record both the version-specific DOI and the concept DOI.

## 5. Update files after DOI minting

- Put the Zenodo badge in `README.md`.
- Replace DOI placeholders in `CITATION.cff`.
- Replace repository and DOI placeholders in the preprint.
- Add the software DOI to the GitHub release notes.
- Optionally add repository and DOI information to the app About dialog.

## 6. Preprint alignment check

- Make sure the title in the preprint matches the software positioning.
- Make sure the abstract and Introduction describe GeneDraft as a workflow-oriented local workbench, not as a novel algorithmic method.
- Keep the scope statement aligned across the preprint, README, and citation file.
- Keep the limitations section explicit and honest.
- Cite both the preprint and the Zenodo release once both exist.

## 7. Recommended minimum repository contents

- `README.md`
- `LICENSE`
- `CITATION.cff`
- `requirements.txt`
- main application file
- at least one screenshot or figure
- release notes for `v1.0.0`

## 8. Metadata still pending in this project

- bioRxiv DOI
- funding statement
- final public release notes
- confirmation of which operating systems have been tested in practice

## 8a. Current identifiers

- GitHub repository: `https://github.com/yaircardenas/GeneDraft`
- Zenodo version DOI (`v1.0.0`): `https://doi.org/10.5281/zenodo.19445034`
- Zenodo concept DOI: `https://doi.org/10.5281/zenodo.19445033`
- bioRxiv DOI: pending

## 9. Suggested citation block

Preprint:

[Authors]. GeneDraft: a local sequence workbench for rapid editing, annotation, and exploratory molecular analysis. bioRxiv, 2026, [doi pending].

Software:

[Authors]. GeneDraft (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.19445034.
