# Windows Packaging Guide

This project can be distributed in three levels of polish. Start with level 1, then move up once the workflow is stable.

## Level 1: Zip Folder With Current-User Installer

Use this while developing with collaborators.

1. Run `powershell -NoProfile -ExecutionPolicy Bypass -File scripts\build_release_zip.ps1`.
2. Send the generated ZIP from the `outputs` folder to testers.
3. Testers extract the zip and double-click `Install_TranscriptoScope.bat`.
4. Testers launch the app from the Desktop or Start Menu shortcut.

Pros: fastest path, easy to debug.

Cons: users still need R installed, and package installation happens on first install.

## Level 2: Windows Installer

Use this when the app workflow is stable.

For a Microsoft Store candidate installer, run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\build_store_installer.ps1
```

This builds `outputs\store_installer\TranscriptoScope_Windows_Store_v<VERSION>.exe` with a bundled R runtime and staged R package dependency library. See `docs\MICROSOFT_STORE_SUBMISSION.md` for Store-specific notes and validation status.

Recommended tools:

- Inno Setup for the installer
- RInno as a helper for Shiny app installer scaffolding
- `renv` to lock R package versions

Recommended release workflow:

```powershell
Rscript scripts\install_packages.R
Rscript scripts\smoke_test.R
```

Then in R:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

After `renv.lock` exists, build an installer that includes:

- The app folder
- A launcher shortcut
- A first-run dependency restore step
- Optional bundled R installer

Notes:

- DESeq2 is a Bioconductor package, so test the installer on a clean Windows machine or VM.
- RInno can include R with `include_R = TRUE`, but Bioconductor dependency handling still needs testing.
- For public downloads, code-sign the installer to reduce Windows SmartScreen warnings.

## Level 3: Polished Desktop App

Use this when you want a product-like app.

Architecture:

- Electron or Tauri front end
- Local R process as the analysis backend
- Same `R/deseq_helpers.R` workflow logic
- Bundled dependency environment

Pros: best user experience.

Cons: more engineering, harder packaging, more installer testing.

## Release Checklist

- Test with bundled example data.
- Test with at least one real dataset.
- Test missing sample names.
- Test duplicated gene IDs.
- Test decimal/normalized input rejection.
- Test no-replicate warnings.
- Test batch design with known non-confounded metadata.
- Test built-in GO enrichment for each bundled organism.
- Test standard and Rosby's Lab-style ORA.
- Confirm release ZIPs do not include `gene_sets/cache/*_kegg.csv` unless explicit KEGG redistribution permission or licensing has been obtained.
- Test optional KEGG enrichment with internet access and confirm the local user cache is created.
- Test ranked GO and KEGG pathway analysis, pathway plots, leading-edge downloads, and sample-level heatmaps.
- Confirm the `annotations/` and `gene_sets/` folders are included in the zip or installer.
- Test on a clean Windows VM.
- Confirm exported CSVs open in Excel.
- Confirm plots export into the result bundle.
- Confirm package licenses and attribution.
