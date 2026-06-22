# Microsoft Store Submission Notes

These notes describe the project-specific path for publishing TranscriptoScope in the Microsoft Store.

## Recommended Store Route

Use the Microsoft Store Win32 MSI/EXE route rather than MSIX for the first Store submission.

Reason: TranscriptoScope is a local R/Shiny desktop workflow that launches a local server, writes user-level result bundles, and depends on an R runtime plus R/Bioconductor packages. The Win32 route lets the Store listing point to a conventional installer, while the app can continue to run as a per-user desktop application.

Microsoft's Win32 Store route requires an installer URL. The installer should be an `.exe` or `.msi`, offline, stable at the submitted URL, and should install only TranscriptoScope.

## Current Project State

Already available:

- Public source release folder.
- MIT license for TranscriptoScope source code.
- Third-party license notes.
- Privacy, support, security, citation, and release notes.
- Windows ZIP installer workflow for testers.
- Clean Windows 11 validation and offline launch validation after dependency setup.
- GitHub repository and Zenodo-ready metadata.
- Store-style Inno Setup `.exe` installer build script.
- Bundled R runtime and frozen R package library staging workflow.

Current Store-installer artifact:

- `outputs/store_installer/TranscriptoScope_Windows_Store_v0.4.3.exe`
- Microsoft Store package URL: `https://tscope743176685.blob.core.windows.net/releases/TranscriptoScope_Windows_Store_v0.4.3.exe`
- Azure Blob resource: storage account `tscope743176685`, container `releases`, resource group `rg-transcriptoscope-store`.
- GitHub Release mirror: `https://github.com/Shettima123/TranscriptoScope/releases/download/v0.4.3/TranscriptoScope_Windows_Store_v0.4.3.exe`
- Size: 181,738,901 bytes.
- SHA256: `80E49ABAF1590AE3735B62B24BFDA8F5BB5317968CCFD7719EBA3FBC46BE73BA`
- Built with Inno Setup 6.7.3.
- Includes bundled R 4.6.0 and 70 staged R package dependencies.
- Local temporary install test passed with no icons, bundled package preflight passed, local Shiny launch probe passed, and temporary uninstall completed.
- Refreshed host validation evidence: `work/store_host_validation_spaces/summary.txt`.
- The Store installer and SHA256 sidecar have been uploaded to the GitHub v0.4.3 release, and the Store installer is also hosted on Azure Blob Storage for a direct non-redirecting package URL.

Still required before Store submission:

- Test the `.exe` installer on a clean Windows 11 machine or VM with no prior R setup.
- Run a true offline installer and launch test.
- Decide whether to code-sign the installer before Store submission.
- Complete Microsoft Partner Center account verification and Store listing fields.

## Store-Ready Installer Requirements

Before submission, build a Store-ready installer that includes:

- TranscriptoScope application files.
- The app icon.
- A Start Menu shortcut and uninstall entry.
- A user-writable app data/result location.
- A pinned R runtime strategy:
  - preferred: bundled R runtime or app-private R distribution;
  - acceptable only if Store certification allows it: clear prerequisite detection with a graceful message when R is missing.
- A frozen R package library or offline package cache for DESeq2, Shiny, plotting, enrichment, and export dependencies.
- Third-party notices covering bundled R, CRAN, Bioconductor, GO/annotation resources, and any other redistributed components.
- No bundled KEGG cache unless redistribution permission is obtained.

## Submission Assets To Prepare

- App name: TranscriptoScope.
- Publisher name: final Partner Center publisher identity.
- Short description.
- Full description.
- Screenshots showing preprocessing, DGE, enrichment, pathway analysis, plots, and downloads.
- Square Store logo and app icon.
- Privacy policy URL.
- Support URL or email.
- Website or repository URL.
- Installer URL from a fixed GitHub Release asset or other stable HTTPS host.
- Version number matching the installer and manuscript/release notes.

A draft listing is available in `docs/MICROSOFT_STORE_LISTING_DRAFT.md`.

## Practical Submission Flow

1. Create or open a Microsoft Store developer account in Partner Center.
2. Reserve the app name `TranscriptoScope`.
3. Build a Store-ready `.exe` or `.msi` offline installer.
4. Test the installer on a clean Windows 11 machine with no prior R setup.
5. Upload the installer to a fixed GitHub Release asset.
6. Create the Store submission and provide the installer URL.
7. Fill listing, privacy, support, category, age rating, pricing, and availability.
8. Submit for Microsoft certification.
9. After approval, update README, manuscript availability statement, and release notes with the Store link.

## Recommended Next Engineering Task

Run the current Store installer on a clean Windows 11 VM with networking disabled after the installer is downloaded. If install, preflight, launch, and uninstall pass there, use the installer hash above as the candidate submitted through Partner Center.
