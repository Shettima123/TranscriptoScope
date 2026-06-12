# Install On Windows

## Quick Install

1. Install R for Windows from https://cran.r-project.org/bin/windows/base/ if R is not already installed.
2. Extract this ZIP folder.
3. Double-click `Install_TranscriptoScope.bat`.
4. Launch `TranscriptoScope` from the Desktop or Start Menu shortcut.

The installer copies the app to `%LOCALAPPDATA%\TranscriptoScope`, installs required R packages, and registers shortcuts. The normal launcher starts the local app quietly and opens it in an Edge or Chrome app-style window when one of those browsers is available.

Creator: Dr. Abubakar Abdulkadir, Southern University A and M.

## If Package Setup Fails

Double-click `Install_Packages.bat` after checking that R is installed and connected to the internet. DESeq2 and fgsea are Bioconductor packages, so the first package setup can take several minutes.

## Run Without Installing

Double-click `Install_Packages.bat`, then double-click `Launch_TranscriptoScope.vbs`.

## Uninstall

Use the Start Menu uninstall shortcut, or double-click `Uninstall_TranscriptoScope.bat` from the installed app folder.
