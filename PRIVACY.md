# Privacy

TranscriptoScope is designed to run locally on the user's Windows computer.

## Data Handling

- Uploaded count matrices, metadata files, normalized expression matrices, fold-change tables, and result files are processed locally.
- The app does not intentionally upload user datasets to an external server.
- Result files are written only when the user downloads or exports them.

## Network Use

The Windows installer may use the internet to install R package dependencies when they are not already available. The app includes bundled Ensembl-derived annotation and GO resources for common workflows. If the user selects KEGG, the app may query the KEGG REST API and write a local user cache for the selected organism.

## Research Use

Users are responsible for following their institution's policies for protected, sensitive, unpublished, or regulated datasets. Do not upload restricted data into any public issue tracker when reporting bugs.

## Contact

For privacy or data-handling questions, open a private support request through the repository owner or institutional contact listed with the public release.
