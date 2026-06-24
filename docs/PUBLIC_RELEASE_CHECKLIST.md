# Public Release Checklist

Use this checklist before publishing TranscriptoScope for free public download.

## 1. Confirm Identity And License

- Confirm the copyright holder in `LICENSE`.
- Confirm whether the app should use MIT, BSD-3-Clause, GPL-3.0, or an institution-required license.
- Ask the institution whether Southern University A and M needs to be listed as copyright holder, funder, or affiliation.

## 2. Verify Third-Party Terms

- Keep `THIRD_PARTY_NOTES.md` in the release.
- Verify exact redistribution terms for bundled Ensembl-derived annotation and GO resources.
- Confirm the public release does not include `gene_sets/cache/*_kegg.csv` unless explicit KEGG redistribution permission or licensing has been obtained.
- Document KEGG as optional user-initiated REST access with local user caching, not as a bundled MIT-licensed resource.
- Preserve package citations for DESeq2, fgsea, R, Shiny, and other major dependencies.

## 3. Create The GitHub Repository

- Create a public repository named `TranscriptoScope`.
- Upload the contents of `public_release/TranscriptoScope`.
- Confirm `CITATION.cff` points to `https://github.com/Shettima123/TranscriptoScope`.
- Add screenshots from `docs/images`.
- Enable GitHub Issues.

## 4. Create The GitHub Release

- Tag: `v0.4.5`
- Title: `TranscriptoScope Windows v0.4.5`
- Upload release asset:
  - `public_release/release_assets/TranscriptoScope_Windows_v0.4.5.zip`
- Paste the release body from:
  - `public_release/release_assets/GITHUB_RELEASE_BODY_v0.4.5.md`
- Include the SHA256 checksum.

## 5. Enable Zenodo DOI

- Log in to Zenodo.
- Enable GitHub integration for the repository.
- Publish a GitHub release so Zenodo archives it.
- Edit Zenodo metadata, author name, affiliation, keywords, license, and description.
- Copy the Zenodo DOI back into `README.md` and `CITATION.cff`.

## 6. Prepare A Research Landing Page

The landing page should include:

- download link;
- installation steps;
- example dataset;
- screenshots;
- citation;
- validation summary;
- privacy statement;
- support link; and
- known limitations.

## 7. Optional Microsoft Store Submission

After GitHub and Zenodo are stable:

- create a Microsoft Partner Center developer account;
- prepare app name, category, privacy URL, support URL, screenshots, package description, and release notes;
- decide whether to submit a packaged installer or link to the GitHub-hosted installer;
- verify code-signing requirements for the chosen package type.
