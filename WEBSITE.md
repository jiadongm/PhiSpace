# Automatic Stereo-seq publishing

The root workflow `.github/workflows/stereoseq-pages.yaml` checks pull requests
and builds and publishes pushes to `main`. It also supports a manual run from
the Actions tab. Only Stereo-seq is rebuilt; the other pages and shared assets
are copied from the committed `docs/` website. Changes to other vignette sources
will not be published until those vignettes are added to this workflow.

## One-time GitHub setup

After these files are pushed, select **Settings → Pages → Build and deployment
→ Source → GitHub Actions**. The workflow uses GitHub's built-in token, with
deployment permissions only in the deployment job. No personal token is needed
as a workflow secret. Pull requests build a downloadable website artifact but
do not deploy. A failed build leaves the previous deployment in place.

The old `pkg/.github/workflows/pkgdown.yaml` was not discoverable by GitHub and
has been replaced by the root workflow.

## Build locally

From the repository root, with R 4.5 or later and Pandoc available, install the
package dependencies plus the R packages listed in the workflow. pkgdown 2.2.0
is installed from its explicit CRAN archive URL to match the existing website's
shared assets. Serialization uses `qs2` from CRAN, not the archived `qs` package.
The vignette reads and writes `.qs2` files with `qs_read()` and `qs_save()`.
Old `.qs` files must be replaced with the migrated files from Dropbox; renaming
the old files does not convert their format. Then run:

```bash
Rscript --vanilla scripts/build-stereoseq.R
```

The Stereo-seq script downloads the public archive, extracts only five required
inputs, and checks their SHA-256 hashes. These include the helper R script and
precomputed bridge annotations. It creates fresh analysis outputs on every run.
The download cache is `Test/stereoseq-inputs`; rendered pages are written to
`Test/site-stereoseq`. Both locations are git-ignored. A session information file
is included in the website artifact. The build refuses to overwrite an existing
site directory; choose a new destination for subsequent runs:

```bash
PHISPACE_SITE_DIR=Test/site-stereoseq-next Rscript --vanilla scripts/build-stereoseq.R
```

Optional variables:

- `PHISPACE_STEREOSEQ_ARCHIVE`: an existing ZIP file to use instead of downloading.
- `PHISPACE_STEREOSEQ_INPUTS`: an alternative verified input-cache directory.
- `PHISPACE_SITE_DIR`: a new directory for the staged website.
- `PHISPACE_STEREOSEQ_DATA`: the extracted data folder when running the vignette
  directly. The automation sets this to a fresh temporary analysis directory.

If upstream inputs change, checksum verification fails. Review those changes
before updating `scripts/stereoseq-inputs.tsv`; changing the manifest also changes
the Actions data-cache key. Computed analysis outputs are never cached by Actions.

## Check deployment

In Actions, confirm both `build` and `deploy` succeed. Open the live Stereo-seq
article and verify it loads `qs2`, uses `.qs2` files, and calls annotation with
`refAssay = "log1p"`. Download the `stereoseq-website` artifact to inspect the
complete staged website if needed.

Dependency installation and resource limits still need validation on the first
GitHub-hosted run; local execution does not test GitHub's runner or permissions.
