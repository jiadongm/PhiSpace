# Automatic vignette publishing

The root workflow `.github/workflows/stereoseq-pages.yaml` checks pull requests
and builds and publishes pushes to `main`. It also supports a manual run from
the Actions tab. Stereo-seq and Visium are rebuilt; the other pages and shared assets
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
is installed from its explicit CRAN archive URL to match the existing website
assets. Both vignettes use `qs2` from the active CRAN repository; the archived
`qs` package is not required. Dependency resolution has been checked locally.
Then run:

```bash
export PHISPACE_SITE_DIR=Test/site-vignettes
Rscript --vanilla scripts/build-stereoseq.R
Rscript --vanilla scripts/build-visium.R
```

The Stereo-seq script downloads the public archive, extracts only five required
inputs, and checks their SHA-256 hashes. These include the helper R script and
precomputed bridge annotations. It creates fresh analysis outputs on every run.
Visium downloads its public folder archive, extracts four verified inputs (about
1.6 GB), and checks their SHA-256 hashes. Both scripts recompute the analyses
shown in their vignettes. The input caches are `Test/stereoseq-inputs` and
`Test/visium-inputs`; the combined website is written to `Test/site-vignettes`.
Each article includes a session information file in the website artifact.
Stereo-seq requires a new site directory; Visium appends to it without replacing
the Stereo-seq page. Separate R processes release memory between articles.
For another full build, choose a new destination:

```bash
export PHISPACE_SITE_DIR=Test/site-vignettes-next
Rscript --vanilla scripts/build-stereoseq.R
Rscript --vanilla scripts/build-visium.R
```

Optional variables:

- `PHISPACE_STEREOSEQ_ARCHIVE`: an existing ZIP file to use instead of downloading.
- `PHISPACE_STEREOSEQ_INPUTS`: an alternative verified input-cache directory.
- `PHISPACE_VISIUM_INPUTS`: an alternative verified Visium input-cache directory.
- `PHISPACE_VISIUM_ARCHIVE`: an existing Visium ZIP file to use instead of downloading.
- `PHISPACE_SITE_DIR`: a new directory for the staged website.
- `PHISPACE_STEREOSEQ_DATA`: the extracted data folder when running the vignette
  directly. The automation sets this to a fresh temporary analysis directory.
- `PHISPACE_VISIUM_DATA`: the data folder when running Visium directly.
  The Visium build script also uses a fresh temporary analysis directory.

For a Visium-only build, run `Rscript --vanilla scripts/build-visium.R`; without
`PHISPACE_SITE_DIR`, it stages the committed website in `Test/site-visium`.

If upstream inputs change, checksum verification fails. Review those changes
before updating the relevant `scripts/*-inputs.tsv`; changing a manifest also
changes the Actions data-cache key. Computed analysis outputs are never cached.

## Check deployment

In Actions, confirm both `build` and `deploy` succeed. Open the live Stereo-seq
and Visium articles and verify that both use `qs2` and `.qs2` data files. Download
the `vignette-website` artifact to inspect the complete staged website if needed.

Dependency installation and resource limits still need validation on the first
GitHub-hosted run; local execution does not test GitHub's runner or permissions.
