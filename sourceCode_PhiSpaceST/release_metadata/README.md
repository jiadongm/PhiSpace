# release_metadata/

This folder holds the Zenodo / citation metadata files that need to live at
the **GitHub repository root** (alongside `pkg/`, `sourceCode_PhiSpaceST/`,
`README.md`, etc.) — not inside `sourceCode_PhiSpaceST/`.

When you sync `sourceCode_PhiSpaceST/` into your GitHub clone of
`jiadongm/PhiSpace`:

```bash
# Move (don't copy) the metadata files to the repo root.
mv release_metadata/CITATION.cff   ../CITATION.cff
mv release_metadata/.zenodo.json   ../.zenodo.json
# Then this folder can be removed.
rmdir release_metadata
```

If both `CITATION.cff` and `.zenodo.json` are present at the repo root,
Zenodo uses **`.zenodo.json`** for the release record metadata and ignores
`CITATION.cff` (`CITATION.cff` is still useful for GitHub's "Cite this
repository" UI and for cff-aware tools).

## Before tagging

1. Verify author list, affiliations, and ORCIDs in both files (currently
   only Jiadong Mao's ORCID is set; Choi and Lê Cao are name-only).
2. Update `version:` and `date-released:` to the actual release tag and
   date (currently `v1.1.0` and a placeholder date).
3. Update `description:` / `abstract:` to the manuscript Summary text once
   the journal version is finalised (currently a generic 3-sentence
   summary).
4. After Zenodo mints the DOI (post-release), add the DOI to:
   - `CITATION.cff` (add `doi: 10.5281/zenodo.XXXXXXX` field)
   - The manuscript Data and Code Availability statement
   - The Key Resources Table PhiSpace row

See also `../RELEASE_CANDIDATE_CHECKLIST.md` for the full pre-release checklist.
