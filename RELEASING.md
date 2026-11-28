# Making an AST release

This document describes how to cut a release of AST: choosing the version number, assembling the release notes, updating the documentation, and tagging.
The build mechanics (bootstrap, configure, `make dist`) are described in `README.md`.

## Overview

A release changes three tracked files:

- `version.number`, the single source of the version number for both CMake and Autotools;
- `ast.news`, the plain-text release notes that Starlink installs into its news directory;
- `sun_master.tex`, the master LaTeX source of SUN/210 and SUN/211, which carries a copy of the notes.

Release notes are not written at release time.
Each pull request adds a fragment to `newsfragments/` (see `newsfragments/README.md`), and `towncrier` assembles the fragments into `ast.news` when the release is made.
The V9.6.0 notes were the last to be written by hand.

## 1. Check the fragments

List the pull requests merged since the last tag and make sure that each one a user would notice has a fragment:

```shell
git log --merges --format='%h %s' v9.6.0..HEAD
ls newsfragments
```

CI checks every pull request for a fragment, but pull requests labelled `no news` skip the check and older ones may predate it.
Write any missing fragments now, as orphan fragments (`+<slug>.<type>`) if the pull request number is not worth recording.

Preview the section that the fragments will produce:

```shell
uvx towncrier build --draft --version "$(cat version.number)"
```

## 2. Choose the version number

AST uses `MAJOR.MINOR.PATCH`.

- Bump `PATCH` for a release that contains only fixes and build changes (`hardening`, `performance`, `serialisation`, `bugfix`, `doc` and `build` fragments).
- Bump `MINOR` if any `feature` or `behaviour` fragment is present, that is if there is new API or if any documented behaviour or numerical result changes.
- Bump `MAJOR` only for a change that breaks source or binary compatibility, such as removing a public function or changing a public struct.
  A `breaking` fragment that merely enforces a rule the documentation already stated (as V9.4.0 did for immutable Mappings) is a `MINOR` bump.

Write the new number to `version.number`.
It is the only place the version is stored; `src/version.h`, the CMake project version and the Autotools package version are all generated from it.

## 3. Build ast.news

```shell
uvx towncrier build --version "$(cat version.number)" --yes
```

This inserts a `Main Changes in VX.Y.Z` section below the preamble of `ast.news` and deletes the fragments.
Read the result and edit it as prose: add an overview paragraph if the release has a theme, merge fragments that describe the same thing, and reorder the items within a section so that the most important come first.
Once built, `ast.news` is the source of truth; the fragments are gone.

Conventions for `ast.news`:

- Headings are underlined with dashes of the same length and text is wrapped at 79 columns.
- Items start with `- `, continuation lines are indented by two spaces, and nested items are indented by two more.
- The preamble above the marker line (`The changes made in each release are listed below, newest first.`) is not touched.
  towncrier looks for that exact line, which is recorded in `towncrier.toml`, and inserts the new section after it.

## 4. Update sun_master.tex

The release-notes section of `sun_master.tex` ends with one `\subsection` per release, followed by a closing paragraph that begins `Programs which are statically linked`.
That sentence appears three times in the file; the release notes go before the last occurrence.

1. Add a new subsection before that paragraph:

   ```latex
   \subsection{\xlabel{changes}\xlabel{list_of_most_recent_changes}Changes Introduced in VX.Y.Z}
   The following describes the most significant changes which have
   occurred in the AST library between versions VP.Q.R and VX.Y.Z (the
   current version):
   ```

2. Remove `\xlabel{changes}\xlabel{list_of_most_recent_changes}` from the previous release's `\subsection`, and remove `(the current version)` from its introductory sentence, so that each label exists exactly once.
3. Copy each `ast.news` section as a `\subsubsection` in title case (`Bug Fixes`, `Build and Test Changes`) holding an `enumerate` list with one `\item` per news item.
   Nested items become an `itemize` list inside the `\item`.
4. Update `\stardocdate` near the top of the file to the release date.

LaTeX conventions:

- Write one sentence per line, as the existing release notes do.
- Escape underscores and percent signs (`AST\_\_MXAXES`, `\%s`) and use ``LaTeX quotes'' rather than straight double quotes.
- Write generic type suffixes as `astMask$<$X$>$` and configure options as `\texttt{-{}-with-yaml}`.
- If no TeX installation is available, at least check that `\begin` and `\end` pairs and braces balance.
  The documentation build in the Autotools workflow is the real check.

## 5. Verify

```shell
cmake -B build-dev -DAST_C_STANDARD=11
grep "current version" build-dev/src/version.h
cmake --build build-dev && ctest --test-dir build-dev --output-on-failure
```

If `starprolog` and `pdflatex` are installed, build the documentation as described in `README.md` under "Standalone Autotools releases" and check that `sun210.pdf` and `sun211.pdf` contain the new section.
Otherwise rely on the Autotools workflow, which builds the PDF-only distribution on every push.

## 6. Commit, tag and publish

Make two commits on a release branch, as previous releases did:

1. `Prepare release notes for VX.Y.Z`, containing `ast.news`, `sun_master.tex` and the deleted fragments.
2. `This will be X.Y.Z`, containing `version.number`.

Open a pull request, wait for CI, and merge it.
Then tag the merge commit on `master` and push the tag:

```shell
git tag -a vX.Y.Z -m "vX.Y.Z release" <merge commit>
git push origin vX.Y.Z
```

The Autotools workflow builds and tests the distribution and publishes the tarball to GitHub Releases once both platform jobs pass.
Check the release page for the tarball and the two PDFs.

Nothing needs to be done to `ast.news` after the release.
The next `towncrier build` inserts its section above the one just written.
