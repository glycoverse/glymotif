## Submission

This is a patch update from glymotif 1.0.0 to 1.0.1.
It fixes source installation on older R versions using C++14.
There are no changes to the public API or dependency requirements.

## Test environments and results

* win-builder, R-devel (2026-10-20 r90728 ucrt), Windows Server 2022:
  0 errors, 0 warnings, 0 notes. PDF and HTML manuals passed.
* GitHub Actions, Ubuntu 24.04, R 4.6.1, with incoming checks and PDF manual:
  0 errors, 0 warnings, 1 note. HTML validation was skipped because the
  runner has no HTML Tidy executable; HTML validation passed on win-builder.
* GitHub Actions, Ubuntu, R 4.1.3 and R 4.2.3, with incoming checks and PDF
  manual: 0 errors, 0 warnings, 1 note each about installed size (7.5 MB,
  primarily native code and vignettes).
* Regular five-platform CI passed on macOS release, Windows release,
  Ubuntu release, Ubuntu devel, and Ubuntu oldrel-1.

All three full GitHub checks passed all 666 test expectations with no
failures, warnings, or skips. They used released dependencies from CRAN/RSPM:
glyrepr 1.1.0, glyparse 0.8.1, glydraw 0.9.0, Rcpp 1.1.2, and BH 1.90.0-1.
URL checking reported all 28 URLs correct. README was rebuilt without changes.

## Reverse dependencies

The CRAN source package index on 2026-10-21 contains no direct reverse
Depends, Imports, LinkingTo, or Suggests dependencies of glymotif.
