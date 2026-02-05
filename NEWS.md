# BicliqueRT 1.0.0

* Initial CRAN submission.
* Replaced `biclust` dependency with embedded Bimax C implementation.
* Added generalized `biclique.decompose()` and `clique_test()` interface for
  arbitrary null hypotheses and exposure mappings.
* Fixed code quality issues (`T`/`F` replaced with `TRUE`/`FALSE`,
  `cat()` replaced with `message()`).
* Added unit tests via `testthat`.
* Added introductory vignette.
* Added `CITATION` file referencing Puelz et al. (2022, JRSSB).
