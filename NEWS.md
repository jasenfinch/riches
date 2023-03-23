# riches 0.3.0

* Added a `NEWS.md` file to track changes to the package.

* Replaced the example data set with a molecular formula assigned [`AnalysisData`](https://jasenfinch.github.io/metabolyseR/reference/AnalysisData-class.html) class object.

* Added a `functionalEnrichment()` method for the [`RandomForest`](https://jasenfinch.github.io/metabolyseR/reference/RandomForest-class.html) S4 class.

* Removed the `EnrichmentParameters` S4 class and associated methods.

* Added example organism [FELLA](https://bioconductor.org/packages/release/bioc/html/FELLA.html) data for *Brachypodium distachyon* (`bdi`).

* Added the `organismData()` function as a convenience wrapper for loading and building KEGG graph data for [FELLA](https://bioconductor.org/packages/release/bioc/html/FELLA.html) enrichment analyses.

* Added accessor methods for the `FunctionalEnrichment` S4 class.

* Added example structural classifications.

* Added the `structuralEnrichment()` method for the [`RandomForest`](https://jasenfinch.github.io/metabolyseR/reference/RandomForest-class.html) S4 class.

* Functional and structural enrichment analyses can now be performed on explanatory features split by feature trends for binary comparisons or regression using the `split` argument of `functionEnrichment()` and `structuralEnrichment()`.
