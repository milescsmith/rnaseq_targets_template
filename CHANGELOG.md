# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.0.0] - 2023-02-17

Not sure what all has changed now.  Check the git diff if you are really curious.

## Feat
  - Replaced use of {drake} and {Rmarkdown} with {targets} and {quarto}
  - Imported changes from milescsmith/rnaseq_targets_pipeline, mixing it with
    the previous CLE analysis

## [2.11.0] - 2022-02-?

### Changed
  - `process_counts.limma` now returns the fit, efit, and differential testing
  results.  `create_results_list.limma` uses those results and only performs
  calculations for parts it needs that are missing 

## [2.10.0] - 2022-02-21

### Fixed
  - Trying to use non-project controls no longer fails to actually add the controls
   due to a grouping column not being present in in the metadata; instead, all
   missing columns are added to the metadata for the controls.  Also added
   options to specify which columns to use in determining what samples are
   controls.

## [2.9.0] - 2021-10-19

### Changed
  - Moved `project_parameters.R` to the project root folder.
  - Overhauled how clustering is performed.  Replaced using random forest
   proximities and the gap statistic and am instead clustering directly
   on a PCA of the variance-stabilized data (with the significant PCs
   selected by either PCAtools::findElbowPoint or PCAtools::parallelPCA),
   selecting the appropriate number of clusters by selecting the
   highest **k** that yields clusters with a Jaccard stability score
   above 0.6.
   
### Fixed
  - Fixed a few stray missing commas or package namespace declarations

## [2.8.0] - 2021-10-13

### Changed
  - Moved PCA of variance-stabilized data for clustering to be performed within
  the `ident_clusters()` function.  Added a "pca_explain_prop" argument to control
  which PCs are used for clustering (essentially, only use those PCs which explain
  a proportion of the variance that is above the given threshold)
  - Moved the volcano plot generation portion of the report to templates, which
  fixed their not being displayed
  - Changed `getRandomPalette` to select from a smaller list of usable palettes instead of from the entire list that {paletteer} contains.  Fixes the issue where two almost indistinguishable colors are chosen to represent a two-factor variable
### Fixed
  - in "report.Rmd", renamed all instances of `res` that were being used to collect 
  markdown text from knit children templates to `template_res` to avoid conflicting 
  with the `res` object respresenting results DEG analysis
  - Volcano plots of differential gene expression display again
  - Fix to keep the metadata a dataframe in the event that only one column of 
  data is used (instead of impliciting casting it to a vector)

## [2.7.0] - 2021-10-13

### Fixed
  - Several fixes for empty enrichment lists

## [2.6.0] - 2021-10-01

### Changed
  - More `%>%` to `|>` replacements
  - Added "extra_controls_metadata_sheet" and "extra_controls_metadata_skip" so that
    it is possible to use an Excel spreadsheet for the extra controls list

### Fixed
  - `generatePalettes` is able to use selecting fuctions from {tidyselect}
  - Deduplicated calculating differential expression and performing shrinkage
  - Correctly added the "degPathwayPlots.Rmd" template
  - Expanded the pattern used by str_split from " - " to " - |_vs_" when trying to
    identify comparison groups
  - `read_md_file` now properly uses the "skip_lines" parameter

## [2.5.1] - 2021-09-29

### Added
  - parameter for including extra columns from the metadata file

### Fixed
  - A few small typos in `01_import_funcs.R/import_metadata()`
  - Pass the sheet name to `01_import_funcs.R/read_md_file()` if the metadata
    file is an Excel spreadsheet

## [2.5.0] - 2021-09-29

### Added
  - Gene ontology and Reactome pathway enrichment analysis for differentially
    expressed genes
    
### Fixed
  -  Volcano plot for differential expression is again properly displayed


## [2.4.0] - 2021-09-27

### Changed
  - Removed {pheatmap}-based heatmaps in favor of {ComplexHeatmap}
  
### Fixed
  - The `primary_report` target now correctly produces a report, with proper formatting and all
  - Fixes for `groupedComplexHeatmap()` and `comparisonComplexHeatmap()`

## [2.3.0] - 2021-09-24

### Added 
  - The function `DimPlotHull` to handle PCA and UMAP plotting
  - The function `groupedComplexHeatmap` to handle heatmap plotting where the data is
    split by a group
  - The function `comparisonHeatMap` to handle heatmap plotting where there is
    a table of upregulated and one of downregulated genes to plot.
  - Add a `conditional_left_join` function to handle instances where joining should
    be conditional on some test (such as "is the right-hand table NULL?")
  - Several child report templates
  - New parameters in `project_parameters.R`:
    "sample_species", "heatmap_row_annotations", "comparison_groups", "row_annotations"
  - Add rule to use `generatePalettes` instead of the old `create_palettes`
  
### Changed
  - Switched many more instances of the {magrittr} pipe `%>%` to base R pipe `|>`
  - `group_pal` target now uses the `generatePalettes` function
  - `groupedHeatMap` no longer tries to use all factor columns to annotate rows, you must
    specify the row annotation columns
  - Calculating stats for within comparison_groups comparisons by iterating over targets
  - Performing random forest modeling by iterating of the targets
  - `module_gsea` now uses `clusterProfiler::enrichGO` instead of `clusterProfiler::enricher`
    and takes an optional `target_species` argument to handle mouse RNAseq datasets
  - Essentially rewrote report.rmd to include use of templates and iterating
  - Replace usage of {pheatmap} with {ComplexHeatmap}.  Pheatmap was causing
    an issue with plotting early (stupid gTable objects)
  - Rewrote the machine learning plotting functions
  - Now pass `corType = "bicor"` to `WGCNA::blockwiseModules` to avoid a bug in {WGCNA}
    where it uses the baseR `cor` instead of its own `cor` (which has different parameters)
  - iterate over `rf_classifier` instead of having one rule per comparison (and allows for
    parallel computation)
  - iterate over `modules_compare_with_stats`
  
### Fixed
  - `annotation_info` target now includes extracting the "cluster" column
  - `pca_results` and `umap_results` no longer try to pass a `cluster_info`
    argument
  - `process_counts.limma` now correctly returns a `res` list member and it
    properly calculates differential gene expression
  - `process_counts.edgeR` now correctly returns a `res` list member
  - `create_deg_tables` no longer shits the bed if you pass it a tibble
  - `getRandomPalette` now actually checks to see if one of the preferred palettes will work,

## [2.2.0] - 2021-09-02

### Added
  - No longer need to create palettes manually thanks to the new `generatePalettes`
    function
  - Functions to generate tables and heatmaps for the report
  
### Changed
  - Added more namespace declarations
  - Changed several anonymous functions to use the new R4.1-style lamba functions
  - Stole `filterThreshold` from {DESeq2} so that it works on {limma}/{edgeR} results
  
### Fixed
  - Currently working through the `analysis/report.rmd` to match the new targets
  - Fix for `read_md_file` so that importing the metadata now accounts for
    `NaN` values
  - Imported data now no longer thinks `initial_concentration_ng_ul` is a character
    column
  - `import_metadata` now actually filters based on the `filter_column` and `filter_value

## [2.1.0] - 2021-08-27

Merging the "generalized" branch

### Changed
  - Separated import of count files from processing of counts
  - removed global import of packages and instead each target loads the packages
    it requires
  - Moved random forest modeling to a new function that handles all steps
    - Also refactored random forest modeling to use {tidymodels}
  - set most targets to not rerun unless dependencies are new
  - Removed many instances where I was unnecessarily re-mutating variables
    to factors
  - Replaced several instances of the {magrittr}-style pipe to the new
    native R pipe
    
### Added
  - Added package namespace declarations in front of 
    all functions
  - Added generics functions for writing results to file
  - Added generics to extract module scores
  - Added generics for differential gene expression testing
  - Added generics for extracting metadata
  
### Fixed
  - Pipeline up to `primary_report` target


## [2.0.0] - 2021-08-13

Pulling changes from updates added during BLAST analysis

### Added
  - ability to manually remove samples
  - custom version of janitor::make_clean_names that optionally allows duplicate values
  
### Changed
  -Replaced some hard coded variables with abstraction.
    - Replaced comparing by disease_class with a variable `comparison_grouping_variable`
  - Palette generation is a little smarter
  - Module data is now split by a custom function (helps with NSE)
  - Update C5 MSigDb file
  
### Removed
  - Eliminated the excessive number of times I was mutating that comparison
    column into a factor.


## [1.2.0] - 2021-08-??

### Changed
  - Separated import of count files from processing of counts
  - Started adding package namespace declarations in front of 
    all functions
    
### Fixed
  - Pipeline up to `sva_graph_data` target

## [1.1.0] - 2021-04-22

### Changed
  - Switched to using a generic metadata template instead of trying to 
    customize the import_metadata function for each new dataset.

### Added
  - Metadata template


## [1.0.0] - 2021-03-08

### Added
  - Started project
  - CHANGELOG.md
  - README.md

### Changed
  - Rearranged directory layout
  - Split analysis plan into parts

[2.10.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.9.0...2.10.0
[2.9.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.8.0...2.9.0
[2.8.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.7.0...2.8.0
[2.7.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.6.0...2.7.0
[2.6.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.5.1...2.6.0
[2.5.1]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.5.0...2.5.1
[2.5.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.4.0...2.5.0
[2.4.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.3.0...2.4.0
[2.3.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.2.0...2.3.0
[2.2.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/2.1.0...2.2.0
[2.1.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/1.0.0...2.1.0
[2.0.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/1.2.0...2.0.0
[1.2.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/milescsmith/rnaseq_targets_pipeline/releases/tag/1.0.0
