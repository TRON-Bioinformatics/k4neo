# Change Log

## v1.0.1

### Added

* mkdocs based annotation

### Changes

* Updated `tronmake-kmer-pipeline` submodule to track new major version v3.0.0

### Fixed


## v1.0.0

### Added

* Refactores and implemented new handling of k-mer indices, aiming to improve maintainability and modularity
  + `KmerIndex` and `KmerMetaIndex` are modular, use a pydantic-based design, support multiple index types, strict validation, and helper methods for querying and mapping indices


### Changed

* Refactored metadata backend to use SQLite3 instead of TinyDB
* Refactored the `Annotator` class to remove direct database file handling from its constructor and moved database interaction to be passed explicitly as a Queries object during annotation, improving separation of concerns and flexibility.
* Consolidated and simplified annotation methods: removed `_annotate_studies`, `_annotate_sample_metadata`, and `_annotate_counts` in favor of a more streamlined annotation flow within `annotate_cts`, and introduced `_count_aggregation` for count aggregation logic.


### Fixed

* Improved logging messages
* Added unit/integration tests
* Improved formatting
