# Changelog
## [0.2.0] - 2025-01-12
### Changed
- [***breaking***] Simplify command line interface (`-I` -> `-i`, `-S` -> `-s`), options related genome copy estimation are removed.
### Added
- Add `--skip-melon` for skipping genome copy estimation.


## [0.1.3] - 2024-11-15
### Changed
- Change default mode (without `--plasmid`) to classifying all reads that carry ARGs.


## [0.1.2] - 2024-11-01
### Changed
- Use total bitscores per cluster for determining final HSPs. This help reduce false positives in certain cases.
- Switch from `dok_matrix` to `coo_matrix` for speed.
### Fixed
- Make database size checking more flexible.
- Reduce peak memory usage.


## [0.1.1] - 2024-09-03
### Fixed
- Fix a bug causing all reads not being classified when `--plasmid` is not given.


## [0.1.0] - 2024-07-21
### Added
- Add option `--plasmid` for controlling whether plasmid reads should be classified.
- Reduce plasmid detection cutoff from 0.9 to 0.5.
- Increase max overlap identity filtering cutoff from 0.01 to 0.05.


## [0.0.1] - 2024-06-30
### Added
- First release.
