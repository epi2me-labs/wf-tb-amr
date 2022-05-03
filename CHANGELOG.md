# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [v1.0.0]
### Changed
- Major restructuring to use bcftools mpileup
- Phasing with whatshap
- VCF files used throughout
- Unit tests
- Parsing of WHO catalogue xlsx file
- Tidying up of reference files
- Added option to set strand bias filter level
- Added versions for primers
### Fixed
- Division by 0 errors for RBK and strand bias
- Sample name given by user can now be an integer
### Added
- Single sample reports

## [v0.0.13]
### Fixed
- Update to primers

## [v0.0.12]
### Fixed
- Typos
- Schema consistency

## [v0.0.11]
### Changed
- Bumping aplanat to >=v0.6.1 for bootstrap 4

## [v0.0.10]
### Fixed
- Schema issues with labslauncher

## [v0.0.9]
### Changed
- Improved report
- Some rationalisation of ancillary files
- Schema changes for labslauncher
### Fixed
- Fastq ingress reverted & unclassified ignored in `main.nf`

## [v0.0.8]
### Changed
- New code to genotype variants with pileup
- No medaka or nanopolish steps
- Updates to schemas to bring in line with other workflows
- Thresholds for negative and positive controls

## [v0.0.7]
### Added
- Fastqingress module for common handling of (possibly
  multiplexed) inputs.
- Optimized container size through removal of various
  conda cruft.
### Changed
- Use mamba by default for building conda environments.
- Cut down README to items specific to workflow.
### Fixed
- Incorrect specification of conda environment file in Nextflow config.

## [v0.0.6]
### Changed
- Explicitely install into base conda env

## [v0.0.5]
### Added
- Software versioning report example.

## [v0.0.4]
### Changed
- Version bump to test CI.

## [v0.0.3]
### Changed
- Moved all CI to templates.
- Use canned aplanat report components.

## [v0.0.2]
### Added
- CI release checks.
- Create pre-releases in CI from dev branch.

## [v0.0.1]

First release.
