# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [v1.0.16]
### Changed
- Updated documentation
- Reconciled workflow with wf-template v5.2.6
- Updated aplanat to v0.6.20

## [v1.0.15]
### Changed
- Updated documentation

### Added
- Memory requirements for each process

## [v1.0.14]
### Fixed
- Updated whole workflow to bring up-to-date with recent template changes

### Added
- Configuration for running demo data in AWS

## [v1.0.13]
### Fixed
- Updated sample sheet to expect a file

## [v1.0.12]
### Changed
- Updated description in manifest

## [v1.0.11]
### Changed
- Documentation
- Nextflow schema

## [v1.0.10]
### Changed
- Report card more intelligent to drop-outs
- Group 2 SNPs now included in reports
- CSV is more straight forward to interpret
- Single sample report uses results from previous reporting step
- Appendix report removed
- Group 3 variants no longer reported
- VariantDB - added ONT reviewed variants
### Added
- Added coverage details to CSV
- Warning to single sample report in the case of target failure
- ONT reviewed group 3 variants are now reported (moved to V3)

## [v1.0.9]
### Fixed
- No reads in NTC causing whatshap failure
- Annotation issue for MNPs
- Phasing issue where variant at pos 1 & 3 but not 2
### Updated
- Base image to v0.2.0

## [v1.0.8]
### Fixed
- Issue with phased variants
- Issue with report when there are unclassified reads
### Changed
- Fastqingress metadata map

## [v1.0.7]
### Added
- Group 2 and group 3 variants to csv
### Changed
- EMB M306I - now ignoring strand bias filter

## [v1.0.6]
### Fixed
- Issue with changes to ingress and reports
### Changed
- Better help text on cli
- Adding variants to csv output
- Keep only mapped reads
### Fixed
- Downgraded pomoxis to prevent issues with empty bams

## [v1.0.5]
### Changed
- Wording on reports
- bcftools mpileup now uses bed file to speed up processing
- `file` to `path` in nextflow process definitions
### Fixed
- Issue where Rv0678 was excluded from variant DB

## [v1.0.4]
### Added
- Crude downsampling step to expedite high coverage samples
### Fixed
- Error in phred strand bias calc if p=0
- Unclassified reads caused failure due to lack of barcode

## [v1.0.3]
### Added
- CSV output of final results for every sample
### Fixed
- Fixed issue with resistance calling
### Changed
- Some wording on reports
- Removed "_blank" from single sample reports for jupyter compatibility

## [v1.0.2]
### Changed
- Report formatting

## [v1.0.1]
### Changed
- New docs format

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
- Sample name given by user can now be an integer
### Added
- Single sample reports
- Additional report on less confident variants

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
