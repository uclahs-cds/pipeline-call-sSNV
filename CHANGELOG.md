# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Mutect2 calls are now scattered by intervals split using GATK SplitIntervals
- Non-canonical regions will be called by default when using Mutect2

### Added
- Ability to specify extra arguments for all GATK commands in the Mutect2 workflow

## [v1.4.0] - 2021-04-27
### Added
- Saved logs for the Strelka2 and Mutect2 workflows
- Created and backfilled CHANGELOG

### Changed
- Mutect2 calls are now scattered by chromosomes

## [v1.3.0] - 2021-04-05
### Added
- Added the Mutect2 algorithm.

## [v1.2.0] - 2021-03-11
### Added
- Added the strelka2 algorithm

### Changed
- Index files are assumed to exist in the same directory as bam/reference files

## [1.1.0] - 2021-02-16
### Changed
- Pipeline rewritten in DSL2
- Docker run as user

## [1.0.0] - 2020-11-24
### Fixed
- Fixed bug to correctly filter pileup file in varFilter step

## [0.0.1-beta] - 2020-11-12
### Added
- The is the first beta release of the call-sSNV pipeline. It implements only 1 SNV caller, somatic sniper. Input and output validation and dynamic resource allocation is implemented.
