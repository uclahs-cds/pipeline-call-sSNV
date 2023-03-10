# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Add NF-test
### Changed
- Update (roll-back) to use `MuSE` v2.0.0
- Update to use `MuSE` v2.0.1 with `MuSE sump` parallelization
- Update `MuSE` retry add memory to 48GB
- Changed `output_dir` to `output_dir_base` (`methods.config` and `main.nf`)

## [6.0.0-rc.1] - 2023-02-08
### Changed
- Update `README`: add Pipeline Steps and Tool descriptions
- Update to use `set_resources_allocation` from pipeline-Nextflow-config repo
- Update SAMtools to v1.16.1
- Switch Docker Hub images to GitHub packages.
- Remove redundant directories in log output directories and intermediate directories.
- Specify `task.index` in log output directories.

### Added
- Add QC output `filteringStats.tsv` from Mutect2's process `run_FilterMutectCalls_GATK`.
- Add `contamination_table` input to `input.yaml`. Contamination estimate table was generated from [CalculateContamination](https://gatk.broadinstitute.org/hc/en-us/articles/9570322332315-CalculateContamination) from GATK.

## [5.0.0] - 2022-10-04
### Changed
- Change index_VCF_tabix process using git submodule from [pipeline-Nextflow-module](https://github.com/uclahs-cds/pipeline-Nextflow-module).
- Standardize output filenames using `generate_standardized_filename` module from [pipeline-Nextflow-module](https://github.com/uclahs-cds/pipeline-Nextflow-module/tree/main/modules/common/generate_standardized_filename).
- Change the `input.yaml` structure to allow separate sample_ids for normal and tumor BAMs.

### Added
- Add retry.config using git submodule from [pipeline-Nextflow-config](https://github.com/uclahs-cds/pipeline-Nextflow-config) to enable processes to retry with more memory.
- Add [MuSE](https://github.com/wwylab/MuSE) workflow.
- Add `pipeline-release.yaml` to [workflow](.github/workflows).

## [4.0.1] - 2022-06-24
### Changed
- Update PR template.

### Fixed
- Fix the bug with M64 detection.
- Fix the bug with interval list input.

## [4.0.0] - 2022-06-13
### Changed
- Update `README.md` to apply the recent changes.
- Update the YAML input structure.

### Fixed
- Add the decoy file for `params.germline_resource_gnomad_vcf`.

## [4.0.0-rc.1] - 2022-05-13
### Changed
- Change the input files to YAML and template.config.
- Change the `sample_name` to `sample_id`.
- Standardize the repository structure.

### Added
- Add Mutect2's multiple samples option into the pipeline.
- Add germline resource to Mutect2.
- Add `--callRegion` option to Stelka2 algorithm.
- Add Manta and Strelka2's intermediate files.
- Add schema.config using git submodule from [pipeline-Nextflow-config](https://github.com/uclahs-cds/pipeline-Nextflow-config) to check if params are valid.

## [3.0.0] - 2022-03-01
### Changed
- Update .gitignore to exclude molecular files.
- Update F72.config to increase the compute efficiency.
### Added
- Add tumor_only_mode in mutect2 options.
- Add Mutect2's orientation bias filter.
### Fixed

## [3.0.0-rc.1] - 2022-01-07
### Changed
- Standardize the output directory.
- Standardize timestamp according to ISO8601.
- Update PipeVal to 2.1.6.
- Rename Docker images and remove Dockerfiles.
- Apply the config file standardization.
### Security
- Update GATK to 4.2.4.1 to address Log4j critical vulnerability [GHSA-jfh8-c2jp-5v3q](https://github.com/advisories/GHSA-jfh8-c2jp-5v3q) and other newly discovered log4j2 vulnerabilities.


## [2.1.1] - 2021-10-15
### Changed
### Added
- Add groups to docker run options.

### Fixed
- Add missing lines of checksum process.

### Deprecated
- Deprecate version 2.1.0.

## [2.1.0] - 2021-10-13
### Changed
- Update the GATK version from 4.2.0 to 4.2.2.
- Standardize process names.
- Change config filenames to F2, F72, and M64.

### Added
- Save logs for the Somaticsniper workflows.
- Add GPL2 license.
- Add pipeline information to the main script.
- Add the checksum file for the final output.
- Add the config file for F32 node.

## [v2.0.0] - 2021-08-19
### Changed
- Allow multiple algorithms to run in one pipeline run.

### Added
- Add pattern to strelka2's filter_vcf_pass publishDir.

### Fixed
- Fix nextflow.config throwing Exception: string interpolation required double quotes.

## [v1.5.0] - 2021-06-18
### Changed
- Mutect2 calls are now scattered by intervals split using GATK SplitIntervals.
- Non-canonical regions will be called by default when using Mutect2.
- Update branch name in CICD-base.yaml from master to main.

### Added
- Specify extra arguments for all GATK commands in the Mutect2 workflow.
- Add steps to index and compress final VCF files.

## [v1.4.0] - 2021-04-27
### Added
- Save logs for the Strelka2 and Mutect2 workflows.
- Create and backfill CHANGELOG.

### Changed
- Mutect2 calls are now scattered by chromosomes.

## [v1.3.0] - 2021-04-05
### Added
- Add the Mutect2 algorithm.

## [v1.2.0] - 2021-03-11
### Added
- Add the strelka2 algorithm.

### Changed
- Index files are assumed to exist in the same directory as bam/reference files.

## [1.1.0] - 2021-02-16
### Changed
- Pipeline rewritten in DSL2.
- Docker run as user.

## [1.0.0] - 2020-11-24
### Fixed
- Fix bug to correctly filter pileup file in varFilter step.

## [0.0.1-beta] - 2020-11-12
### Added
- The is the first beta release of the call-sSNV pipeline. It implements only 1 SNV caller, somatic sniper. Input and output validation and dynamic resource allocation is implemented.
