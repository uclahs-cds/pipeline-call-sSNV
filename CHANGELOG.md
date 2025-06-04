# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [8.2.0] - 2025-05-01

### Added

- Use `methods.setup_process_afterscript()` for process logs
- Output pipeline parameters to log directory using `store_object_as_json`
- Add `panel_of_normals_vcf` for MuTect2.
- Add `a_mini-all-tools-vcf-input` to `nftest`
- Add option to input VCFs for intersection
- Add VAF stripplot for variant overlaps across tools

### Changed

- Key changes to Strelka2 and SomaticSniper VAF calculations
- Switch to generalized resource handling
- Update NFTest paths
- Fix single tool run logic
- Update PipeVal 4.0.0-rc.2 -> 5.1.0
- Update GATK 4.5.0.0 -> 4.6.1.0
- Update bam_readcount 0.8.0 -> 1.0.1
- Update MuSE 2.0.4 -> MuSE 2.1.2
- Update BCFtools 1.17 -> 1.21
- Update NFtest for new tool versions

### Fixed

- Avoid warning when run without panel of normals

## [8.1.0] - 2024-05-15

### Added

- Add workflow to build and publish documentation to GitHub Pages
- Add workflow to run Nextflow configuration regression tests
- Add one regression test
- Add workflow to respond to "/fix-tests" comments

### Changed

- Update M64 resource allocations
- Update resource allocations

### Fixed

- Grant explicit permissions for Nextflow configuration test workflow
- Update Nextflow configuration test workflows

## [8.0.0] - 2024-01-29

### Changed

- `Strelka2` retry triggered by error code `1`
- Pass reference index to `SomaticSniper` processes
- Use external `indexFile` function
- Update submodules
- Update GATK `v4.4.0.0` to `v4.5.0.0`
- add `F8.config`

## [8.0.0-rc.1] - 2023-12-13

### Changed

- Sample names sanitized for all output
- Sample names parsed from input BAMs
- Update `MuSE` to `v2.0.4`
- Resource limit check now from submodule
- Add BigDecimal to `check_if_number` validation

## [7.0.0] - 2023-10-18

### Changed

- Use `bzip2` directly for compression

## [7.0.0-rc.2] - 2023-10-05

### Added

- Add .github/CODEOWNERS
- Add check for MuSE or Mutect2 on F2 node

### Changed

- Resource allocations changed for F32 and F72
- Update `MuSE` to `v2.0.3`
- Reorder all VCFs before intersection
- Move `filter_VCF_BCFtools` to `common.nf`
- Fix blarchive compression log output directory
- Delay readcount compression until original file is no longer needed

## [7.0.0-rc.1] - 2023-08-28

### Added

- Custom resource allocation updates through configuration parameters
- Add assertions to `nftest`
- Add compression of `SomaticSniper` `bam-readcount` output and move to `intermediate` directory
- Add `ncbi_build` parameter
- Add conversion of concatenated VCF to MAF
- Add concatenation of consensus variants to one VCF
- Add variant intersection Venn diagram
- Add regions filter to variant intersections
- Add second BCFtools step to create full presence/absence variant table (including private)
- Add workflow to create a `consensus.vcf` that includes SNVs found by two or more variant callers
- Add `fix_sample_names_VCF`, tumor and normal sample IDs from input BAMs used in output VCFs
- Add `split_VCF_bcftools` to `Mutect2` workflow, separating SNVs, MNVs and Indels

### Changed

- Update plot-venn.R to work with all numbers of algorithms greater than two
- Fix CPU allocation behavior with Docker
- Remove redundant directories in Intersect log output directories
- Change compression of intersect MAF file to bzip2
- Update `README.md`
- Use `set_env` from `pipeline-Nextflow-config`
- Update resource allocation to include new processes
- Reconfigure `intersect_regions` to use all contigs except `decoy`
- Reconfigure `call_regions` to `intersect_regions`
- Update to BCFtools v1.17
- Keep `bam-readcount` output in `SomaticSniper` QC folder
- Update `MuSE` to `v2.0.2`
- Update to use sample ID from input BAM files (single tumor/normal BAM input only)
- Use BCFtools to compress PASS variants instead of bgzip
- Use BCFtools to extract PASS variants instead of awk
- Update to use external `run_validate_PipeVal`

## [6.0.0] - 2023-04-05

### Added

- Add Mutect2 flow chart
- Add plantUML action and MuSE flow chart
- Add NF-test

### Changed

- Update LearnReadOrientationModel allocated memory and cpus
- Update to GATK v4.4.0.0
- Update `MuSE` retry add memory to 48GB
- Changed `output_dir` to `output_dir_base` (`methods.config` and `main.nf`)

### Fixed

- Specify empty string as default for bgzip and tabix extra args

## [6.0.0-rc.1] - 2023-02-08

### Added

- Add QC output `filteringStats.tsv` from Mutect2's process `run_FilterMutectCalls_GATK`.
- Add `contamination_table` input to `input.yaml`. Contamination estimate table was generated from [CalculateContamination](https://gatk.broadinstitute.org/hc/en-us/articles/9570322332315-CalculateContamination) from GATK.

### Changed

- Update `README`: add Pipeline Steps and Tool descriptions
- Update to use `set_resources_allocation` from pipeline-Nextflow-config repo
- Update SAMtools to v1.16.1
- Switch Docker Hub images to GitHub packages.
- Remove redundant directories in log output directories and intermediate directories.
- Specify `task.index` in log output directories.

## [5.0.0] - 2022-10-04

### Added

- Add retry.config using git submodule from [pipeline-Nextflow-config](https://github.com/uclahs-cds/pipeline-Nextflow-config) to enable processes to retry with more memory.
- Add [MuSE](https://github.com/wwylab/MuSE) workflow.
- Add `pipeline-release.yaml` to [workflow](.github/workflows).

### Changed

- Change index_VCF_tabix process using git submodule from [pipeline-Nextflow-module](https://github.com/uclahs-cds/pipeline-Nextflow-module).
- Standardize output filenames using `generate_standardized_filename` module from [pipeline-Nextflow-module](https://github.com/uclahs-cds/pipeline-Nextflow-module/tree/main/modules/common/generate_standardized_filename).
- Change the `input.yaml` structure to allow separate sample_ids for normal and tumor BAMs.

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

### Added

- Add Mutect2's multiple samples option into the pipeline.
- Add germline resource to Mutect2.
- Add `--callRegion` option to Stelka2 algorithm.
- Add Manta and Strelka2's intermediate files.
- Add schema.config using git submodule from [pipeline-Nextflow-config](https://github.com/uclahs-cds/pipeline-Nextflow-config) to check if params are valid.

### Changed

- Change the input files to YAML and template.config.
- Change the `sample_name` to `sample_id`.
- Standardize the repository structure.

## [3.0.0] - 2022-03-01

### Added

- Add tumor_only_mode in mutect2 options.
- Add Mutect2's orientation bias filter.

### Changed

- Update .gitignore to exclude molecular files.
- Update F72.config to increase the compute efficiency.

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

### Added

- Add groups to docker run options.

### Deprecated

- Deprecate version 2.1.0.

### Fixed

- Add missing lines of checksum process.

## [2.1.0] - 2021-10-13

### Added

- Save logs for the Somaticsniper workflows.
- Add GPL2 license.
- Add pipeline information to the main script.
- Add the checksum file for the final output.
- Add the config file for F32 node.

### Changed

- Update the GATK version from 4.2.0 to 4.2.2.
- Standardize process names.
- Change config filenames to F2, F72, and M64.

## [2.0.0] - 2021-08-19

### Added

- Add pattern to strelka2's filter_vcf_pass publishDir.

### Changed

- Allow multiple algorithms to run in one pipeline run.

### Fixed

- Fix nextflow.config throwing Exception: string interpolation required double quotes.

## [1.5.0] - 2021-06-18

### Added

- Specify extra arguments for all GATK commands in the Mutect2 workflow.
- Add steps to index and compress final VCF files.

### Changed

- Mutect2 calls are now scattered by intervals split using GATK SplitIntervals.
- Non-canonical regions will be called by default when using Mutect2.
- Update branch name in CICD-base.yaml from master to main.

## [1.4.0] - 2021-04-27

### Added

- Save logs for the Strelka2 and Mutect2 workflows.
- Create and backfill CHANGELOG.

### Changed

- Mutect2 calls are now scattered by chromosomes.

## [1.3.0] - 2021-04-05

### Added

- Add the Mutect2 algorithm.

## [1.2.0] - 2021-03-11

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

[0.0.1-beta]: https://github.com/uclahs-cds/pipeline-call-sSNV/releases/tag/v0.0.1-beta
[1.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v0.0.1-beta...v1.0.0
[1.1.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v1.0.0...v1.1.0
[1.2.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v1.1.0...v1.2.0
[1.3.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v1.2.0...v1.3.0
[1.4.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v1.3.0...v1.4.0
[1.5.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v1.4.0...v1.5.0
[2.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v1.5.0...v2.0.0
[2.1.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v2.0.0...v2.1.0
[2.1.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v2.1.0...v2.1.1
[3.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v3.0.0-rc.1...v3.0.0
[3.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v2.1.1...v3.0.0-rc.1
[4.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v4.0.0-rc.1...v4.0.0
[4.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v3.0.0...v4.0.0-rc.1
[4.0.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v4.0.0...v4.0.1
[5.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v4.0.1...v5.0.0
[6.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v6.0.0-rc.1...v6.0.0
[6.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v5.0.0...v6.0.0-rc.1
[7.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v7.0.0-rc.2...v7.0.0
[7.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v6.0.0...v7.0.0-rc.1
[7.0.0-rc.2]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v7.0.0-rc.1...v7.0.0-rc.2
[8.0.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v8.0.0-rc.1...v8.0.0
[8.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v7.0.0...v8.0.0-rc.1
[8.1.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v8.0.0...v8.1.0
[8.2.0]: https://github.com/uclahs-cds/pipeline-call-sSNV/compare/v8.1.0...v8.2.0
