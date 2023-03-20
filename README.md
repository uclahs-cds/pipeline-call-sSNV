# pipeline-call-sSNV

- [call-sSNV](#pipeline-call-ssnv)
  - [Overview](#overview)
  - [How To Run](#how-to-run)
  - [Flow Diagrams](#flow-diagrams)
  - [Pipeline Steps](#pipeline-steps)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Testing and Validation](#testing-and-validation)
    - [Test Data Set](#test-data-set)
    - [Performance Validation](#performance-validation)
  - [References](#references)
  - [License](#license)

## Overview
The call-sSNV nextflow pipeline performs somatic SNV calling given a pair of tumor/normal BAM files. Four somatic SNV callers are available: SomaticSniper, Strelka2, Mutect2 and MuSE. The user may request one or more callers, and each caller produces an independently generated filtered VCF file.  

SomaticSniper, Strelka2, and MuSE require there to be **exactly one pair of input tumor/normal** BAM files, but Mutect2 will take tumor-only input (no paired normal), as well as tumor/normal BAM pairs from multiple samples from the same individual.

### Somatic SNV callers:
* [SomaticSniper](https://github.com/genome/somatic-sniper) is an older tool yielding high specificity single nucleotide somatic variants.

* [Strelka2](https://github.com/Illumina/strelka) here uses candidate indels from `Manta` and calls somatic short mutations (single nucleotide and small indel) filtered with a random forest model.

* [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) calls somatic short mutations via local assembly of haplotypes.

* [MuSE](https://github.com/wwylab/MuSE) accounts for tumor heterogeneity and calls single nucleotide somatic variants.

## How To Run
Below is a summary of how to run the pipeline.  See [here](https://confluence.mednet.ucla.edu/pages/viewpage.action?spaceKey=BOUTROSLAB&title=How+to+run+a+nextflow+pipeline) for more information on running Nextflow pipelines.

> **Note**: Because this pipeline uses an image stored in the GitHub Container Registry, you must follow the steps listed in the [Docker Introduction](https://confluence.mednet.ucla.edu/display/BOUTROSLAB/Docker+Introduction#DockerIntroduction-GitHubContainerRegistryGitHubContainerRegistry|Setup) on Confluence to set up a PAT for your GitHub account and log into the registry on the cluster before running this pipeline.
> 
1. The recommended way of running the pipeline is to directly use the source code located here: `/hot/software/pipeline/pipeline-call-sSNV/Nextflow/release/`, rather than cloning a copy of the pipeline.

    * The source code should never be modified when running our pipelines

2. Copy and edit the [input config file](config/template.config)
> Make sure the reference .fa file in config file matches the reference genome in the input BAM files.
3. Copy and edit the [input YAML](input/call-sSNV-template.yaml)
4. The pipeline can be executed locally using the command below:

```bash
nextflow run path/to/main.nf -config path/to/input.config -params-file input.yaml`
```

For example, 
* `path/to/main.nf` could be: `/hot/software/pipeline/pipeline-call-sSNV/Nextflow/release/5.0.0/main.nf`
* `path/to/input.config` is the path to where you saved your project-specific copy of [template.config](config/template.config) 
* `path/to/input.yaml` is the path to where you saved your project-specific copy of [template.yaml](input/call-sSNV-template.yaml) 

To submit to UCLAHS-CDS's Azure cloud, use the submission script [here](https://github.com/uclahs-cds/tool-submit-nf) with the command below:

```bash
python path/to/submit_nextflow_pipeline.py \
    --nextflow_script path/to/main.nf \
    --nextflow_config path/to/input.config\
    --nextflow_yaml path/to/input.yaml \
    --pipeline_run_name <run_name> \
    --partition_type F72 \
    --email jdoe@ucla.edu
```
> **Note**: Although --partition_type F2 is an available option for small data sets, Mutect2 and Muse will fail due to lack of memory.


---

## Flow Diagrams

### SomaticSniper
![Diagram](image/somatic-sniper.svg)
#### Tools
##### SomaticSniper
SomaticSniper source: https://github.com/genome/somatic-sniper
Version: SomaticSniper v1.0.5.0 (Released on Jul 16, 2015)
GitHub Package: ghcr.io/uclahs-cds/somaticsniper:1.0.5.0
##### bam-readcount
bam-readcount source: https://github.com/genome/bam-readcount
Version: v0.8.0 Release (Released on Oct 21, 2016)
GitHub Package: ghcr.io/uclahs-cds/bam-readcount:0.8.0

### Strelka2
![Diagram](image/strelka2.svg)
#### Tools
##### Manta
Manta source: https://github.com/Illumina/manta
Version: v1.6.0 (Released on Jul 9, 2019)
GitHub Package: ghcr.io/uclahs-cds/manta:1.6.0
##### Strelka2
Strelka2 source: https://github.com/Illumina/strelka
Version: v2.9.10 (Released on Nov 7, 2018)
GitHub Package: ghcr.io/uclahs-cds/strelka2:2.9.10

### Mutect 2
![alt text](docs/mutect2_chart.svg)
#### Tools
##### GATK
GATK source: https://github.com/broadinstitute/gatk
Version: 4.2.4.1 (Released on Jan 4, 2022)
Docker Image: broadinstitute/gatk:4.2.4.1

### MuSE
![alt text](docs/muse_chart.svg?raw=true)
#### Tools
##### MuSE
MuSE source: https://github.com/wwylab/MuSE
Version: 2.0 (Released on Aug 25, 2021)
GitHub Package: https://github.com/uclahs-cds/docker-MuSE/pkgs/container/muse


---

## Pipeline Steps

### SomaticSniper
#### 1. `SomaticSniper` v1.0.5.0
Compare a pair of tumor and normal bam files and output an unfiltered list of single nucleotide positions that are different between tumor and normal, in VCF format.
#### 2. Filter out ambiguous positions.
This takes several steps, listed below, and starts with the same input files given to `SomaticSniper`.
##### a. Get pileup summaries
Summarize counts of reads that support reference, alternate and other alleles for given sites.  This is done for both of the input bam files and the results are used in the next step.
##### b. Filter pileup outputs
Use `samtools.pl varFilter` to filter each pileup output (tumor and normal), then further filters each to keep only indels with QUAL > 20. `samtools.pl` is packaged with `SomaticSniper`. 
##### c. Filter SomaticSniper vcf
Use `snpfilter.pl` (packaged with `SomaticSniper`):
i. filter vcf using normal indel pileup (from step `b`).
ii. filter vcf output from step `i` using tumor indel pileup (from step `b`).
##### d. Summarize alignment information for retained variant positions
Extract positions from filtered vcf file and use with `bam-readcount` to generate a summary of read alignment metrics for each position.
##### e. Final filtering of variants using metrics summarized above
Use `fpfilter.pl` and `highconfidence.pl` (packaged with SomaticSniper), resulting in a final high confidence vcf file.

### Strelka2
#### 1. `Manta` v1.6.0
The input pair of tumor/normal bam files are used by Manta to produce candidate small indels via the `Manta` somatic configuration protocol. *Note, larger (structural) variants are also produced and can be retrieved from the intermediate files directory if save intermediate files is enabled.* 
#### 2. `Strelka2` v2.9.10
The input pair of tumor/normal bam files, along with the candidate small indel file produced by `Manta` are used by `Strelka2` to create lists of somatic single nucleotide and small indel variants, both in vcf format.  Lower quality variants that did not pass filtering are subsequently removed, yielding `somatic_snvs_pass.vcf` and `somatic_indels_pass.vcf` files.


### GATK Mutect 2

#### 1. Intervals not provided
  ##### a. Call non-canonical
  Call somatic variants in non-canonical chromosomes with `Mutect2`.
  ##### b. Split canonical
  Split the set of canonical chromosomes into x intervals for parallelization, where x is defined by the input `params.scatter_count`.
  ##### c. Call canonical
  Call somatic variant in canonical chromosomes with `Mutect2`.
  ##### d. Merge
  Merge scattered canonical and non-canonical chromosome outputs (vcfs, statistics).
  ##### e. Learn read orientations
  Create artifact prior table based on read orientations with GATK's `LearnReadOrientationModel`.
  ##### f. Filter
  Filter variants with GATK's `FilterMutectCalls`, using read orientation prior table and contamination table as well as standard filters.

#### 2. Intervals provided
  ##### a. Split
  Split the set of provided intervals into x intervals for parallelization, where x is defined by the input `params.scatter_count`. 
  ##### b. Call
  Call somatic variants for the provided intervals with `Mutect2`.
  ##### c. Merge
  Merge scattered outputs (vcfs, statistics).
  ##### d. Learn read orientations
  Create artifact prior table based on read orientations with GATK's `LearnReadOrientationModel`.
  ##### e. Filter
  Filter variants with GATK's `FilterMutectCalls`, using read orientation prior table as well as standard filters.


### MuSE
#### 1.`MuSE call`
This step carries out pre-filtering and calculating position-specific summary statistics using the Markov substitution model.
#### 2.`MuSE sump`
This step computes tier-based cutoffs from a sample-specific error model.
#### 3.Filter vcf
`MuSE` output has variants labeled as `PASS` or one of `Tier 1-5` for the lower confidence calls (`Tier 5` is lowest). This step keeps only variants labeled `PASS`.


## Inputs
To run the pipeline, one `input.yaml` and one `input.config` are needed, as follows.

### input.yaml. ([see template](input/call-sSNV-template.yaml))

| Input       | Type   | Description                               |
|-------------|--------|-------------------------------------------|
| patient_id | string | The name/ID of the patient
| tumor_BAM | string | The path to the tumor .bam file (.bai file must exist in same directory) |
| tumor_id | string | The name/ID of the tumor sample    |
| normal_BAM | string | The path to the normal .bam file (.bai file must exist in same directory) |
| normal_id | string | The name/ID of the normal sample      |
| contamination_table | path | Optional, but only for tumor samples. The path of the `contamination.table`, which is generated from the GATK's `CalculateContamination` in `pipeline-call-gSNP`. The contamination.table path can be found under `pipeline-call-gSNP`'s output `QC` folder.

* `input.yaml` should follow the standardized structure:
```
patient_id: 'patient_id'
input:
  normal:
    - id: normal_id
      BAM: /path/to/normal.bam
  tumor:
    - id: tumor_id
      BAM: /path/to/tumor.bam
      contamination_table: /path/to/contamination.table
```

* `Mutect2` can take other inputs: tumor-only sample and one patient's multiple samples. The pipeline will define `params.tumor_only_mode`, `params.multi_tumor_sample`, and `params.multi_normal_sample`. For tumor-only samples, remove the normal input in `input.yaml`, e.g. [template_tumor_only.yaml](input/example-test-tumor-only.yaml). For multiple samples, put all the input BAMs in the `input.yaml`, e.g. [template_multi_sample.yaml](input/example-test-multi-sample.yaml). Note, for these non-standard inputs, the configuration file must have 'mutect2' listed as the only algorithm. 


### input.config ([see template](config/template.config))
| Input | Required | Type   | Description                               |
|--------|---|--------|-------------------------------------------|
| algorithm   | yes | list   | List containing a combination of somaticsniper, strelka2, mutect2 and muse |
| reference   | yes | string | The reference .fa file (.fai and .dict file must exist in same directory) |
| output_dir  | yes | string | The location where outputs will be saved  |
| dataset_id | yes | string | The name/ID of the dataset    |
| exome       | yes | boolean | The option will be used by `Strelka2` and `MuSE`. When `true`, it will add the `--exome` option  to Manta and Strelka2, and `-E` option to MuSE. |
| save_intermediate_files | yes | boolean | Whether to save intermediate files |
| work_dir | no | string | The path of working directory for Nextflow, storing intermediate files and logs. The default is `/scratch` with `ucla_cds` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. |
| docker_container_registry | no | string | Registry containing tool Docker images, optional. Default: `ghcr.io/uclahs-cds` |

#### Module Specific Configuration
| Input       | Required | Type   | Description                               |
|-------------|----|--------|-------------------------------------------|
| bgzip_extra_args       | no | string | The extra option used for compressing VCFs |
| tabix_extra_args       | no | string | The extra option used for indexing VCFs |

#### Strelka2 Specific Configuration
| Input       | Required | Type   | Description                               |
|-------------|----|--------|-------------------------------------------|
| call_region | no | string | Adds '--callRegions' option when running manta and strelka2 |
* Manta and Strelka2 call the entire genome by default, however variant calling may be restricted to an arbitrary subset of the genome by providing a region file in BED format with the `--callRegions` configuration option. See the `--callRegions` documentations here: [Strelka2](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#call-regions), [Manta](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#call-regions). `--callRegions` is optional for Strelka2, but can be used to specify canonical regions to save the running time. An example of call region's bed.gz can be found and used here: `/hot/ref/tool-specific-input/Strelka2/GRCh38/strelka2_call_region.bed.gz`.

* The BED file's index file `bed.gz.tbi` needs to be stored in the same folder.
* In particular, as noted in Strelka's [User Guide](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#call-regions):
> Even when `--callRegions` is specified, the `--exome` flag is still required for exome or targeted data to get appropriate depth filtration behavior for non-WGS cases.

#### Mutect2 Specific Configuration
| Input       | Required | Type | Description                               |
|-------------|--|----|-------------------------------------------|
| split_intervals_extra_args | no | string | Additional arguments for the SplitIntervals command |
| mutect2_extra_args | no | string | Additional arguments for the Mutect2 command |
| filter_mutect_calls_extra_args | no | string | Additional arguments for the FilterMutectCalls command |
| gatk_command_mem_diff | yes | nextflow.util.MemoryUnit | How much to subtract from the task's allocated memory where the remainder is the Java heap max. (should not be changed unless task fails for memory related reasons) |
| scatter_count | yes | int | Number of intervals to split the desired interval into. Mutect2 will call each interval seperately. |
| intervals   | no | string | A GATK accepted interval list file containing intervals to search for somatic mutations. <br/> If empty or missing, will optimally partition canonical genome based on scatter_count and process non-canonical regions separately. This is the default use case. <br/> If specified and evaluates to a valid path, will pass that path to GATK to restrict the genomic regions searched. |
| germline_resource_gnomad_vcf | no | path | A stripped down version of the [gnomAD VCF](https://gnomad.broadinstitute.org/) stripped of all unneeded INFO fields, keeping only AF, currently available for GRCh38:`/hot/ref/tool-specific-input/GATK/GRCh38/af-only-gnomad.hg38.vcf.gz` and GRCh37: `/hot/ref/tool-specific-input/GATK/GRCh37/af-only-gnomad.raw.sites.vcf`. |


#### MuSE Specific Configuration
| Input       | Required | Type   | Description                               |
|-------------|----|--------|-------------------------------------------|
| dbSNP | yes | path | The path to dbSNP database's `*.vcf.gz` |

## Outputs
| Output                                         | Type         | Description                   |
|------------------------------------------------|--------------|-------------------------------|
| SomaticSniper-{version}_{sample_id}_hc.vcf.gz             | .vcf.gz         | Filterd SNV VCF (somaticsniper)|
| Strelka2-{version}_{sample_id}_somatic-snvs-pass.vcf.gz   | .vcf.gz         | Filterd SNV VCF(strelka2)     |
| Strelka2-{version}_{sample_id}_somatic-indels-pass.vcf.gz | .vcf.gz         | Filterd Indel VCF (strelka2)     |
| Mutect2-{version}_{sample_id}_filtered-pass.vcf.gz        | .vcf.gz         | Filterd SNV VCF (mutect2)      |
| MuSE-{version}_{sample_id}_filtered-pass.vcf.gz        | .vcf.gz         | Filterd SNV VCF (MuSE)   |
| report.html, timeline.html, trace.txt          | .html, .txt | Nextflow logs                 |

## Testing and Validation

Testing was performed primarily in the Boutros Lab SLURM Development cluster using F72 node. Metrics below will be updated where relevant with additional testing and tuning outputs.

### Test Data Set

| Data Set | Run Configuration | Output Dir | Control Sample | Tumor Sample |
| ------ | ------ | ------- | ------ | ------- |
| A-full-P2 |/hot/software/pipeline/pipeline-call-sSNV/Nextflow/development/unreleased/maotian-update-README/analysis/all/A-full/nextflow.config | /hot/software/pipeline/pipeline-call-sSNV/Nextflow/development/unreleased/maotian-update-README/analysis/all/A-full/output | HG002.N | P2 |

### Performance Validation
Testing was performed in the Boutros Lab SLURM Development cluster. Metrics below will be updated where relevant with additional testing and tuning outputs. Pipeline version used here is v4.0.0-rc.1

#### Mutect2
Duration: 3h 25m 24s
* Process `call_sSNVInAssembledChromosomes_Mutect2` has been split into 50 intervals, so the following table shows one of those processes:

|process_name                                 |max_duration     |max_cpu |max_peak_vmem |
|:--------------------------------------------|:----------------|:-------|:-------------|
|call_sSNVInNonAssembledChromosomes_Mutect2   | 32m 44s         | 142.0% |33.1 GB       |
|call_sSNVInAssembledChromosomes_Mutect2      |1h 20m 12s       | 123.8% |7.8 GB        |
|run_LearnReadOrientationModel_GATK           |31m 5s         |106.8%  |10.2 GB       |

#### SomaticSniper
Duration: 9h 21m 23s

|process_name                           |max_duration           |max_cpu |max_peak_vmem |
|:--------------------------------------|:----------------------|:-------|:-------------|
|convert_BAM2Pileup_SAMtools            |4h 18m 29s             | 98.2%  | 1.9 GB       |
|call_sSNV_SomaticSniper                |8h 48m 45s             |98.7%   | 511.6 MB     |
|generate_ReadCount_bam_readcount       |29m 33s          |75.9%   | 261.5 MB     |


#### Strelka2
Strelka2's runtime will be significantly improved when using `--callRegions` option to exclude the non-canoincal regions of the genome, here is the results of CPCG0196:
Sample: CPCG0196
Normal BAM: `/hot/software/pipeline/pipeline-align-DNA/Nextflow/development/outputs/bwa-mem2_and_hisat2-2.2.1/bwa-mem2/bams/a-full-CPCG0196-B1/align-DNA-20210424-024139/pipeline-alignDNA.inputs.CPCG0196-B1.bam`
Tumor BAM: `/hot/resource/pipeline_testing_set/WGS/GRCh38/A/full/CPCG0000000196-T001-P01-F.bam`

##### without `--callRegions`:

|process_name             |max_duration        |max_cpu |max_peak_vmem |
|:------------------------|:-------------------|:-------|:-------------|
|call_sIndel_Manta        |1h 24m 26s      |2724.2% |23.2 GB       |
|call_sSNV_Strelka2       |22h 32m 24s      |511.3%  |17.4 GB       |

##### with `--callRegions`:

|process_name             |max_duration        |max_cpu |max_peak_vmem |
|:------------------------|:-------------------|:-------|:-------------|
|call_sIndel_Manta        |1h 35m 25s         |1848.6% |11.7 GB        |
|call_sSNV_Strelka2       |59m 19s        |3234.0%  |8.2 GB       |

Therefore, we strongly suggest to use the `--callRegions` if the non-canonical region is unnecessary. `-callRegions`'s input `bed.gz` file can be found here: `/hot/ref/tool-specific-input/Strelka2/GRCh38/strelka2_call_region.bed.gz`. For other genome version, you can use [UCSC Liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert.

#### MuSE v2.0
MuSE v2.0 was tested with a normal/tumor paired CPCG0196 WGS sample on a F32 slurm-dev node.
Duration: 1d 11h 6m 54s

|process_name             |max_duration        |max_cpu |max_peak_vmem |
|:------------------------|:-------------------|:-------|:-------------|
|call_sSNV_MuSE        | 3h 44m 15s   | 3181.7% | 60.4 GB   |
|run_sump_MuSE         | 1d 7h 22m 2s | 100.0%  | 41.6 GB   |

---

## References
1.	Larson, D. E. et al. SomaticSniper: identification of somatic point mutations in whole genome sequencing data. Bioinformatics 28, 311–317 (2012).
2.	Kim, S. et al. Strelka2: fast and accurate calling of germline and somatic variants. Nat. Methods 15, 591–594 (2018).
3.	McKenna, A. et al. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20, 1297–1303 (2010).
4.	Fan, Y. et al. MuSE: accounting for tumor heterogeneity using a sample-specific error model improves sensitivity and specificity in mutation calling from sequencing data. Genome Biol. 17, 178 (2016).

## License

Authors: Mao Tian (maotian@mednet.ucla.edu), Bugh Caden, Helena Winata (HWinata@mednet.ucla.edu), Sorel Fitz-Gibbon (sfitzgibbon@mednet.ucla.edu).

pipeline-call-sSNV is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

This pipeline performs somatic SNV calling on a pair of normal/tumor BAMs, utilizing SomaticSniper, Strelka2, Mutect2 and MuSE.

Copyright (C) 2020-2023 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
