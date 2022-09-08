# pipeline-call-sSNV

## Overview
This pipeline performs somatic SNV calling given a pair of tumor/normal BAM. 3 somatic SNV callers are available and described below. Each caller will run independently of each other.
The mutect2 algorithm can also take multiple samples and tumor only samples.

## Somatic SNV callers:
* Somatic Sniper
* Strelka2
* Mutect2


### SomaticSniper
![Diagram](image/diagram.svg)
#### Tools
##### SomaticSniper
SomaticSniper source: https://github.com/genome/somatic-sniper
Version: SomaticSniper v1.0.5.0 (Released on Jul 16, 2015)
Docker image: blcdsdockerregistry/call-ssnv:somaticsniper-v1.0.5.0
##### bam-readcount
bam-readcount source: https://github.com/genome/bam-readcount
Version: v0.8.0 Release (Released on Oct 21, 2016)
Docker image: blcdsdockerregistry/call-ssnv:bam-readcount-v0.8.0

### Strelka2
![Diagram](image/strelka2.svg)
#### Tools
##### Manta
Manta source: https://github.com/Illumina/manta
Version: v1.6.0 (Released on Jul 9, 2019)
Docker image: blcdsdockerregistry/call-ssnv:manta-v1.6.0
##### Strelka2
Strelka2 source: https://github.com/Illumina/strelka
Version: v2.9.10 (Released on Nov 7, 2018)
Docker image: blcdsdockerregistry/call-ssnv:strelka2-v2.9.10

### Mutect 2
#### Tools
##### GATK
GATK source: https://github.com/broadinstitute/gatk
Version: 4.2.4.1 (Released on Jan 4, 2022)
Docker image: broadinstitute/gatk:4.2.4.1

## Inputs
To run the pipeline, one `input.yaml` and one `template.config` are needed. When running a batch of samples, `template.config` can be shared, while `input` is unique for each sample.

### Input YAML

| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| patient_id | string | The name/ID of the patient    | YAML File |
| tumor_BAM | string | The path to the tumor .bam file (.bai file must exist in same directory) | YAML File |
| tumor_id | string | The name/ID of the tumor sample    | YAML File |
| normal_BAM | string | The path to the normal .bam file (.bai file must exist in same directory) | YAML File |
| normal_id | string | The name/ID of the normal sample      | YAML File |

* `input.YAML` should follow the standardized structure:
```
patient_id: 'patient_id'
input:
  normal:
    - id: normal_id
      BAM: /path/to/normal.bam
  tumor:
    - id: tumor_id
      BAM: /path/to/tumor.bam
```
* A template of `input.YAML` can be found [here](./input/call-sSNV-template.yaml).

### Input Config
| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| dataset_id | string | The name/ID of the dataset    | Config File |
| algorithm   | list   | List containing a combination of somaticsniper, strelka2 or mutect2 | Config File |
| reference   | string | The reference .fa file (.fai and .dict file must exist in same directory) | Config File |
| output_dir  | string | The location where outputs will be saved  | Config File |
| output_log_dir | string | The location where log files (.command.\*) will be saved |
Config File |
| save_intermediate_files | boolean | Whether to save intermediate files | Config File |
| work_dir | string | The path of working directory for Nextflow, storing intermediate files and logs. The default is `/scratch` with `ucla_cds` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. | Config File |

#### Strelka2 Specific Configuration
| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| exome       | string | Adds the '--exome' option when running manta and strelka2 | Config File |
| call_region | string | Adds '--callRegions' option when running manta and strelka2 | Config File |
* Manta and Strelka2 call the entire genome by default, however variant calling may be restricted to an arbitrary subset of the genome by providing a region file in BED format with the `--callRegions` configuration option. See the `--callRegions` documentations here: [Strelka2](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#call-regions), [Manta](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#call-regions). `--callRegions` is optional for Strelka2, but can be used to specify canonical regions to save the running time. An example of call region's bed.gz can be found and used here: `/hot/ref/tool-specific-input/Strelka2/GRCh38/strelka2_call_region.bed.gz`.

* The BED file's index file `bed.gz.tbi` needs to be stored in the same folder.
* In particular, as noted in Strelka's [User Guide](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#call-regions):
> Even when `--callRegions` is specified, the `--exome` flag is still required for exome or targeted data to get appropriate depth filtration behavior for non-WGS cases.


#### Mutect2 Specific Configuration
| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| split_intervals_extra_args | string | Additional arguments for the SplitIntervals command | Config File |
| mutect2_extra_args | string | Additional arguments for the Mutect2 command | Config File |
| filter_mutect_calls_extra_args | string | Additional arguments for the FilterMutectCalls command | Config File |
| gatk_command_mem_diff | nextflow.util.MemoryUnit | How much to subtract from the task's allocated memory where the remainder is the Java heap max. (should not be changed unless task fails for memory related reasons) | Config File |
| scatter_count | int | Number of intervals to split the desired interval into. Mutect2 will call each interval seperately. | Config File |
| intervals   | string | A GATK accepted interval list file containing intervals to search for somatic mutations. <br/> If empty or missing, will optimally partition canonical genome based on scatter_count and process non-canonical regions separately. This is the default use case. <br/> If specified and evaluates to a valid path, will pass that path to GATK to restrict the genomic regions searched. | Config File |
| germline_resource_gnomad_vcf | path | A copy of the gnomAD VCF only kept AF but stripped of all unnecessary INFO fields, currently available for GRCh38:`/hot/ref/tool-specific-input/GATK/GRCh38/af-only-gnomad.hg38.vcf.gz` and GRCh37: `/hot/ref/tool-specific-input/GATK/GRCh37/af-only-gnomad.raw.sites.vcf`. | Config File |

For special input, such as tumor-only sample and one patient's multiple samples, the pipeline will define `params.tumor_only_mode`, `params.multi_tumor_sample`, and `params.multi_normal_sample`. For tumor-only samples, leave the normal input in `input.YAML` empty, as [template_tumor_only.yaml](input/example-test-tumor-only.yaml). For multiple samples, put all the input bams in the `input.YAML`, as [template_multi_sample.yaml](input/example-test-multi-sample.yaml).

## Outputs
| Output                                         | Type         | Description                   |
|------------------------------------------------|--------------|-------------------------------|
| {SomaticSniper-version}_{sample_id}_hc.vcf             | .vcf         | Final VCF file (somaticsniper)|
| {Strelka2-version}_{sample_id}_somatic_snvs_pass.vcf   | .vcf         | Final VCF file (strelka2)     |
| {Strelka2-version}_{sample_id}_somatic_indels_pass.vcf | .vcf         | Indel VCF file (strelka2)     |
| {Mutect2-version}_{sample_id}_filtered_pass.vcf        | .vcf         | Final VCF file (mutect2)      |
| report.html, timeline.html, trace.txt          | .html & .txt | Nextflow logs                 |

#### How to run the pipeline
1. Using the [stable release](https://github.com/uclahs-cds/pipeline-call-sSNV/releases) stored under `/hot/software/pipeline/pipeline-call-sSNV/Nextflow/release/` or the development version by cloning the GitHub repository to your machine.
2. Fill in the params section of the [config file](config/template.config) and [input YAML](input/call-sSNV-template.yaml)
3. Run the pipeline using the [Nextflow submission script](https://github.com/uclahs-cds/tool-submit-nf) with the command below:
```bash
python path/to/submit_nextflow_pipeline.py \
    --nextflow_script path/to/call-sSNV.nf \
    --nextflow_config path/to/nextflow\
    --nextflow_yaml path/to/input.yaml \
    --pipeline_run_name <sample_id> \
    --partition_type F72 \
    --email jdoe@mednet.ucla.edu
```
<b><i>Notes:</i></b>
> The reference .fa file in config file should be the same with the reference genome that genereates the input bam files.

---


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
* Process `call_sSNVInAssembledChromosomes_Mutect2` has been splited into 50 intervals, so the following table shows one of those processes:

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


---

## License

Authors: Mao Tian (maotian@mednet.ucla.edu), Bugh Caden, Helena Winata (HWinata@mednet.ucla.edu).

Call-sSNV is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

This pipeline performs somatic SNV calling on a pair of normal/tumor BAMs. Mutect2, SomaticSniper, and Strelka2 are currently available in this pipeline.

Copyright (C) 2020-2022 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

