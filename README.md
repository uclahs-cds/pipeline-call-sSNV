# pipeline-call-sSNV

## Overview
This pipeline performs somatic SNV calling given a pair of tumor/normal BAM. 3 somatic SNV callers are available and described below. Each caller will run independently of each other.
The mutect2 algorithm can also take multiple samples and tumor only samples.

## Somatic SNV callers:
* Somatic Sniper
* Strelka2
* Mutect2


### SomaticSniper
![Diagram](docs/diagram.svg)
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
![Diagram](docs/strelka2.svg)
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
To run the pipeline, one input.yaml and one template.config are needed. When running a batch of samples, template.config can be shared, while input.yaml is unique for each sample.

| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| sample_id | string | The name/ID of the sample                 | YAML File |
| algorithm   | list   | List containing a combination of somaticsniper, strelka2 or mutect2 | Config File |
| tumor       | string | The path to the tumor .bam file (.bai file must exist in same directory) | YAML File |
| normal      | string | The path to the normal .bam file (.bai file must exist in same directory) | YAML File |
| reference   | string | The reference .fa file (.fai and .dict file must exist in same directory) | Config File |
| output_dir  | string | The location where outputs will be saved  | Config File |
| output_log_dir | string | The location where log files (.command.\*) will be saved |
Config File |
| save_intermediate_files | boolean | Whether to save intermediate files | Config File |
| workdir | string | The path of working directory for Nextflow, storing intermediate files and logs. The default is `/scratch` with `ucla_cds` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. | Config File |

## Strelka2 Specific Configuration
| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| exome       | string | Adds the '--exome' option when running manta and strelka2 | Config File |
| call_region | string | Adds '--callRegions' option when running manta and strelka2 | Config File |
* Manta and Strelka2 call the entire genome by default, however variant calling may be restricted to an arbitrary subset of the genome by providing a region file in BED format with the `--callRegions` configuration option. See the `--callRegions` documentations here: [Strelka2](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#call-regions), [Manta](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#call-regions). `--callRegions` is optional for Strelka2, but can be used to specify canonical regions to save the running time. An example of call region's bed.gz can be found and used here: `/hot/ref/tool-specific-input/Strelka2/GRCh38/strelka2_call_region.bed.gz`.

* The BED file's index file `bed.gz.tbi` needs to be stored in the same folder.
* In particular, as noted in Strelka's [User Guide](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#call-regions):
> Even when `--callRegions` is specified, the `--exome` flag is still required for exome or targeted data to get appropriate depth filtration behavior for non-WGS cases.


## Mutect2 Specific Configuration
| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| split_intervals_extra_args | string | Additional arguments for the SplitIntervals command | Config File |
| mutect2_extra_args | string | Additional arguments for the Mutect2 command | Config File |
| filter_mutect_calls_extra_args | string | Additional arguments for the FilterMutectCalls command | Config File |
| gatk_command_mem_diff | nextflow.util.MemoryUnit | How much to subtract from the task's allocated memory where the remainder is the Java heap max. (should not be changed unless task fails for memory related reasons) | Config File |
| scatter_count | int | Number of intervals to split the desired interval into. Mutect2 will call each interval seperately. | Config File |
| intervals   | string | A GATK accepted interval list file containing intervals to search for somatic mutations. <br/> If empty or missing, will optimally partition canonical genome based on scatter_count and process non-canonical regions separately. This is the default use case. <br/> If specified and evaluates to a valid path, will pass that path to GATK to restrict the genomic regions searched. | Config File |
| germline_resource_gnomad_vcf | path | A copy of the gnomAD VCF only kept AF but stripped of all unnecessary INFO fields, currently available for GRCh38:`/hot/ref/tool-specific-input/GATK/GRCh38/af-only-gnomad.hg38.vcf.gz` and GRCh37: `/hot/ref/tool-specific-input/GATK/GRCh37/af-only-gnomad.raw.sites.vcf`. | Config File |

For special input, such as tumor-only sample and one patient's multiple samples, the pipeline will define `params.tumor_only_mode`, `params.multi_tumor_sample`, and `params.multi_normal_sample`. For tumor-only samples, leave the normal input in input.YAML empty, as [template_tumor_only.yaml](pipeline/inputs/template_tumor_only.yaml). For multiple samples, put all the input bams in the input.YAML, as [template_multi_sample.yaml](pipeline/inputs/template_multi_sample.yaml).

## Outputs
| Output                                         | Type         | Description                   |
|------------------------------------------------|--------------|-------------------------------|
| somaticsniper_{sample_id}_hc.vcf             | .vcf         | Final VCF file (somaticsniper)|
| strelka2_{sample_id}_somatic_snvs_pass.vcf   | .vcf         | Final VCF file (strelka2)     |
| strelka2_{sample_id}_somatic_indels_pass.vcf | .vcf         | Indel VCF file (strelka2)     |
| mutect2_{sample_id}_filtered_pass.vcf        | .vcf         | Final VCF file (mutect2)      |
| report.html, timeline.html, trace.txt          | .html & .txt | Nextflow logs                 |

#### How to run the pipeline
1. Download the stable release [here](https://github.com/uclahs-cds/pipeline-call-sSNV/releases) or the development version by cloning the GitHub repository to your machine.
2. Fill in the params section of the [config file](https://github.com/uclahs-cds/pipeline-call-sSNV/blob/main/pipeline/nextflow.config).
3. Run the pipeline using the [Nextflow submission script](https://github.com/uclahs-cds/tool-submit-nf) with the command below:
```bash
python path/to/submit_nextflow_pipeline.py \
    --nextflow_script path/to/call-sSNV.nf \
    --nextflow_config path/to/nextflow.config \
    --nextflow_yaml path/to/input.yaml \
    --pipeline_run_name <sample_id> \
    --partition_type F72 \
    --email jdoe@mednet.ucla.edu
```
<b><i>Notes:</i></b>
> The reference .fa file in config file should be the same with the reference genome that genereates the input bam files.

---

## License

Authors: Bugh Caden (YBugh@mednet.ucla.edu), Helena Winata (HWinata@mednet.ucla.edu), Mao Tian (maotian@mednet.ucla.edu).

Call-sSNV is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

This pipeline performs somatic SNV calling on a pair of normal/tumor BAMs. Mutect2, SomaticSniper, and Strelka2 are currently available in this pipeline.

Copyright (C) 2021 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

