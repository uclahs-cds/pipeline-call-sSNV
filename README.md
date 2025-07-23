# pipeline-call-sSNV

[![GitHub release](https://img.shields.io/github/v/release/uclahs-cds/pipeline-call-sSNV)](https://github.com/uclahs-cds/pipeline-call-sSNV/actions/workflows/prepare-release.yaml)

- [call-sSNV](#pipeline-call-ssnv)
  - [Overview](#overview)
  - [How To Run](#how-to-run)
  - [Flow Diagrams](#flow-diagrams---variant-calling)
    - [Variant Calling](#flow-diagrams---variant-calling)
    - [Variant Intersection](#flow-diagrams---variant-intersection)
  - [Pipeline Steps](#pipeline-steps---variant-calling)
    - [Variant Calling](#pipeline-steps---variant-calling)
       - [SomaticSniper](#somaticsniper-1)
       - [Strelka2](#strelka2-1)
       - [Mutect2](#gatk-mutect2)
       - [MuSE](#muse-1)
    - [Variant Intersection](#pipeline-steps---variant-intersection)
       - [BCFtools and VennDiagram](#pipeline-steps---variant-intersection)
       - [vcf2maf](#vcf2maf)
       - [VAF Plotting](#vaf-plotting)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Performance Validation and Resource Requirements](#performance-validation)
  - [References](#references)
  - [Discussions](https://github.com/uclahs-cds/pipeline-call-sSNV/discussions)
  - [Contributors](https://github.com/uclahs-cds/template-NextflowPipeline/graphs/contributors)
  - [License](#license)

## Overview
The call-sSNV nextflow pipeline performs somatic SNV calling given a pair of tumor/normal BAM files. Four somatic SNV callers are available: SomaticSniper, Strelka2, Mutect2 and MuSE. The user may request one or more callers, and each caller produces an independently generated filtered VCF file.

If two or more callers are requested, additional output includes both a VCF and an MAF file with the set of SNVs shared by two or more callers, and a Venn Diagram showing counts of shared and private SNVs.

SomaticSniper, Strelka2, and MuSE require there to be **exactly one pair of input tumor/normal** BAM files, but Mutect2 will take tumor-only input (no paired normal), as well as tumor/normal BAM pairs for multiple samples from the same individual.

### Somatic SNV callers:
* [SomaticSniper](https://github.com/genome/somatic-sniper) is an older tool yielding high specificity single nucleotide somatic variants.

* [Strelka2](https://github.com/Illumina/strelka) uses candidate indels from `Manta` and calls somatic short mutations (single nucleotide and small indel), filtered with a random forest model.

* [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) calls somatic short mutations via local assembly of haplotypes.

* [MuSE](https://github.com/wwylab/MuSE) accounts for tumor heterogeneity and calls single nucleotide somatic variants.

## How To Run
Below is a summary of how to run the pipeline.  See [here](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3197004/How+to+run+a+nextflow+pipeline) for more information on running Nextflow pipelines.
>
1. Copy and edit the [input config file](config/template.config)
> Make sure the reference .fa file in config file matches the reference genome in the input BAM files.
1. Copy and edit the [input YAML](input/call-sSNV-template.yaml)
1. The pipeline can be executed locally using the command below:

```bash
nextflow run path/to/main.nf -config path/to/input.config -params-file input.yaml
```

---

## Flow Diagrams - Variant Calling
### [SomaticSniper](docs/flowcharts.md#somaticsniper)
### [Strelka2](docs/flowcharts.md#strelka2)
### [Mutect2](docs/flowcharts.md#mutect2)
### [MuSE](docs/flowcharts.md#muse)
## Flow Diagrams - Variant Intersection
### [BCFtools intersection and VennDiagram visualization using R](docs/flowcharts.md#intersect)

---

## Pipeline Steps - Variant Calling

### SomaticSniper
#### 1. `SomaticSniper` v1.0.5.0
Compare a pair of tumor and normal BAM files and output an unfiltered list of single nucleotide positions that are different between tumor and normal, in VCF format.
#### 2. Filter out ambiguous positions.
This takes several steps, listed below, and starts with the same input files given to `SomaticSniper`. These are used to generate a list of high confidence indels to assist SNV filtering.
##### a. Get indel pileup summaries
Summarize counts of reads that support reference, alternate and other alleles for given sites.  This is done for both of the input BAM files and the results are used in the next step.
##### b. Filter indel pileup outputs
Use `samtools.pl varFilter` to filter each pileup output (tumor and normal), then further filter each to keep only indels with QUAL > 20. `samtools.pl` is packaged with `SomaticSniper`.
##### c. Filter SomaticSniper VCF
Use `snpfilter.pl` (packaged with `SomaticSniper`):
i. filter VCF using normal indel pileup (from step `b`).
ii. filter VCF output from step `i` using tumor indel pileup (from step `b`).
##### d. Summarize alignment information for retained SNV positions
Extract positions from filtered VCF file and use with `bam-readcount` to generate a summary of read alignment metrics for each position.
##### e. Final filtering of SNVs using metrics summarized above
Use `fpfilter.pl` and `highconfidence.pl` (packaged with SomaticSniper), resulting in a final high confidence VCF file.

### Strelka2
#### 1. `Manta` v1.6.0
The input pair of tumor/normal BAM files are used by Manta to produce candidate small indels via the `Manta` somatic configuration protocol. *Note, larger (structural) variants are also produced and can be retrieved from the intermediate files directory if save intermediate files is enabled.*
#### 2. `Strelka2` v2.9.10
The input pair of tumor/normal BAM files, along with the candidate small indel file produced by `Manta` are used by `Strelka2` to create lists of somatic single nucleotide and small indel variants, both in VCF format.  Lower quality variants that did not pass filtering are subsequently removed, yielding `.SNV-pass.vcf.gz` and `.Indel-pass.vcf.gz` files.

### GATK Mutect2

#### 1. Define intervals for scattering
The `params.intersect_regions` of the reference genome are split into x intervals for parallelization, where x is defined by `params.scatter_count`.
#### 2. Call small somatic variants
Call somatic variants with `Mutect2`.
#### 3. Merge
Merge scattered outputs (VCFs, statistics).
#### 4. Learn read orientations
Create artifact prior table based on read orientations with GATK's `LearnReadOrientationModel`.
#### 5. Filter
Filter variants with GATK's `FilterMutectCalls`, using read orientation prior table and contamination table as well as standard filters.
#### 6. Split VCF
Split filtered VCF into separate files for each variant type: SNVs, MNVs and INDELs.

### MuSE
#### 1. `MuSE call`
Pre-filtering and calculating position-specific summary statistics using the Markov substitution model.
#### 2. `MuSE sump`
Computes tier-based cutoffs from a sample-specific error model.
#### 3.Filter VCF
`MuSE` output has SNVs labeled as `PASS` or one of `Tier 1-5` for the lower confidence calls (`Tier 5` is lowest). This step keeps only SNVs labeled `PASS`.

## Pipeline Steps - Variant Intersection
If two or more algorithms were selected the Intersect workflow will run. Currently the resulting VCF and MAF files include any SNVs found by two or more algorithms.
### BCFtools isec -n +1; VennDiagram
Determines presence/absence of each SNV within each algorithm's set of filtered SNVs. Results are listed in the output files: `isec-1-or-more/README.txt` and `isec-1-or-more/sites.txt`, and are summarized in a Venn Diagram plot (TIFF format).
### BCFtools isec -n +2
Determines presence/absence of SNVs found in two or more of each algorithm's set of filtered SNVs, and outputs a `consensus` VCF for each algorithm containing SNVs found by that algorithm plus at least one other algorithm. Results are also listed in the output files: `isec-2-or-more/README.txt` and `isec-2-or-more/sites.txt`.
### BCFtools concat
Concatenates the 2+ algorithm `consensus` SNVs into one VCF (SNV-concat.vcf.gz).  The output header is a uniquified concatenation of all input VCF headers.  The output fields `INFO`, `FORMAT`, `NORMAL` and `TUMOR` are from the first listed VCF that has the SNV. Input VCFs are sorted alphanumerically by the algorithm name.
### vcf2maf
Converts SNV-concat.vcf.gz from step 3 into [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/).  Output includes allele counts and flanking basepairs, but most fields are blank.  Details can be found [here](https://github.com/uclahs-cds/pipeline-call-sSNV/discussions/222#discussion-5512332).

### VAF Plotting
A stripplot is generated to display the distribution of allele frequencies for sets of SNVs categorized by the number of callers that include the variant.

## Inputs
To run the pipeline, one `input.yaml` and one `input.config` are needed, as follows.

### input.yaml. ([see template](input/call-sSNV-template.yaml))

| Input       | Type   | Description                               |
|-------------|--------|-------------------------------------------|
| patient_id | string | The name/ID of the patient
| tumor_BAM | path | The path to the tumor .bam file (.bai file must exist in same directory) |
| normal_BAM | path | The path to the normal .bam file (.bai file must exist in same directory) |
| contamination_table | path | Optional, but only for tumor samples. The path of the `contamination.table`, which is generated from the GATK's `CalculateContamination` in `pipeline-call-gSNP`. The contamination.table path can be found under `pipeline-call-gSNP`'s output `QC` folder

* `input.yaml` should follow the standardized structure:
```
patient_id: 'patient_id'
input:
  normal:
    - BAM: /path/to/normal.bam
  tumor:
    - BAM: /path/to/tumor.bam
      contamination_table: /path/to/contamination.table
```

* `Mutect2` can take other inputs: tumor-only sample and one patient's multiple samples. For tumor-only samples, remove the normal input in `input.yaml`, e.g. [template_tumor_only.yaml](input/example-test-tumor-only.yaml). For multiple samples, put all the input BAMs in the `input.yaml`, e.g. [template_multi_sample.yaml](input/example-test-multi-sample.yaml). Note, for these non-standard inputs, the configuration file must have 'mutect2' listed as the only algorithm.


### input.config ([see template](config/template.config))
| Input | Required | Type   | Description                               |
|--------|---|--------|-------------------------------------------|
| `algorithm`   | yes | list   | List containing a combination of somaticsniper, strelka2, mutect2 and muse |
| `reference`   | yes | string | The reference .fa file (.fai and .dict file must exist in same directory) |
| `intersect_regions`* | yes | string | A bed file listing the genomic regions for variant calling. Excluding `decoy` regions is HIGHLY recommended *
| `output_dir`  | yes | string | The location where outputs will be saved  |
| `dataset_id` | yes | string | The name/ID of the dataset    |
| `exome`       | yes | boolean | The option will be used by `Strelka2` and `MuSE`. When `true`, it will add the `--exome` option  to Manta and Strelka2, and `-E` option to MuSE |
| `save_intermediate_files` | yes | boolean | Whether to save intermediate files |
| `work_dir` | no | string | The path of working directory for Nextflow, storing intermediate files and logs. The path to a temporary working directory for Nextflow, storing intermediate files and logs. It is recommended to use fast, local storage with high I/O performance. |
| `docker_container_registry` | no | string | Registry containing tool Docker images, optional. Default: `ghcr.io/uclahs-cds` |
| `base_resource_update` | optional | namespace | Namespace of parameters to update base resource allocations in the pipeline. Usage and structure are detailed in `template.config` and below. |

 *Providing `intersect_regions` is required and will limit the final output to just those regions.  All regions of the reference genome could be provided as a `bed` file with all contigs, however it is HIGHLY recommended to remove `decoy` contigs from the human reference genome. Including these thousands of small contigs will require the user to increase available memory for `Mutect2` and will cause a very long runtime for `Strelka2`. See [Discussion here](https://github.com/uclahs-cds/pipeline-call-sSNV/discussions/216).

 ### Base resource allocation updaters
To optionally update the base resource (cpus or memory) allocations for processes, use the following structure and add the necessary parts to the [input.config](config/template.config) file. The default allocations can be found in `config/resources.json`. If available resources have matched cpus and memory within `90% - 1GB` of one of the pre-specified configurations, that configuration will be used.  Otherwise the default configuration will be used. A spreadsheet view of the resource configuration as of Dec 2024 is [here](https://github.com/uclahs-cds/pipeline-call-sSNV/discussions/328).  For very large or challanging input samples, we suggest using the `m64` configuration or similar.


```Nextflow
base_resource_update {
    memory = [
        [['process_name', 'process_name2'], <multiplier for resource>],
        [['process_name3', 'process_name4'], <different multiplier for resource>]
    ]
    cpus = [
        [['process_name', 'process_name2'], <multiplier for resource>],
        [['process_name3', 'process_name4'], <different multiplier for resource>]
    ]
}
```
> **Note** Resource updates will be applied in the order they're provided so if a process is included twice in the memory list, it will be updated twice in the order it's given.

Examples:

- To double memory of all processes:
```Nextflow
base_resource_update {
    memory = [
        [[], 2]
    ]
}
```
- To double memory for `call_sSNV_Mutect2` and triple memory for `run_validate_PipeVal` and `run_sump_MuSE`:
```Nextflow
base_resource_update {
    memory = [
        ['call_sSNV_Mutect2', 2],
        [['run_validate_PipeVal', 'run_sump_MuSE'], 3]
    ]
}
```
- To double CPUs and memory for `run_sump_MuSE` and double memory for `run_validate_PipeVal`:
```Nextflow
base_resource_update {
    cpus = [
        ['run_sump_MuSE', 2]
    ]
    memory = [
        [['run_sump_MuSE', 'run_validate_PipeVal'], 2]
    ]
}
```

#### Module Specific Configuration
| Input       | Required | Type   | Description                               |
|-------------|----|--------|-------------------------------------------|
| bgzip_extra_args       | no | string | The extra option used for compressing VCFs |
| tabix_extra_args       | no | string | The extra option used for indexing VCFs |

#### Mutect2 Specific Configuration
| Input       | Required | Type | Description                               |
|-------------|--|----|-------------------------------------------|
| split_intervals_extra_args | no | string | Additional arguments for the SplitIntervals command |
| mutect2_extra_args | no | string | Additional arguments for the Mutect2 command |
| filter_mutect_calls_extra_args | no | string | Additional arguments for the FilterMutectCalls command |
| gatk_command_mem_diff | yes | nextflow.util.MemoryUnit | How much to subtract from the task's allocated memory where the remainder is the Java heap max. (should not be changed unless task fails for memory related reasons) |
| scatter_count | yes | int | Number of intervals to split the desired interval into. Mutect2 will call each interval seperately. |
| germline_resource_gnomad_vcf | no | path | A stripped down version of the [gnomAD VCF](https://gnomad.broadinstitute.org/) stripped of all unneeded INFO fields, keeping only AF. |
| panel_of_normals_vcf | no | path | VCF file of sites observed in normal. This could be useful for tumor only mode. |

#### MuSE Specific Configuration
| Input       | Required | Type   | Description                               |
|-------------|----|--------|-------------------------------------------|
| dbSNP | yes | path | The path to [NCBI's dbSNP database](https://www.ncbi.nlm.nih.gov/snp/) of known SNPs in VCF format, e.g. `GCF_000001405.40.gz` |

#### Variant Intersection Specific Configuration
| Input       | Required | Type   | Description                               |
|-------------|----|--------|-------------------------------------------|
| ncbi_build | yes | string | vcf2maf requires the reference genome build ID, e.g. GRCh38 |
| vcf2maf_extra_args | no | string | additional arguments for the vcf2maf command|

## Outputs
| Tool Outputs                                         | Type         | Description                   |
|------------------------------------------------|--------------|-------------------------------|
| SomaticSniper-{version}_{sample_id}_SNV.vcf.gz             | .vcf.gz         | Filtered SNV VCF (somaticsniper)|
| Strelka2-{version}_{sample_id}_SNV.vcf.gz   | .vcf.gz         | Filtered SNV VCF(strelka2)     |
| Strelka2-{version}_{sample_id}_Indel.vcf.gz | .vcf.gz         | Filtered Indel VCF (strelka2)     |
| Mutect2-{version}_{sample_id}_SNV.vcf.gz        | .vcf.gz         | Filtered SNV VCF (mutect2)      |
| Mutect2-{version}_{sample_id}_Indel.vcf.gz        | .vcf.gz         | Filtered Indel VCF (mutect2)      |
| Mutect2-{version}_{sample_id}_MNV.vcf.gz        | .vcf.gz         | Filtered MNV VCF (mutect2)      |
| Mutect2-{version}_{sample_id}_filteringStats.tsv        | .tsv         | FilterMutectCalls output (mutect2 QC)      |
| MuSE-{version}_{sample_id}_SNV.vcf.gz        | .vcf.gz         | Filtered SNV VCF (MuSE)   |
| report.html, timeline.html, trace.txt          | .html, .txt | Nextflow logs                 |

| Intersect Outputs                                         | Type         | Description                   |
|------------------------------------------------|--------------|-------------------------------|
| isec-1-or-more | directory | BCFtools isec README.txt and sites.txt, all variants |
| isec-2-or-more | directory | BCFtools isec README.txt and sites.txt, variants shared by 2 or more tools |
| SomaticSniper-{version}_{sample_id}_consensus-variants.vcf.gz             | .vcf.gz         | `2-or-more` SNV VCF|
| Strelka2-{version}_{sample_id}_consensus-variants.vcf.gz   | .vcf.gz         | `2-or-more` SNV VCF     |
| Mutect2-{version}_{sample_id}_consensus-variants.vcf.gz        | .vcf.gz         | `2-or-more` SNV VCF      |
| MuSE-{version}_{sample_id}_consensus-variants.vcf.gz        | .vcf.gz         | `2-or-more` SNV VCF   |
| BCFtools-{version}_{sample_id}_Venn-diagram.tiff | .tiff | Venn Diagram with intersection counts for all variants (`1-or-more`) |
| BCFtools-{version}_{sample_id}_SNV-concat.vcf.gz | .vcf.gz | Single SNV VCF with all `2-or-more` variants and mixed annotation |
| BCFtools-{version}_{sample_id}_SNV-concat.maf.gz | .maf.gz | Single SNV MAF with all `2-or-more` variants and mixed annotation |
| BPG-{version}_{dataset_id}_{sample_id}_adjVAF.png | .png | Stripplot of adjusted VAFs with combinations of callers |

### Performance Validation
Testing was performed in the Boutros Lab SLURM Development cluster. Metrics below will be updated where relevant with additional testing and tuning outputs. Pipeline version used here is v4.0.0-rc.1

#### Whole Exomes
General estimates, with wide variation, are that whole exome sequences require 16 cpus and 32 GB of memory to run all of the pipeline algorithms.  If MuSE is excluded 8 cpus and 16 GB of memory are sufficient.  Run time for a test pair of exome tumor/normal input BAMs of 4 GB/5 GB was in both cases 1 to 2 hours.


#### Whole Genomes
General estimates, with wide variation, are that whole genome sequences require 72 cpus and 144 GB of memory to run all of the pipeline algorithms. If MuSE is excluded 8 cpus and 16 GB of memory are sufficient, but run time could be very long. Run time for a test pair of WGS tumor/normal input BAMs of 400 GB/200 GB was 15 hours for 72 cpus/144 GB, and 52 hours for 8 cpus/16 GB excluding MuSE.

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
Strelka2's runtime will be significantly improved when using `--callRegions` option to exclude the non-canonincal regions of the genome. Here are results from a typical BAM pair

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

Therefore, we strongly suggest to use the `--callRegions` if the non-canonical region is unnecessary.

#### MuSE v2.0
MuSE v2.0 was tested with a normal/tumor paired CPCG0196 WGS sample on an exclusive node with 32 cpus and 64 GB memory.
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

Copyright (C) 2020-2024 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
