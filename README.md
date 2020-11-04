# pipeline-call-sSNV

## Somatic SNV callers:
* Somatic Sniper


### SomaticSniper
#### Tools
##### SomaticSniper
SomaticSniper source: https://github.com/genome/somatic-sniper  
Version: SomaticSniper v1.0.5.0 (Released on Jul 16, 2015)  
Docker image: blcdsdockerregistry/call-ssnv:somaticsniper-v1.0.5.0
##### bam-readcount
bam-readcount source: https://github.com/genome/bam-readcount  
Version: v0.8.0 Release (Released on Oct 21, 2016)  
Docker image: blcdsdockerregistry/call-ssnv:bam-readcount-v0.8.0
#### Inputs
| Input       | Type   | Description                               | Location    |
|-------------|--------|-------------------------------------------|-------------|
| sample_name | string | The name/ID of the sample                 | Config File |
| tumor       | string | The path to the tumor .bam file           | Config File |
| tumor_index | string | The path to the tumor index .bam.bai file | Config File |
| normal      | string | The path to the normal .bam file          | Config File |
| reference   | string | The reference .fa file                    | Config File |
| output_dir  | string | The location where outputs will be saved  | Config File |

#### Outputs
| Output                                | Type         | Description    |
|---------------------------------------|--------------|----------------|
| somaticsniper_{sample_name}_hc.vcf    | .vcf         | Final VCF file |
| report.html, timeline.html, trace.txt | .html & .txt | Nextflow logs  |

#### How to run the pipeline
1. Fill in the params section of the config file
2. Run the pipeline using the [submission script](https://github.com/uclahs-cds/tool-submit-nf)

![Diagram](docs/diagram.svg)