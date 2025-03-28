log.info """\
====================================
         I N T E R S E C T
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_r_VennDiagram: ${params.docker_image_r_VennDiagram}
Intersect Options:
- ncbi_build:             ${params.ncbi_build}
- vcf2maf_extra_args:  ${params.vcf2maf_extra_args}
====================================
"""
process reorder_samples_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf.gz",
        enabled: params.save_intermediate_files
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}-${algorithm}" }

    input:
    tuple val(algorithm), path(gzvcf)
    path indices
    val tumor_id
    val normal_id

    output:
    path "*-reorder.vcf.gz", emit: gzvcf

    script:
    """
    set -euo pipefail
    infile=\$(basename ${gzvcf} .vcf.gz)
    outfile="\${infile}-reorder.vcf.gz"
    bcftools view -s ${tumor_id},${normal_id} --output \${outfile} ${gzvcf}
    """
    }

process intersect_VCFs_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.vcf.gz*"
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "isec-2-or-more/*.txt",
        saveAs: { "${file(it).getParent().getName()}/${params.output_filename}_${file(it).getName()}" }
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "isec-1-or-more/*.txt",
        saveAs: { "${file(it).getParent().getName()}/${params.output_filename}_${file(it).getName()}" }
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}" }

    input:
    path gzvcf
    path indices
    path intersect_regions
    path intersect_regions_index

    output:
    path "*.vcf.gz", emit: gzvcf
    path "*.vcf.gz.tbi", emit: idx
    path "isec-2-or-more/*.txt"
    path "isec-1-or-more/*.txt", emit: isec

    script:
    vcf_list = gzvcf.join(' ')
    regions_command = params.use_intersect_regions ? "--regions-file ${intersect_regions}" : ""
    """
    set -euo pipefail
    # intersect keeping only variants that are present in at least 2 VCFs
    # Use README.txt to rename output files to include sample names
    bcftools isec --nfiles +2 \
        --output-type z \
        --prefix isec-2-or-more \
        ${regions_command} \
        ${vcf_list}
    awk '/Using the following file names:/{x=1;next} x' isec-2-or-more/README.txt  \
        | sed 's/-reorder.vcf.gz\$/-intersect.vcf.gz/' \
        | while read a b c d; do
            mv \$a \$d
            mv \$a.tbi \$d.tbi
            done
    # intersect, keeping all variants, to create presence/absence list of variants in each VCF
    bcftools isec --nfiles +1\
        --output-type z \
        --prefix isec-1-or-more \
        ${regions_command} \
        ${vcf_list}
    """
    }

process plot_VennDiagram_R {
    container params.docker_image_r_VennDiagram
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.tiff"
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}" }

    input:
    path script_dir
    path isec

    output:
    path "*.tiff"

    script:
    """
    set -euo pipefail
    Rscript ${script_dir}/plot-venn.R --isec_readme README.txt --isec_sites sites.txt --outfile ${params.output_filename}_Venn-diagram.tiff
    """
    }

process concat_VCFs_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*concat.vcf",
        enabled: params.save_intermediate_files
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}" }

    input:
    path vcfs
    path indices

    output:
    path "*concat.vcf", emit: vcf

    script:
    vcf_list = vcfs.join(' ')
    """
    set -euo pipefail
    # BCFtools concat to create a single VCF with all nfiles +2 variants
    # output header is a uniquified concatenation of all headers
    # output `INFO` `FORMAT` `NORMAL` and `TUMOR` fields are from the first listed VCF that has the variant
    bcftools concat \
        --output-type v \
        --output ${params.output_filename}_SNV-concat.vcf \
        --allow-overlaps \
        --rm-dups all \
        ${vcf_list}
    """
    }

process convert_VCF_vcf2maf {
    container params.docker_image_vcf2maf
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.maf",
        enabled: params.save_intermediate_files
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}" }

    input:
    path vcf
    path reference
    val normal_id
    val tumor_id

    output:
    path "*.maf", emit: maf

    script:
    """
    set -euo pipefail
    perl /opt/vcf2maf.pl --inhibit-vep \
        --filter-vcf 0 \
        --ncbi-build ${params.ncbi_build} \
        --input-vcf ${vcf} \
        --normal-id ${normal_id} \
        --tumor-id ${tumor_id} \
        --output-maf ${params.output_filename}_SNV-concat.maf \
        --ref-fasta ${reference} \
        ${params.vcf2maf_extra_args}
    """
    }
