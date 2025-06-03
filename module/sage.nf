include { run_REDUX_SAGE as run_REDUX_SAGE_tumor; run_REDUX_SAGE as run_REDUX_SAGE_normal; call_sSNV_SAGE } from './sage-processes'
include { filter_VCF_BCFtools; split_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common' addParams(
    log_dir_prefix: "SAGE-${params.sage_version}"
    )
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'  addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: "${params.log_output_dir}/process-log/SAGE-${params.sage_version}",
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])

workflow sage {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        // Run REDUX on tumor and normal BAMs separately
        run_REDUX_SAGE_tumor(
            tumor_bam.map { [params.tumor_id, 'tumor', it] }
                .combine(tumor_index),
            params.reference,
            params.reference_index,
            params.reference_dict,
            params.redux_unmap_regions ?: "${params.work_dir}/NO_FILE.tsv",
            params.redux_ref_genome_msi_file
        )
        
        run_REDUX_SAGE_normal(
            normal_bam.map { [params.normal_id, 'normal', it] }
                .combine(normal_index),
            params.reference,
            params.reference_index,
            params.reference_dict,
            params.redux_unmap_regions ?: "${params.work_dir}/NO_FILE.tsv",
            params.redux_ref_genome_msi_file
        )
        
        // Get REDUX results for tumor and normal
        tumor_redux_results = run_REDUX_SAGE_tumor.out.redux_results
            .filter { it[1] == 'tumor' }
        normal_redux_results = run_REDUX_SAGE_normal.out.redux_results  
            .filter { it[1] == 'normal' }
        
        // Determine which BAMs to use for SAGE
        if (params.redux_jitter_msi_only) {
            // Use original BAMs but REDUX-generated jitter params and ms tables
            sage_tumor_bam = tumor_bam
            sage_normal_bam = normal_bam
        } else {
            // Use REDUX-processed BAMs
            sage_tumor_bam = tumor_redux_results.map { it[2] }
            sage_normal_bam = normal_redux_results.map { it[2] }
        }
        
        call_sSNV_SAGE(
            sage_tumor_bam,
            tumor_redux_results.map { it[3] }, // jitter_params
            tumor_redux_results.map { it[4] }, // ms_table
            sage_normal_bam,
            normal_redux_results.map { it[3] }, // jitter_params
            normal_redux_results.map { it[4] }, // ms_table
            params.reference,
            params.reference_index,
            params.reference_dict,
            params.sage_hotspots,
            params.sage_panel_bed,
            params.sage_ensembl_data_dir,
            params.sage_high_confidence_bed
        )
        
        // Filter VCF to keep only PASS variants and compress
        filter_VCF_BCFtools(call_sSNV_SAGE.out.sage_vcf
            .map{ it -> ['all', it] }
        )
        
        // Split VCF by variant type
        split_VCF_BCFtools(filter_VCF_BCFtools.out.gzvcf
            .map{ it -> it[1] },
            ['snps', 'mnps', 'indels']
        )
        
        // Rename samples in VCF to match expected format
        // Create ID mapping channel for renaming samples
        Channel
            .from([
                [orig_id: 'tumor', id: params.tumor_id],
                [orig_id: 'normal', id: params.normal_id]
            ])
            .set { id_ch }
        
        rename_samples_BCFtools(
            // combine with split_VCF_BCFtools output to duplicate the id input for each file.
            id_ch
                .collect()
                .combine(split_VCF_BCFtools.out.gzvcf)
                .map { it.take(it.size() -2) } //remove the split_VCF_BCFtools files
            ,
            split_VCF_BCFtools.out.gzvcf
        )
        
        // Compress and index the final VCF
        compress_index_VCF(
            rename_samples_BCFtools.out.gzvcf
        )
        
        // Generate checksums
        generate_sha512sum(
            compress_index_VCF.out.index_out
                .map { it -> ["sage-${it[0]}-vcf", it[1]] }
                .mix(compress_index_VCF.out.index_out
                    .map { it -> ["sage-${it[0]}-index", it[2]] })
        )

    emit:
        gzvcf = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[2]}"] }
} 