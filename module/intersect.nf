include { generate_sha512sum } from './common'
include { intersect_VCFs_BCFtools; plot_venn_R } from './intersect-processes.nf'

workflow intersect {
    // pass bin directory in project folder as channel into docker
    script_dir_ch = Channel.fromPath("$projectDir/r-scripts", checkIfExists: true)

    take:
    tool_vcfs
    tool_indices

    main:
        intersect_VCFs_BCFtools(
            tool_vcfs,
            tool_indices,
            params.call_region,
            params.call_region_index
            )
        file_for_sha512 = intersect_VCFs_BCFtools.out.consensus_vcf
            .flatten()
            .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-vcf", it]}
            .mix(intersect_VCFs_BCFtools.out.consensus_idx
                .flatten()
                .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-idx", it]}
                )
        generate_sha512sum(file_for_sha512)
        plot_venn_R(
            script_dir_ch,
            intersect_VCFs_BCFtools.out.isec_dir,
            params.call_region
        )
    }
