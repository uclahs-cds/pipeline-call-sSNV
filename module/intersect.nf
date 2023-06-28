include { generate_sha512sum } from './common'
include { intersect_VCFs_BCFtools } from './intersect-processes.nf'

workflow intersect {
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
    }
