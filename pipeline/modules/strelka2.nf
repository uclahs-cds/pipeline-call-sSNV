include { strelka2_somatic; manta } from './strelka2-processes'

workflow strelka2 {
    main:
        manta(
            params.tumor,
            "${params.tumor_index}.bai",
            params.normal,
            "${params.normal_index}.bai",
            params.reference,
            "${params.reference_index}.fai"
        )
        strelka2_somatic(
            params.tumor,
            "${params.tumor_index}.bai",
            params.normal,
            "${params.normal_index}.bai",
            params.reference,
            "${params.normal_index}.bai",
            manta.out[0],
            manta.out[1]
        )
    emit:
        strelka2_somatic.out
}