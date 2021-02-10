include { strelka2_somatic; manta } from './strelka2-processes'

workflow strelka2 {
    main:
        manta(
            params.tumor,
            params.tumor_index,
            params.normal,
            params.normal_index,
            params.reference,
            params.reference_index
        )
        strelka2_somatic(
            params.tumor,
            params.tumor_index,
            params.normal,
            params.normal_index,
            params.reference,
            params.reference_index,
            manta.out[0],
            manta.out[1]
        )
    emit:
        strelka2_somatic.out
}