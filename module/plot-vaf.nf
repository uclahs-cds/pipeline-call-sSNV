include { calculate_adjVAF_Python; plot_adjVAF_R } from './plot-vaf-processes'

workflow plot_vaf {
    take:
    identified_gzvcfs
    all_files

    main:
    calculate_adjVAF_Python(
        identified_gzvcfs,
        all_files
        )

    plot_adjVAF_R(
        calculate_adjVAF_Python.out.adjusted_vafs
        )
}
