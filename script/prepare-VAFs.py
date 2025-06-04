# pylint: disable=invalid-name
""" Calculate and format VAFs for stripplot """
import os
import argparse
from typing import Tuple
from collections import defaultdict
from itertools import combinations
from statistics import mean
import vcf

ALLOWED_TOOLS = [
    'Mutect2',
    'Strelka2',
    'SomaticSniper',
    'MuSE'
]

def parse_args() -> argparse.Namespace:
    """ Parse arguments """
    parser = argparse.ArgumentParser(
        description="Format VAFs from SNV calls"
    )

    parser.add_argument(
        '--vcf-data',
        action='append',
        nargs=3,
        help='Input VCF file from a caller',
        required=True,
        metavar=('TOOL', 'PATH', 'PURITY')
    )

    parser.add_argument(
        '--sample',
        help='Sample ID being processed',
        required=True
    )

    parser.add_argument(
        '--output-dir',
        help='Output directory for writing formatted file',
        default='./'
    )

    return parser.parse_args()

def validate_tool(tool: str) -> None:
    """ Ensure given tool is valid """
    if tool not in ALLOWED_TOOLS:
        raise argparse.ArgumentTypeError(
            f"`{tool}` is not supported. Please select from {ALLOWED_TOOLS}")

def validate_purity(purity: str) -> None:
    """ Ensure purity is a number in (0, 1] """
    try:
        purity_float = float(purity)
    except ValueError as ve:
        raise argparse.ArgumentTypeError(f"Invalid purity `{purity}`. Must be a float.") from ve

    if (purity_float <= 0 or purity_float > 1):
        raise ValueError(f"`{purity_float}` must be in (0, 1].")

def validate_vcf(vcf_path: str) -> None:
    """ Validate path exists and is readable """
    if not os.path.isfile(vcf_path):
        raise argparse.ArgumentTypeError(
            f"Invalid file path: {vcf_path}, Please enter a valid path to VCF"
        )

    if not os.access(vcf_path, os.R_OK):
        raise PermissionError(f"`{vcf_path}` is not readable. Please ensure it is readable.")

def validate_output_dir(output_dir: str) -> None:
    """ Ensure output dir is writeable """
    if not os.path.isdir(output_dir):
        raise argparse.ArgumentTypeError(
            f"Invalid file path: {output_dir}, Please enter a valid path to VCF"
        )

    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"`{output_dir}` is not writeable. Please ensure it is readable.")

def validate_args(opts: argparse.Namespace) -> None:
    """ Ensure proper set of inputs is provided per VCF input """
    for (tool, vcf_path, purity) in opts.vcf_data:
        validate_tool(tool)
        validate_purity(purity)
        validate_vcf(vcf_path)

    # If a tool is repeated, the latest value will be used
    # Not erroring out in this case in case we want to support
    # multiple tools/samples in the future

    validate_output_dir(opts.output_dir)

def get_vaf_strelka2(variant: dict, sample: str) -> Tuple[int, int]:
    """ Calculate read count for Strelka2 format """
    sample_index = variant.samples.index(sample)
    ref_base = variant.REF.strip(' ') + 'U'
    ref_reads = int(variant.samples[sample_index][ref_base][0])
    # pylint: disable=R1728
    # Take alternate allele with most number of reads as the ALT
    alt_reads = max([int(variant.samples[sample_index][allele.sequence + 'U'][0]) \
        for allele in variant.ALT])
    total_reads = ref_reads + alt_reads
    if total_reads == 0:
        return 0
    return alt_reads/total_reads

def get_vaf_somaticsniper(variant: dict, sample: str) -> Tuple[int, int]:
    """ Calculate read count for SomaticSniper format """
    # using DP4
    sample_index = variant.samples.index(sample)
    ref_reads = int(variant.samples[sample_index]["DP4"][0]) \
        + int(variant.samples[sample_index]["DP4"][1])
    alt_reads = int(variant.samples[sample_index]["DP4"][2]) \
        + int(variant.samples[sample_index]["DP4"][3])
    total_reads = ref_reads + alt_reads
    if total_reads == 0:
        return 0
    return alt_reads/total_reads

def get_vaf_mutect2(variant: dict, sample: str) -> Tuple[int, int]:
    """ Calculate read count for Mutect2 format """
    sample_index = variant.samples.index(sample)
    ref_reads = int(variant.samples[sample_index]["AD"][0])
    alt_reads = int(variant.samples[sample_index]["AD"][1])
    total_reads = ref_reads + alt_reads
    if total_reads == 0:
        return 0
    return alt_reads/total_reads

def does_variant_pass(variant: dict) -> bool:
    """ Does variant pass filter checks """
    return (variant.FILTER is None or (len(variant.FILTER) == 0))

def calculate_adjusted_vaf(sample_id: str, variant: dict, purity: float) -> float:
    """ Determine function to calculate VAF and return adjusted VAF """
    sample = None
    for a_sample in variant.samples:
        if a_sample.sample == sample_id:
            sample = a_sample
            break

    if not sample:
        raise ValueError(f"{sample_id} not found.")

    sample_index = variant.samples.index(sample)
    variant_data_keys = variant.samples[sample_index].data._asdict().keys()

    raw_vaf = 0

    if 'AD' in variant_data_keys:
        raw_vaf = get_vaf_mutect2(variant, sample)
    elif 'DP4' in variant_data_keys:
        raw_vaf = get_vaf_somaticsniper(variant, sample)
    elif all(base+'U' in variant_data_keys for base in ['A', 'C', 'G', 'T']):
        raw_vaf = get_vaf_strelka2(variant, sample)
    else:
        raise ValueError(
            f'Unable to identify variant source for: {variant} with {variant_data_keys}'
        )

    return raw_vaf * purity

def parse_variants(sample: str, vcf_path: str, purity: float) -> dict:
    """ Parse variants and convert to dictionary of SNV ID and calculated VAF """
    vcf_rd = vcf.Reader(filename=vcf_path)
    variants = {}
    for variant in vcf_rd:
        if not does_variant_pass(variant):
            continue

        variants[f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"] = \
            calculate_adjusted_vaf(sample, variant, purity)

    return variants

def calculate_variant_vafs(raw_data: list, sample: str) -> dict:
    """ Take raw data and calculate the variant VAF per caller """
    variant_vafs = defaultdict(lambda: {})
    tool_vafs = {}

    for (tool, vcf_path, purity) in raw_data:
        tool_vafs = parse_variants(sample, vcf_path, float(purity))

        for variant_key, variant_vaf in tool_vafs.items():
            variant_vafs[variant_key][tool] = variant_vaf

    return variant_vafs

def get_combination_key(tools: list) -> str:
    """ Generate a key for the combination of tools """
    tool_id = '.'.join(sorted(tools))
    return f"{len(tools)}algorithm-{tool_id}"

def generate_averaged_vafs(variant_vafs: dict, all_tools: list) -> dict:
    """ Average VAFs for all combinations of tools """
    averaged_vafs = {}
    for num_tools in range(1, len(all_tools) + 1):
        all_tool_combinations = list(combinations(all_tools, num_tools))

        for tool_combination in all_tool_combinations:
            combination_id = get_combination_key(list(tool_combination))
            averaged_vafs[combination_id] = []

    for _, variant_data in variant_vafs.items():
        combination_key = get_combination_key(list(variant_data.keys()))
        averaged_vafs[combination_key].append(mean(variant_data.values()))

    for combination_key, vafs in averaged_vafs.items():
        averaged_vafs[combination_key] = mean(vafs) if vafs else 0

    return averaged_vafs

def write_vafs(data: dict, sample: str, outfile: str) -> None:
    """ Write data to output file """
    with open(outfile, 'wt', encoding='utf-8') as wr:
        wr.write("sample\tcombination\tadjVAF\n")
        for combination, adjvaf in data.items():
            wr.write(f"{sample}\t{combination}\t{adjvaf}\n")

def main() -> None:
    """ Main entrypoint function """
    opts = parse_args()
    validate_args(opts)

    variant_vafs = calculate_variant_vafs(opts.vcf_data, opts.sample)

    all_tools = list({x[0] for x in opts.vcf_data})
    all_tools.sort()

    melted_data = generate_averaged_vafs(variant_vafs, all_tools)

    write_vafs(melted_data, opts.sample, os.path.join(opts.output_dir, 'adjusted_vafs.tsv'))

if __name__=='__main__':
    main()