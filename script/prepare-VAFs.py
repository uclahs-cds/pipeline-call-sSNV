""" Calculate and format VAFs for stripplot """
import os
import argparse

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
        '--vcf',
        action='append',
        nargs=3,
        help='Input VCF file from a caller',
        required=True,
        metavar=('TOOL', 'PATH', 'PURITY')
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

def validate_vcf(vcf: str) -> None:
    """ Validate path exists and is readable """
    if not os.path.isfile(vcf):
        raise argparse.ArgumentTypeError(
            f"Invalid file path: {vcf}, Please enter a valid path to VCF"
        )

    if not os.access(vcf, os.R_OK):
        raise PermissionError(f"`{vcf}` is not readable. Please ensure it is readable.")

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
    for (tool, vcf_path, purity) in opts.vcf:
        validate_tool(tool)
        validate_purity(purity)
        validate_vcf(vcf_path)

    validate_output_dir(opts.output_dir)

def main() -> None:
    opts = parse_args()
    print(opts)
    validate_args(opts)


if __name__=='__main__':
    main()
