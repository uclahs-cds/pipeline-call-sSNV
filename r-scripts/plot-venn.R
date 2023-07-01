# Script to plot the Venn diagram of the intersection of the VCF files
# Initial commit: Sorel Fitz-Gibbon 2023-06-29
# Input:
#  -i, --isec_dir: The directory containing the output from BCFtools intersect
#  -d, --dataset: The dataset ID passed from nextflow
# Output:
# - A Venn diagram of the intersection of the VCF files

## 1. Setup the environment ########################################################################
library('argparse');
library('BoutrosLab.utilities');
library('VennDiagram');

# tmp: for testing
#args <- list();
#args$path <- '/hot/software/pipeline/pipeline-call-sSNV/Nextflow/development/unreleased/sfitz-intersect-vcfs/call-sSNV-6.0.0/CPCG0000000196-T001-P01'
#args$valid_algorithm <- 'NULL';
#args$dataset.id <- 'CPCG';
#args$output.dir <- paste(getwd(), paste(gsub(' |:', '-', Sys.time()), 'VCF-Comparison-Result', sep = '_'), sep = '/');

## 2. Parse the arguments ##########################################################################
parser <- ArgumentParser();
parser$add_argument('-i', '--isec_dir', help = 'The directory containing the output from BCFtools intersect', type = 'character');
parser$add_argument('-d', '--dataset', help = 'The dataset ID passed from nextflow', type = 'character');
args <- parser$parse_args();

### Main ##########################################################
algorithms <- readLines(paste0(args$isec_dir,'/README.txt'));
algorithms <- algorithms[grep(paste0('^', args$isec_dir), algorithms)];
algorithms <- gsub(paste0(args$isec_dir,'.*\t'), '', algorithms);
algorithms <- gsub('-.*', '', algorithms);
sites <- read.table(paste0(args$isec_dir,'/sites.txt'), header = FALSE, colClasses = 'character');
split.col <- strsplit(as.character(sites$V5), '');
sites$col1 <- sapply(split.col, '[', 1);
sites$col2 <- sapply(split.col, '[', 2);
sites$col3 <- sapply(split.col, '[', 3);
sites$col4 <- sapply(split.col, '[', 4);
sites$V5 <- NULL;
header <- c('chrom', 'pos', 'ref', 'alt', algorithms);
colnames(sites) <- header
variants <- paste(sites$chrom, sites$pos, sep = '_');
tool.variants <- lapply(sites[, algorithms], function(x) variants[x == 1])
tool.variants.ordered <- tool.variants[order(lengths(tool.variants), decreasing = TRUE)]

VennDiagram::venn.diagram(
    tool.variants.ordered,
    filename = generate.filename(args$dataset, 'Venn-diagram', 'tiff' ),
    fill = c('orange', 'red', 'green', 'blue'),
    lty = 'dashed',
    cex = 1,
    cat.cex = 0.8,
    cat.col = 'black'
    );
