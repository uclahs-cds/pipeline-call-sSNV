# Script to plot a Venn diagram of shared variants from the different SNV calling algorithms, using the output of BCFtools isec
# Initial commit: Sorel Fitz-Gibbon 2023-06-29
# Input:
#  -i, --isec_dir: The directory containing the output from BCFtools intersect
#  -d, --dataset: The dataset ID passed from nextflow
# Output:
# - A Venn diagram of shared variant counts from the BCFtools intersection of the VCF files

## Setup the environment ###########################################################################
library('argparse');
library('VennDiagram');

## Parse the arguments #############################################################################
parser <- ArgumentParser();
parser$add_argument('-i', '--isec_dir', help = 'The directory containing the output from BCFtools intersect', type = 'character');
parser$add_argument('-o', '--outfile', help = 'Output filename', type = 'character');
args <- parser$parse_args();

## Function: plot venn diagram #####################################################################
plot.venn <- function(tool.variants, outfile) {
    VennDiagram::venn.diagram(
        tool.variants.ordered,
        filename = outfile,
        fill = c('orange', 'red', 'green', 'blue'),
        lty = 'dashed',
        cex = 1,
        cat.cex = 0.8,
        cat.col = 'black'
        );
    }

### Main ###########################################################################################
# Get intersection counts from BCFtools isec output and format for plotting
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
tool.variants.ordered <- tool.variants[order(lengths(tool.variants), decreasing = TRUE)];
plot.venn(tool.variants.ordered, args$outfile);
