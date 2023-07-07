# Script to plot the Venn diagram of the intersection of the VCF files
# Initial commit: Sorel Fitz-Gibbon 2023-06-29
# Input:
#  -i, --isec_dir: The directory containing the output from BCFtools intersect
#  -d, --dataset: The dataset ID passed from nextflow
#  -r, --regions: call regions (for single chromosome QC plotting)
# Output:
# - A Venn diagram of the intersection of the VCF files

## Setup the environment ###########################################################################
library('argparse');
library('BoutrosLab.utilities');
library('VennDiagram');

# tmp: for testing
#args <- list();
#args$isec_dir <- 'isec-1-or-more';
#args$dataset<- 'CPCG';
#args$regions_file <- '/hot/ref/tool-specific-input/Strelka2/GRCh38/strelka2_call_region.bed.gz';

## Parse the arguments #############################################################################
parser <- ArgumentParser();
parser$add_argument('-i', '--isec_dir', help = 'The directory containing the output from BCFtools intersect', type = 'character');
parser$add_argument('-d', '--dataset', help = 'The dataset ID passed from nextflow', type = 'character');
parser$add_argument('-r', '--regions', help = 'call regions (for single chromosome QC plotting)', type = 'character');
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
plot.venn(tool.variants.ordered, generate.filename(args$dataset, 'Venn-diagram', 'tiff'));

## Single chromosome QC ############################################################################
regions <- read.table(args$regions, header = FALSE, colClasses = c('character', 'NULL', 'NULL'));
# keep tools ordered by number of variants in full set
tool.order <- names(tool.variants.ordered);
for chr in unique(regions$V1) {
    chr.sites <- sites[sites$chrom == chr,];
    chr.variants <- paste(chr.sites$chrom, chr.sites$pos, sep = '_');
    chr.tool.variants <- lapply(chr.sites[, tool.order], function(x) chr.variants[x == 1])
    #... heatmap
    }
