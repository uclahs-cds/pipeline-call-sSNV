### PREAMBLE #######################################################################################
library(BoutrosLab.plotting.general);
library(argparse);
options(stringsAsFactors = FALSE);

### PREPROCESSING ##################################################################################
parser <- ArgumentParser();
parser$add_argument(
    '-i',
    '--input',
    type = 'character',
    required = TRUE,
    help = 'Melted data ready for plotting'
    );

parser$add_argument(
    '-o',
    '--output-dir',
    type = 'character',
    default = './',
    help = 'Plot output directory'
    );
args <- parser$parse_args();

## Edit stripplot data -----------------------------------------------------------------------------
stripplot.data <- read.table(
    args$input,
    sep = '\t',
    header = TRUE
    );

# Indicate whether the combination has One/Two/Three/Four sSNV callers
snvcaller.combination <- stripplot.data$combination;
stripplot.data$num_of_callers <- 'Four';
stripplot.data$num_of_callers[grep('1algorithm',snvcaller.combination)] <- 'One';
stripplot.data$num_of_callers[grep('2algorithm',snvcaller.combination)] <- 'Two';
stripplot.data$num_of_callers[grep('3algorithm',snvcaller.combination)] <- 'Three';

# Remove groups in which no sSNVs were called
stripplot.data <- stripplot.data[stripplot.data$adjVAF != 0, ];

# Convert `num_of_callers` to factors and sort
stripplot.data$num_of_callers <- factor(
    stripplot.data$num_of_callers,
    levels = c(
        'One',
        'Two',
        'Three',
        'Four'
        )
    );

stripplot.data <- stripplot.data[order(stripplot.data$num_of_callers),];

### MAIN ###########################################################################################
# Prepare a legend showing sample size
total.samples <- length(unique(stripplot.data$sample));
total.samples.label <- as.character(paste0('n = ', total.samples));
total.samples.label <- as.expression(substitute(bold(x), list(x = total.samples.label)));
legend.total.samples <- legend.grob(
    list(
        legend = list(
            labels = total.samples.label,
            size = 0
            )
        ),
    label.cex = 1.5
    );

output.filename <- 'stripplot.png';

total.samples <- unique(stripplot.data$sample);
sample.pchs <- sapply(
    total.samples,
    FUN = function(x) which(total.samples == x)[1]
    );

stripplot.data$pch <- sapply(
    stripplot.data['sample'],
    FUN = function(x) which(total.samples == x)[1]
    );

stripplot.data$color <- sapply(
    stripplot.data['sample'],
    FUN = function(x) 'black'
    );

stripplot.data;

create.stripplot(
    formula = adjVAF ~ num_of_callers,
    data = stripplot.data,
    filename = paste0(args$output_dir,output.filename),
    jitter.data = TRUE,
    jitter.factor = 0.6,
    xlab.label = 'Number of sSNV Callers',
    ylab.label = 'Adjusted VAF',
    yaxis.rot = 0,
    xaxis.tck = 0,
    yaxis.tck = 0.5,
    ylab.axis.padding = 2.5,
    pch = stripplot.data$pch,
    cex = 1,
    col = stripplot.data$color,
    key = list(
        space = 'right',
        text = list(
            lab = total.samples,
            cex = 1.25,
            col = 'black'
            ),
        points = list(
            pch = sample.pchs,
            col = 'black',
            cex = 1.25
            ),
        padding.text = 3
        ),
    legend = list(
        inside = list(
            fun = legend.total.samples,
            x = 0.99,
            y = 0.99,
            corner = c(1,1)
            )
        ),
    top.padding = 2,
    left.padding = 2,
    right.padding = 2,
    height = 8,
    width = 10,
    resolution = 200
    );
