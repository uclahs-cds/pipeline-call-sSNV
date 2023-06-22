# This is the script to compare the vcf results ################
# Initial commit: Mao Tian 09-28-2022
## 1. Setup the environment ########################
library('argparse');
library('BoutrosLab.plotting.general');
library('BoutrosLab.utilities');
library('VennDiagram');

# tmp: for testing
#args <- list();
##args$path <- 'NULL';
##args$file <- '/hot/software/pipeline/pipeline-call-sSNV/Nextflow/development/unreleased/sfitz-intersect-vcfs/maotian06-compare-vcf/input-vcfs.txt';
#args$path <- '/hot/software/pipeline/pipeline-call-sSNV/Nextflow/development/unreleased/sfitz-intersect-vcfs/call-sSNV-6.0.0/CPCG0000000196-T001-P01'
#args$file <- 'NULL';
#args$valid_algorithm <- 'NULL';
#args$dataset.id <- 'CPCG';
#args$output.dir <- paste(getwd(), paste(gsub(' |:', '-', Sys.time()), 'VCF-Comparison-Result', sep = '_'), sep = '/');
#args$pairwise_comparsion <- 'FALSE';
#args$global_comparsion <- 'TRUE';

### 2. Dataset Argument Parser ########################################
parse.args <- function() {
    data.parser <- ArgumentParser();
    data.parser$add_argument('-p', '--path', help = 'directory that stores the vcf.gz files', default = 'NULL');
    data.parser$add_argument('-f', '--file', help = 'an input file that lists all the vcf.gz files', default = 'NULL');
    data.parser$add_argument('-v', '--valid_algorithm', help = 'valid algorithms, default are the four call-sSNV algorithms', default = 'NULL');
    data.parser$add_argument('-d', '--dataset.id', help = 'dataset id will be used to generate the output directory that will store both intermediate and output files', default = 'test');
    data.parser$add_argument('-o', '--output.dir', help = 'the output directory, default is current working directory', default = 'NULL');
    data.parser$add_argument('--pairwise_comparsion', help = 'the option to do pairwise comparsion, deflaut is FALSE.', default = 'FALSE');
    data.parser$add_argument('--global_comparsion', help = 'the output directory, default is current working directory, deflaut is TRUE', default = 'TRUE');
    ### Set the arg_parase
    args <- data.parser$parse_args();
    if ( args$output.dir == 'NULL') {
        args$output.dir <- paste(getwd(), paste(gsub(' |:', '-', Sys.time()), 'VCF-Comparison-Result', sep = '_'), sep = '/');
    }
    return(args);
    };
### 3. Read Input VCF Lists ##################################
read.inputs <- function(args) {
    if (args$path != 'NULL') {
        vcf.list <- list.files(path = args$path,
                                    pattern = '*snvs.vcf.gz$',
                                    all.files = TRUE,
                                    full.names = TRUE,
                                    recursive = TRUE);
        } else if (args$file != 'NULL') {
            vcf.list <- as.list(t(read.table(args$file, header = FALSE)));
        } else {
            print( 'No input file or input directory!');
        };
    # set the default valid algorithm, the default are four algorithms used in call-sSNV
    if (args$valid_algorithm == 'NULL') {
        algorithm.list <- c('MuSE', 'Mutect2', 'SomaticSniper', 'Strelka2');
#        args$valid_algorithm <- algorithm.list;
        } else {
            algorithm.list <- args$valid_algorithm
        };
    # this function is used to extract the algorithm of the vcf.gz file
    # need to pre-define the algorithm.list
    extract.algorithm <- function(x) {
        for (i in 1:length(algorithm.list)) {
            if (grepl(algorithm.list[i], x)) {
                return(algorithm.list[i]);
            };
        };
    };
    # Strelka2 has an output with indel variants, removed it here
#    vcf.list <- vcf.list[ - grep('indel', vcf.list)];
    vcf.list <- vcf.list[ - grep('intermediate', vcf.list)];
    vcf.list.name <- unlist(lapply(vcf.list, extract.algorithm));
    names(vcf.list) <- vcf.list.name;
    vcf.count <- c();
    # count each vcf file variant
    for (file in vcf.list) {
        count <- system(paste('zcat', file, '| grep -v "^#" | wc -l', sep = ' '), intern = TRUE);
        vcf.count <- c(vcf.count, as.numeric(count));
    }
    vcf.summary <- data.frame(vcf.list, vcf.list.name, vcf.count);
    colnames(vcf.summary) <- c('filename', 'algorithm', 'variant.count');
    # rank vcf by variant.count
    vcf.summary <- vcf.summary[order(vcf.summary$variant.count, decreasing = TRUE),];
    vcf.list <- as.character(vcf.summary$filename);
    names(vcf.list) <- as.character(vcf.summary$algorithm);
    system(paste('mkdir', args$output.dir));
    system(paste('mkdir', paste(args$output.dir, 'intermediate', sep = '/')));
    setwd(args$output.dir);
    write.table(vcf.summary, 'vcf_summary.tsv', sep = '\t', row.names = FALSE, quote = FALSE);
    return(vcf.list);
    };

### 4. Run VCF compare with system() ##################################
run.vcf.compare <- function(vcf.list, args) {
    docker.command <- 'docker run --rm -v /hot/:/hot/ -u `id -u` biocontainers/vcftools:v0.1.16-1-deb_cv1 vcf-compare';
    if (args$global_comparsion == 'TRUE') {
        vcf.list.text <- paste(vcf.list, collapse = ' ');
        vcf.compare.command <- paste(docker.command,  vcf.list.text, '> global_comparison.txt', sep = ' ');
        futile.logger::flog.info('Running the following command to check VCF overlap:');
        print(vcf.compare.command);
        system(vcf.compare.command, intern = TRUE, ignore.stderr = TRUE);

    };
    if (args$pairwise_comparsion == 'TRUE') {
        for (i in 1:(length(vcf.list) - 1)) {
            for (j in (i + 1) : length(vcf.list)) {
                file1 <- vcf.list[i];
                file2 <- vcf.list[j];
                name <- paste(names(vcf.list)[i], names(vcf.list)[j], 'pairwise', sep = '-' );
                system(paste('bash vcf_check.sh', file1, file2, name, args$output.dir, sep = ' '), intern = TRUE, ignore.stderr = TRUE);
            }
        };
    };
    system('mv *.txt intermediate/');
};

### 5. Gather VCF compare results ##################################
gather.vcf.compare.result <- function( args ){
    setwd(args$output.dir);

    if (args$global_comparsion == 'TRUE') {
        futile.logger::flog.info('Reading the following file:');
        file <- list.files(path = getwd(),
                                    pattern = 'global_comparison.txt$',
                                    full.names = TRUE,
                                    recursive = TRUE);
        system(paste('cat', file, '| grep "^VN" | cut -f 2- >temp_overlap_count.txt', sep = ' '));
        # start read temp_overlap_count.txt by each line
        con <- file('temp_overlap_count.txt', "r");
        overlap.count <- c();
        category <- c();
        vcf.list.name <- names(vcf.list);
        while ( TRUE ) {
            line = readLines(con, n = 1);
            if ( length(line) == 0 ) {
            break
            }
            line.content <- unlist(strsplit(line, '\t'));
            overlap.count <- c(overlap.count, line.content[1]);
            line.content <- gsub('\\s+[(]+[0-9]+.[0-9]%[)]', '', line.content);
            if (length(line.content) == 2) {
                filename <- line.content[2];
                algorithm <- vcf.list.name[match(filename, vcf.list)];
            } else {
                algorithm <- c()
                for (i in (2 : length(line.content))){
                    filename <- line.content[i];
                    algorithm <- c(algorithm, vcf.list.name[match(filename, vcf.list)]);
                }
            }
            category <- c(category, paste(algorithm, collapse = '-'));
        }
        close(con)
        overlap.count <- data.frame(overlap.count, category);
        colnames(overlap.count) <- c('count', 'category');
        overlap.count$algorithm.count <- (1 + lengths(regmatches(overlap.count$category, gregexpr("-", overlap.count$category))));
        overlap.count <- overlap.count[order(overlap.count$algorithm.count),];
        vcf.list.name.level <- as.numeric(factor(vcf.list.name));
        # generate compare group, such as 1, 2, 3, 12, 13, 23
        for (i in 1:nrow(overlap.count)) {
            if ( overlap.count$algorithm.count[i] == 1) {
                overlap.count$compare.group[i] <- match(overlap.count$category[i], vcf.list.name);
            } else {
                category.list <- unlist(strsplit(overlap.count$category[i], '-'));
                category.list.number <- as.character(sapply(category.list, function(x) {
                    return(match(x, vcf.list.name))
                    }));
                category.list.number <- category.list.number[order(category.list.number)];
                overlap.count$compare.group[i] <- paste(category.list.number, collapse = '');
                rm(category.list.number);
            }
        };
        overlap.count <- overlap.count[order(overlap.count$compare.group),];
        futile.logger::flog.info('Here is the result will be stored in overlap_count.tsv');

        print(overlap.count);
        write.table(overlap.count,
                    'overlap_count.tsv',
                    sep = '\t',
                    row.names = FALSE,
                    quote = FALSE);
        return(overlap.count);
        }
    if (args$pairwise_comparsion == 'TRUE') {
        file.list <- list.files(path = getwd(),
                                pattern = 'intermediate/*_vcf_comparsion.txt$',
                                full.names = TRUE);
        for (i in 1:length(file.list)) {
            futile.logger::flog.info('Reading the following file:');
            system(paste('cat', file.list[i], '| grep "^VN" | cut -f 2- >temp_overlap_count.txt', sep = ' '));
            overlap.df <- read.table('temp_overlap_count.txt', header = FALSE, allowEscapes = FALSE);
            system(paste('cat', file.list[i], '| grep "^VN" | cut -f 3 | head -2 >filename.txt', sep = ' '));
            filename.df <- read.table('filename.txt', header = FALSE, sep = ' ');
            algo1 <- unlist(strsplit(basename(filename.df$V1[1]), '-'))[1];
            algo2 <- unlist(strsplit(basename(filename.df$V1[2]), '-'))[1];
            system(paste('cat', file.list[i], '| grep "^VN" | cut -f 2 | head -3 >temp_number.txt', sep = ' '));
            number.df <- read.table('temp_number.txt', header = FALSE);
            number.df$V1 <- as.numeric(number.df$V1);
            cross.area <- number.df$V1[3];
            area1 <- number.df$V1[1] + cross.area;
            area2 <- number.df$V1[2] + cross.area;
            venn.plot <- draw.pairwise.venn(
                            area1, area2, cross.area,
                            category = c(algo1, algo2),
                            cat.pos = c(0, 180),
                            euler.d = TRUE,
                            sep.dist = 0.03,
                            rotation.degree = 45
                        );
            tiff(filename = paste(algo1, algo2, 'Pairwise-Venn-diagram.tiff', sep = '_'),
                    compression = "lzw");
            grid.draw(venn.plot);
            dev.off();
            };
        };
    }
### 5. Plot Global Venn Plot ######################################
plot.quad.venn.plot <- function (overlap.count, vcf.list) {
    device <- pdf(file = NULL);
    futile.logger::flog.info('Reading vcf_summary.tsv');
    vcf.summary <- read.table('vcf_summary.tsv', sep = '\t', header = TRUE);
    if (all(vcf.summary$filename == vcf.list)) {
        futile.logger::flog.info('algorithm sequence matched, run variant count check');
    } else {
        futile.logger::flog.info('algorithm sequence not matched, reordering');
        vcf.summary <- vcf.summary[match(vcf.list, vcf.summary$filename),];
        if (all(vcf.summary$algorithm == vcf.list)) {
        futile.logger::flog.info('reordered, run variant count check');
        }
    }
    # check if the sum of the overlap.count equals to the summary of vcf
    futile.logger::flog.info('check if the sum of the overlap.count equals to the summary of vcf');
    overlap.count$count <- as.numeric(overlap.count$count);
    for (algo in vcf.summary$algorithm) {
        count.from.summary <- vcf.summary$variant.count[vcf.summary$algorithm == algo];
        count.from.overlap.sum <- sum(
            overlap.count$count[grep(algo, overlap.count$category)]
            );
        if (count.from.summary == count.from.overlap.sum) {
            futile.logger::flog.info(paste(algo, 'Match!', sep = ' '));
        } else {
            futile.logger::flog.info(paste(algo, 'Not Match! The sum of the vcf_compare tool does not equal to the summary of vcf', sep = ' '));
        }
    }
#    u1 = vcf.summary$variant.count[1];
#    u2 = vcf.summary$variant.count[2];
#    u3 = vcf.summary$variant.count[3];
#    u4 = vcf.summary$variant.count[4];
#
#    combination2 <- combn(1:4, 2);
#    combination3 <- combn(1:4, 3);
#    combination4 <- combn(1:4, 4);
#    comb.list <- c(
#        apply(combination2, 2, paste0, collapse=""),
#        apply(combination3, 2, paste0, collapse=""),
#        apply(combination4, 2, paste0, collapse="")
#        );

    count <- function(x) {
        cnt <- overlap.count$count[match(x, overlap.count$compare.group)];
        if (is.na(cnt)) {
            return(0);
        } else {
            return(cnt);
        }
    };
#    u12 = overlap.count$count[match('12', overlap.count$compare.group)];
#    u13 = overlap.count$count[match('13', overlap.count$compare.group)];
#    u14 = overlap.count$count[match('14', overlap.count$compare.group)];
#    u23 = overlap.count$count[match('23', overlap.count$compare.group)];
#    u24 = overlap.count$count[match('24', overlap.count$compare.group)];
#    u34 = overlap.count$count[match('34', overlap.count$compare.group)];
#    u123 = overlap.count$count[match('123', overlap.count$compare.group)];
#    u124 = overlap.count$count[match('124', overlap.count$compare.group)];
#    u134 = overlap.count$count[match('134', overlap.count$compare.group)];
#    u234 = overlap.count$count[match('234', overlap.count$compare.group)];
#    u1234 = overlap.count$count[match('1234', overlap.count$compare.group)];
#    area.vector <- c(u3, u34, u4, u13, u134, u1234, u234, u24, u1, u14,  u124, u123, u23, u2, u12);
#    print(area.vector);
    # plot qual.venn.plot
    venn.plot <- draw.quad.venn(
        area1 = vcf.summary$variant.count[1],
        area2 = vcf.summary$variant.count[2],
        area3 = vcf.summary$variant.count[3],
        area4 = vcf.summary$variant.count[4],
        n12 = (count('12') + count('123') + count('124') + count('1234')),
        n13 = (count('13') + count('123') + count('134') + count('1234')),
        n14 = (count('14') + count('124') + count('134') + count('1234')),
        n23 = (count('23') + count('123') + count('234') + count('1234')),
        n24 = (count('24') + count('124') + count('234') + count('1234')),
        n34 = (count('34') + count('134') + count('234') + count('1234')),
        n123 = (count('123') + count('1234')),
        n124 = (count('124') + count('1234')),
        n134 = (count('134') + count('1234')),
        n234 = (count('234') + count('1234')),
        n1234 = count('1234'),
        category = vcf.summary$algorithm,
        fill = c("orange", "red", "green", "blue"),
        lty = "dashed",
        cex = 1,
        cat.cex = 0.8,
        cat.col = 'black');
    tiff(
        filename = generate.filename(args$dataset.id, 'Quad-Venn-diagram', '.tiff' ),
        compression = "lzw" );
    grid.rect(x = unit(1, "npc"), y = unit(1, "npc"),
          width = unit(5, "npc"), height = unit(5, "npc"),
          just = "centre", default.units = "npc", name = NULL,
          gp=gpar(), draw = TRUE, vp = NULL);
    grid.draw(venn.plot);
    dev.off();
}



### Main ##########################################################

args <- parse.args();
vcf.list <- read.inputs(args);
run.vcf.compare(vcf.list, args);
if (args$global_comparsion == 'TRUE') {
    overlap.count <- gather.vcf.compare.result(args);
} else {
    gather.vcf.compare.result(args);
};
plot.quad.venn.plot(overlap.count, vcf.list);
