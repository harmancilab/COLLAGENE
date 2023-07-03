#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2)
{
	stop(sprintf("USAGE: %s [2 column S/T stats] [Output file]", args[2]));
}


ST_stats_file <- args[1];
op_file <- args[2];

print(sprintf("Calculating 1-DOF chi-square p-values using S/T statistics in %s and saving to %s", ST_stats_file, op_file));

#ST_stats_file <- "../Secure_Client_KeyMaker_Keys/PER_CLIENT_LOCAL_DIRS/SITE_1/POOLED_PVALS.txt"
all_ST_stats <- read.delim(ST_stats_file, header=FALSE);
ST_stats <- all_ST_stats[, 1] / all_ST_stats[, 2];
pvals <- pchisq(ST_stats, 1, lower.tail=FALSE);
pvals_w_stats <- cbind(all_ST_stats[, 1], all_ST_stats[, 2], pvals);
write.table(pvals_w_stats, sep = "\t", file = op_file, row.names = FALSE, col.names=FALSE);


