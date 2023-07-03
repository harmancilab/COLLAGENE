#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2)
{
	stop(sprintf("USAGE: %s [mat1 fp] [mat2 fp]", args[2]));
}

mat1_fp <- args[1];
mat2_fp <- args[2];

raw_mat1 <- read.delim(mat1_fp, header=FALSE);
raw_mat2 <- read.delim(mat2_fp, header=FALSE);

#print(raw_mat1[,-1])
#print(raw_mat2[,-1])

mat1 <- as.matrix(raw_mat1[, -1]);
mat2 <- as.matrix(raw_mat2[, -1]);

n_diff_row <- min(nrow(mat1), nrow(mat2));
n_diff_col <- min(ncol(mat1), ncol(mat2));

diff_mat <- mat1[1:n_diff_row, 1:n_diff_col]-mat2[1:n_diff_row, 1:n_diff_col];

print("Difference dimension:");
dim(diff_mat);
print(sprintf("Mean difference: %f", mean(abs(diff_mat))));
print(sprintf("Max difference: %f", max(abs(diff_mat))));

