#!/usr/bin/env Rscript

### This script is used to compare results of the secure federated weighted k-means protocol.
### Note that there are several hardcoded parameters.

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("USAGE: [exec] [# iterations]");
} 

library(ggplot2);

graphics.off();

N_iters <- as.numeric(args[1]);

pooled_data_df <- read.delim("DATA/per_cluster_data.txt", header=TRUE, sep = "\t");
per_cl_known_centroids <- read.delim("DATA/known_centroids.txt", header=TRUE, sep = "\t");

pooled_data_df$cluster_i <- factor(pooled_data_df$cluster_i, levels=sort(unique(pooled_data_df$cluster_i)))
pooled_data_df$site_i <- factor(pooled_data_df$site_i, levels=sort(unique(pooled_data_df$site_i)))

# Plot the data first.
pdf('original_data.pdf')
p<-ggplot(pooled_data_df, aes(X1, X2, col=cluster_i))+geom_point();
print(p);
p<-ggplot(pooled_data_df, aes(X1, X2, col=site_i))+geom_point();
print(p);
dev.off()

################################################################################################

cur_cl_centroids <- data.matrix(read.delim('DATA/Kclust_centroids.txt', header=TRUE, sep="\t"));
K_clust <- nrow(cur_cl_centroids);
n_dim <- ncol(cur_cl_centroids);
N_sites <- length(unique(pooled_data_df$site_i));

cl_iter <- 0;
print(sprintf("Iter %d:", cl_iter));
print(cur_cl_centroids);

for(cl_iter in 1:N_iters)
{
	print(sprintf("@ iteration %d", cl_iter));

	# We keep track of the total weighted coordinates for each cluster's centroid.
	all_site_update_per_cl <- matrix(0, nrow=K_clust, ncol=n_dim);

	# We also keep track of the total weight that is used for cluster's centroid.
	all_site_total_prob_per_cl <- matrix(0, nrow=K_clust, ncol=1);
	  
	cur_site_data_df <- data.matrix(pooled_data_df[, 1:n_dim]);

	# Update the means and variances for each cluster using soft k-means:
	per_subj_per_cl_dists <- matrix(nrow=nrow(cur_site_data_df), ncol=K_clust);

	################################################################
	# Calculate the distance of each sample to every K-cluster's centroid.
	for(cl_i in 1:K_clust)
	{
		# Following creates a matrix of replicated rows for the centroid so that we can subtract it.
		cur_clust_centroids_mat <- matrix(rep(cur_cl_centroids[cl_i, ], nrow(cur_site_data_df)), nrow=nrow(cur_site_data_df), byrow = TRUE);

		# Calculate the distance squared.
		cur_cl_dists <- rowSums((cur_site_data_df - cur_clust_centroids_mat)^2);

		# Set the distances to a new matrix.
		per_subj_per_cl_dists[, cl_i] <- (cur_cl_dists);
	} # cl_i loop.

	################################################################
	# For each row (subject), assign a centroid contribution score based on distance:
	################################################################
	# This is the exact exponential.
	cur_cl_probs <- exp(-1*per_subj_per_cl_dists/5);

	################################################################
	# This is the exponential approximation.
	# cur_cl_probs <- (1+(per_subj_per_cl_dists)/5);

	################################################################
	# # This is the original k-means that uses max.
	# cur_cl_probs <- matrix(0, nrow=nrow(per_subj_per_cl_dists), ncol=K_clust);
	# for(subj_i in 1:nrow(per_subj_per_cl_dists))
	# {
	#   max_k_clust_i <- which(per_subj_per_cl_dists[subj_i, ] == min(per_subj_per_cl_dists[subj_i, ]));
	#   cur_cl_probs[subj_i, ] <- 0;
	#   cur_cl_probs[subj_i, max_k_clust_i] <- 1;
	# }
	################################################################

	# ################################################################
	# # Normalize each score row-wise:
	# for(subj_i in 1:nrow(cur_cl_probs))
	# {
	#   cur_cl_probs[subj_i, ] <- cur_cl_probs[subj_i, ] / sum(cur_cl_probs[subj_i, ]);
	# }
	# ################################################################

	# Now calculate a weighted update for all subjects.
	for(cl_i in 1:K_clust)
	{
		rep_per_subj_cl_probs <- matrix(rep(cur_cl_probs[, cl_i], n_dim), nrow=nrow(cur_site_data_df));
		all_site_update_per_cl[cl_i, ] <- all_site_update_per_cl[cl_i, ] + colSums(cur_site_data_df * rep_per_subj_cl_probs);
	}

	################################################################
	# Calculate, for each cluster, the total score that was contributed to it, this is used to normalize its coordinates.
	for(cl_i in 1:K_clust)
	{
		# Add all the probabilities from this site assigned to this cluster.
		all_site_total_prob_per_cl[cl_i] <- all_site_total_prob_per_cl[cl_i] + sum(cur_cl_probs[, cl_i]);
	}
	################################################################
	  
	# Division should be by all points that contribute to this center by the cluster's total probability from all points.
	# print(all_site_update_per_cl)
	for(cl_i in 1:K_clust)
	{
		cur_cl_centroids[cl_i, ] <- (all_site_update_per_cl[cl_i, ] / all_site_total_prob_per_cl[cl_i]);
	}

	print(sprintf("Iter %d:", cl_iter));
	print(cur_cl_centroids);
}

exp_cl_centroids <- cur_cl_centroids;

# ###########################################################################################
data_w_cents <- pooled_data_df[, 1:n_dim];
data_w_cents$p_type <- 'data';
data_w_cents$alpha <- 0.5;
data_w_cents$p_size <- 1.6;

cur_cl_centroids_df <- data.frame(exp_cl_centroids);
cur_cl_centroids_df$p_type <- 'plain-centroid';
cur_cl_centroids_df$alpha <- 1.0;
cur_cl_centroids_df$p_size <- 2.2;

cur_sec_cl_centroids_df <- read.delim("secure_centroids.txt", header=FALSE, sep="\t");
cur_sec_cl_centroids_df$p_type <- 'secure-centroid';
cur_sec_cl_centroids_df$alpha <- 1.0;
cur_sec_cl_centroids_df$size <- 1.4;

colnames(cur_cl_centroids_df) <- colnames(data_w_cents);
colnames(cur_sec_cl_centroids_df) <- colnames(data_w_cents);

data_w_cents <- rbind(data_w_cents, cur_cl_centroids_df, cur_sec_cl_centroids_df);

pdf('data_plot.pdf');
p<-ggplot(data_w_cents, aes(x=X1, y=X2, col=p_type, alpha=alpha, shape=p_type))+geom_point(size=data_w_cents$p_size);
print(p);
dev.off();
# ###########################################################################################

