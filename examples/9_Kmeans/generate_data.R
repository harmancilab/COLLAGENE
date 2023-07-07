#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("USAGE: [exec] [Subject-2-cluster mapping file] [K-cluster centroids to generate] [Data dimension] [Centroid mean sigma] [Centroid sigma sigma]");
} 

library('pracma');
library(MASS)

subj_2_cluster_matrix_file <- args[1];
K_clust <- as.numeric(args[2]);
n_dim <- as.numeric(args[3]);
centroid_mean_sigma <- as.numeric(args[4]);
centroid_sigma_sigma <- as.numeric(args[5]);

# Read the number of subjects per cluster.
per_site_per_cl_n_subj <- data.matrix(read.delim(subj_2_cluster_matrix_file, header=FALSE, sep="\t"));

print(per_site_per_cl_n_subj);

N_sites <- nrow(per_site_per_cl_n_subj);

# This has to be hardcoded.
if(N_sites != 3)
{
	stop("MUST HAVE 3 SITES");
}

n_clusters <- ncol(per_site_per_cl_n_subj);

print(sprintf("Simulating mapping file: %s ;; n_clusters: %d ;; N_sites: %d ;; K_clust: %d ;; n_dim: %d ;; sigma(mean): %f ;; sigma(sigma): %f", 
		subj_2_cluster_matrix_file, 
		n_clusters,
		N_sites,
		K_clust, 
		n_dim, 
		centroid_mean_sigma, 
		centroid_sigma_sigma));

WRITE_NEW_DATA <- 1;

pooled_data_df <- NULL;
per_cl_known_centroids <- matrix(0, n_clusters,n_dim);

if(WRITE_NEW_DATA == 1)
{
  print("Writing new data..");
  pooled_data_df <- data.frame();
  
  for(cl_i in 1:n_clusters)
  {
    cur_cl_cent <- rnorm(n_dim) * centroid_mean_sigma;
    cur_cl_sigma <- abs(rnorm(n_dim) * centroid_sigma_sigma);
    
    # Store the current cluster's centroid.
    per_cl_known_centroids[cl_i, ] <- cur_cl_cent;
    
    cur_cl_sigma_mat <- zeros(n_dim, n_dim);
    
    # Set the sigma to a diagonal matrix.
    cur_cl_sigma_mat[cbind(1:n_dim, 1:n_dim)] = cur_cl_sigma;
    
    # Get the total # of subjects on this cluster.
    cur_cl_n_total_subj <- colSums(per_site_per_cl_n_subj)[cl_i];
    
    print(sprintf("Sampling %d subjects for cluster-%d", cur_cl_n_total_subj, cl_i));
    
    # Sample current cluster's subjects.
    cur_cl_data <- mvrnorm(cur_cl_n_total_subj, cur_cl_cent, cur_cl_sigma_mat);
    
    # Allocate data frame.
    cur_cl_data_df <- data.frame(cur_cl_data);
    cur_cl_data_df$cluster_i <- cl_i;
    cur_cl_data_df$site_i <- -1;
    cumsum_n_subj <- cumsum(per_site_per_cl_n_subj[, cl_i])
    for(site_i in 1:N_sites)
    {
      if(per_site_per_cl_n_subj[site_i, cl_i] > 0)
      {
        cur_ind <- which(cur_cl_data_df$site_i[1:cumsum_n_subj[site_i]] == -1);
        print(sprintf("cl_i: %d ; site_i: %d: %d subjects.", 
                      cl_i, site_i, length(cur_ind)));
        cur_cl_data_df$site_i[cur_ind] = site_i;
      }
    }
    
    pooled_data_df <- rbind(pooled_data_df, cur_cl_data_df);
  } # cluster index loop.
  
  pooled_data_df$cluster_i <- factor(pooled_data_df$cluster_i, levels=sort(unique(pooled_data_df$cluster_i)));

  # Save data.
  write.table(file="DATA/per_cluster_data.txt", x=pooled_data_df, sep = "\t", quote = FALSE, row.names = FALSE);
  
  write.table(file="DATA/known_centroids.txt", x=per_cl_known_centroids, sep = "\t", quote = FALSE, row.names = FALSE);
  
  # Write the data for each site.
  for(site_i in 1:N_sites)
  {
    cur_site_data_matrix <- pooled_data_df[pooled_data_df$site_i == site_i, 1:n_dim];
    write.table(file=sprintf("DATA/site_%d_data.txt", site_i), x=cur_site_data_matrix, sep = "\t", quote = FALSE, row.names = FALSE);
    
    cur_site_feats_matrix <- pooled_data_df[pooled_data_df$site_i == site_i, (n_dim+1):ncol(pooled_data_df)];
    write.table(file=sprintf("DATA/site_%d_feats.txt", site_i), x=cur_site_feats_matrix, sep = "\t", quote = FALSE, row.names = FALSE);  
  }

  # Write the initial k-means clusters.
	cur_cl_centroids <- matrix(rnorm(n_dim*K_clust)*10, nrow=K_clust, ncol=n_dim);

	write.table(file="DATA/Kclust_centroids.txt", x=cur_cl_centroids, sep="\t", quote = FALSE, row.names = FALSE);

  print("Wrote the data, exiting after data generation..");  
} 
###########################################################################################

