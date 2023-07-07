#!/bin/bash

if [[ $# -lt 1 ]]
then
	echo "USAGE: $0 [option]
Options:
	-generate_data
	-setup_directories
	-run_protocol
	-plot_results"
	exit 1
fi

PT_DATA_DIR=DATA

N_SITES=3
N_ITERS=5

# Check the keys directory.
for((i_site=0; i_site<${N_SITES}; i_site++))
do
	if [[ ! -d "SITE_${i_site}" ]]
	then
		echo "Could not find the keys directory SITE_${i_site}"
		exit 1
	fi
done

if [[ ! -f "COLLAGENE.sh" ]]
then
	echo "Could not fine COLLAGENE.sh script, copy it from under scripts directory."
	exit 1
fi


if [[ $1 == "-generate_data" ]]
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [K-clusters for initing centroids] [Data dimension]"
		exit 1	
	fi

	K_clust=$2
	n_dim=$3

	# These are hardcoded here: Note that this is about data generation and not about the number of clusters we track.
	# The number of generated clusters is set to 4. You can change the number of inherent clusters by changing this table.
	# Rows: The number of sites, 
	# Columns: The number of inherent clusters in the data.
	# The number of subjects at each site must be a power of 2, i.e., summation of rows in this table.
	echo -e "64\t32\t0\t32" >  site_2_cluster_mapper.list
        echo -e "32\t16\t16\t0" >> site_2_cluster_mapper.list
        echo -e "0\t64\t32\t32" >> site_2_cluster_mapper.list	

	rm -f -r DATA
	mkdir DATA
	mean_sigma=6
	sigma_sigma=2
	Rscript generate_data.R site_2_cluster_mapper.list ${K_clust} ${n_dim} ${mean_sigma} ${sigma_sigma}

	site_iter=`seq 1 $N_SITES`
	for i_site in ${site_iter[@]}
	do
		./COLLAGENE.sh -convert_plaintext_matrix_text_2_bin DATA/site_${i_site}_data.txt 0 1 DATA/site_${i_site}_data.bin
	done

	./COLLAGENE.sh -convert_plaintext_matrix_text_2_bin DATA/Kclust_centroids.txt 0 1 DATA/Kclust_centroids.bin

	# Row expand the Kclust's by 1.
	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_0 DATA/Kclust_centroids.bin DATA/Kclust_centroids.bin.enc

	./COLLAGENE.sh -write_encrypted_matrix_dimensions SITE_0 DATA/Kclust_centroids.bin.enc DATA/CENTROID_DIMS.txt

	echo "==============="
	echo "Initial K-means cluster dimensions:"
	cat DATA/CENTROID_DIMS.txt
	echo ""
	echo "==============="

	mkdir DATA/Kcluster_centroids
	./COLLAGENE.sh -row_expand_plaintext_matrix SITE_0 DATA/Kclust_centroids.bin 1 DATA/Kcluster_centroids

	echo "Renaming cluster centroids.."
	n_clusters=`awk {'print $1'} DATA/CENTROID_DIMS.txt`
	cl_iter=`seq 1 ${n_clusters}`
	for cl_i in ${cl_iter[@]}
	do
		cl_i_min_one=`echo $cl_i | awk {'print $1-1'}`
		mv DATA/Kcluster_centroids/reprow_${cl_i_min_one}.bin.enc DATA/Kcluster_centroids/centroid_${cl_i_min_one}.bin.enc
	done

	exit 0
fi

if [[ $1 == "-setup_directories" ]]
then
	# Setting up SITE data:
	echo "Setting up site-local directories for ${N_SITES} sites and ${N_ITERS} iterations."
	site_iters=`seq 1 ${N_SITES}`
	iter_iters=`seq 1 ${N_ITERS}`

	rm -f *.OP

	rm -f -r SITE_DATA
	mkdir SITE_DATA
	for i_iter in ${iter_iters[@]}
	do
		mkdir SITE_DATA/ITER_${i_iter}
		for i_site in ${site_iters[@]}
		do
			mkdir SITE_DATA/ITER_${i_iter}/SITE_${i_site}
		done
	done

	exit 0
fi

if [[ $1 == "-run_protocol" ]]
then
	CUR_CLUSTER_CENTROIDS_DIR=ITER_0_CENTROIDS

	site_iters=`seq 1 ${N_SITES}`
	iter_iters=`seq 1 ${N_ITERS}`

	# Copy the current centroids.
	rm -f -r ${CUR_CLUSTER_CENTROIDS_DIR}
	cp -r ${PT_DATA_DIR}/Kcluster_centroids ${CUR_CLUSTER_CENTROIDS_DIR}

	# Start the iterations.
	for iter_i in ${iter_iters[@]}
	do
		CUR_ITER_WORKDIR=SITE_DATA/ITER_${iter_i}
	
		if [[ ! -d ${CUR_ITER_WORKDIR} ]]
		then
			echo "Could not find iteration directory @ ${CUR_ITER_WORKDIR}"
			exit 1
		fi

		echo "ITER-${iter_i}: Centroids=>${CUR_CLUSTER_CENTROIDS_DIR}"

		# Process all sites for this iteration.
		for i_site in ${site_iters[@]}
		do
			CUR_SITE_ITER_WORKDIR=${CUR_ITER_WORKDIR}/SITE_${i_site}
		
			{
			if [[ ! -d ${CUR_SITE_ITER_WORKDIR} ]]
			then
				echo "Could not find site directory @ ${CUR_SITE_ITER_WORKDIR}"
				exit 1
			fi

			site_i_min_one=`echo ${i_site} | awk {'print $1-1'}`

			echo "@ Iteration-${iter_i}::Site-${i_site}: -get_X_2_centroids_dist"

			./KMEANS_UTILITIES.sh -get_X_2_centroids_dist SITE_${site_i_min_one} ${PT_DATA_DIR}/site_${i_site}_data.bin ${CUR_CLUSTER_CENTROIDS_DIR} ${CUR_SITE_ITER_WORKDIR} >& get_X_2_centroids_dist_${iter_i}_${i_site}.OP

					echo "@ ${iter_i}::Site-${i_site}: -collab_exponentiate_cluster_sq_dists"

			./KMEANS_UTILITIES.sh -collab_exponentiate_cluster_sq_dists SITE_${site_i_min_one} ${CUR_SITE_ITER_WORKDIR}/X_min_per_clust_centroids_dir ${CUR_SITE_ITER_WORKDIR} >& collab_exponentiate_cluster_sq_dists_${iter_i}_${i_site}.OP

					echo "@ ${iter_i}::Site-${i_site}: -calculate_centroid_updates"

			./KMEANS_UTILITIES.sh -calculate_centroid_updates SITE_${site_i_min_one} ${PT_DATA_DIR}/site_${i_site}_data.bin ${CUR_SITE_ITER_WORKDIR} >& calculate_centroid_updates_${iter_i}_${i_site}.OP
			} &
		done # site iteration.

		# Wait for all sites to finish.
		wait

		echo "@ ${iter_i}:: pool_normalize_centroid_updates"

		find ${CUR_ITER_WORKDIR} -name 'Updated_Coords_Weights' > updated_coords_weights_dir.list

		echo "Number of updated coordinates:"
		wc -l updated_coords_weights_dir.list

		rm -f -r ITER_${iter_i}_CENTROIDS
		mkdir ITER_${iter_i}_CENTROIDS
		./KMEANS_UTILITIES.sh -pool_normalize_centroid_updates SITE_0 updated_coords_weights_dir.list updated_coords_weights_dir.list SITE_DATA/ITER_${iter_i}/SITE_1 ITER_${iter_i}_CENTROIDS >& pool_normalize_centroid_updates_${iter_i}.OP

		# Update the cluster centroids directory.
		CUR_CLUSTER_CENTROIDS_DIR=ITER_${iter_i}_CENTROIDS
	done 
fi

if [[ $1 == "-plot_results" ]]
then
	if [[ ! -d ITER_${N_ITERS}_CENTROIDS ]]
	then
		echo "Could not find the final centroids @ ITER_${N_ITERS}_CENTROIDS"
		exit 1
	fi

	echo "Decrypting the final centroids."
	ls ITER_${N_ITERS}_CENTROIDS/*.enc | awk -v N_SITES=${N_SITES} {'print "./collaborative_decrypt_matrix.sh "$0" "N_SITES" >& /dev/null"'} > temp_dec.sh
	chmod 755 temp_dec.sh
	./temp_dec.sh

	# Parse the coordinates, only.
	cat ITER_${N_ITERS}_CENTROIDS/centroid_*.txt | cut -f2- > secure_centroids.txt

	# Plot results.
	Rscript weighted_Kmeans.R ${N_ITERS}
fi # analysis check.

