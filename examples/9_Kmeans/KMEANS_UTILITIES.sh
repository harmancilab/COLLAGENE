#!/bin/bash

# This is the backend script that runs mean and weight updates using the current site's data.

if [[ $# -lt 1 ]]
then
	echo "$0 [option] [arguments]
Pipeline:
	-get_X_2_centroids_dist
	-collab_exponentiate_cluster_sq_dists
	-calculate_centroid_updates
	-pool_normalize_centroid_updates
"
	exit 1
fi

N_SITES=3
inv_D2_sigma_sq="-0.2"
weight_refresh_add_noise_sigma="0.06"

if [[ $1 == "-get_X_2_centroids_dist" ]]
then
	if [[ $# -ne 5 ]]
	then
		echo "USAGE: $0 $1 [params dir] [plaintext data matrix file] [Encrypted centroids directory] [Output directory]"
		exit 1
	fi

	params_dir=$2
	pt_data_matrix_file=$3
	enc_centroids_dir=$4
	op_dir=$5

	./COLLAGENE.sh -write_plaintext_matrix_dimensions ${pt_data_matrix_file} ${op_dir}/data_dims.txt

	n_subj=`cat ${op_dir}/data_dims.txt | cut -f1,1`
	n_dim=`cat ${op_dir}/data_dims.txt | cut -f2,2`

        echo "Data is [${n_subj}x${n_dim}]"
	
	echo "Encrypting ${pt_data_matrix_file}"
	./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${pt_data_matrix_file} ${op_dir}/X.enc

	rm -f -r ${op_dir}/X_min_per_clust_centroids_dir
	mkdir ${op_dir}/X_min_per_clust_centroids_dir
	n_k_clust=`cut -f1,1 DATA/CENTROID_DIMS.txt`
	clust_iter=`seq 1 ${n_k_clust}`
	for clust_i in ${clust_iter[@]}
	do
		echo "Calculating enc(X-K_cluster[${clust_i}])"

		# Each centroid is processed separately.
		clust_i_min_one=`echo ${clust_i} | awk {'print $1-1'}`
	        echo "Row expanding encrypted centroids."
		enc_centroids_file=${enc_centroids_dir}/centroid_${clust_i_min_one}.bin.enc
		rm -f -r ${op_dir}/centroid_row_exps
	        mkdir ${op_dir}/centroid_row_exps
	        ./COLLAGENE.sh -secure_row_expand_encrypted_matrix ${params_dir} ${enc_centroids_file} ${n_subj} ${op_dir}/centroid_row_exps

		./COLLAGENE.sh -secure_sub_matrices ${params_dir} ${op_dir}/X.enc ${op_dir}/centroid_row_exps/reprow_0.bin.enc ${op_dir}/X_min_per_clust_centroids_dir/${clust_i}_diff.bin.enc

		echo "Calculating Eucledian dist.."
		./COLLAGENE.sh -secure_row2row_multiply_matrices ${params_dir} ${op_dir}/X_min_per_clust_centroids_dir/${clust_i}_diff.bin.enc ${op_dir}/X_min_per_clust_centroids_dir/${clust_i}_diff.bin.enc ${op_dir}/X_min_per_clust_centroids_dir/sq_${clust_i}_diff.bin.enc

		echo "Negating the distances.."
		./COLLAGENE.sh -generate_constant_plaintext_matrix ${n_subj} 1 ${inv_D2_sigma_sq} ${op_dir}/neg_vector.bin
		./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${op_dir}/neg_vector.bin ${op_dir}/neg_vector.bin.enc

		./COLLAGENE.sh -secure_elementwise_multiply_matrices ${params_dir} ${op_dir}/X_min_per_clust_centroids_dir/sq_${clust_i}_diff.bin.enc ${op_dir}/neg_vector.bin.enc ${op_dir}/X_min_per_clust_centroids_dir/neg_sq_${clust_i}_diff.bin.enc

	done

	# We have the distance of each subject to every current K-means cluster.

	exit 0
fi

if [[ $1 == "-collab_exponentiate_cluster_sq_dists" ]]
then
	if [[ $# -ne 4 ]]
	then
		echo "USAGE: $0 $1 [params dir] [per cluster sq dists directory] [Output directory]"
		exit 1
	fi

	params_dir=$2
	per_cluster_sq_dists_dir=$3
	op_dir=$4

        n_subj=`cat ${op_dir}/data_dims.txt | cut -f1,1`
        n_dim=`cat ${op_dir}/data_dims.txt | cut -f2,2`

        echo "Data is [${n_subj}x${n_dim}]"

	rm -f -r ${op_dir}/per_cluster_weights_dir
	mkdir ${op_dir}/per_cluster_weights_dir

        n_k_clust=`cut -f1,1 DATA/CENTROID_DIMS.txt`
        clust_iter=`seq 1 ${n_k_clust}`
        for clust_i in ${clust_iter[@]}
        do
		# Do collaborative exponentiation on the current sq-distance matrix.
		cur_sq_dists_file=${per_cluster_sq_dists_dir}/neg_sq_${clust_i}_diff.bin.enc

		if [[ ! -f ${cur_sq_dists_file} ]]
		then
			echo "Could not find the cluster cluster-${clust_i} D2 file @ ${cur_sq_dists_file}"
			exit 1
		fi

		echo "Collectively exponentiating ${cur_sq_dists_file}..."

		./function_approximation.sh ${params_dir} ${cur_sq_dists_file} ${op_dir} >& ${op_dir}/CLUST_${clust_i}_D2_EXPONENTIATE.OP

		# Check the eixstence of the output file.
		if [[ ! -f  "${op_dir}/final_exponentiated_random_matrix.bin.enc" ]]
		then
			echo "Could not find collectively exponentiated output: final_exponentiated_random_matrix.bin.enc"
			exit 1
		fi

		# We have the result for this cluster, move to the next.
		mv ${op_dir}/final_exponentiated_random_matrix.bin.enc ${op_dir}/per_cluster_weights_dir/K_cluster_${clust_i}_exp_D2.bin.enc
	done

	exit 0
fi

if [[ $1 == "-calculate_centroid_updates" ]]
then
	if [[ $# -ne 4 ]]
	then
		echo "USAGE: $0 $1 [params dir] [Plaintext data matrix file] [output dir.]"
		exit 1
	fi

	params_dir=$2
	pt_data_matrix=$3
	op_dir=$4

	per_cluster_weights_dir=${op_dir}/per_cluster_weights_dir

	n_subj=`cut -f1,1 ${op_dir}/data_dims.txt`
	n_dim=`cut -f2,2 ${op_dir}/data_dims.txt`

	echo "Transposing data matrix ${pt_data_matrix}"
	./COLLAGENE.sh -transpose_plaintext_matrix ${pt_data_matrix} ${op_dir}/X_trans.bin

	# Following code works but it is too slow, rather, we assume both dimensions are 2^n and use row2row multiplications.
	echo "Encrypting transposed coordinates."
	./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${op_dir}/X_trans.bin ${op_dir}/X_trans.bin.enc

	mkdir ${op_dir}/Updated_Coords_Weights

        n_k_clust=`cut -f1,1 DATA/CENTROID_DIMS.txt`
        clust_iter=`seq 1 ${n_k_clust}`
        for clust_i in ${clust_iter[@]}
        do
		echo "Calculating weighted average coordinates for cluster-${clust_i}"

                if [[ ! -f "${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc" ]]
                then
                        echo "Could not find the weights in ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc"
                        exit 1
                fi

		echo "Transposing weights."
	
		# Transpose the current cluster's weights.
		./COLLAGENE.sh -transpose_encrypted_vector ${params_dir} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc

		if [[ ! -f "${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc" ]] 
		then
			echo "Could not find transposed weights @ ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc"
			exit 1
		fi

                ./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc ${op_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc_vitals.txt

		# Additively refreshing weights.
		# We need to refresh the weights here additively since they lost a whole bunch of levels.
		./additive_refresh_ct.sh ${params_dir} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc ${weight_refresh_add_noise_sigma} ${N_SITES} ${op_dir}/reenc_weights_${clust_i}.bin.enc

		echo "Row expanding transposed weights to ${n_dim}."

		# Row expand to n_dims.
		rm -f -r ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps
                mkdir ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps
		#./COLLAGENE.sh -secure_row_expand_encrypted_matrix ${params_dir} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc ${n_dim} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps
		./COLLAGENE.sh -secure_row_expand_encrypted_matrix ${params_dir} ${op_dir}/reenc_weights_${clust_i}.bin.enc ${n_dim} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps

		echo "======================================="
		echo "Files in expansion (Should be 1 file.):"
		ls -sort ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps
		echo "======================================="

		# This should have generated only one row expansion since the vector has one row.
		echo "Row-row multiplying.."

		# Report vital stats.
		./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${op_dir}/X_trans.bin.enc ${op_dir}/X_trans.bin.enc_vital_stats.txt
                ./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps/reprow_0.bin.enc ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps/reprow_0.bin.enc_vital_stats.txt

		# Do row-2-row multiplication.
		./COLLAGENE.sh -secure_row2row_multiply_matrices ${params_dir} ${op_dir}/X_trans.bin.enc ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_row_exps/reprow_0.bin.enc ${op_dir}/Updated_Coords_Weights/weighted_K_cluster_${clust_i}_coords_update.bin.enc

                if [[ ! -f "${op_dir}/Updated_Coords_Weights/weighted_K_cluster_${clust_i}_coords_update.bin.enc" ]]
                then
                        echo "Could not find the weighted coordinates @ weighted_K_cluster_${clust_i}_coords_update.bin.enc"
                        exit 1
                fi

                echo "Calculating total weight for cluster-${clust_i}."

                # Calculate the total weight for this cluster.
                ./COLLAGENE.sh -generate_constant_plaintext_matrix 1 ${n_subj} 1 ${op_dir}/ones_vec.bin
		./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${op_dir}/ones_vec.bin ${op_dir}/ones_vec.bin.enc

		#./COLLAGENE.sh -secure_row2row_multiply_matrices ${params_dir} ones_vec.bin.enc ${per_cluster_weights_dir}/K_cluster_${clust_i}_exp_D2.bin.enc_trans.bin.enc weighted_K_cluster_${clust_i}_weights_update.bin.enc
		./COLLAGENE.sh -secure_row2row_multiply_matrices ${params_dir} ${op_dir}/ones_vec.bin.enc ${op_dir}/reenc_weights_${clust_i}.bin.enc ${op_dir}/Updated_Coords_Weights/weighted_K_cluster_${clust_i}_weights_update.bin.enc

		if [[ ! -f "${op_dir}/Updated_Coords_Weights/weighted_K_cluster_${clust_i}_weights_update.bin.enc" ]]
		then
			echo "Could not find the total weight update @ weighted_K_cluster_${clust_i}_weights_update.bin.enc"
		fi
	done

	exit 0
fi

if [[ $1 == "-pool_normalize_centroid_updates" ]]
then
	if [[ $# -ne 6 ]]
	then
		echo "USAGE: $0 $1 [params dir] [List of weighted coordinate updates] [Lost of total weights] [Output directory] [New centroids directory]"
		exit 1
	fi

	params_dir=$2
	coord_update_dirs_list=$3
	total_weight_dirs_list=$4
	op_dir=$5
	new_centroids_dir=$6

        n_k_clust=`cut -f1,1 DATA/CENTROID_DIMS.txt`
        clust_iter=`seq 1 ${n_k_clust}`

	mkdir ${op_dir}/Pooled_Normalized_Centroid_Coords

        for clust_i in ${clust_iter[@]}
        do
		echo "Pooling coordinates and total weights for Cluster-${clust_i}."

		awk -v clust_i=${clust_i} {'print $0"/weighted_K_cluster_"clust_i"_coords_update.bin.enc"'} ${coord_update_dirs_list} > ${op_dir}/temp_cur_clust_coord_updates_list.txt
                awk -v clust_i=${clust_i} {'print $0"/weighted_K_cluster_"clust_i"_weights_update.bin.enc"'} ${total_weight_dirs_list} > ${op_dir}/temp_cur_clust_total_weight_updates_list.txt

		echo "Number of coord list:"
		wc -l ${op_dir}/temp_cur_clust_coord_updates_list.txt
		echo "Number of total weights update:"
		wc -l ${op_dir}/temp_cur_clust_total_weight_updates_list.txt

		echo "Pooling coordinates and weights.."
		./COLLAGENE.sh -secure_add_matrix_list ${params_dir} ${op_dir}/temp_cur_clust_coord_updates_list.txt ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc
		./COLLAGENE.sh -secure_add_matrix_list ${params_dir} ${op_dir}/temp_cur_clust_total_weight_updates_list.txt ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_weights.bin.enc

		if [[ ! -f "${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc" ]]
		then
			echo "Could not find ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc"
			exit 1
		fi

		if [[ ! -f "${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_weights.bin.enc" ]]
		then
			echo "Could not find ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_weights.bin.enc"
			exit 1
		fi

		echo "Collaborative normalizing coordinates.."
		./collaborative_normalize_coords.sh ${params_dir} ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_weights.bin.enc ${N_SITES} ${op_dir}	

		if [[ ! -f ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc_normalized.bin.enc_refreshed.bin.enc ]]
		then
			echo "Could not find normalized coordinates @ ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc_normalized.bin.enc_refreshed.bin.enc"
			exit 1
		fi

		# Transpose and save to current centroids; also change the indexing to match here.
		clust_i_min_one=`echo ${clust_i} | awk {'print $1-1'}`
		./COLLAGENE.sh -transpose_encrypted_vector ${params_dir} ${op_dir}/Pooled_Normalized_Centroid_Coords/K_cluster_${clust_i}_pooled_coords.bin.enc_normalized.bin.enc_refreshed.bin.enc ${new_centroids_dir}/centroid_${clust_i_min_one}.bin.enc
	done 
	
fi

