#!/bin/bash 

if [[ $# -ne 5 ]]
then
	echo "USAGE: $0 [params dir] [Pooled coords file] [Pooled weight file] [# sites] [output dir]"
	exit 1
fi

params_dir=$1
pooled_coords_enc_mat=$2
pooled_weights_enc_mat=$3
N_SITES=$4
op_dir=$5

if [[ ! -f ${pooled_coords_enc_mat} ]]
then
	echo "Could not find ${pooled_coords_enc_mat}"
	exit 1
fi

if [[ ! -f ${pooled_weights_enc_mat} ]]
then
	echo "Could not find ${pooled_weights_enc_mat}"
	exit 1
fi

site_iter=`seq 1 ${N_SITES}`

./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${pooled_coords_enc_mat} ${pooled_coords_enc_mat}_coords_vital.txt

./COLLAGENE.sh -write_encrypted_matrix_dimensions ${params_dir} ${pooled_coords_enc_mat} ${pooled_coords_enc_mat}_coord_dims.txt

n_rows=`awk {'print $1'} ${pooled_coords_enc_mat}_coord_dims.txt`
n_cols=`awk {'print $2'} ${pooled_coords_enc_mat}_coord_dims.txt`

echo "Expanding pooled weights to ${n_rows} rows."

rm -f -r ${pooled_weights_enc_mat}_rowexps
mkdir ${pooled_weights_enc_mat}_rowexps
./COLLAGENE.sh -secure_row_expand_encrypted_matrix ${params_dir} ${pooled_weights_enc_mat} ${n_rows} ${pooled_weights_enc_mat}_rowexps

# Collab invert the total weights.
rm -f ${pooled_weights_enc_mat}_noise_matrices.txt

for i_site in ${site_iter[@]}
do
        i_site_min_one=`echo ${i_site} | awk {'print $1-1'}`

        ./COLLAGENE.sh -generate_random_plaintext_matrix ${n_rows} ${n_cols} ${pooled_weights_enc_mat}_SITE_${i_site}_noise.bin

        ./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site_min_one} ${pooled_weights_enc_mat}_SITE_${i_site}_noise.bin ${pooled_weights_enc_mat}_SITE_${i_site}_noise.bin.enc

        echo ${pooled_weights_enc_mat}_SITE_${i_site}_noise.bin.enc >> ${pooled_weights_enc_mat}_noise_matrices.txt
done

./COLLAGENE.sh -secure_add_matrix_list ${params_dir} ${pooled_weights_enc_mat}_noise_matrices.txt ${pooled_weights_enc_mat}_pooled_noise.bin.enc

./COLLAGENE.sh -secure_elementwise_multiply_matrices ${params_dir} ${pooled_weights_enc_mat}_pooled_noise.bin.enc ${pooled_coords_enc_mat} ${pooled_coords_enc_mat}_noisy.bin.enc
./COLLAGENE.sh -secure_elementwise_multiply_matrices ${params_dir} ${pooled_weights_enc_mat}_pooled_noise.bin.enc ${pooled_weights_enc_mat}_rowexps/reprow_0.bin.enc ${pooled_weights_enc_mat}_noisy.bin.enc

./collaborative_decrypt_matrix.sh ${pooled_weights_enc_mat}_noisy.bin.enc ${N_SITES}

./COLLAGENE.sh -transform_plaintext_elementwise_per_callback ${pooled_weights_enc_mat}_noisy.bin.enc_collaborative_dec.bin inv  ${pooled_weights_enc_mat}_noisy.bin.enc_collaborative_dec.bin_inv.bin

./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${pooled_weights_enc_mat}_noisy.bin.enc_collaborative_dec.bin_inv.bin ${pooled_weights_enc_mat}_noisy.bin.enc_collaborative_dec.bin_inv.bin.enc

# Final normalization.
./COLLAGENE.sh -secure_elementwise_multiply_matrices ${params_dir} ${pooled_weights_enc_mat}_noisy.bin.enc_collaborative_dec.bin_inv.bin.enc ${pooled_coords_enc_mat}_noisy.bin.enc ${pooled_coords_enc_mat}_normalized.bin.enc

./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${pooled_coords_enc_mat}_normalized.bin.enc ${pooled_coords_enc_mat}_normalized.bin.enc_vitals.txt

# Transpose the coordinates and get them ready for the next iteration.
refresh_sigma=10
./additive_refresh_ct.sh ${params_dir} ${pooled_coords_enc_mat}_normalized.bin.enc ${refresh_sigma} ${N_SITES} ${pooled_coords_enc_mat}_normalized.bin.enc_refreshed.bin.enc

# This should be a fresh ciphertext for the centroid coordinates.
./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${pooled_coords_enc_mat}_normalized.bin.enc_refreshed.bin.enc ${pooled_coords_enc_mat}_normalized.bin.enc_refreshed.bin.enc_vitals.txt

exit 0

