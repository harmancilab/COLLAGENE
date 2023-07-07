#!/bin/bash

if [[ $# != 5 ]]
then
	echo "USAGE: $0 [params directory] [Encrypted matrix path] [noise var] [# sites] [Output enc matrix file]"
	exit 1
fi

params_dir=$1
enc_matrix_file=$2
noise_sigma=$3
N_sites=$4
op_enc_mat_file=$5

./COLLAGENE.sh -write_encrypted_matrix_dimensions ${params_dir} ${enc_matrix_file} ${enc_matrix_file}_dims.txt

nrow=`awk {'print $1'} ${enc_matrix_file}_dims.txt`
ncol=`awk {'print $2'} ${enc_matrix_file}_dims.txt`

site_iter=`seq 1 ${N_sites}`

rm -f ${enc_matrix_file}_noise_matrices.txt

for i_site in ${site_iter[@]}
do
	i_site_min_one=`echo ${i_site} | awk {'print $1-1'}`

	./COLLAGENE.sh -generate_random_plaintext_matrix ${nrow} ${ncol} ${enc_matrix_file}_SITE_${i_site}_noise.bin

	./COLLAGENE.sh -scalar_multiply_plaintext_matrix ${enc_matrix_file}_SITE_${i_site}_noise.bin ${noise_sigma} ${enc_matrix_file}_SITE_${i_site}_noise.bin_normalized.bin

	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site_min_one} ${enc_matrix_file}_SITE_${i_site}_noise.bin_normalized.bin ${enc_matrix_file}_SITE_${i_site}_noise.bin.enc

	echo ${enc_matrix_file}_SITE_${i_site}_noise.bin.enc >> ${enc_matrix_file}_noise_matrices.txt
done

./COLLAGENE.sh -secure_add_matrix_list ${params_dir} ${enc_matrix_file}_noise_matrices.txt ${enc_matrix_file}_pooled_noise.bin.enc

./COLLAGENE.sh -secure_add_matrices ${params_dir} ${enc_matrix_file}_pooled_noise.bin.enc ${enc_matrix_file} ${enc_matrix_file}_noisy_mat.bin.enc

./collaborative_decrypt_matrix.sh ${enc_matrix_file}_noisy_mat.bin.enc ${N_sites}

./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${enc_matrix_file}_noisy_mat.bin.enc_collaborative_dec.bin ${enc_matrix_file}_noisy_mat.bin.enc_collaborative_dec.bin.enc

./COLLAGENE.sh -secure_sub_matrices ${params_dir} ${enc_matrix_file}_noisy_mat.bin.enc_collaborative_dec.bin.enc ${enc_matrix_file}_pooled_noise.bin.enc ${enc_matrix_file}_refreshed_mat.bin.enc

./COLLAGENE.sh -write_encrypted_matrix_vital_stats ${params_dir} ${enc_matrix_file}_refreshed_mat.bin.enc ${enc_matrix_file}_refreshed_mat.bin.enc_vitals.txt

cp ${enc_matrix_file}_refreshed_mat.bin.enc ${op_enc_mat_file}

./collaborative_decrypt_matrix.sh ${enc_matrix_file} ${N_sites}
./collaborative_decrypt_matrix.sh ${op_enc_mat_file} ${N_sites}



