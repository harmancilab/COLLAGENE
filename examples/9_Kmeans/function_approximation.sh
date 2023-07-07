#!/bin/bash

if [[ $# != 3 ]]
then
	echo "USAGE: $0 [params dir] [encrypted matrix to exponentiate]"
	exit 1
fi

enc_input_matrix=$2
params_dir=$1
op_dir=$3

if [[ ! -f ${enc_input_matrix} ]]
then
	echo "Could not find encrypted matrix ${enc_input_matrix} to exponentiate."
	exit 1
fi

# Approximate exponential function evaluation to high accuracy using collaboratively predicted weights.
# We focus on exponential as it underlies many functions such as sigmoid.

./COLLAGENE.sh -write_encrypted_matrix_dimensions ${params_dir} ${enc_input_matrix} ${op_dir}/INPUT_MATRIX_DIMs.txt
n_row=`awk {'print $1'} ${op_dir}/INPUT_MATRIX_DIMs.txt`
n_col=`awk {'print $2'} ${op_dir}/INPUT_MATRIX_DIMs.txt`

N_SITES=3

# We have a -- the sensitive value that noone knows.
echo "Read ${N_SITES} from data config file."
N_SITES_MIN_ONE=`echo ${N_SITES} | awk {'print $1-1'}`
site_iters=`seq 0 ${N_SITES_MIN_ONE}`

# We demonstrate generation of a random matrix, processing it, then collectively calculating exp(A) using all sites.
# We assume the data is owned by site 0 and is multiplied by some random matrices.
#i_site=0

#./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix.bin
#./COLLAGENE.sh -scalar_multiply_plaintext_matrix random_matrix.bin 0.5 random_matrix.bin_scaled.bin
#./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix.bin_scaled.bin random_matrix.bin.enc

# Do partial decryptions of one matrix, for site 0
echo "Generating site specific noises"
#echo random_matrix.bin.enc > masking_matrices.list
echo ${enc_input_matrix} > ${op_dir}/masking_matrices.list
for cur_masking_site_i in ${site_iters[@]}
do
	./COLLAGENE.sh -generate_plaintext_mask_matrix SITE_${cur_masking_site_i} ${enc_input_matrix} 1 ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin

	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${cur_masking_site_i} ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin.enc
	echo ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin.enc >> ${op_dir}/masking_matrices.list
done

# Calculate enc(noisy_a)=enc(a + n1+n2+n3) -- HE manner.
./COLLAGENE.sh -secure_add_matrix_list ${params_dir} ${op_dir}/masking_matrices.list ${op_dir}/collective_masked_random_matrix.bin.enc

# Collective decrypt to get enc(noisy_a)=(a+n123) -- This is now noisy value of the noisy_a
./collaborative_decrypt_matrix.sh ${op_dir}/collective_masked_random_matrix.bin.enc ${N_SITES}

# Calculate exp(noisy_a)
./COLLAGENE.sh -transform_plaintext_elementwise_per_callback ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin exp ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin

# Scale s*exp(noisy_a)
noisy_a_scaler=1
./COLLAGENE.sh -scalar_multiply_plaintext_matrix ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin ${noisy_a_scaler} ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin 

# Re-encrypt exp(noisy_a)
./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin.enc

# Now calculate exp(-n1)xexp(-n2)xexp(-n3)
./COLLAGENE.sh -generate_constant_plaintext_matrix ${n_row} ${n_col} 1 ${op_dir}/cumul_unmasking_matrix.bin
./COLLAGENE.sh -encrypt_plaintext_matrix ${params_dir} ${op_dir}/cumul_unmasking_matrix.bin ${op_dir}/cumul_unmasking_matrix.bin.enc
for cur_masking_site_i in ${site_iters[@]}
do
	# Negate the matrix.
	./COLLAGENE.sh -scalar_multiply_plaintext_matrix ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin -1 ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin_neg.bin

	# Take exp(-n_s)
	./COLLAGENE.sh -transform_plaintext_elementwise_per_callback ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin_neg.bin exp ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin
	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${cur_masking_site_i} ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin.enc

	./COLLAGENE.sh -secure_elementwise_multiply_matrices SITE_${cur_masking_site_i} ${op_dir}/cumul_unmasking_matrix.bin.enc ${op_dir}/site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin.enc ${op_dir}/cumul_unmasking_matrix.bin.enc
done

# Finally, remove the final noise from the exponentiated matrix.
# Multiply with exp(-n123).exp(noisy_a) elementwise.
./COLLAGENE.sh -secure_elementwise_multiply_matrices ${params_dir} ${op_dir}/cumul_unmasking_matrix.bin.enc ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin.enc ${op_dir}/final_exponentiated_random_matrix.bin.enc

# Do final collaborative decryption.
./collaborative_decrypt_matrix.sh ${op_dir}/final_exponentiated_random_matrix.bin.enc ${N_SITES}

# This is the masked matrix that is collaboratively decrypted.
./COLLAGENE.sh -save_matrix_text ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin ${op_dir}/collective_masked_random_matrix.bin.enc_collaborative_dec.bin.txt



