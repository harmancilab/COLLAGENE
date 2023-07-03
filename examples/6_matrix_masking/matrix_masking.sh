if [[ ! -d "SITE_0" ]]
then
	echo "Could not find SITE_0 directory, copy it from previous exercise with setup keys."
	exit 1
fi

if [[ ! -f "./COLLAGENE.sh" ]]
then
	echo "Could not find COLLAGENE script, copy it from source code scripts/ directory."
	exit 1
fi

N_SITES=3

rm -f *.list
rm -f *.enc
rm -f *.partdec
rm -f *.bin 

# Generate random matrices of differnt scales and showcase the importance of data scaling especially for padded matrices.
mat_A_nrows=10
mat_A_ncols=10

n_sites_min_one=`echo ${N_SITES} | awk {'print $1-1}'`
site_iters=`seq 0 ${n_sites_min_one}`

rm -f site_matrices.list
for i_site in ${site_iters[@]}
do
	# Generate a random rescaled matrix.
	./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix_${i_site}.bin

	./COLLAGENE.sh -scalar_multiply_plaintext_matrix random_matrix_${i_site}.bin 5 random_matrix_${i_site}.bin_scaled.bin

	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix_${i_site}.bin_scaled.bin random_matrix_${i_site}.bin_scaled.bin.enc

	echo random_matrix_${i_site}.bin_scaled.bin.enc >> site_matrices.list
done

echo "Pooling matrices from sites."
./COLLAGENE.sh -secure_add_matrix_list SITE_0 site_matrices.list pooled_random_matrix.bin.enc

# Exhaust the random data matrix as we did in previous examples.
ENC_FILE=pooled_random_matrix.bin.enc

./COLLAGENE.sh -write_encrypted_matrix_vital_stats SITE_0 ${ENC_FILE} orig_vital_stats.txt

ct_scale=`cat orig_vital_stats.txt | awk {'print $3'}`

# Now we process the matrix until it is not possible to multiply it any more.
echo "Exhausting the total matrix.."
cp ${ENC_FILE} cur_exhausted_ct.enc
scale_iter=`seq 1 ${ct_scale}`
for scale_i in ${scale_iter[@]}
do
        echo "In exhaustion loop @ ${scale_i}"
        ./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 cur_random_matrix.bin

        ./COLLAGENE.sh -scalar_multiply_plaintext_matrix cur_random_matrix.bin 2 cur_rand_matrix.bin_norm.bin

        ./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} cur_rand_matrix.bin_norm.bin cur_rand_matrix.bin_norm.bin.enc >& /dev/null

        ./COLLAGENE.sh -secure_elementwise_multiply_matrices SITE_${i_site} cur_rand_matrix.bin_norm.bin.enc cur_exhausted_ct.enc temp_cur_exhausted_ct.enc >& /dev/null

        # Update the current exhausted ciphertext.
        cp temp_cur_exhausted_ct.enc cur_exhausted_ct.enc
done

# Write the vital stats for the exhausted matrix, this should be on the first level.
./COLLAGENE.sh -write_encrypted_matrix_vital_stats SITE_0 cur_exhausted_ct.enc exhausted_vital_stats.txt

# At this point, we cannot process this matrix any more, so we do masking and refresh it by decypting and re-encrypting. This is an ad-hoc way to refresh the matrix by using an error term. In principle the error can be aedded using DP-type approaches, too. As long as the CKKS parameters are appropriately set, the error should not dominate.
rm -f enc_mask_matrix_list.txt
for i_site in ${site_iters[@]}
do
	# Generate a masking matrix from this site, this is simple a random matrix with a large variance.
	./COLLAGENE.sh -generate_plaintext_mask_matrix SITE_${i_site} cur_exhausted_ct.enc 30 site_${i_site}_mask_matrix.bin

	# Encrypt the mask matrix.
	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} site_${i_site}_mask_matrix.bin site_${i_site}_mask_matrix.bin.enc
	
	echo site_${i_site}_mask_matrix.bin.enc >> enc_mask_matrix_list.txt
done

# Add all the noise levels from all sites.
./COLLAGENE.sh -secure_add_matrix_list SITE_0 enc_mask_matrix_list.txt pooled_random_mask_matrix.bin.enc

# Mask the encrypted matrix:
./COLLAGENE.sh -mask_encrypted_matrix SITE_0 cur_exhausted_ct.enc pooled_random_mask_matrix.bin.enc masked_exhausted_ct.enc

# Now we can decrypt the echausted matrix as it contains noise from all sites.
./collaborative_decrypt_matrix.sh masked_exhausted_ct.enc ${N_SITES}

# Re-encrypt the data matrix.
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_0 masked_exhausted_ct.enc_collaborative_dec.bin refreshed_exhausted_ct.enc

# Unmask noise.
./COLLAGENE.sh -unmask_encrypted_matrix SITE_0 refreshed_exhausted_ct.enc pooled_random_mask_matrix.bin.enc unmasked_refreshed_exhausted_ct.enc

# Write the vital statistics of the refreshed ciphertext to ensure we are at the top level. 
./COLLAGENE.sh -write_encrypted_matrix_vital_stats SITE_0 unmasked_refreshed_exhausted_ct.enc refreshed_vital_stats.txt

# Decrypt for comparison of the exhausted and refreshed ciphertexts:
./collaborative_decrypt_matrix.sh masked_exhausted_ct.enc ${N_SITES}
./collaborative_decrypt_matrix.sh unmasked_refreshed_exhausted_ct.enc ${N_SITES}
./collaborative_decrypt_matrix.sh cur_exhausted_ct.enc ${N_SITES}

echo "You can compare cur_exhausted_ct.enc_collaborative_dec.bin.txt and unmasked_refreshed_exhausted_ct.enc_collaborative_dec.bin.txt"




