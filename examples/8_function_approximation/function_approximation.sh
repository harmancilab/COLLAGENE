# Approximate exponential function evaluation to high accuracy using collaboratively predicted weights.
# We focus on exponential as it underlies many functions such as sigmoid.
# It should also be noted that the accuracy in this example may not satisfy different protocols. It is therefore necessary to tune the noise parameters accordingly.

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

# Make all scripts executable.
#dos2unix *.sh
chmod 755 *.sh

N_SITES=3

# We have a -- the sensitive value that noone knows.
echo "Read ${N_SITES} from data config file."
N_SITES_MIN_ONE=`echo ${N_SITES} | awk {'print $1-1'}`
site_iters=`seq 0 ${N_SITES_MIN_ONE}`

# We demonstrate generation of a random matrix, processing it, then collectively calculating exp(A) using all sites.
# We assume the data is owned by site 0 and is multiplied by some random matrices.
i_site=0

./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix.bin
./COLLAGENE.sh -scalar_multiply_plaintext_matrix random_matrix.bin 0.5 random_matrix.bin_scaled.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix.bin_scaled.bin random_matrix.bin.enc

# Do partial decryptions of one matrix, for site 0
echo "Generating site specific noises"
echo random_matrix.bin.enc > masking_matrices.list
for cur_masking_site_i in ${site_iters[@]}
do
	# We use uniform noise in this example. Note that exponentials can easily overflow or underflow so the noise variance should not be set too high.
	./COLLAGENE.sh -generate_plaintext_mask_matrix SITE_${cur_masking_site_i} random_matrix.bin.enc 3 site_${cur_masking_site_i}_mask_matrix.bin

	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${cur_masking_site_i} site_${cur_masking_site_i}_mask_matrix.bin site_${cur_masking_site_i}_mask_matrix.bin.enc
	echo site_${cur_masking_site_i}_mask_matrix.bin.enc >> masking_matrices.list
done

# Calculate enc(noisy_a)=enc(a + n1+n2+n3) -- HE manner.
./COLLAGENE.sh -secure_add_matrix_list SITE_0 masking_matrices.list collective_masked_random_matrix.bin.enc

# Collective decrypt to get enc(noisy_a)=(a+n123) -- This is now noisy value of the noisy_a
./collaborative_decrypt_matrix.sh collective_masked_random_matrix.bin.enc ${N_SITES}

# Calculate exp(noisy_a)
./COLLAGENE.sh -transform_plaintext_elementwise_per_callback collective_masked_random_matrix.bin.enc_collaborative_dec.bin exp collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin

# Re-encrypt exp(noisy_a)
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_0 collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin.enc

# Now calculate exp(-n1)xexp(-n2)xexp(-n3)
./COLLAGENE.sh -generate_constant_plaintext_matrix 10 10 1 cumul_unmasking_matrix.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_0 cumul_unmasking_matrix.bin cumul_unmasking_matrix.bin.enc
for cur_masking_site_i in ${site_iters[@]}
do
	# Negate the matrix.
	./COLLAGENE.sh -scalar_multiply_plaintext_matrix site_${cur_masking_site_i}_mask_matrix.bin -1 site_${cur_masking_site_i}_mask_matrix.bin_neg.bin

	# Take exp(-n_s)
	./COLLAGENE.sh -transform_plaintext_elementwise_per_callback site_${cur_masking_site_i}_mask_matrix.bin_neg.bin exp site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin
	./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${cur_masking_site_i} site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin.enc

	./COLLAGENE.sh -secure_elementwise_multiply_matrices SITE_${cur_masking_site_i} cumul_unmasking_matrix.bin.enc site_${cur_masking_site_i}_mask_matrix.bin_neg_exp.bin.enc cumul_unmasking_matrix.bin.enc
done

# Finally, remove the final noise from the exponentiated matrix.
# Multiply with exp(-n123).exp(noisy_a) elementwise.
./COLLAGENE.sh -secure_elementwise_multiply_matrices SITE_0 cumul_unmasking_matrix.bin.enc collective_masked_random_matrix.bin.enc_collaborative_dec.bin_exp.bin.enc final_exponentiated_random_matrix.bin.enc

# Do final collaborative decryption.
./collaborative_decrypt_matrix.sh final_exponentiated_random_matrix.bin.enc ${N_SITES}

# This is the masked matrix that is collaboratively decrypted.
./COLLAGENE.sh -save_matrix_text collective_masked_random_matrix.bin.enc_collaborative_dec.bin collective_masked_random_matrix.bin.enc_collaborative_dec.bin.txt

# Generate the plaintext result for comparison.
# Write the original scaled matrix to check the noise swamping.
./COLLAGENE.sh -save_matrix_text random_matrix.bin_scaled.bin random_matrix.bin_scaled.bin.txt

# Calculate the exponential on the original plaintext matrix.
./COLLAGENE.sh -transform_plaintext_elementwise_per_callback random_matrix.bin_scaled.bin exp random_matrix.bin_scaled.bin_exp.bin

# Save as a text readable matrix.
./COLLAGENE.sh -save_matrix_text random_matrix.bin_scaled.bin_exp.bin random_matrix.bin_scaled.bin_exp.bin.txt

echo "Compare the exponentiated matrices: random_matrix.bin_scaled.bin_exp.bin.txt and final_exponentiated_random_matrix.bin.enc_collaborative_dec.bin.txt"

echo "You can also compare the original matrix (random_matrix.bin_scaled.bin.txt), and the collectively decrypted masked matrix (collective_masked_random_matrix.bin.enc_collaborative_dec.bin.txt)"
