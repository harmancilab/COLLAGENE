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

# Site indexing is always done with 0-based indices. This is important for key generation step and for collective decryption.
# You also do not have to keep N_SITES parameter in data_config.params file. This file is not required by COLLAGENE.sh but it is needed by file I/O script to describe the common sharing storage options (AWs, SCP, etc).
# Note that the site indexing is not necessary, it is just here for generating an iterator to loop through partial decryptions.
# In principle, the sites only need to know their index to make sure the files they upload can be tracked by other sites.
echo "Read ${N_SITES} from data config file."
N_SITES_MIN_ONE=`echo ${N_SITES} | awk {'print $1-1'}`
site_iters=`seq 0 ${N_SITES_MIN_ONE}`

# We demonstrate generation of a random matrix, processing it, then collectively decrypting using all sites.
# We assume the data is owned by site 0 and is multiplied by some random matrices.
i_site=0

./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix.bin random_matrix.bin.enc

# Multiply the matrix with itself to process it.
./COLLAGENE.sh -secure_elementwise_multiply_matrices SITE_${i_site} random_matrix.bin.enc random_matrix.bin.enc mult_random_matrix.bin.enc

# KEep track of the plaintext result for comparison later.
./COLLAGENE.sh -multiply_elementwise_plaintext_matrices random_matrix.bin random_matrix.bin mult_rand_matrix.bin

ENC_FILE=mult_random_matrix.bin.enc
PT_FILE=mult_rand_matrix.bin

./COLLAGENE.sh -write_encrypted_matrix_vital_stats SITE_0 ${ENC_FILE} orig_vital_stats.txt

ct_scale=`cat orig_vital_stats.txt | awk {'print $3'}`

# Now we process the matrix until it is not possible to multiply it any more.
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

	# Update pt results in parallel.
	./COLLAGENE.sh -multiply_elementwise_plaintext_matrices cur_rand_matrix.bin_norm.bin ${PT_FILE} ${PT_FILE}
done

./COLLAGENE.sh -save_matrix_text ${PT_FILE} ${PT_FILE}.txt

# Replace the exhausted ciphertext.
ENC_FILE=cur_exhausted_ct.enc

# Check the vital's of the exhausted ciphertext, it should be at the bottom of the chain index.
./COLLAGENE.sh -write_encrypted_matrix_vital_stats SITE_0 ${ENC_FILE} exhausted_vital_stats.txt

# Partially decrypt matrix at all sites.
echo "Partially decrypting matrix at all sites.."
rm -f *.partdec
smdging_noise_var_bit_size=40
for decrypting_i_site in ${site_iters[@]}
do
	echo "Partially decrypting ${i_site}.'s data matrix by site-${decrypting_i_site}.."
	./COLLAGENE.sh -partial_decrypt_matrix SITE_${decrypting_i_site} ${ENC_FILE} ${decrypting_i_site} ${smdging_noise_var_bit_size} ${ENC_FILE}_partdec_by_${decrypting_i_site}.partdec
done

# Final step: Pool the partial decryptions on site 0, in principle, all sites would do this locally.
echo "Pooling partial decryptions.."

# save the partdec file names in another file.
ls *.partdec > PARTDECS.list
n_partdecs=`wc -l PARTDECS.list | awk {'print $1'}`
echo "Pooling ${n_partdec} parts.."
./COLLAGENE.sh -pool_partially_decrypted_matrices SITE_${i_site} PARTDECS.list random_matrix_${i_site}.bin.enc_collaborative_dec.bin

echo "Converting decrypted matrix to text.."
./COLLAGENE.sh -save_matrix_text random_matrix_${i_site}.bin.enc_collaborative_dec.bin random_matrix_${i_site}.bin.enc_collaborative_dec.bin.txt

echo "Converting the plaintext overall result to text.."
./COLLAGENE.sh -save_matrix_text ${PT_FILE} overall_plaintext_results.txt

echo "You can compare random_matrix_${i_site}.bin.enc_collaborative_dec.bin.txt and overall_plaintext_results.txt to make sure that they are matching from plaintext calculations and collaborative decryption."


