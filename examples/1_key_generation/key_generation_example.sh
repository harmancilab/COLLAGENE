if [[ $# != 1 ]]
then
	echo "USAGE: $0 [option]
Options:
	-example_generate_DSK_encryption_keys
	-example_decrypt_set_ckks_keys
	-example_test_ckks_keys"

	exit 1
fi

N_SITES=3
cmd_option=$1

# Key generation relies only on COLLAGENE script.
if [[ ! -f "./COLLAGENE.sh" ]]
then
	echo "Could not find COLLAGENE script, copy it from source code scripts/ directory."
	exit 1
fi

# Make all scripts executable.
#dos2unix *.sh
chmod 755 *.sh

# Site indexintg is always done with 0-based indices. This is important only for downloading data from other sites but protocol decides on the file names
# We use 0-based indexing here for naming the local site key folders, which is not necessary. Local sites can setup their keys in any folder. 
# It is, however, necessary to assign a 0-based index to all sites while protocols are executed.
# You also do not have to keep N_SITES in data_config.params file. This file is not required by COLLAGENE.sh but it is needed by FILE I/O module. 
# Site indexing is used to keep track of the sites. The sites need to know their index to make sure the files they upload can be tracked by other sites. 
# Also, importantly, site indices are needed while collectively decrypting the data.
# The keys that we generate and setup in this exercise can be re-used in all remaining exercises.
echo "Using ${N_SITES} sites."
N_SITES_MIN_ONE=`echo ${N_SITES} | awk {'print $1-1'}`
site_iters=`seq 0 ${N_SITES_MIN_ONE}`

echo "N_SITES: ${N_SITES} (${N_SITES_MIN_ONE})"

if [[ ${cmd_option} == "-example_generate_DSK_encryption_keys" ]]
then
	echo "Cleaning directory.."
	rm -f *.txt *.list *.OP *.bin *.enc *.partdec *.dec *_key

	for i_site in ${site_iters[@]}
	do
		./COLLAGENE.sh -generate_DSK_encryption_key ${i_site}
	done

	if [[ ! -f ckks.params ]]
	then
		echo "Could not find ckks.params"
		exit
	fi

	# Pool all public keys into the tar file:
	DSK_ENCRYPTION_KEYS_DIR=CLIENT_DSK_ENCRYPTION_KEYS
	rm -f -r ${DSK_ENCRYPTION_KEYS_DIR}*
	mkdir ${DSK_ENCRYPTION_KEYS_DIR}
	cp *.dsk_enc_public_key ${DSK_ENCRYPTION_KEYS_DIR}
	cp ckks.params ${DSK_ENCRYPTION_KEYS_DIR}

	# Tar the directory
	tar -cvf ${DSK_ENCRYPTION_KEYS_DIR}.tar ${DSK_ENCRYPTION_KEYS_DIR}

	echo "###############################################"
	echo "SUBMIT ${DSK_ENCRYPTION_KEYS_DIR}.tar to KeyMaker @ https://www.secureomics.org/KeyMaker"
	echo "###############################################"

	exit 0
fi

########################################################################################################################
## SERVER SIDE: Generate keys, encrypt them, send back to the clients.
########################################################################################################################

if [[ ${cmd_option} == "-example_decrypt_set_ckks_keys" ]]
then
	if [[ ! -f "DSK_KEYS_FROM_SERVER.tar" ]]
	then
		echo "Could not find the keys from server @ \"DSK_KEYS_FROM_SERVER.tar\""
		exit
	fi

	echo "UnTARring the received keys.."
	rm -f -r RECEIVED_KEYS
	mkdir RECEIVED_KEYS
	tar -xvf DSK_KEYS_FROM_SERVER.tar --strip-components 2 -C RECEIVED_KEYS

	echo "Decrypting keys and testing using directory RECEIVED_KEYS/"
	for i_site in ${site_iters[@]}
	do
		echo "Decrypting the keys for site ${i_site}"
		./COLLAGENE.sh -decrypt_site_DSK ${i_site} RECEIVED_KEYS

		echo "Setting up the keys for site ${i_site}"
		rm -f -r SITE_${i_site}
		mkdir SITE_${i_site}
		./COLLAGENE.sh -set_params ckks.params RECEIVED_KEYS/pooled.public_key RECEIVED_KEYS/pooled.relin_keys RECEIVED_KEYS/pooled.galois_keys RECEIVED_KEYS/site_${i_site}.secret_key SITE_${i_site}

		# Copy the partially decrypted data encryption key, which is a symmetric key.
		cp RECEIVED_KEYS/partdec_data_enc_hash.symmetric_key SITE_${i_site}
	done

	exit 0
fi

if [[ ${cmd_option} == "-example_test_ckks_keys" ]]
then
	for i_site in ${site_iters[@]}
        do
		./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix_${i_site}.bin
		sleep 1

		./COLLAGENE.sh -save_matrix_text random_matrix_${i_site}.bin random_matrix_${i_site}.bin.txt

		# Encrypt.
		./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix_${i_site}.bin random_matrix_${i_site}.bin.enc
	done

	# Do partial decryptions of all sites.
	# Note that there are two loop since we have the site from whose matrix we are decoding and also the "decrypting_i_site", who is currently decrypting the data.
	for i_site in ${site_iters[@]}
	do
		echo "Partially decrypting matrix for site 0.."
		rm -f *.partdec
		for decrypting_i_site in ${site_iters[@]}
		do
			echo "Partially decrypting ${i_site}.'s data matrix by site-${decrypting_i_site}.."
			smdging_noise_var_bit_size=40
			./COLLAGENE.sh -partial_decrypt_matrix SITE_${decrypting_i_site} random_matrix_${i_site}.bin.enc ${decrypting_i_site} ${smdging_noise_var_bit_size} random_matrix_${i_site}.bin.enc_partdec_by_${decrypting_i_site}.partdec

			################################################################################################################################
			## We skip the partdec data encryption/decryption here but you can use partdec encryption key to perform these as:
			#./COLLAGENE.sh -symmetric_encrypt_partdec_data SITE_${decrypting_i_site} random_matrix_${i_site}.bin.enc_partdec_by_${decrypting_i_site}.partdec SITE_${decrypting_i_site}/partdec_data_enc_hash.symmetric_key random_matrix_${i_site}.bin.enc_partdec_by_${decrypting_i_site}.partdec.enc
			#./COLLAGENE.sh -symmetric_decrypt_partdec_data SITE_${decrypting_i_site} random_matrix_${i_site}.bin.enc_partdec_by_${decrypting_i_site}.partdec.enc SITE_${decrypting_i_site}/partdec_data_enc_hash.symmetric_key random_matrix_${i_site}.bin.enc_partdec_by_${decrypting_i_site}.partdec.enc.dec
			# Note that the symmetric partdec key is used for encryption and also for decryption in the above commands.
			################################################################################################################################
		done

		# for this site, decrypt the data using the partial decryptions.
		# Final step: Pool the matrix at all sites.
		echo "Pooling partial decryption for site-${i_site}'s partial decryptions.."
		ls *.partdec > PARTDECS.list
		n_partdecs=`wc -l PARTDECS.list | awk {'print $1'}`
		echo "Pooling ${n_partdec} parts.."
		./COLLAGENE.sh -pool_partially_decrypted_matrices SITE_${i_site} PARTDECS.list random_matrix_${i_site}.bin.enc_collaborative_dec.bin

		echo "Converting decrypted matrix to text.."
		./COLLAGENE.sh -save_matrix_text random_matrix_${i_site}.bin.enc_collaborative_dec.bin random_matrix_${i_site}.bin.enc_collaborative_dec.bin.txt
	done
	
	echo "You can compare random_matrix_*.bin.enc_collaborative_dec.bin.txt and random_matrix_*.bin.txt to make sure that the collaborative decryption yielded the correct matrix."

	exit 0

fi



