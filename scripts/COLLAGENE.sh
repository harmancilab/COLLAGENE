#!/bin/bash

# These are the hardcoded key file names, they are here to make it easier to set them up.
PUBLIC_KEY_FILE=pooled.public_key
RELIN_KEY_FILE=pooled.relin_keys
GALOIS_KEY_FILE=pooled.galois_keys
PRIVATE_KEY_FILE=local_site.secret_key
PARTDEC_SYMKEY_FILE=partdec_data_enc_hash.symmetric_key
TEXT_PARAMS_PATH=ckks.params
COLLAGENE_SECURE_EXEC=COLLAGENE_Release
MAX_N_WORKERS=5

# This is the wrapper script for COLLAGENE's utilities. You can add this script to your path and use it to process data.
if [[ $# -lt 2 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
Options:
	Key Generation:
		-generate_DSK_encryption_key [Site Index (0-based)]
		-decrypt_site_DSK [Site Index (0-based)] [Encrypted DSK keys directory]
		-generate_openssl_hash [Site Index (0-based)]
	Initialization:
		-set_params [Text parameters path] [Public key path] [Relin. key path] [Galois key path] [Site's secret key share path] [Parameter directory]
		-validate_ckks_text_params [Parameter directory] [Validity output file]
	Secure Matrix Operations:
		-encrypt_plaintext_matrix [Parameter directory] [Plaintext matrix path] [Output matrix path]
		-full_decrypt_encrypted_matrix [Parameter directory] [Encrypted matrix path] [Output matrix path]
		-transpose_encrypted_vector [Parameter directory] [Encrypted matrix path (1xn or nx1)] [Output matrix path]
		-secure_multiply_rowcol_expansion [Parameter directory] [Encrypted matrix A column expansion Directory] [Encrypted matrix B row expansion Directory] [Output matrix path]
		-row_expand_plaintext_matrix [Parameter directory] [Plaintext matrix path] [# rows per expansion] [Encrypted row expansion directory]
		-secure_row_expand_encrypted_matrix [Parameter directory] [Encrypted matrix A path] [# rows per expansion] [Encrypted row expansion directory]
		-col_expand_plaintext_matrix [Parameter directory] [Plaintext matrix path] [# columns per expansion] [Encrypted column expansion directory]
		-pool_col_expanded_matrices [Parameter directory] [Column expanded matrix directories list file] [Output matrix file]
		-pool_row_expanded_matrices [Parameter directory] [Row expanded matrix directories list file] [Output matrix file]
		-secure_row2row_multiply_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]
		-secure_elementwise_multiply_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]
		-secure_add_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]
		-secure_add_matrix_list [Parameter directory] [Encrypted matrix list path] [Output matrix path]
		-secure_sub_matrices [Parameter directory] [Encrypted matrix A path] [matrix B path] [Output matrix path]
		-write_encrypted_matrix_dimensions [Parameter directory] [Encrypted matrix path] [Dimensions file path]
		-write_encrypted_matrix_vital_stats [Parameter directory] [Encrypted matrix path] [Vital stats file path]
	Collaborative Decryption:
		-partial_decrypt_matrix [Parameter directory] [Encrypted matrix path] [Site index (0-based)] [Decryption noise variance bit size (40-bits by default)] [Output matrix path]
		-pool_partially_decrypted_matrices [Parameter directory] [Partial decryptions list path] [Output matrix path]
	Symmetric Encryption/Decryption (openssl):
		-symmetric_encrypt_partdec_data [Parameter directory] [Partdec file] [openssl hash file] [Output file]
		-symmetric_decrypt_partdec_data [Parameter directory] [Encrypted partdec file] [openssl hash file] [Output file]
	Secure Matrix Masking:
		-generate_plaintext_mask_matrix [Parameter directory] [Encrypted matrix file] [Mask noise standard deviation] [Output plaintext matrix file]
		-mask_encrypted_matrix [Parameter directory] [Encrypted matrix path] [Encrypted mask matrix path] [Output matrix path]
		-unmask_encrypted_matrix [Parameter directory] [Encrypted matrix path] [Encrypted mask matrix path] [Output matrix path]
	Plaintext Matrix Operations (Does not require keys or parameter setup):
		-generate_random_plaintext_matrix [# rows] [# columns] [Output matrix path]
		-generate_constant_plaintext_matrix [# rows] [# columns] [Scalar value] [Output matrix path]
		-generate_diagonal_plaintext_matrix [# rows] [# columns] [Scalar value] [Output matrix path]
		-scalar_multiply_plaintext_matrix [Plaintext matrix path] [Scalar value] [Output matrix path]
		-invert_plaintext_matrix [Plaintext matrix path] [Output matrix path]
		-transpose_plaintext_matrix [Plaintext matrix path] [Output matrix path]
		-multiply_plaintext_matrix [Matrix A file] [Matrix B file] [Output plaintext matrix file]
		-multiply_elementwise_plaintext_matrices [Matrix A file] [Matrix B file] [Output plaintext matrix file]
		-row2row_multiply_plaintext_matrices [Matrix A file] [Matrix B file] [Output plaintext matrix file]
		-add_plaintext_matrix [Matrix A file] [Matrix B file] [Output plaintext matrix file]
		-sub_plaintext_matrix [Matrix A file] [Matrix B file] [Output plaintext matrix file]
		-add_plaintext_matrix_list [Plaintext matrix list file] [Output plaintext matrix file]
		-pad_plaintext_matrix_2_to_n [Plaintext matrix path] [Output matrix path]
		-pad_plaintext_matrix_row_2_to_n [Plaintext matrix path] [Output matrix path]
		-pad_plaintext_matrix_col_2_to_n [Plaintext matrix path] [Output matrix path]
		-pad_plaintext_matrix_custom [Plaintext matrix path] [Padded # rows] [Padded # columns] [Output matrix path]
		-save_matrix_text [Plaintext matrix path] [Output text matrix path]
		-transform_plaintext_elementwise_per_callback [Plaintext matrix path] [Function: \"log\"/\"exp\"/\"sigmoid\"/\"inv\"] [Output matrix path]
		-convert_plaintext_matrix_text_2_bin [Matrix A file] [Row ids in input (first col.)? (0/1)] [Col. id's in input (first row)? (0/1)] [Output plaintext matrix file]
		-write_plaintext_matrix_dimensions [Parameter directory] [Plaintext matrix path] [Dimensions file path]"
		
	exit
fi

cmd_option=$1

collagene_check=`type -P ${COLLAGENE_SECURE_EXEC}`
if [[ "${collagene_check}" == "" ]]
then
	echo "Could not find COLLAGENE executable."
	exit 1
fi

openssl_check=`type -P openssl`
if [[ "${openssl_check}" == "" ]]
then
	echo "Could not find \"openssl\"."
	exit 1
fi

if [[ "${cmd_option}" == "-validate_params_dir" ]]
then
	if [[ $# != 2 ]]
	then
		echo "$0 $1 [Parameters directory]"
		exit 1
	fi

	params_dir=$2
	
	if [[ ! -d ${params_dir} ]]
	then
		echo "Could not find ${params_dir}"
		exit 1
	fi

	if [[ ! -f ${params_dir}/${TEXT_PARAMS_PATH} ]]
	then
		echo "Could not find text params under ${params_dir}/"
		exit 1
	fi

	if [[ ! -f ${params_dir}/${PUBLIC_KEY_FILE} ]]
	then
		echo "Could not find public key under ${params_dir}/"
		exit 1
	fi
	
	if [[ ! -f ${params_dir}/${RELIN_KEY_FILE} ]]
	then
		echo "Could not find relinearization keys under ${params_dir}/"
		exit 1
	fi

	if [[ ! -f ${params_dir}/${GALOIS_KEY_FILE} ]]
	then
		echo "Could not find Galois keys under ${params_dir}/"
		exit 1
	fi

	if [[ ! -f ${params_dir}/${PRIVATE_KEY_FILE} ]]
	then 
		echo "Could not find site private key under ${params_dir}/"
		exit 1
	fi

	exit 0
fi


if [[ "${cmd_option}" == "-generate_DSK_encryption_key" ]]
then
	if [[ $# != 2 ]]
	then
		echo "USAGE: $0 $1 [Site Index (0-based)]"
		exit 1
	fi

	site_i=$2

	# generate the private key.
	openssl genrsa -aes128 -out site_${site_i}.dsk_enc_secret_key 1024
	res=$?

	if [[ ${res} != 0 ]]
	then
		exit ${res} 
	fi

	# Extract the public key.
	openssl rsa -in site_${site_i}.dsk_enc_secret_key -pubout > site_${site_i}.dsk_enc_public_key
	res=$?

	exit ${res}
fi

if [[ "${cmd_option}" == "-decrypt_site_DSK" ]] 
then
        if [[ $# != 3 ]]
        then
                echo "USAGE: $0 $1 [Site Index (0-based)] [Encrypted DSK keys directory]"
                exit 1
        fi

	i_site=$2
	RECEIVED_KEYS_DIR=$3

	if [[ ! -f site_${i_site}.dsk_enc_secret_key ]]
	then
		echo "Could not find the DSK encryption secret key @ site_${i_site}.dsk_enc_secret_key"
		exit
	fi

	# Decrypt the hash for this site.
	echo "Decrypting hash for site ${i_site}"
	if [[ ! -f ${RECEIVED_KEYS_DIR}/site_${i_site}_dsk_rand_hash.bin.enc ]]
	then
		echo "Could not find the hash file @ ${RECEIVED_KEYS_DIR}/site_${i_site}_dsk_rand_hash.bin.enc"
		exit
	fi

	openssl rsautl -decrypt -inkey site_${i_site}.dsk_enc_secret_key -in ${RECEIVED_KEYS_DIR}/site_${i_site}_dsk_rand_hash.bin.enc -out ${RECEIVED_KEYS_DIR}/site_${i_site}_dsk_rand_hash.bin
	if [[ $? != 0 ]]
	then
		echo "Could not decrypt the hash for site ${i_site}"
		exit 1
	fi

	# Decrypt the private key share for this site.
	echo "Decrypting secret key share for site ${i_site}"
	if [[ ! -f ${RECEIVED_KEYS_DIR}/site_${i_site}.secret_key.enc ]]
	then
		echo "Could not find the encrypted secret key share @ ${RECEIVED_KEYS_DIR}/site_${i_site}.secret_key.enc"
		exit
	fi

	openssl enc -d -aes-256-cbc -in ${RECEIVED_KEYS_DIR}/site_${i_site}.secret_key.enc -out ${RECEIVED_KEYS_DIR}/site_${i_site}.secret_key -pass file:./${RECEIVED_KEYS_DIR}/site_${i_site}_dsk_rand_hash.bin
	if [[ $? != 0 ]]
	then
		echo "Could not decrypt the secret key for site ${i_site}"
		exit 1
	fi

	# Decrypt the partdec data symmetric has for this site.
	echo "Decrypting partdec data encryption symmetric key (hash) for site ${i_site}"
	if [[ ! -f ${RECEIVED_KEYS_DIR}/site_${i_site}_partdec_data_enc_hash.symmetric_key.enc ]]
	then
		echo "Could not find the encrypted partdec data encrypted hash @ ${RECEIVED_KEYS_DIR}/site_${i_site}_partdec_data_enc_hash.symmetric_key.enc"
		exit
	fi

	openssl enc -d -aes-256-cbc -in ${RECEIVED_KEYS_DIR}/site_${i_site}_partdec_data_enc_hash.symmetric_key.enc -out ${RECEIVED_KEYS_DIR}/partdec_data_enc_hash.symmetric_key -pass file:./${RECEIVED_KEYS_DIR}/site_${i_site}_dsk_rand_hash.bin
	if [[ $? != 0 ]]
	then
		echo "Could not decrypt the symmetric key."
		exit 1
	fi
fi

if [[ "${cmd_option}" == "-generate_openssl_hash" ]]
then
	if [[ $# != 2 ]]
	then
		echo "USAGE: $0 $1 [Site Index (0-based)]"
		exit 1
	fi

	site_i=$2

	# Create a simple hash
	openssl rand -base64 32 > SITE_${site_i}.hash
	res=$?

	exit ${res}
fi

if [[ "${cmd_option}" == "-set_params" ]]
then
	#-set_params [Text parameters path] [Public key path] [Relin. key path] [Galois key path] [Master secret key path (optional)] [Parameter directory]
	if [[ $# != 7 ]]
	then
		echo "USAGE: $0 $1 [Text parameters path] [Public key path] [Relin. key path] [Galois key path] [Site's secret key share path] [Parameter directory]"
		exit
	fi

	text_params_file=$2
	public_key_path=$3
	relin_keys_path=$4
	galois_keys_path=$5
	private_key_share_path=$6
	params_dir=$7

	# Make sure keys/params exist.
	if [[ ! -f ${text_params_file} ]]
	then
		echo "Could not find the text params file @ ${text_params_file}"
		exit
	elif [[ ! -f ${public_key_path} ]]
	then
		echo "Could not find the public key file @ ${public_key_path}"
		exit
	elif [[ ! -f ${relin_keys_path} ]]
	then
		echo "Could not find the relinearization keys file @ ${relin_keys_path}"
		exit
	elif [[ ! -f ${galois_keys_path} ]]
	then
		echo "Could not find the Galois keys file @ ${galois_keys_path}"
		exit
	elif [[ ! -f ${private_key_share_path} ]]
	then
		echo "Could not find the site's secret key share file @ ${private_key_share_path}"
		exit
	fi

	# Make sure parameters directory exist.
	if [[ ! -d ${params_dir} ]]
	then
		echo "Creating parameters directory @ ${params_dir}"
		mkdir ${params_dir}
	fi

	if [[ ! -d ${params_dir} ]]
	then
		echo "Could not create parameters directory @ ${params_dir}"
		exit
	fi

	# Copy the keys.
	echo "Copying keys to ${params_dir}"
	cp ${text_params_file} ${params_dir}/${TEXT_PARAMS_PATH}
	cp ${public_key_path} ${params_dir}/${PUBLIC_KEY_FILE}
	cp ${relin_keys_path} ${params_dir}/${RELIN_KEY_FILE}
	cp ${galois_keys_path} ${params_dir}/${GALOIS_KEY_FILE}
	cp ${private_key_share_path} ${params_dir}/${PRIVATE_KEY_FILE}

	echo "Done, you can use params_dir=${PWD}/${params_dir} for calculations."
fi


if [[ "${cmd_option}" == "-encrypt_plaintext_matrix" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Plaintext matrix path] [Output matrix path]"
		exit 1
	fi

	params_dir=$2
	plaintext_matrix_file=$3
	enc_matric_file=$4

	${COLLAGENE_SECURE_EXEC} -continuous_encrypt_data_matrix ${plaintext_matrix_file} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${enc_matric_file}

	exit $?
fi

if [[ "${cmd_option}" == "-full_decrypt_encrypted_matrix" ]] 
then
if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Master secret key file] [Encrypted matrix path] [Output plaintext matrix path]"
		exit 1
	fi

	params_dir=$2
	master_key_file=$3
	enc_matric_file=$4
	plaintext_matrix_file=$5

	#-full_decrypt_encrypted_matrix [Parameter directory] [Encrypted matrix path] [Output matrix path]
	${COLLAGENE_SECURE_EXEC} -fully_decrypt_continuous_encrypted_matrix ${enc_matric_file} ${params_dir}/${TEXT_PARAMS_PATH} ${master_key_file} ${plaintext_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-transpose_encrypted_vector" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Encrypted matrix path] [Output encrypted matrix path]"
		exit
	fi

	params_dir=$2
	enc_matric_file=$3
	enc_trans_matrix_file=$4

	${COLLAGENE_SECURE_EXEC} -transpose_continuous_encrypted_vector ${enc_matric_file} ${params_dir}/${TEXT_PARAMS_PATH} ${enc_trans_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-secure_multiply_rowcol_expansion" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [A col. expansion directory] [B row expansion directory] [Output encrypted matrix path]"
		exit
	fi

	params_dir=$2
	matA_col_exp_dir=$3
	matB_row_exp_dir=$4
	enc_mult_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -secure_multiply_matrices_Acol_Brow_expansions ${matA_col_exp_dir} ${matB_row_exp_dir} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} ${params_dir}/${PRIVATE_KEY_FILE} ${enc_mult_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-row_expand_plaintext_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Plaintext matrix path] [# rows per expansion] [Encrypted row expansion directory]"
		exit
	fi

	params_dir=$2
	plaintext_matrix_file=$3
	n_rows_in_expansions=$4
	enc_row_expansion_dir=$5

	${COLLAGENE_SECURE_EXEC} -row_expand_dense_encrypt_matrix ${plaintext_matrix_file} ${n_rows_in_expansions} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${enc_row_expansion_dir}

	exit $?
fi

if [[ "${cmd_option}" == "-secure_row_expand_encrypted_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [encrypted matrix path] [# rows per expansion] [Encrypted row expansion directory]"
		exit
	fi

	params_dir=$2
	encrypted_matrix_file=$3
	n_rows_in_expansions=$4
	enc_row_expansion_dir=$5

	${COLLAGENE_SECURE_EXEC} -row_expand_continuous_encrypted_matrix ${encrypted_matrix_file} ${n_rows_in_expansions} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} ${params_dir}/${PRIVATE_KEY_FILE} ${enc_row_expansion_dir}

	exit $?
fi

if [[ "${cmd_option}" == "-col_expand_plaintext_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Plaintext matrix path] [# columns per expansion] [Encrypted column expansion directory]"
		exit
	fi

	params_dir=$2
	plaintext_matrix_file=$3
	n_cols_in_expansions=$4
	enc_row_expansion_dir=$5

	${COLLAGENE_SECURE_EXEC} -col_expand_dense_encrypt_matrix ${plaintext_matrix_file} ${n_cols_in_expansions} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${enc_row_expansion_dir}

	exit $?
fi

if [[ "${cmd_option}" == "-pool_col_expanded_matrices" ]]
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Column expanded matrix directories list file] [Pooled column expanded matrix directory]"
		exit
	fi

	params_dir=$2
	col_exp_dirs_list_file=$3
	pooled_col_exp_dir=$4

	if [[ ! -f ${col_exp_dirs_list_file} ]]
	then
		echo "Could not find the list of col-expansion directories @ \"${col_exp_dirs_list_file}\""
		exit 1
	fi

	if [[ ! -d ${pooled_col_exp_dir} ]]
	then
		echo "Could not find the list of output directory @ \"${pooled_col_exp_dir}\""
		exit 1
	fi

	n_expansions=0
	exp_i=0
	n_workers=0
	while [[ 1 == 1 ]]
	do
		# Pool current expansion.
		temp_colexp_list_file=`mktemp`
		awk -v exp_i=${exp_i} {'print $1"/repcol_"exp_i".bin.enc"'} ${col_exp_dirs_list_file} > ${temp_colexp_list_file}

		cur_files=`cat ${temp_colexp_list_file}`
		found_all=1
		for cur_file in ${cur_files[@]}
		do
			if [[ ! -f ${cur_file} ]]
			then
				found_all=0
				break
			fi
		done

		if [[ ${found_all} == 0 ]]
		then
			break
		fi

		n_files=`wc -l ${temp_colexp_list_file} | awk {'print $1'}`

		if [[ ${n_expansions} == 0 ]]
		then
			n_expansions=${n_files}
		elif [[ ${n_expansions} != ${n_files} ]]
		then
			echo "Mismatching expansion fize @ expansion index ${exp_i}"
			exit 1
		fi

		echo "Pooling ${n_files} expansions for ${exp_i}.."
		$0 -secure_add_matrix_list ${params_dir} ${temp_colexp_list_file} ${pooled_col_exp_dir}/repcol_${exp_i}.bin.enc &
		
		# Update the number of workers.
		n_workers=`awk -v n_workers=${n_workers} 'BEGIN{print n_workers+1}'`

		# If we have started maximum number, wait for them to finish.
		if [[ $n_workers == ${MAX_N_WORKERS} ]]
		then
			echo "Waiting 5 col workers.."
			n_workers=0
			wait
		fi

		# Update the expansion id.
		exp_i=`echo ${exp_i} | awk {'print $1+1'}`	
	done

	# Wait for all expansions to pool.
	wait
	
	exit 0
fi

if [[ "${cmd_option}" == "-pool_row_expanded_matrices" ]]
then
    if [[ $# != 4 ]]
    then
		echo "USAGE: $0 $1 [Parameter directory] [Row expanded matrix directories list file] [Pooled row expanded matrix directory]"
		exit
    fi

    params_dir=$2
    row_exp_dirs_list_file=$3
    pooled_row_exp_dir=$4

	if [[ ! -f ${row_exp_dirs_list_file} ]]
	then
		echo "Could not find the list of row-expansion directories @ \"${row_exp_dirs_list_file}\""
		exit 1
	fi

	if [[ ! -d ${pooled_row_exp_dir} ]]
	then
		echo "Could not find the list of output directory @ \"${pooled_row_exp_dir}\""
		exit 1
	fi

    n_expansions=0
    exp_i=0
	n_workers=0
    while [[ 1 == 1 ]]
    do
		# Pool current expansion.
		temp_rowexp_list_file=`mktemp`
		awk -v exp_i=${exp_i} {'print $1"/reprow_"exp_i".bin.enc"'} ${row_exp_dirs_list_file} > ${temp_rowexp_list_file}

		cur_files=`cat ${temp_rowexp_list_file}`
		found_all=1
		for cur_file in ${cur_files[@]}
		do
		if [[ ! -f ${cur_file} ]]
		then
			found_all=0
			break
		fi
		done

		if [[ ${found_all} == 0 ]]
		then
			break
		fi

		n_files=`wc -l ${temp_rowexp_list_file} | awk {'print $1'}`

		if [[ ${n_expansions} == 0 ]]
		then
			n_expansions=${n_files}
		elif [[ ${n_expansions} != ${n_files} ]]
		then
			echo "Mismatching expansion fize @ expansion index ${exp_i}"
			exit 1
		fi

		echo "Pooling ${n_files} expansions for ${exp_i}.."
		$0 -secure_add_matrix_list ${params_dir} ${temp_rowexp_list_file} ${pooled_row_exp_dir}/reprow_${exp_i}.bin.enc &

		# Update the number of workers.
		n_workers=`awk -v n_workers=${n_workers} 'BEGIN{print n_workers+1}'`

		# If we have started maximum number, wait for them to finish.
		if [[ $n_workers == ${MAX_N_WORKERS} ]]
		then
			echo "Waiting 5 row workers.."
			n_workers=0
			wait
		fi

		# Update the expansion id.
		exp_i=`echo ${exp_i} | awk {'print $1+1'}`
    done

	# Wait for all expansions to pool.
	wait

    exit 0
fi

if [[ "${cmd_option}" == "-secure_row2row_multiply_matrices" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [matrix A file] [matrix B file] [Output matrix file]"
		exit
	fi

	params_dir=$2
	enc_matA_file=$3
	enc_matB_file=$4
	res_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -secure_row2row_inner_prod_continuous_encrypted_matrices ${enc_matA_file} ${enc_matB_file} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} NO_POOLED_KEY ${res_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-secure_elementwise_multiply_matrices" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [matrix A file] [matrix B file] [Output matrix file]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	enc_matB_file=$4
	res_matrix_file=$5

	$0 -validate_params_dir ${params_dir}
	validate_res=$?
	if [[ ${validate_res} != 0 ]]
	then
		echo "Could not validate parameters directory."
		exit 1
	fi

	if [[ ! -f ${enc_matA_file} ]]
	then
		echo "Could not find matrix ${enc_matA_file}"
		exit 1
	fi

	if [[ ! -f ${enc_matB_file} ]]
	then
		echo "Could not find matrix ${enc_matB_file}"
		exit 1
	fi

	${COLLAGENE_SECURE_EXEC} -secure_elementwise_mul_cont_ct_matrices ${enc_matA_file} ${enc_matB_file} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} ${params_dir}/${PRIVATE_KEY_FILE} ${res_matrix_file}

	exit $?
fi
if [[ "${cmd_option}" == "-secure_add_matrices" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [matrix A file] [matrix B file] [Output matrix file]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	enc_matB_file=$4
	res_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices ${enc_matA_file} ${enc_matB_file} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} ${params_dir}/${PRIVATE_KEY_FILE} ${res_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-secure_add_matrix_list" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [matrix list file] [Output matrix file]"
		exit 1
	fi

	params_dir=$2
	mat_list_file=$3
	res_matrix_file=$4

	${COLLAGENE_SECURE_EXEC} -secure_add_cont_ct_matrices_per_list ${mat_list_file} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} ${params_dir}/${PRIVATE_KEY_FILE} ${res_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-secure_sub_matrices" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix file]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	enc_matB_file=$4
	res_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -secure_sub_cont_ct_matrices ${enc_matA_file} ${enc_matB_file} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PUBLIC_KEY_FILE} ${params_dir}/${RELIN_KEY_FILE} ${params_dir}/${GALOIS_KEY_FILE} ${params_dir}/${PRIVATE_KEY_FILE} ${res_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-write_encrypted_matrix_dimensions" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [matrix A file] [Dimensions file]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	matrix_dim_file=$4

	${COLLAGENE_SECURE_EXEC} -write_enc_matrix_dimensions ${enc_matA_file} ${matrix_dim_file}

	exit $?
fi

if [[ "${cmd_option}" == "-validate_ckks_text_params" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Validity output file]"
		exit 1
	fi

	params_dir=$2
	param_valid_op_fp=$3

	# First, validate the parameters directory.
	$0 -validate_params_dir ${params_dir}
	validate_res=$?
	if [[ ${validate_res} != 0 ]]
	then
		echo "Could not validate parameters directory."
		exit 1
	fi

	${COLLAGENE_SECURE_EXEC} -validate_ckks_text_params ${params_dir}/${TEXT_PARAMS_PATH} ${param_valid_op_fp}

	exit $?
fi

if [[ "${cmd_option}" == "-write_encrypted_matrix_vital_stats" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [matrix A file] [Vital statistics file]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	matrix_vital_stats_file=$4

	if [[ ! -f ${enc_matA_file} ]]
	then
		echo "Could not find encrypted matrix file @ \"${enc_matA_file}\""
		exit 1
	fi

	$0 -validate_params_dir ${params_dir}
	validate_res=$?
	if [[ ${validate_res} != 0 ]]
	then
		echo "Could not validate parameters directory."
		exit 1
	fi

	${COLLAGENE_SECURE_EXEC} -write_continuous_encrypted_ciphertext_vital_stats ${enc_matA_file} ${params_dir}/${TEXT_PARAMS_PATH} ${matrix_vital_stats_file}

	exit $?
fi

if [[ "${cmd_option}" == "-partial_decrypt_matrix" ]] 
then
	if [[ $# != 6 ]]
	then
		echo "USAGE: $0 $1 -partial_decrypt_matrix [Parameter directory] [Encrypted matrix path] [Site index (0-based)] [Decryption noise variance (set to 0 for default value of 40-bits)] [Output matrix path]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	site_i=$4
	smdging_noise_variance_bit_size=$5
	partdec_file=$6

	if [[ ! -f ${enc_matA_file} ]]
	then
		echo "Could not find encrypted matrix file @ \"${enc_matA_file}\""
		exit 1
	fi

	$0 -validate_params_dir ${params_dir}
	validate_res=$?
	if [[ ${validate_res} != 0 ]]
	then
		echo "Could not validate parameters directory."
		exit 1
	fi

	# Set site0 flag.
	is_site0=0
	if [[ $site_i == 0 ]]
	then
		is_site0=1
	fi

	echo "Decrypting for site ${site_i} (${is_site0})"
	${COLLAGENE_SECURE_EXEC} -partial_decrypt_continuous_enc_per_noisy_secretkey ${enc_matA_file} ${is_site0} ${params_dir}/${TEXT_PARAMS_PATH} ${params_dir}/${PRIVATE_KEY_FILE} ${smdging_noise_variance_bit_size} ${partdec_file}

	exit $?
fi


if [[ "${cmd_option}" == "-pool_partially_decrypted_matrices" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 -pool_partially_decrypted_matrices [Parameter directory] [Partially decrypted matrix list file] [Output matrix path]"
		exit 1
	fi

	params_dir=$2
	partdec_files_list=$3
	fulldec_file=$4

	# First, validate the params directory.
	$0 -validate_params_dir ${params_dir}
	validate_res=$?
	if [[ ${validate_res} != 0 ]]
	then
		echo "Could not validate parameters directory."
		exit 1
	fi

	if [[ ! -f ${partdec_files_list} ]]
	then
		echo "Could not find partially decrypted files list @ \"${partdec_files_list}\""
		exit 1
	fi

	${COLLAGENE_SECURE_EXEC} -pool_partial_decrypted_continuous_enc_data ${partdec_files_list} ${params_dir}/${TEXT_PARAMS_PATH} ${fulldec_file}

	exit $?
fi

if [[ ${cmd_option} == "-symmetric_encrypt_partdec_data" ]]
then
    if [[ $# != 5 ]]
    then
            echo "USAGE: $0 $1 -symmetric_encrypt_partdec_data [Data params dir] [partdec file] [openssl hash file] [Output file]"
            exit
    fi

    params_dir=$2
    partdec_file=$3
	openssl_hash_file=$4
    output_file=$5

	openssl enc -aes-256-cbc -in ${partdec_file} -out ${output_file} -pass file:${openssl_hash_file}
	exit $?
fi

if [[ ${cmd_option} == "-symmetric_decrypt_partdec_data" ]]
then
    if [[ $# != 5 ]]
    then
            echo "USAGE: $0 $1 -symmetric_decrypt_partdec_data [Data params dir] [encrypted partdec file] [openssl hash file] [Output file]"
            exit
    fi

    params_dir=$2
    enc_partdec_file=$3
	openssl_hash_file=$4
    output_file=$5

	openssl enc -d -aes-256-cbc -in ${enc_partdec_file} -out ${output_file} -pass file:${openssl_hash_file}
	exit $?
fi


if [[ "${cmd_option}" == "-generate_plaintext_mask_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Encryptd matrix file] [Mask noise standard deviation] [Output plaintext matrix file]"
		exit 1
	fi

	params_dir=$2
	matA_file=$3
	mask_sigma=$4
	plaintext_mask_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -generate_plaintext_mask_per_continuous_encrypted_data ${matA_file} ${mask_sigma} ${plaintext_mask_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-mask_encrypted_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Encrypted matrix file] [Encrypted mask matrix file] [Output plaintext matrix file]"
		exit 1
	fi

	params_dir=$2
	matA_file=$3
	enc_mask_matrix_file=$4
	masked_enc_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -additive_mask_continuous_encrypted_data ${matA_file} ${enc_mask_matrix_file} ${params_dir}/${TEXT_PARAMS_PATH} ${masked_enc_matrix_file}

	exit $?
fi

if [[ "${cmd_option}" == "-unmask_encrypted_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Encrypted matrix path] [Encrypted mask matrix path] [Output matrix path]"
		exit 1
	fi

	params_dir=$2
	enc_matA_file=$3
	enc_mask_matrix_file=$4
	masked_enc_matrix_file=$5

	${COLLAGENE_SECURE_EXEC} -additive_unmask_continuous_encrypted_data ${enc_matA_file} ${enc_mask_matrix_file} ${params_dir}/${TEXT_PARAMS_PATH} ${masked_enc_matrix_file}

	exit $?
fi

###########################################################################
# Following are the plaintext matrix interface.
###########################################################################
if [[ "${cmd_option}" == "-generate_random_plaintext_matrix" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [# rows] [# cols] [Output plaintext matrix file]"
		exit 1
	fi

	n_rows=$2
	n_cols=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -generate_mult_full_noise_matrix ${n_rows} ${n_cols} ${op_pt_matrix_fp}

	exit $?
fi

if [[ "${cmd_option}" == "-generate_diagonal_plaintext_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [# rows] [# columns] [Scalar value] [Output matrix path]"
		exit 1
	fi

	n_rows=$2
	n_cols=$3
	diag_val=$4
	op_pt_matrix_fp=$5

	${COLLAGENE_SECURE_EXEC} -generate_constant_diagonal_matrix ${n_rows} ${n_cols} ${diag_val} ${op_pt_matrix_fp}

	exit $?
fi

if [[ "${cmd_option}" == "-generate_constant_plaintext_matrix" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [# rows] [# columns] [Scalar value] [Output matrix path]"
		exit 1
	fi

	n_rows=$2
	n_cols=$3
	init_val=$4
	op_pt_matrix_fp=$5

	${COLLAGENE_SECURE_EXEC} -generate_constant_full_matrix ${n_rows} ${n_cols} ${init_val} ${op_pt_matrix_fp}

	exit $?
fi

if [[ "${cmd_option}" == "-scalar_multiply_plaintext_matrix" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Plaintext matrix A file] [Scalar value] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	scalar_multiplier_value=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -scalar_multiply_matrix_plain ${plaintext_matA_file} ${scalar_multiplier_value} ${op_pt_matrix_fp}

	exit $?
fi

if [[ "${cmd_option}" == "-invert_plaintext_matrix" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Matrix A file] [Inverted plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	op_pt_matrix_fp=$3

	${COLLAGENE_SECURE_EXEC} -plain_invert_matrix ${plaintext_matA_file} ${op_pt_matrix_fp}

	exit $?
fi

if [[ "${cmd_option}" == "-transpose_plaintext_matrix" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Transposed plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	op_pt_matrix_fp=$3

	${COLLAGENE_SECURE_EXEC} -transpose_pt_matrix ${plaintext_matA_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-multiply_plaintext_matrix" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Matrix B file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	plaintext_matB_file=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -multiply_matrices_pt ${plaintext_matA_file} ${plaintext_matB_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-multiply_elementwise_plaintext_matrices" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Matrix B file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	plaintext_matB_file=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -multiply_matrices_elementwise_pt ${plaintext_matA_file} ${plaintext_matB_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-row2row_multiply_plaintext_matrices" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 -row2row_multiply_plaintext_matrices [Matrix A file] [Matrix B file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	plaintext_matB_file=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -row2row_multiply_pt ${plaintext_matA_file} ${plaintext_matB_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-sub_plaintext_matrix" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Matrix B file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	plaintext_matB_file=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -plain_sub_matrices ${plaintext_matA_file} ${plaintext_matB_file} ${op_pt_matrix_fp}
	
	exit $?
fi



if [[ "${cmd_option}" == "-add_plaintext_matrix" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Matrix B file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	plaintext_matB_file=$3
	op_pt_matrix_fp=$4

	${COLLAGENE_SECURE_EXEC} -plain_add_matrices ${plaintext_matA_file} ${plaintext_matB_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-add_plaintext_matrix_list" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Plaintext matrix list file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matrix_list_file=$2
	op_pt_matrix_fp=$3

	${COLLAGENE_SECURE_EXEC} -plain_add_matrices_per_list ${plaintext_matrix_list_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-pad_plaintext_matrix_2_to_n" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	op_pt_matrix_fp=$3

	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_to_next_power_of_2 ${plaintext_matA_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-pad_plaintext_matrix_row_2_to_n" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	op_pt_matrix_fp=$3

	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_rows_to_next_power_of_2 ${plaintext_matA_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-pad_plaintext_matrix_col_2_to_n" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	op_pt_matrix_fp=$3

	${COLLAGENE_SECURE_EXEC} -plain_pad_matrix_cols_to_next_power_of_2 ${plaintext_matA_file} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-pad_plaintext_matrix_custom" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [# rows] [# cols] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	n_rows=$3
	n_cols=$4
	op_pt_matrix_fp=$5

	${COLLAGENE_SECURE_EXEC} -plain_unpad_matrix_to_size ${plaintext_matA_file} ${n_rows} ${n_cols} ${op_pt_matrix_fp}
	
	exit $?
fi

if [[ "${cmd_option}" == "-save_matrix_text" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	op_file=$3
	
	${COLLAGENE_SECURE_EXEC} -dump_matrix_plain ${plaintext_matA_file} ${op_file}

	exit $?
fi

if [[ "${cmd_option}" == "-transform_plaintext_elementwise_per_callback" ]] 
then
	if [[ $# != 4 ]]
	then
		echo "USAGE: $0 $1 [Plaintext matrix path] [Function: \"log\"/\"exp\"/\"sigmoid\"] [Output matrix path]"
		exit 1
	fi

	plaintext_matA_file=$2
	function_type=$3
	op_file=$4
	
	${COLLAGENE_SECURE_EXEC} -transform_elementwise_per_callback_pt ${plaintext_matA_file} ${function_type} ${op_file}

	exit $?
fi

if [[ "${cmd_option}" == "-convert_plaintext_matrix_text_2_bin" ]] 
then
	if [[ $# != 5 ]]
	then
		echo "USAGE: $0 $1 [Matrix A file] [Row ids in input (first col.)? (0/1)] [Col. id's in input (first row)? (0/1)] [Output plaintext matrix file]"
		exit 1
	fi

	plaintext_matA_file=$2
	row_ids_in_input=$3
	col_ids_in_input=$4
	op_file=$5
	
	${COLLAGENE_SECURE_EXEC} -save_matrix_plain_2_bin ${plaintext_matA_file} ${row_ids_in_input} ${col_ids_in_input} ${op_file}

	exit $?
fi

if [[ "${cmd_option}" == "-write_plaintext_matrix_dimensions" ]] 
then
	if [[ $# != 3 ]]
	then
		echo "USAGE: $0 $1 [Parameter directory] [Matrix A file] [Output dimensions file]"
		exit 1
	fi

	plaintext_matA_file=$2
	plain_mat_dims_file=$3

	${COLLAGENE_SECURE_EXEC} -write_plain_matrix_dimensions ${plaintext_matA_file} ${plain_mat_dims_file}

	exit $?
fi

