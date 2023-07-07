# All sites perform the same set of steps in this protocol except for the different site indices below.
if [[ $# != 1 ]]
then
	echo "Networked partial decryption example using COLLAGENE and file-based communication example."
	echo "USAGE: $0 [site index for the current running site]"
	exit 1
fi

# i_site is always necessary for COLLAGENE. It must be 0-based counting from 0 to N_sites-1. 
i_site=$1
N_sites=3

if [[ ${i_site} -ge ${N_sites} ]]
then
	echo "Cannot have ${i_site} >= ${N_sites}"
	exit 1
fi

if [[ ! -f "./COLLAGENE.sh" ]]
then
	echo "Could not find COLLAGENE script, copy it from source code scripts/ directory."
	exit 1
fi

if [[ ! -f "./FILE_IO_UTILITIES.sh" ]]
then
	echo "Could not find FILE_IO_UTILITIES.sh script, copy it from source code scripts/ directory."
	exit 1
fi

# We need N_sites_min_one to run the loops below; it is only necessary since all sites are running in the same directory and site keys directories are named SITE_0, SITE_1, SITE_2.
N_sites_min_one=`echo ${N_sites} | awk {'print $1-1'}`

# The random matrix can represent an intermediate result (e.g., masked feature correlation matrix) that will be used by sites to perform an operation.
echo "SITE${i_site}: Generating matrix and encrypting.."
./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 matrix_from_site${i_site}.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} matrix_from_site${i_site}.bin matrix_from_site${i_site}.bin.enc

# Site uploads the matrix to be decrypted so other sites can download and partially decrypt it.
echo "SITE${i_site}: Uploading matrix.."
echo "matrix_from_site${i_site}.bin.enc" > files.list 
./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params files.list 

if [[ $? == 1 ]]
then
	echo "Upload failed."
	exit 1
fi

# Partially decrypt all sites:
for ((i_source_site=0; i_source_site<${N_sites}; i_source_site++))
do
	# Site partially decrypts and uploads the data to shared directory.
	enc_matrix_file="matrix_from_site${i_source_site}.bin.enc"

	# Do not download if the same site's data.
	echo "SITE${i_site}: Downloading the matrix of ${i_source_site}"
	echo "${enc_matrix_file}" > files.list

	./FILE_IO_UTILITIES.sh -wait_for_files_in_shared data_config.params files.list 
	./FILE_IO_UTILITIES.sh -download_files_from_shared data_config.params files.list . 

	if [[ $? == 1 ]]
	then
		echo "SITE${i_site}: Could not download matrix_from_site${i_source_site}.enc"
		exit 
	fi

	if [[ ! -f ${enc_matrix_file} ]]
	then
		echo "Error downloading ${enc_matrix_file}"
		exit 1
	fi

	# Do partial decryption on this site's matrix.
	./COLLAGENE.sh -partial_decrypt_matrix SITE_${i_site} ${enc_matrix_file} ${i_site} 40 ${enc_matrix_file}_partdec_by_${i_site}.partdec
	./COLLAGENE.sh -symmetric_encrypt_partdec_data SITE_${i_site} ${enc_matrix_file}_partdec_by_${i_site}.partdec SITE_${i_site}/partdec_data_enc_hash.symmetric_key ${enc_matrix_file}_partdec_by_${i_site}.partdec.enc

	# Upload the partial decryption.
	# Note that we do not upload the partial decryption if 
	echo "matrix_from_site${i_source_site}.bin.enc_partdec_by_${i_site}.partdec.enc" > files.list
	./FILE_IO_UTILITIES.sh -upload_files_to_shared data_config.params files.list
done

# At this point, site download the partial decryptions for its matrix:
for ((i_decrypting_site=0; i_decrypting_site<${N_sites}; i_decrypting_site++))
do
	# The name of the current partdec file for this site.
	partdec_file="matrix_from_site${i_site}.bin.enc_partdec_by_${i_decrypting_site}.partdec"

	# Wait the file to be ready in the shared folder.
	echo "${partdec_file}.enc" > files.list
	./FILE_IO_UTILITIES.sh -wait_for_files_in_shared data_config.params files.list

	# Download the file.
	./FILE_IO_UTILITIES.sh -download_files_from_shared data_config.params files.list .
	
	# Error check.
	if [[ $? == 1 ]]
	then
		echo "SITE${i_site}: Could not download partial decryptions from other sites."
		exit 1
	fi

	# Remove symmetric encryption from the sites using partdec encryption key.
	./COLLAGENE.sh -symmetric_decrypt_partdec_data SITE_${i_site} ${partdec_file}.enc SITE_${i_site}/partdec_data_enc_hash.symmetric_key ${partdec_file}
done

# Pool partial decryptions:
echo "SITE${i_site}: Aggregating partial decryptions.."
ls matrix_from_site${i_site}.bin.enc_partdec_by_*.partdec > site_${i_site}.partdecs
./COLLAGENE.sh -pool_partially_decrypted_matrices SITE_${i_site} site_${i_site}.partdecs full_decrypted_matrix_from_site${i_site}.bin

# Write the matrix to text.
./COLLAGENE.sh -save_matrix_text full_decrypted_matrix_from_site${i_site}.bin full_decrypted_matrix_from_site${i_site}.bin.txt
./COLLAGENE.sh -save_matrix_text matrix_from_site${i_site}.bin matrix_from_site${i_site}.bin.txt

echo "Original matrix:"
cat matrix_from_site${i_site}.bin.txt

echo "Decrypted matrix:"
cat full_decrypted_matrix_from_site${i_site}.bin.txt
