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

# Check the validity of the original parameters:
./COLLAGENE.sh -validate_ckks_text_params SITE_0 is_valid.txt

# Update the ckks.params file and check the validity of parameter:
# Use following prime sizes: 60 30 30 30 30 30 60
# Don't forget to replace scale (3rd line in SITE_0/ckks.params) to 30 
./COLLAGENE.sh -validate_ckks_text_params SITE_0 is_valid.txt

# Use following prime sizes: 60 20 20 20 20 20 60
# This should be invalid.
./COLLAGENE.sh -validate_ckks_text_params SITE_0 is_valid.txt

# Write vitals of a ciphertext. We use a ciphertext that was generated in the previous example. 
# The vitals indicate how many operations we can perform on a ciphertext.
# The protocols should keep track of the vitals to ensure that the protocol can be correctly executed.
./COLLAGENE.sh -write_encrypted_matrix_vital_stats SITE_0 random_matrix_0.bin.enc vital_stats.txt
