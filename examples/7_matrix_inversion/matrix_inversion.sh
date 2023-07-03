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

rm -f *.list
rm -f *.enc
rm -f *.partdec
rm -f *.bin 

# In this exercise, the sites generate random matrices, row/col expand them and encrypt them. These are shared with the other sites.
# The sites download the row/col expansions of other sites and pool them locally.
# The sites generate multiplicative noise matrices, encrypt their expansions and exchange with other sites.
# Finally, the sites add the multiplicative noise to the pooled matrices by row/col multiplication.
# The sites locally invert the matrix, do row/col expansion to it, and multiply with the noise matrix on the right to obtain final inverted encrypted matrix.
# The results are decrypted and compared.
# Note that this example writes a lot of messages on the screen, which can be suppressed by piping into /dev/null or to a file.
mat_A_nrows=20
mat_A_ncols=20

n_sites_min_one=`echo ${N_SITES} | awk {'print $1-1}'`
site_iters=`seq 0 ${n_sites_min_one}`

rm -f site_rowexp_dirs.list
rm -f site_colexp_dirs.list
rm -f data_matrices.list
for i_site in ${site_iters[@]}
do
	# Generate a random rescaled matrix.
	./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix_${i_site}.bin

	./COLLAGENE.sh -scalar_multiply_plaintext_matrix random_matrix_${i_site}.bin 5 random_matrix_${i_site}.bin_scaled.bin

	echo random_matrix_${i_site}.bin_scaled.bin >> data_matrices.list

        # Row/col expand the matrix.
        rm -f -r site_${i_site}_data_colexp
        rm -f -r site_${i_site}_data_rowexp

        # Following row expands the padded second matrix.
        rm -f -r site_${i_site}_data_rowexp
        mkdir site_${i_site}_data_rowexp
        ./COLLAGENE.sh -row_expand_plaintext_matrix SITE_${i_site} random_matrix_${i_site}.bin_scaled.bin ${mat_A_nrows} site_${i_site}_data_rowexp &

        rm -f -r site_${i_site}_data_colexp
        mkdir site_${i_site}_data_colexp
        ./COLLAGENE.sh -col_expand_plaintext_matrix SITE_${i_site} random_matrix_${i_site}.bin_scaled.bin ${mat_A_ncols} site_${i_site}_data_colexp &

	echo site_${i_site}_data_rowexp >> site_rowexp_dirs.list
	echo site_${i_site}_data_colexp >> site_colexp_dirs.list

	wait
done

# This is the matrix that we are protecting, we save it here for comparison of final results.
./COLLAGENE.sh -add_plaintext_matrix_list data_matrices.list pooled_pt_data_matrices.bin
./COLLAGENE.sh -save_matrix_text pooled_pt_data_matrices.bin pooled_pt_data_matrices.bin.txt

echo "POOLING COL/ROW EXPANSIONS."
rm -f -r pooled_data_colexp_dir
mkdir pooled_data_colexp_dir
./COLLAGENE.sh -pool_col_expanded_matrices SITE_0 site_colexp_dirs.list pooled_data_colexp_dir >& /dev/null

rm -f -r pooled_data_rowexp_dir
mkdir pooled_data_rowexp_dir
./COLLAGENE.sh -pool_row_expanded_matrices SITE_0 site_rowexp_dirs.list pooled_data_rowexp_dir >& /dev/null

res_n_rows=${mat_A_nrows}
res_n_cols=${mat_A_ncols}

# We would like to invert this matrix using multiplicative noise from all sites.
rm -f col_exp_dirs_list.txt
rm -f row_exp_dirs_list.txt
for i_site in ${site_iters[@]}
do
	# Generate a masking matrix from this site, this is simple a random matrix with a very large variance.
	./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} site_${i_site}_mask_matrix.bin

	# Row/col expand the matrix.
	rm -f -r site_${i_site}_mask_colexp
	rm -f -r site_${i_site}_mask_rowexp

	# Following row expands the padded second matrix.
	rm -f -r site_${i_site}_mask_rowexp
	mkdir site_${i_site}_mask_rowexp
	./COLLAGENE.sh -row_expand_plaintext_matrix SITE_${i_site} site_${i_site}_mask_matrix.bin ${mat_A_nrows} site_${i_site}_mask_rowexp &

        rm -f -r site_${i_site}_mask_colexp
        mkdir site_${i_site}_mask_colexp
        ./COLLAGENE.sh -col_expand_plaintext_matrix SITE_${i_site} site_${i_site}_mask_matrix.bin ${mat_A_ncols} site_${i_site}_mask_colexp &

	# Wait for the expansions to be generated.
	wait

        echo site_${i_site}_mask_matrix.bin.enc >> enc_mask_matrix_list.txt
	echo site_${i_site}_mask_colexp >> col_exp_dirs_list.txt
	echo site_${i_site}_mask_rowexp >> row_exp_dirs_list.txt
done

echo "POOLING COL/ROW EXPANSIONS OF NOISE MATRIX"
rm -f -r pooled_noise_colexp
mkdir pooled_noise_colexp 
./COLLAGENE.sh -pool_col_expanded_matrices SITE_0 col_exp_dirs_list.txt pooled_noise_colexp >& /dev/null

rm -f -r pooled_noise_rowexp
mkdir pooled_noise_rowexp
./COLLAGENE.sh -pool_row_expanded_matrices SITE_0 row_exp_dirs_list.txt pooled_noise_rowexp >& /dev/null

# Multiplicatively add the noise on the left side.
./COLLAGENE.sh -secure_multiply_rowcol_expansion SITE_0 pooled_data_colexp_dir pooled_noise_rowexp mult_noise_data.enc

# Sites decrypt the noisy matrices.
./collaborative_decrypt_matrix.sh mult_noise_data.enc ${N_SITES}

echo "INVERTING MULTIPLICATIVELY MASKED DATA MATRIX.."
./COLLAGENE.sh -invert_plaintext_matrix mult_noise_data.enc_collaborative_dec.bin mult_noise_data.enc_collaborative_dec.bin_inv.bin

echo "ROW EXPANDING MASKED INVERTED MATRIX.."
rm -f -r inv_matrix_rowexp_dir
mkdir inv_matrix_rowexp_dir
./COLLAGENE.sh -row_expand_plaintext_matrix SITE_${i_site} mult_noise_data.enc_collaborative_dec.bin_inv.bin ${mat_A_nrows} inv_matrix_rowexp_dir

# Finally, multiply with the noise on the right side to remove the noise in secure domain.
./COLLAGENE.sh -secure_multiply_rowcol_expansion SITE_0 pooled_noise_colexp inv_matrix_rowexp_dir unmasked_inv_matrix_data.enc

# Decrypt to make sure what we get from secure inversion is the same as plain inversion.
./collaborative_decrypt_matrix.sh unmasked_inv_matrix_data.enc ${N_SITES}
./COLLAGENE.sh -invert_plaintext_matrix pooled_pt_data_matrices.bin inv_pooled_pt_data_matrices.bin
./COLLAGENE.sh -save_matrix_text inv_pooled_pt_data_matrices.bin inv_pooled_pt_data_matrices.bin.txt

echo "You can compare inverted matrices are inv_pooled_pt_data_matrices.bin.txt and unmasked_inv_matrix_data.enc_collaborative_dec.bin.txt"

echo "You can also compare the original pooled matrix (pooled_pt_data_matrices.bin.txt) and multiplicatively masked pooled matrix (mult_noise_data.enc_collaborative_dec.bin.txt)."


