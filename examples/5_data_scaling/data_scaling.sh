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

# This is the data scaling exercise of COLLAGENE. The basic idea is to demonstrate that when the data range is large across two matrices, scaling can help 
# alleviate error accumulation in multiplication of padded matrices.

# This is the scaler value that we use to scale matrix2 below.
# It serves to match the data scales between multiplied matrices.
# 
SCALER=1

# of sites, used only for collaborative decryption in this case.
N_SITES=3

# Generate random matrices of differnt scales and showcase the importance of data scaling especially for padded matrices. 
mat_A_nrows=3
mat_A_ncols=15

mat_B_nrows=15
mat_B_ncols=3

# Generate a matrix with large values.
MATRIX1_DATA_RANGE=10000
./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix1.bin
./COLLAGENE.sh -scalar_multiply_plaintext_matrix random_matrix1.bin ${MATRIX1_DATA_RANGE} random_matrix1.bin
./COLLAGENE.sh -pad_plaintext_matrix_2_to_n random_matrix1.bin random_matrix1_padded.bin

# Generate a second matrix with relatively small values.
./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_B_nrows} ${mat_B_ncols} random_matrix2.bin
cp -r random_matrix2.bin random_matrix2.bin_orig.bin

# Now scale matrix2, this will be used to test the effect of scaling.
./COLLAGENE.sh -scalar_multiply_plaintext_matrix random_matrix2.bin ${SCALER} random_matrix2_scaled.bin
./COLLAGENE.sh -pad_plaintext_matrix_2_to_n random_matrix2_scaled.bin random_matrix2_scaled_padded.bin

# Multiply them to shwocase the impact of added zeros.
./COLLAGENE.sh -write_plaintext_matrix_dimensions random_matrix1_padded.bin matrix1_dims.txt
./COLLAGENE.sh -write_plaintext_matrix_dimensions random_matrix2_scaled_padded.bin matrix2_dims.txt
padded_res_n_rows=`cat matrix1_dims.txt | awk {'print $1'}`
padded_res_n_cols=`cat matrix2_dims.txt | awk {'print $2'}`

# Column expand the first matrix as we did before.
# Following row expands the padded second matrix.
rm -f -r padded_colexp_dir
mkdir padded_colexp_dir
./COLLAGENE.sh -col_expand_plaintext_matrix SITE_0 random_matrix1_padded.bin ${padded_res_n_cols} padded_colexp_dir

rm -f -r padded_rowexp_dir
mkdir padded_rowexp_dir
./COLLAGENE.sh -row_expand_plaintext_matrix SITE_0 random_matrix2_scaled_padded.bin ${padded_res_n_rows} padded_rowexp_dir

./COLLAGENE.sh -secure_multiply_rowcol_expansion SITE_0 padded_colexp_dir padded_rowexp_dir ct_mat_mul_random_matrix12.bin.enc

# Decrypt then un-scale the result.
UNSCALER=`echo ${SCALER} | awk {'print 1/$1'}`
./collaborative_decrypt_matrix.sh ct_mat_mul_random_matrix12.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -scalar_multiply_plaintext_matrix ct_mat_mul_random_matrix12.bin.enc_collaborative_dec.bin ${UNSCALER} ct_mat_mul_random_matrix12.bin.enc_collaborative_dec.bin_rescaled.bin
./COLLAGENE.sh -save_matrix_text ct_mat_mul_random_matrix12.bin.enc_collaborative_dec.bin_rescaled.bin ct_mat_mul_random_matrix12.bin.enc_collaborative_dec.bin_rescaled.bin.txt

# Now test the plaintext values with each other.
./COLLAGENE.sh -multiply_plaintext_matrix random_matrix1.bin random_matrix2.bin pt_mat_mul_randon_matrix12.bin
./COLLAGENE.sh -save_matrix_text pt_mat_mul_randon_matrix12.bin pt_mat_mul_randon_matrix12.bin.txt

echo "You can compare ct_mat_mul_random_matrix12.bin.enc_collaborative_dec.bin_rescaled.bin.txt and pt_mat_mul_randon_matrix12.bin.txt."

echo "The error accumulation in the padded portion of ct_mat_mul_random_matrix12.bin.enc_collaborative_dec.bin_rescaled.bin.txt is important when padded matrices are processed in a pipeline."

echo "Note that if SCALER parameter above is increased, the noise in the padded portion of product matrix should decrease."

