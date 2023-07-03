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

echo "Using ${N_SITES} sites.."
N_SITES_MIN_ONE=`echo ${N_SITES} | awk {'print $1-1'}`
site_iters=`seq 0 ${N_SITES_MIN_ONE}`

i_site=0

################################################################################################################################
echo "TESTING SUMMATIONS.."

./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix1.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix1.bin random_matrix1.bin.enc

./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix2.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix2.bin random_matrix2.bin.enc

# Add/sub matrices.
./COLLAGENE.sh -secure_add_matrices SITE_${i_site} random_matrix1.bin.enc random_matrix2.bin.enc ct_added_random_matrix12.bin.enc
./COLLAGENE.sh -add_plaintext_matrix random_matrix1.bin random_matrix2.bin pt_added_random_matrix12.bin

echo "Decrypting summation."
./collaborative_decrypt_matrix.sh ct_added_random_matrix12.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -save_matrix_text pt_added_random_matrix12.bin pt_added_random_matrix12.bin.txt

# Element-element multiplication:
./COLLAGENE.sh -secure_elementwise_multiply_matrices SITE_${i_site} random_matrix1.bin.enc random_matrix2.bin.enc ct_elm_mul_random_matrix12.bin.enc
./COLLAGENE.sh -multiply_elementwise_plaintext_matrices random_matrix1.bin random_matrix2.bin pt_elm_mul_random_matrix12.bin

echo "Decrypting elementwise multiplications."
./collaborative_decrypt_matrix.sh ct_elm_mul_random_matrix12.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -save_matrix_text pt_elm_mul_random_matrix12.bin pt_elm_mul_random_matrix12.bin.txt

# Add a list of encrypted matrices:
./COLLAGENE.sh -generate_random_plaintext_matrix 10 10 random_matrix3.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix3.bin random_matrix3.bin.enc
echo "random_matrix1.bin.enc" > ct_matrices.list
echo "random_matrix2.bin.enc" >> ct_matrices.list
echo "random_matrix3.bin.enc" >> ct_matrices.list
./COLLAGENE.sh -secure_add_matrix_list SITE_${i_site} ct_matrices.list ct_list_add_matrices.bin.enc

echo "random_matrix1.bin" > pt_matrices.list
echo "random_matrix2.bin" >> pt_matrices.list
echo "random_matrix3.bin" >> pt_matrices.list
./COLLAGENE.sh -add_plaintext_matrix_list pt_matrices.list pt_list_add_matrices.bin

echo "Decrypting matrix list additions."
./collaborative_decrypt_matrix.sh ct_list_add_matrices.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -save_matrix_text pt_list_add_matrices.bin pt_list_add_matrices.bin.txt

################################################################################################################################
echo "TESTING ROW/COL EXPANSIONS AND MATRIX MULTIPLICATIONS.."

# Row and col expand the matrices: We generate two matrices to be multiplied and save them.
mat_A_nrows=10
mat_A_ncols=30
mat_B_nrows=${mat_A_ncols}  # necessary for valid multiplication.
mat_B_ncols=20
./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix1.bin
./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_B_nrows} ${mat_B_ncols} random_matrix2.bin

# Col-expand first matrix and row expand second one.
# col expansion of the first matrix generates resulting size matrices.
res_n_cols=${mat_B_ncols}
res_n_rows=${mat_A_nrows}
rm -f -r colexp_dir
mkdir colexp_dir
./COLLAGENE.sh -col_expand_plaintext_matrix SITE_0 random_matrix1.bin ${res_n_cols} colexp_dir
rm -f -r rowexp_dir
mkdir rowexp_dir
./COLLAGENE.sh -row_expand_plaintext_matrix SITE_0 random_matrix2.bin ${res_n_rows} rowexp_dir

# Now we can multiply the expanded matrices.
./COLLAGENE.sh -secure_multiply_rowcol_expansion SITE_0 colexp_dir rowexp_dir ct_mat_mul_random_matrix12.bin.enc

# Do plaintext multiplication to compare.
./COLLAGENE.sh -multiply_plaintext_matrix random_matrix1.bin random_matrix2.bin pt_mat_mul_random_matrix12.bin

# Decrypt and test. 
./collaborative_decrypt_matrix.sh ct_mat_mul_random_matrix12.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -save_matrix_text pt_mat_mul_random_matrix12.bin pt_mat_mul_random_matrix12.bin.txt

#################################################################################################################################
# Pad the matrices.
echo "TESTING MATRIX PADDING"
mat_A_nrows=10
mat_A_ncols=30
mat_B_nrows=${mat_A_ncols}  # necessary for valid multiplication.
mat_B_ncols=20

./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix1.bin
./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_B_nrows} ${mat_B_ncols} random_matrix2.bin

# We just pad the first matrix and column expand it.
./COLLAGENE.sh -pad_plaintext_matrix_2_to_n random_matrix1.bin random_matrix1.bin_padded.bin

# We pad the second matrix and encrypt it. The idea is to test the row expansion in encrypted domain.
./COLLAGENE.sh -pad_plaintext_matrix_2_to_n random_matrix2.bin random_matrix2.bin_padded.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix2.bin_padded.bin random_matrix2.bin_padded.bin.enc

# Get the padded matrix dimensions
./COLLAGENE.sh -write_plaintext_matrix_dimensions random_matrix1.bin_padded.bin matrix1_dims.txt
./COLLAGENE.sh -write_encrypted_matrix_dimensions SITE_0 random_matrix2.bin_padded.bin.enc matrix2_dims.txt

padded_res_n_rows=`cat matrix1_dims.txt | awk {'print $1'}`
padded_res_n_cols=`cat matrix2_dims.txt | awk {'print $2'}`

# Column expand the first matrix as we did before.
rm -f -r padded_colexp_dir
mkdir padded_colexp_dir
./COLLAGENE.sh -col_expand_plaintext_matrix SITE_0 random_matrix1.bin_padded.bin ${padded_res_n_cols} padded_colexp_dir

# Following row expands the padded second matrix.
rm -f -r padded_rowexp_dir
mkdir padded_rowexp_dir
./COLLAGENE.sh -secure_row_expand_encrypted_matrix SITE_${i_site} random_matrix2.bin_padded.bin.enc ${padded_res_n_rows} padded_rowexp_dir

# Perform multiplication as we did before.
./COLLAGENE.sh -secure_multiply_rowcol_expansion SITE_0 padded_colexp_dir padded_rowexp_dir ct_mat_mul_random_matrix12.bin.enc

# Do plaintext multiplication to compare, note that this multiplication does not use padded matrices, it will return a smaller result below.
./COLLAGENE.sh -multiply_plaintext_matrix random_matrix1.bin random_matrix2.bin pt_mat_mul_random_matrix12.bin

# Decrypt and compare: Note that the encrypted matrix is padded, which means there will be a large number of 0's at the ends of the columns and rows.
./collaborative_decrypt_matrix.sh ct_mat_mul_random_matrix12.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -save_matrix_text pt_mat_mul_random_matrix12.bin pt_mat_mul_random_matrix12.bin.txt

#################################################################################################################################
# Do row-2-row multiplication: This makes vector-2-matrix multiplications faster.
mat_A_nrows=10
mat_A_ncols=30

./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix1.bin
./COLLAGENE.sh -generate_random_plaintext_matrix ${mat_A_nrows} ${mat_A_ncols} random_matrix2.bin

# We need to pad the columns before row-2-row multiplication.
./COLLAGENE.sh -pad_plaintext_matrix_col_2_to_n random_matrix1.bin random_matrix1.bin_col_padded.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix1.bin_col_padded.bin random_matrix1.bin_col_padded.bin.enc

./COLLAGENE.sh -pad_plaintext_matrix_col_2_to_n random_matrix2.bin random_matrix2.bin_col_padded.bin
./COLLAGENE.sh -encrypt_plaintext_matrix SITE_${i_site} random_matrix2.bin_col_padded.bin random_matrix2.bin_col_padded.bin.enc

# Note that we are wasting a lot of space in this small matrix as this operation processes a ciphertext that is much larger than the matrix size. This is ok here for demonstration purposes but the parameters can be optimized if this is the only operation that is being performed.
./COLLAGENE.sh -secure_row2row_multiply_matrices SITE_${i_site} random_matrix1.bin_col_padded.bin.enc random_matrix2.bin_col_padded.bin.enc ct_row2row_mul_random_matrix12.bin.enc
./COLLAGENE.sh -row2row_multiply_plaintext_matrices random_matrix1.bin random_matrix2.bin pt_row2row_mul_random_matrix12.bin

./collaborative_decrypt_matrix.sh ct_row2row_mul_random_matrix12.bin.enc ${N_SITES} >& /dev/null
./COLLAGENE.sh -save_matrix_text pt_row2row_mul_random_matrix12.bin pt_row2row_mul_random_matrix12.bin.txt



