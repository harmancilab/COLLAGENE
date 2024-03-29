# Command Line Interface for COLLAGENE
The command line options of COLLAGENE and network file I/O are described here.

## COLLAGENE.sh:

##### >> COLLAGENE.sh -generate_DSK_encryption_key [Site Index (0-based)]
This option generates a DSK encryption key to be submitted to KeyMaker. 

*__Inputs:__* 0-based index for this site. Each site should have a unique site index.

*__Outputs:__* A private/public key pair is written with names *site_0.dsk_enc_secret_key* and *site_0.dsk_enc_public_key*. 

*__Remarks:__* This option uses openssl to generate DSK encryption keys. A passphrase will be prompted by openssl three times while key pair is generated. The public key file is used for encrypting the DSK's by KeyMaker. Private key is used by the client to decrypt its DSK.

***

##### >> COLLAGENE.sh -download_DSK_archive [KeyMaker folder identifier] [Key archive tar file name]
This option probes the completion of DSK key generation directly from *KeyMaker* and downloads the tar archive if the key generation task is completed.

*__Inputs:__* Folder identifier in *KeyMaker*. This is the identifier that is assigned while DSK encryption archive is uploaded to *KeyMaker*. It is shown at the top of the page.

*__Outputs:__* The CKKS DSK archive (tar formatted) is saved if the key generation is completed. Otherwise a message is printed and option returns.

*__Remarks:__* This option is currently experimental and relies on an external library. This option directly interacts with *KeyMaker* and exposes site's IP to *KeyMaker*. If regulations do not allow exposing the IP address, this option should not be used to obtain the key archive. Users should obtain the link to download the DSK archive directly from *KeyMaker* website.

***

##### >> COLLAGENE.sh -decrypt_site_DSK [Site Index (0-based)] [Encrypted DSK keys directory]
This option decrypts a DSK file for the site after the keys are generated by KeyMaker. 

*__Inputs:__* 0-based index for this site and the encrypted DSK keys directory received from KeyMaker.

*__Outputs:__* The secret key is decrypted and stored inside encrypted DSK keys directory with names *site_0.secret_key*.

*__Remarks:__* After the keys are generated, they are downloaded by each site. Each site must first locally uncompress the keys file using '*tar -xvf*'. After decompressing, the sites can use this option to decrypt their own DSK. The passphrase that is used while generating DSK encryption key should be used while decrypting DSKs.

Note that the uncompressed directory contains one set of pooled keys (*pooled.galois_keys*, *pooled.relin_keys*, *pooled.public_key*), a partial data encryption key (symmetric) and a DSK for each site. After decryption, the partial data encryption key (*partdec_data_enc_hash.symmetric_key*) and DSKs are decrypted.

***

##### >> COLLAGENE.sh -generate_openssl_hash [Site Index (0-based)]
This option can be used to generate a symmetric hash for encrypting/decrypting datasets. 

*__Inputs:__* 0-based index for the site.

*__Outputs:__* A stream of random bytes will be stored in a file named, e.g., *SITE_0.hash*

*__Remarks:__* This option is not explicitly requires since KeyMaker generates a common partial data encryption key to be used at each site. The sites, however, may wish to generate a separate symmetric key for other purposes.

***

###  Initialization:
##### >> COLLAGENE.sh -set_params [Text parameters path] [Public key path] [Relin. key path] [Galois key path] [Site's secret key share path] [Parameter directory]
This option copies the received keys and the text parameters to a local directory.

*__Inputs:__* Following options are used:

<ul> 

<li> [Text parameters path]: The path to the text parameters (e.g., *ckks.params*) that will be used for collaborative task</li> 

<li> [Public key path]: Path to the pooled public key, *pooled.public_key*.</li> 

<li> [Relin. key path]: Path to the pooled public key, *pooled.relin_keys*.</li>

<li> [Galois key path]: Path to the pooled public key, *pooled.galois_keys*.</li>

<li> [Site's secret key share path]: Path to the pooled public key, *site_0.secret_key*.</li>

<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 
</ul>

*__Outputs:__* None.

*__Remarks:__* This option is used to store the keys and parameter settings in one directory. *COLLAGENE.sh* takes the path to this directory to read the necessary keys and the parameters.

***

##### >> COLLAGENE.sh -validate_ckks_text_params [Parameter directory] [Validity output file]
This option checks the validity of parameters under 128-bit security.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Validity output file]: This file contains one string "VALID" or "INVALID".</li> 
</ul>

*__Outputs:__* Validity output file is written.

*__Remarks:__* This option checks for the validity of the parameter settings in ckks file under 128-bit security as defined by homomorphicencryption.org's white paper.

***

### Secure Matrix Operations:

##### >> COLLAGENE.sh -encrypt_plaintext_matrix [Parameter directory] [Plaintext matrix path] [Output matrix path]

This option encrypts a plaintext matrix file.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li>

<li> [Plaintext matrix path]: The path to the binary plaintext matrix file. </li>

<li> [Output matrix path]: Path to the encrypted matrix file.</li>
</ul>

*__Outputs:__* Encrypted matrix file.

*__Remarks:__* This option uses a simple binary formatted matrix file. This file can be generated from a tab-delimited text file using option *-convert_plaintext_matrix_text_2_bin* below.

***

##### >> COLLAGENE.sh -full_decrypt_encrypted_matrix [Parameter directory] [Master secret key file] [Encrypted matrix path] [Output matrix path]

This option encrypts a plaintext matrix file.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.

<li> [Master secret key file]: This is the master secret key file.</li>

<li> [Encrypted matrix path]: Path to the encrypted matrix file.</li>

<li> [Output matrix path]: The decrypted output binary matrix file.</li> 
</ul>

*__Outputs:__* Decrypted binary matrix file.

*__Remarks:__* This option should be only available in debugging environment where a pooled master key is available. Generally, the users will not have direct access to this key unless they are running their own server.

***

##### >> COLLAGENE.sh -transpose_encrypted_vector [Parameter directory] [Encrypted matrix path (1xn or nx1)] [Output matrix path]

Transpose an encrypted vector.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li>

<li> [Encrypted matrix path (1xn or nx1)]: The path to the encrypted matrix file that stored a vector.</li>

<li> [Output matrix path]: Transposed vector file.
</ul>

*__Outputs:__* Transposed vector stored in the output file.

*__Remarks:__* The input matrix must be a vector, i.e., number of columns or rows is 1.

***

##### >> COLLAGENE.sh -secure_multiply_rowcol_expansion [Parameter directory] [Encrypted matrix A column expansion Directory] [Encrypted matrix B row expansion Directory] [Output matrix path]

Multiply row and column expanded matrices. 

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix A column expansion Directory]: The directory that stores the column-expansion of left matrix.</li> 

<li> [Encrypted matrix B row expansion Directory]: The directory that stores the row-expansion of right matrix.</li> 

<li> [Output matrix path]: Encrypted result matrix.</li>
</ul>

*__Outputs:__* Product matrix.

*__Remarks:__* The sizes of the input matrices must be conformant for matrix multiplication, i.e., the number of columns of column-expanded (left) matrix must be same as the number of rows in row-expanded (right) matrix.

***

##### >> COLLAGENE.sh -row_expand_plaintext_matrix [Parameter directory] [Plaintext matrix path] [# rows per expansion] [Encrypted row expansion directory]

Row expand a plaintext binary formatted matrix.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Plaintext matrix path]: Input binary formatted plaintext matrix.</li> 

<li> [Encrypted row expansion directory]: Directory that stores the encrypted row expanded matrix. </li> 
</ul>

*__Outputs:__* Matrix row-expansion.

*__Remarks:__* The directory must be created before this option can be used.

***

##### >> COLLAGENE.sh -secure_row_expand_encrypted_matrix [Parameter directory] [Encrypted matrix A path] [# rows per expansion] [Encrypted row expansion directory]

Row expand an encrypted matrix. Matrix must have dimensions that are power of 2.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix A path]: Encrypted input matrix.</li> 

<li> [# rows per expansion]: Number of rows in the expansion.</li> 

<li> [Encrypted row expansion directory]: Directory that stores the encrypted row-expanded matrix.</li> 
</ul>

*__Outputs:__* Matrix row-expansions.

*__Remarks:__* The number of rows and columns of the encrypted matrix must be powers of 2. The directory must be created before this option can be used.

***

##### >> COLLAGENE.sh -col_expand_plaintext_matrix [Parameter directory] [Plaintext matrix path] [# columns per expansion] [Encrypted column expansion directory]

Column expand an encrypted matrix. Matrix must have dimensions that are power of 2.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Plaintext matrix path]: Encrypted input matrix.</li> 

<li> [# columns per expansion]: Number of columns in the expansion.</li> 

<li> [Encrypted columns expansion directory]: Directory that stores the encrypted column-expanded matrix.</li> 
</ul> 

*__Outputs:__* Matrix column-expansions.

*__Remarks:__* The directory must be created before this option can be used.

***

##### >> COLLAGENE.sh -pool_col_expanded_matrices [Parameter directory] [Column expanded matrix directories list file] [Pooled column expanded matrix directory]

Pool a list of column expanded matrices by summation.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.

<li> [Column expanded matrix directories list file]: A file that contains the paths to the column-expanded matrix directories, i.e., each row in this file should point to a column-expansion directory.</li>

<li> [Pooled column expanded matrix directory]: The directory in which the pooled columns expansions will be stored.</li>
</ul> 

*__Outputs:__* Column expansion of the pooled matrices.

*__Remarks:__* This option makes it convenient to pool expansions. The dimensions of the expansions must be the same. The output directory must be created before this option can be used.

***

##### >> COLLAGENE.sh -pool_row_expanded_matrices [Parameter directory] [Row expanded matrix directories list file] [Pooled row expanded matrix directory]

Pool a list of row expanded matrices by summation.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li>

<li> [Row expanded matrix directories list file]: A file that contains the paths to the row-expanded matrix directories, i.e., each row in this file should point to a row-expansion directory.</li>

<li> [Pooled row expanded matrix directory]: Directory where pooled expansions will be stored.</li>
</ul>

*__Outputs:__* Row expansion of the pooled matrices.

*__Remarks:__* The dimensions of the expansions must be the same. The directory must be created before this option can be used.

***

##### >> COLLAGENE.sh -secure_row2row_multiply_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]

Calculate row-by-row inner products of two matrices, namely A and B.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix A path]: The path to the first encrypted matrix.</li> 

<li> [Encrypted matrix B path]: The path to the second encrypted matrix.</li> 

<li> [Output matrix path]: An encrypted output vector with the row-row inner product values. The length of the vector is equal to the number of rows of the matrices.</li> 
</ul>

*__Outputs:__* Encrypted vector that holds the row-row inner product values.

*__Remarks:__* The dimensions of the matrices must be the same. The number of columns must be a power of 2. This can be accomplished by padding the columns to next power of two.

***

##### >> COLLAGENE.sh -secure_elementwise_multiply_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]

Calculate element-by-element multiplication of two matrices, namely A and B.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li>

<li> [Encrypted matrix A path]: The path to the first encrypted matrix.</li>

<li> [Encrypted matrix B path]: The path to the second encrypted matrix.</li>

<li> [Output matrix path]: An encrypted output matrix containing the elementwise multiplication of A and B.</li>
</ul> 

*__Outputs:__* Encrypted vector that holds the product matrix.

*__Remarks:__* The dimensions of the matrices must be the same. The size of output matrix is the same as input matrices.

***

##### >> COLLAGENE.sh -secure_add_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]

Calculate summation of two matrices, namely A and B.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li>

<li> [Encrypted matrix A path]: The path to the first encrypted matrix.</li>

<li> [Encrypted matrix B path]: The path to the second encrypted matrix.</li>
</ul>

<li> [Output matrix path]: An encrypted output matrix that contains sum of A and B.</li>

*__Outputs:__* Encrypted matrix that holds the sum matrix.

*__Remarks:__* The dimensions of the matrices must be the same. The size of output matrix is the same as input matrices.

***

##### >> COLLAGENE.sh -secure_add_matrix_list [Parameter directory] [Encrypted matrix list path] [Output matrix path]

Pool a list of matrices by summing.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix list path]: The list of encrypted matrix files that will be summed (pooled). Each row in this list is an encrypted matrix file.</li> 

<li> [Output matrix path]: Output encrypted matrix of the sum matrix.</li> 
</ul>

*__Outputs:__* Encrypted matrix that holds the sum matrix.

*__Remarks:__* The dimensions of the matrices must be the same. The size of output matrix is the same as input matrices.

***

##### >> COLLAGENE.sh -secure_sub_matrices [Parameter directory] [Encrypted matrix A path] [Encrypted matrix B path] [Output matrix path]

Subtract right matrix (B) elementwise from left matrix (A).

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li>

<li> [Encrypted matrix A path]: The path to the first encrypted matrix.</li>

<li> [Encrypted matrix B path]: The path to the second encrypted matrix.</li>

<li> [Output matrix path]: An encrypted output matrix that contains A-B.</li>
</ul>

*__Outputs:__* Encrypted matrix that holds the subtracted matrix.

*__Remarks:__* The dimensions of the matrices must be the same. The size of output matrix is the same as input matrices.

***

##### >> COLLAGENE.sh -write_encrypted_matrix_dimensions [Parameter directory] [Encrypted matrix path] [Dimensions file path]

Write the dimensions of a matrix into a text file.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix path]: The path to the encrypted matrix.</li> 

<li> [Dimensions file path]: 2-column output file that contains the dimensions of the encrypted matrix.</li> 
</ul>

*__Outputs:__* Matrix dimensions.

*__Remarks:__* This function is useful for programmatically track matrix dimensions in pipelines, e.g., matrix padding is used.

***

##### >> COLLAGENE.sh -write_encrypted_matrix_vital_stats [Parameter directory] [Encrypted matrix path] [Vital stats file path]

Write the scale, size, and chain index for the matrix.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix path]: The path to the encrypted matrix.</li> 

<li> [Dimensions file path]: Multicolumn file that contains the size, chain index, modulus degree, and scale of the ciphertexts in the matrix.</li> 
</ul>

*__Outputs:__* Encrypted matrix ciphertext information.

*__Remarks:__* This function is useful for programmatically track ciphertext information in pipelines such as the chain index, which relates directly to the number of multiplications.

***

### Collaborative Decryption:
##### >> COLLAGENE.sh -partial_decrypt_matrix [Parameter directory] [Encrypted matrix path] [Site index (0-based)] [Decryption noise variance bit size (40-bits by default)] [Output matrix path]

Partially decrypts a matrix locally using the site's DSK.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix path]: The path to the encrypted matrix to be partially decrypted.</li> 

<li> [Site index (0-based)]: This is the unique index of the local site starting from *0, 1 ,..., n_sites-1*.</li> 

<li> [Decryption noise variance bit size (40-bits by default)]: The variance of the noise in number of bits, i.e., smudging noise. If this is set to 0, 40-bits of noise is added to the partially decrypted data.</li> 

<li> [Output matrix path]: Partially decrypted matrix file. This file is not immediately readable by the site.</li> 
</ul>

*__Outputs:__* Partially decrypted matrix file.

*__Remarks:__* The site index has a special meaning in the decryption and should be correctly set.

***

##### >> COLLAGENE.sh -pool_partially_decrypted_matrices [Parameter directory] [Partial decryptions list path] [Output matrix path]

Pool partially decrypted matrices from all sites.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Partial decryptions list path]: List of partially decrypted files obtained from other sites. Each file's path in in a separate line, i.e., if there are 3 sites, there must be 3 lines in the list file pointing to the partially decrypted matrices from the corresponding sites.</li> 

<li> [Output matrix path]: Fully decrypted matrix file.</li> 
</ul>

*__Outputs:__* Fully decrypted matrix file.

*__Remarks:__* Each site should have obtained the partially decrypted matrices before running this option, e.g., downloaded from the shared space. The final decrypted matrix is stored in a plaintext binary matrix.

***

### Symmetric Encryption/Decryption (openssl):

##### >> COLLAGENE.sh -symmetric_encrypt_partdec_data [Parameter directory] [Partdec file] [openssl hash file] [Output file]

Symmetrically encrypt partially decrypted matrices.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Partdec file]: Path to the partially decrypted, "partdec", matrix file generated by partial decryption (i.e., *-partial_decrypt_matrix* option).</li> 

<li> [openssl hash file]: This is the encryption key that is generated by KeyMaker (i.e., *partdec_data_enc_hash.symmetric_key*)</li> 

<li> [Output file]: Symmetrically encrypted "partdec" matrix.</li> 
</ul> 

*__Outputs:__* Symmetrically encrypted partdec matrix file.

*__Remarks:__* This option should be used for protecting the partial decryptions when the shared space is not trusted. Every time a site partially decrypts a matrix, it can be encrypted using this option.

***

##### >> COLLAGENE.sh -symmetric_decrypt_partdec_data [Parameter directory] [Encrypted partdec file] [openssl hash file] [Output file]

Symmetrically decrypt the partially decrypted matrices that were encrypted by the symmetric partdec key.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted partdec file]: The path to the encrypted partdec file.</li> 

<li> [openssl hash file]: Partdec encryption key that is generated by KeyMaker (i.e., *partdec_data_enc_hash.symmetric_key*)</li> 

<li> [Output file]: Decrypted "partdec" matrix.</li> 
</ul>

*__Outputs:__* Fully decrypted matrix file.

*__Remarks:__* This option is used when the partdec matrices are symmetrically encrypted.

***

### Secure Matrix Masking:

##### >> COLLAGENE.sh -generate_plaintext_mask_matrix [Parameter directory] [Encrypted matrix file] [Mask noise standard deviation] [Output plaintext matrix file]

Generate an plaintext mask matrix for an encrypted matrix file.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix file]: The path to the encrypted matrix file.</li> 

<li> [Mask noise standard deviation]: The range for the uniformly distributed noise at each entry of the mask file.</li> 

<li> [Output plaintext matrix file]: Output plaintext matrix file.</li> 
</ul>

*__Outputs:__* Plaintext mask file.

*__Remarks:__* This option generates a mask file that can be used for masking an encrypted matrix file.

***

##### >> COLLAGENE.sh -mask_encrypted_matrix [Encrypted matrix file] [Encrypted mask matrix path] [Output matrix file]

Additively mask an encrypted matrix.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix file]: The path to the encrypted matrix file.</li> 

<li> [Encrypted mask matrix path]: Matrix to be used for masking.</li> 

<li> [Output plaintext matrix file]: Encrypted masked matrix.</li> 
</ul>

*__Outputs:__* Encrypted masked matrix file.

*__Remarks:__* This option uses a generated mask file and adds it to the encrypted matrix.

***

##### >> COLLAGENE.sh -unmask_encrypted_matrix [Parameter directory] [Encrypted matrix file] [Encrypted mask matrix file] [Output matrix path]

Removes additive mask from an encrypted matrix.

*__Inputs:__*

<ul>
<li> [Parameter directory]: The directory that stores the copies of keys and parameter file.</li> 

<li> [Encrypted matrix file]: The path to the encrypted matrix file.</li> 

<li> [Encrypted mask matrix path]: Matrix to be used for unmasking.</li> 

<li> [Output plaintext matrix file]: Encrypted unmasked matrix.</li> 
</ul> 

*__Outputs:__* Encrypted unmasked matrix file.

*__Remarks:__* This option removes a previous mask from the encrypted matrix. The mask used for unmasking should match the matrix that was used for masking.

***

### Plaintext Matrix Operations (Does not require keys or parameter setup):
##### >> COLLAGENE.sh -generate_random_plaintext_matrix [# rows] [# columns] [Output matrix path]

Generate a random matrix with specific size using unit variance Gaussian distribution.

*__Inputs:__*

<ul>
<li> [# rows]: Number of rows.</li> 

<li> [# columns]: Number of columns.</li> 

<li> [Output matrix path]: Output plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Plaintext matrix file with random entries.

*__Remarks:__* This option generates Gaussian distributed random matrix that can be used as a multiplicative mask matrix.

***

##### >> COLLAGENE.sh -generate_constant_plaintext_matrix [# rows] [# columns] [Scalar value] [Output matrix path]

Generate a plaintext matrix with constant value.

*__Inputs:__*

<ul>
<li> [# rows]: Number of rows.</li> 

<li> [# columns]: Number of columns.</li> 

<li> [Scalar value]: Constant scalar value of the entries in the matrix.</li> 

<li> [Output matrix path]: Output plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Plaintext matrix file with constant entries.

***

##### >> COLLAGENE.sh -generate_diagonal_plaintext_matrix [# rows] [# columns] [Scalar value] [Output matrix path]

Generate a plaintext matrix with constant value only at the diagonals.

*__Inputs:__*

<ul>
<li> [# rows]: Number of rows.</li> 

<li> [# columns]: Number of columns.</li> 

<li> [Scalar value]: Constant scalar value of the diagonal entries in the matrix.</li> 

<li> [Output matrix path]: Output plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Plaintext matrix file with constant entries at the diagonal.

***

##### >> COLLAGENE.sh -scalar_multiply_plaintext_matrix [Plaintext matrix path] [Scalar value] [Output matrix path]

Multiple a plaintext matrix with a scalar value.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Scalar value]: Scalar value.</li> 

<li> [Output matrix path]: Output plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Plaintext matrix file with random entries.

***

<!-- ##### >> COLLAGENE.sh -invert_plaintext_matrix [Plaintext matrix path] [Output matrix path]
##### >> COLLAGENE.sh -transpose_plaintext_matrix [Plaintext matrix path] [Output matrix path]
##### >> COLLAGENE.sh -multiply_plaintext_matrix [Matrix A file] [Matrix B file] [Output plaintext matrix file]
##### >> COLLAGENE.sh -multiply_elementwise_plaintext_matrices [Matrix A file] [Matrix B file] [Output plaintext matrix file]
##### >> COLLAGENE.sh -row2row_multiply_plaintext_matrices [Matrix A file] [Matrix B file] [Output plaintext matrix file]
##### >> COLLAGENE.sh -add_plaintext_matrix [Matrix A file] [Matrix B file] [Output plaintext matrix file]
##### >> COLLAGENE.sh -add_plaintext_matrix_list [Plaintext matrix list file] [Output plaintext matrix file]
##### >> COLLAGENE.sh -write_plaintext_matrix_dimensions [Parameter directory] [Plaintext matrix path] [Dimensions file path]--->

##### >> COLLAGENE.sh -pad_plaintext_matrix_2_to_n [Plaintext matrix path] [Output matrix path]

Pad a plaintext matrix to the nearest 2^n dimensions by 0 padding in rows and columns.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Output matrix path]: Output padded plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Padded plaintext matrix file with 0 padding.

*__Remarks:__* This option identifies the next closest power of 2 for rows and columns and pads the matrix with zeroes to match this size.

***

##### >> COLLAGENE.sh -pad_plaintext_matrix_row_2_to_n [Plaintext matrix path] [Output matrix path]

Pad a plaintext matrix to the nearest 2^n dimensions by 0 padding in the rows.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Output matrix path]: Output row padded plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Padded plaintext matrix file with 0 padding in the rows.

*__Remarks:__* This option identifies the next closest power of 2 for rows and pads the matrix with zeroes to match this size.

***

##### >> COLLAGENE.sh -pad_plaintext_matrix_col_2_to_n [Plaintext matrix path] [Output matrix path]

Pad a plaintext matrix to the nearest 2^n dimensions by 0 padding in the columns.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Output matrix path]: Output column padded plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Padded plaintext matrix file with 0 padding in the columns.

*__Remarks:__* This option identifies the next closest power of 2 for columns and pads the matrix with zeroes to match this size.

***

##### >> COLLAGENE.sh -pad_plaintext_matrix_custom [Plaintext matrix path] [Padded # rows] [Padded # columns] [Output matrix path]

Pad a plaintext matrix to custom dimensions by 0 padding in the rows and columns.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Padded # rows]: Padded row number.</li> 

<li> [Padded # columns]: Padded column number.</li> 

<li> [Output matrix path]: Output column padded plaintext binary matrix file.</li> 
</ul>

*__Outputs:__* Padded plaintext matrix file with 0 padding in the columns.

*__Remarks:__* This option pads the matrix with 0's to match the requested size. If the requested number of rows and columns are smaller, the output matrix uses first set of rows and columns, i.e., rows/columns are trimmed at the ends.

***

##### >> COLLAGENE.sh -save_matrix_text [Plaintext matrix path] [Output text matrix path]

Save a plaintext binary matrix file in text format.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Output matrix path]: Tab-delimited matrix file.</li> 
</ul>

*__Outputs:__* Text file for the plaintext matrix.

*__Remarks:__* This option is useful to convert decrypted matrices into text format.

***

##### >> COLLAGENE.sh -transform_plaintext_elementwise_per_callback [Plaintext matrix path] [Function: "log"/"exp"/"sigmoid"/"inv"] [Output matrix path]

Elementwise transform all entries in a plaintext matrix using a callback function.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Plaintext matrix file.</li> 

<li> [Function: "log"/"exp"/"sigmoid"/"inv"]: The function to be used for transforming.</li> 

<li> [Output matrix path]: Transformed plaintext output matrix.</li> 
</ul>

*__Outputs:__* Plaintext matrix.

*__Remarks:__* This option provides convenience when it is being processed elementwise using the callback functions.

***

##### >> COLLAGENE.sh -convert_plaintext_matrix_text_2_bin [Matrix A file] [Row ids in input (first col.)? (0/1)] [Col. id's in input (first row)? (0/1)] [Output plaintext matrix file]

Convert a tab-delimited plaintext matrix into binary formatted matrix.

*__Inputs:__*

<ul>
<li> [Plaintext matrix path]: Tab-delimited text file containing the matrix.</li> 

<li> [Row ids in input (first col.)? (0/1)]: Indicator for existence of row identifiers. </li>

<li> [Col. id's in input (first row)? (0/1)]: Indicator for existence of column identifiers. </li> 

<li> [Output matrix path]: Binary formatted matrix file.</li> 
</ul> 

*__Outputs:__* Text file for the plaintext matrix.

*__Remarks:__* This option is useful to convert tab-delimited matrix into binary format that can be encrypted.

***

##### Binary Executable Options
There are similar options to the above commands that can be used with COLLAGENE binary executable, i.e., *bin/COLLAGENE_Release*. The users can have more control over the options provided by the executable. For example, it is not explicitly necessary to setup the parameters directory when working with the executable file. Using executable file options may be more appropriate when optimizing workflows.

***

## FILE_IO_UTILS.sh:
File Input/Output is used for sharing encrypted datasets such as partdec files, mask files, encrypted summary matrices.

##### ./FILE_IO_UTILS.sh -clean_shared_directory [Data config file]

Remotely delete files in shared directory.

*__Inputs:__*

<ul>
<li> [Data config file]: Configuration file for I/O operations, it should describe the type of file I/O used for sharing data files (i.e., *data_config.params*).</li>
</ul>

*__Outputs:__* None.

*__Remarks:__* This option deletes all of the encrypted intermediate files in the shared directory. This is necessary whenever new collaborative calculations will be executed between sites. Users must be careful to ensure that these directories do not include important files.

***

##### ./FILE_IO_UTILS.sh -test_file_IO [Data config file]

Test file I/O to shared directory.

*__Inputs:__*

<ul>
<li> [Data config file]: Configuration file for I/O operations, it should describe the type of file I/O used for sharing data files (i.e., *data_config.params*).</li>
</ul

*__Outputs:__* None. Script returns 0 upon successful completion, this can be used for error checking in scripts:
```
./FILE_IO_UTILS.sh -test_file_IO data_config.params

if [[ $? != 0 ]]
then
	echo "Network I/O failed!"
	exit 1
fi

...
```

*__Remarks:__* This option uploads a random file, probes it, downloads, and compares the download for testing. Each site should run a test before protocols are executed to ensure that file I/O works correctly.

***

##### ./FILE_IO_UTILS.sh -probe_files_in_shared [Data config file] [File/Directory list]

Probe a list of files/directories in shared directory. 

*__Inputs:__*

<ul>
<li> [Data config file]: Configuration file for I/O operations, it should describe the type of file I/O used for sharing data files (i.e., *data_config.params*).</li> 

<li> [File/Directory list]: A file with the list of files and directories to be probed in the shared space.</li> 
</ul>

*__Outputs:__* The return code indicates whether all of the files and directories are in the shared space.

*__Remarks:__* This option probes for a specific file or directory in the shared space. It returns immediately without blocking. This option works only on base level directories and files.

***

##### ./FILE_IO_UTILS.sh -wait_for_files_in_shared [Data config file] [File/Directory list]

Waits for a list of files/directories in shared directory. 

*__Inputs:__*

<ul>
<li> [Data config file]: Configuration file for I/O operations, it should describe the type of file I/O used for sharing data files (i.e., *data_config.params*).</li> 

<li> [File/Directory list]: A file with the list of files and directories to be waited on within the shared space. This option works only on base level directories and files.</li> 
</ul>

*__Outputs:__* None.

*__Remarks:__* This option probes for a list of specific files or directories in the shared space before they can be downloaded and processed. It blocks until all files are ready to be downloaded. Non-zero return code indicates error (See error message).

***

##### ./FILE_IO_UTILS.sh -download_files_from_shared [Data config file] [File/Directory list]

Download a list of files/directories from the shared directory. 

*__Inputs:__*

<ul>
<li> [Data config file]: Configuration file for I/O operations, it should describe the type of file I/O used for sharing data files (i.e., *data_config.params*).</li> 

<li> [File/Directory list]: A file with the list of files and directories to be downloaded from the shared space. This option works only on base level directories and files.</li> 
</ul>

*__Outputs:__* None.

*__Remarks:__* This option downloads a list of specific files or directories from the shared space. The files that do not exist are skipped. To avoid this, the files should be waited on. Non-zero return code indicates error (See error message).

***

##### ./FILE_IO_UTILS.sh -upload_files_to_shared [Data config file] [File/Directory list]

Upload a list of files/directories to the shared directory. 

*__Inputs:__*

<ul>
<li> [Data config file]: Configuration file for I/O operations, it should describe the type of file I/O used for sharing data files (i.e., *data_config.params*).</li> 

<li> [File/Directory list]: A file with the list of files and directories to be uploaded to the shared space. This option works only on base level directories and files.</li> 
</ul>

*__Outputs:__* None.

*__Remarks:__* This option uploads a list of specific files or directories to the shared space. The operations are blocked until all uploads are completed. Non-zero return code indicates error (See error message).
