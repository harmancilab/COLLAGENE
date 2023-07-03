# Matrix Arithmetic (matrix_arithmetic.sh)
COLLAGENE includes numerous options for matrix arithmetic including elementwise addition, multiplication and matrix multiplication using row/col expansions. This folder contains the script to perform these operations.

---

## Elementwise Operations
Elementwise operations operate at the element-level without any inner product operations. 

Additions can be performed on pairwise or list of encrypted matrices.

## Inner-Product Operations
Opeartions such as matrix multiplications and row-row products are implemented by using special representations.

## Row and Column Expansions
COLLAGENE's current matrix multiplications rely on row and column "expansions" of matrices to simplify multiplications. 

Given a plaintext matrix, COLLAGENE row-expands (column-expands) it into multiple encrypted matrices by repeating its rows (columns) for a pre-determined number of times.

To multiple an axb matrix M1 and bxc matrix M2, M1 is column expanded into b-many axc matrices, where ith column expansion contains the concatenation of M1's ith column c times. Similarly, M2 is row-expanded into b-many axc matrices where ith expansion consists of ith row of M2 repeated a times. COLLAGENE stores column and row expansions directories.

Next, the multiplication of M1 and M2 is computed as the summation of the elementwise multiplication of the corresponding matrices in colexp(M1) and rowexp(M2).

Row and column expansions are distributive over addition operation. To make it simpler to process expansions, COLLAGENE includes options for pooling multiple expansions in parallel.

## Matrix Padding
It is occasionally necessary to change the matrix size by padding zeros at the end of it. COLLAGENE includes options to process plaintext matrices and pad them to custom sizes before encryption or expansions.

## Row-Row inner product
Given two encrypted matrices whose rows are equal to a power of 2, row-row multiplications is implemented as an option in COLLAGENE. 

## Row-expansion of Encrypted Matrices
Given a matrix whose column size is equal to a power of 2, COLLAGENE can generate the row expansion of the encrypted matrix for further processing. This can enables some of the consecutive encrypted matrix multiplications.

It should be noted that these operations may not be sufficient to perform some of the necessary matrix operations such as more than 3 consecutive multiplications. In these cases, it is necessary to reformulate multiplications with respect to the final decrypted result.

## Matrix Masking Operations
Matrix masking provides flexibility for performing some of the more complex operations such as matrix inversion in collaborative setting using multiparty computation-like protocols by adding additive or multiplicative masks at each site.

COLLAGENE includes mask matrix generation options to generate mask matrices that can be encrypted and used for masking matrices.

---

*matrix_arithmetic.sh* script implements numerous options that are discussed above using calls to COLLAGENE.sh script.