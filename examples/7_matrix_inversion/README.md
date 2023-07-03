# Matrix Inversion (matrix_inversion.sh)

*matrix_inversion.sh* demonstrates an example of multiplicative matrix masking to perform matrix inversion. 

The matrix inversion is accomplished by: 

<ol>
<li> Adding multiplicative noise in encrypted domain,</li>
<li> Collectively decrypting the noisy matrix,</li>
<li> Inverting the decrypted noisy matrix,</li>
<li> Removing the noise by multiplying with the noise matrix.</li>
</ol>

One of the important factors is the usage of row and column expansions in different steps. 

The noise is added on the right side and the masking uses the row expansions of the noise matrix.

Later, while unmasking the matrix, the noise must be multiplied on the left side with the inverted masked matrix. This is when the column expansion of the noise matrix is used with the row expansion of the inverted matrix.

Overall, the row and column expansions should be managed appropriately to make sure multiplicate masks are removed correctly.



