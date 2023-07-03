# Function Approximation (function_approximation.sh)
 *function_approximation.sh* demonstrates a collaborative protocol to approximate exponential value of a matrix using numerous options in *COLLAGENE.sh*. 

 In this example, the sites would like to perform element-wise exponentiation of an encrypted matrix, denote by enc(A). This is accomplished by additively masking the matrix with site-specific additive masks, then removing the mask in *inverse exponentiated form*:

<ol>
<li>The sites first generate a uniform noise locally, encrypt it and share with other sites, enc(n1), enc(n2), enc(n3)</li>
<li>The noise matrices are pooled among all sites: enc(n1+n2+n3)</li>
<li>The encrypted noise matrix is added to A: enc(A+n1+n2+n3)</li>
<li>The noisy matrix is collectively decrypted: coldec(enc(A+n1+n2+n3))=A+n1+n2+n3</li>
<li>This noisy matrix is exponentiated elementwise: exp(A+n1+n2+n3)</li>
<li>The exponentiated noisy matrix is encrypted: enc(exp(A+n1+n2+n3))</li>
<li>Next, the sites *multiplicatively* pool the exponentiated-inverse-noise-levels: n_inv=enc(exp(-n1)) x enc(exp(-n2)) x enc(exp(-n3))=enc(exp(-n1-n2-n3))</li>
<li>Finally, the inverted noise is multiplied (elementwise) with the noisy exponentiated matrix: enc(exp(A+n1+n2+n3))xenc(exp(-n1-n2-n3))=enc(exp(A))</li>
</ol>

The main protection in this example stems from pooled masked term (n1+n2+n3) that is added to A. Note that since exponentials can easily overflow/underflow, this mask should not be too strong to avoid numerical issues. 