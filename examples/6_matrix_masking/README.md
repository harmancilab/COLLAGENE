# Matrix Masking (matrix_masking.sh)

Numerous operations such as matrix inversion are challenging for conversion into HE. While these opeartions can be performed using a procedure called *bootstrapping*, it is computationally very intensive an is not practical.

COLLAGENE implements a commonly used approach that is based on masking of the encrypted data using site-specific random noise patterns to alleviate computational requirements. 

The masked-plaintext data procedure makes use of an MPC-type approach and takes advantage of the collaborative or federated nature of analysis. It integrates a collective decryption step into analysis after data is hidden under site-specific noise patterns. After decryption, the masked (noisy) data can be processed in plaintext to perform heavy calculations.

The basic idea is as following:
1. Mask the encrypted data matrix using site-specific random noise, 
2. Collectively decrypting the masked data matrix,
3. Process the decrypted masked data in plaintext domain,
4. Re-encrypt the data matrix using the public key,
5. Remove the masking noise from the encrypted matrix to generate the final encrypted result.

This procedure serves two purposes: 
1. Firstly, appropriately selected noise levels protect the data in plaintext calculations and alleviates challenges around performing a full HE calculation. 
2. Secondly, the final encrypted results is a fresh (or almost fresh) encryption that can be further operated on by the sites.

COLLAGENE provides a matrix masking procedure that can be used to generate additive or multiplicative noise into data in encrypted domain that can be removed after processing.

A similar masking-based procedure enables refreshing of ciphertext.

We use matrix masking in the encrypted domain before collective decryption for simplifying numerous operations in collaborative setings.

For instance, matrix inversion can be performed after the matrix is masked with multiplicative noise from all sites.

---

*matrix_masking.sh* script implements several examples to demonstrate usage of matrix masking.