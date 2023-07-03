# Data Scaling (data_scaling.sh)
Numeric precision of secure calculations are impacted by the dynamic range of the data matrices that are being processed.

For instance, assume we have a matrix with large values and we re-size the matrix to another dimension by padding with 0's. We encrypt this matrix.

We would like to multiply this matrix with another matrix that contains small values. We pad this matrix to a concordant dimension to be multiplied with the first matrix.

When the padded matrices are multiplied, the product matrix should have exactly 0 entries in the padded portions. 

This, however, occasionally does not hold since the encryptions introduce small errors, which accumulate to non-zero entries in the padded portions of the product matrix.

When the padded product matrix is furhter multiplied by other matrices, the error may accumulate to high levels and may impact results adversely.

In these scenarios, we can scale the second matrix with an appropriate value before padding and encryption. We next perform matrix product. We finally unscale the product matrix. 

This effectively removes the very small noise terms in the padded portions of the product matrix and substantially amplify the signal component of the product matrix.

---

*data_scaling.sh* script implements a padded matrix multiplication to demonstrate usage of data scaling.