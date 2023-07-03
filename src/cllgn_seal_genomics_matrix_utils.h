#ifndef __SEAL_GENOMICS_MATRIX_UTILS__
#define __SEAL_GENOMICS_MATRIX_UTILS__


#include <stdio.h>	
#include <stdlib.h>
#include <seal/seal.h>
#include "SEAL-main/native/examples/examples.h"
#include <seal/util/common.h>
#include "seal/keygenerator.h"
#include "seal/randomtostd.h"
#include "seal/util/common.h"
#include "seal/util/galois.h"
#include "seal/util/ntt.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include "seal/util/rlwe.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintcore.h"

//#include <chrono>
//#include <sys/resource.h>   // check the memory
//#include <random>
//#include <iomanip>          // std::setprecision

#include <vector>

using namespace std;
using namespace seal;
using namespace seal::util;

void write_enc_matrix_dimensions(char* enc_mat_fp, char* op_fp);
void write_plain_matrix_dimensions(char* plain_mat_fp, char* op_fp);

void encrypt_encoded_pt_matrix(char* encoded_matrix_pt_path, 
	char* text_params_path,
	char* pooled_public_key_path,
	char* encrypted_data_path);

void transpose_continuous_encrypted_vector(char* enc_mat_fp, char* text_params_path, char* op_fp);

void plain_unpad_matrix_to_size(char* matrix_fp, int new_n_row, int new_n_col, char* op_fp);

void pad_matrix_cols_to_next_power_of_2(char* mat_fp, char* op_fp);
void pad_matrix_rows_to_next_power_of_2(char* mat_fp, char* op_fp);
void pad_matrix_sizee_to_next_power_of_2(char* mat_fp, char* op_fp);

void plain_invert_matrix(char* A_mat_fp, char* op_fp);

//void plain_add_matrices_per_list(char* matrix_path_list_fp, char* op_fp);

void encrypt_plaintext_matrix_continuous_ct(double** matrix, int nrow, int ncol,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_fp);

void encrypt_plaintext_matrix(double** matrix, int nrow, int ncol,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_fp);

void save_enc_matrix(vector<Ciphertext>* matrix_cts, int nrow, int ncol,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_fp);

void save_continuous_enc_matrix(vector<Ciphertext>* matrix_cts, int nrow, int ncol,
	char* text_params_path,
	char* op_fp);

vector<Ciphertext>* load_continuous_enc_matrix(char* enc_matrix_fp, int& loaded_nrow, int& loaded_ncol, char* text_params_path);

vector<Ciphertext>* load_enc_matrix(char* enc_matrix_fp, int& loaded_nrow, int& loaded_ncol, int& loaded_n_cts_per_row, char* text_params_path);

// This is used for partially decrypting other user's data matrices.
void partial_decrypt_matrix(char* encA_mat_fp, bool is_site0, char* text_params_path, char* noisy_sk_path, double smdg_noise_variance, char* partial_decrypt_data_path);

// This is used for pooling other site's decryptions into a final result.
void pool_partially_decrypted_matrix(char* partial_decrypted_matrix_path_list_fp, char* text_params_path, char* full_decrypted_matrix_fp);

void fully_decrypt_matrix(char* enc_matrix_path, char* text_params_path, char* pooled_secret_key_path, char* full_decrypted_matrix_fp);

void fully_decrypt_continuous_encrypted_matrix(char* enc_matrix_fp, char* text_params_path,
	char* pooled_private_key_path,
	char* full_decrypted_matrix_fp);

// This function multiplies the matrices AxB -- Note that B is given as Bt here.
// The problem with this approach is that each inner product provides one entry and we need to do a lot of shifts to move them into the result matrix.
// Therefore this function will not be efficient for large matrices.
void secure_multiply_matrices(char* encA_mat_fp, char* encBt_mat_fp, char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // TB Removed..
	char* op_fp);

void secure_multiply_matrices_Acol_Brow_expansions(char* A_dir, char* B_dir,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp);

void col_expand_dense_encrypt_matrix(char* A_fp, int n_cols_per_expanded_matrix,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_dir);

void row_expand_dense_encrypt_matrix(char* A_fp, int n_rows_per_expanded_matrix,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_dir);

void row_expand_continuous_encrypted_matrix(char* enc_A_fp, int n_rows_per_expanded_matrix, char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_dir);

void encrypt_plaintext_matrix_continuous_ct(double** matrix, int nrow, int ncol,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_fp);

void secure_add_continuous_encrypted_matrices_per_list(char* enc_mat_list_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp);

void secure_add_continuous_encrypted_matrices(char* enc_A_fp, char* enc_B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp);

void secure_subtract_continuous_encrypted_matrices(char* A_fp, char* B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp);

void secure_multiply_elementwise_continuous_encrypted_matrices(char* A_fp, char* B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp);

void secure_row2row_inner_prod_continuous_encrypted_matrices(char* A_fp, char* B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp);

void write_vital_stats_per_continuous_encrypted_matrix(char* enc_mat_fp,
	char* text_params_path,
	char* vital_stats_fp);

void add_remove_encrypted_mask_2_continuous_encrypted_data_matrix(char* encrypted_data_matrix_path, char* encrypted_mask_matrix_fp,
	bool add_mask_flag,
	char* text_params_path,
	char* masked_encrypted_data_matrix_fp);

void generate_mask_matrix_per_data_matrix(char* encrypted_data_matrix_path,
	double mask_variance,
	char* plain_mask_matrix_op_fp);

void collaborative_pool_partial_decrypted_continuous_enc_plaintext_results(char* partial_decrypted_paths_list_path, char* text_params_path, char* op_fp);
void partial_decrypt_continuous_enc_per_noisy_secretkey_w_smdgng_noise(char* encrypted_data_path, bool is_site0, char* text_params_path, char* noisy_key_path, int n_bits_per_max_smdging_noise, char* partial_decrypt_data_path);

#endif // __SEAL_GENOMICS_MATRIX_UTILS__