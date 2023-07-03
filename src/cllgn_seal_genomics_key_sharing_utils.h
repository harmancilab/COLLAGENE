#ifndef __SEAL_GENOMICS_KEY_SHARING_UTILS__
#define __SEAL_GENOMICS_KEY_SHARING_UTILS__

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

struct t_text_params
{
	size_t poly_modulus_degree;
	vector<int>* decomp_bit_sizes_ptr;
	int scale_n_bits;
	int per_site_key_noise_var_bit_size;
};

bool validate_coeff_modulus_bit_size(char* text_params_path);

vector<double*>* load_data_matrix(char* data_matrix_fp, vector<char*>* row_ids, vector<char*>* col_ids,
	bool no_header, bool no_row_names,
	int& n_cols);


void generate_mask_matrix_per_data_matrix(char* encrypted_data_matrix_path,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	double mask_variance,
	char* plain_mask_matrix_op_fp);

void add_remove_encrypted_mask_2_data_matrix(char* encrypted_data_matrix_path, char* encrypted_mask_matrix_fp,
	bool add_mask_flag,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* masked_encrypted_data_matrix_fp);

void write_ciphertext_vital_stats(char* encrypted_data_matrix_path,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* vital_stats_fp);

//
// These functions implement KeyMaker, Client, and Processor functions. All functionalities are included in one executable and can be run at different entities using the same executable.
//

void NTT_transform_coefficients_array_in_place(uint64_t* non_NTT_coeffs_array,
	size_t coeff_modulus_size, size_t coeff_count,
	const NTTTables* small_ntt_tables);

// The per site secret-key masking arrays.
uint64_t** SG_generate_per_site_sk_coeff_selector(int n_sites, int poly_modulus_degree, int coeff_modulus_size, const vector<Modulus> coeff_modulus);

// Generate smudging noise for one site.
uint64_t* SG_generate_smudging_noise_SEAL_RNG(EncryptionParameters context_parms, double smdging_noise_variance);
uint64_t* SG_generate_smudging_noise(int poly_modulus_degree, int coeff_modulus_size, double smdging_noise_variance, const vector<Modulus> coeff_modulus);

// The per site complementary noise arrays.
uint64_t** SG_generate_per_site_complementary_noise_random_complementation(int n_sites, int poly_modulus_degree, int coeff_modulus_size, double noise_variance, const vector<Modulus> coeff_modulus);

t_text_params* load_text_params(char* text_params_path);

vector<double*>* load_data_processing_vectors_for_encryption(char* data_file_path, int& l_processing_vectors);

vector<double*>* load_data_matrix(char* data_matrix_fp, vector<char*>* row_ids, vector<char*>* col_ids,
	bool no_header, bool no_row_names,
	int& n_cols);

uint64_t* SEAL_Genomics_Get_Normal_Random(EncryptionParameters parms);
uint64_t* SEAL_Genomics_Get_Normal_Random_Custom(EncryptionParameters parms, double std_dev, double max_dev);
uint64_t* SEAL_Genomics_Get_Unif_Random(EncryptionParameters parms);
uint64_t* SEAL_Genomics_Get_Unif_Random_Custom(EncryptionParameters parms, const int n_bits_per_max_val);

// Generates the keys for each site; this is the keymaker function. This function also saves the parameters that are used for keys so that they can be used in encryption, evaluation, and decryption steps.
void generate_share_private_keys(int n_sites, char* text_params_path,
	bool perform_key_testing,
	char* op_dir);

// This is called on each site using the pooled public key.
void encrypt_data_per_pooled_public_key(char* data_file_path, char* params_path, char* pooled_pk_path, char* op_fp);

// Secure evaluation code is under seal_genomics_custom_evaluation.cpp

// site0 flag is necessary for keeping track of c0 term in pooling step.
void partial_decrypt_data_per_noisy_key(char* encrypted_data_path, bool is_site0, char* text_params_path, char* noisy_key_path, double smdg_noise_variance, char* partial_decrypt_data_path);

// This function does not do anything other than summing the plaintext's up. The c0 issue is handled in the per site decryption function above.
void collaborative_pool_partial_decrypted_plaintext_results(char* partial_decrypted_paths_list_path, char* params_path, char* op_fp);

// Perform full decryption.
void full_decrypt_encrypted_data_per_pooled_key(char* enc_data_matrix_path, char* text_params_path, char* pooled_secret_key_path, char* op_fp);

#endif // __SEAL_GENOMICS_KEY_SHARING_UTILS__