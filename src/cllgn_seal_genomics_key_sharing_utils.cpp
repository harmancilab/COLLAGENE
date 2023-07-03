#include <stdio.h>
#include <stdlib.h>
#include <seal/seal.h>
#include "cllgn_ansi_string.h"
#include "cllgn_vector_macros.h"
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_matrix_linalg_utils.h"
#include "cllgn_file_utils.h"
#include "SEAL-main/native/examples/examples.h"
#include <seal/util/common.h>
#include "seal/keygenerator.h"
#include "seal/randomtostd.h"
#include "seal/util/common.h"
#include "seal/util/clipnormal.h"
#include "seal/util/galois.h"
#include "seal/util/ntt.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include "seal/util/rlwe.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintcore.h"
#include "cllgn_seal_genomics_key_sharing_utils.h"

//#include <chrono>
//#include <sys/resource.h>   // check the memory
//#include <random>
//#include <iomanip>          // std::setprecision

#include <vector>

using namespace std;
using namespace seal;
using namespace seal::util;

bool __DUMP_KEY_SHARING_MSGS__ = false;

bool validate_coeff_modulus_bit_size(char* text_params_path)
{
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
	}

	// Catch exceptions.
	try
	{
		t_text_params* text_params = load_text_params(text_params_path);

		if (text_params->decomp_bit_sizes_ptr->size() < 3)
		{
			fprintf(stderr, "CKKS needs at least 1 inner level modulus.\n");
			return(false);
		}

		int total_mod_bit_cnt = 0;
		vector<int> coeff_modulus_bit_sizes;
		for (int i_dec = 0; i_dec < (int)text_params->decomp_bit_sizes_ptr->size(); i_dec++)
		{
			coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
			total_mod_bit_cnt += text_params->decomp_bit_sizes_ptr->at(i_dec);

			// Check to make sure if scale is close to the inner bit modulus sizes.
			if (i_dec > 0 &&
				i_dec < (int)text_params->decomp_bit_sizes_ptr->size() - 1)
			{
				if (text_params->scale_n_bits != text_params->decomp_bit_sizes_ptr->at(i_dec))
				{
					fprintf(stderr, "The inner modulus sizes should be the same as CKKS scale parameter: %d vs %d @ %d. modulus.\n",
						text_params->scale_n_bits,
						text_params->decomp_bit_sizes_ptr->at(i_dec),
						i_dec);

					return(false);
				}
			}
		} // i_dec loop.

		EncryptionParameters parms(scheme_type::ckks);

		size_t poly_modulus_degree = text_params->poly_modulus_degree;
		parms.set_poly_modulus_degree(poly_modulus_degree);
		parms.set_coeff_modulus(CoeffModulus::Create(
			poly_modulus_degree, coeff_modulus_bit_sizes));

		if (__DUMP_KEY_SHARING_MSGS__)
		{
			fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
		}

		if (total_mod_bit_cnt >= CoeffModulus::MaxBitCount(poly_modulus_degree))
		{
			return(false);
		}
		else
		{
			return(true);
		}
	}
	catch (std::logic_error& err)
	{
		fprintf(stderr, "Fatal logic error: %s\n", err.what());
		return(false);
	}
	catch (std::runtime_error& err)
	{
		fprintf(stderr, "Fatal runtime error: %s\n", err.what());
		return(false);
	}
	catch (std::exception& exc)
	{
		fprintf(stderr, "Fatal exception: %s\n", exc.what());
		return(false);
	}
}


// Get the data sizes and generate masking matrix: For each subject's value, generate a masking value and add it.
void generate_mask_matrix_per_data_matrix(char* encrypted_data_matrix_path,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	double mask_variance,
	char* plain_mask_matrix_op_fp)
{
	fprintf(stderr, "Generating the plaintext mask matrix:\n\
Encrypted matrix file : %s\n\
Text params file : %s\n\
pooled public key path : %s\n\
pooled relin. key path : %s\n\
pooled Galois key path : %s\n\
Output file : %s\n", encrypted_data_matrix_path, text_params_path, pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path, plain_mask_matrix_op_fp);

	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
	}

	long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
	fprintf(stderr, "Mask seed: %lu\n", rand_seed);
	t_rng* rng = new t_rng(rand_seed);

	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	//double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_public_key_path) ||
		!check_file(pooled_relin_key_path) ||
		!check_file(pooled_galois_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s, %s, %s\n", pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path);
		exit(0);
	}

	// Instantiate and load the public key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}
	
	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	// Load galois key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	}

	ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	GaloisKeys pooled_galois_key;
	pooled_galois_key.load(context, ifs_galois_key);
	ifs_galois_key.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded Galois key.\n");
	}

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	//int n_samples_per_ct = encoder.slot_count();

	// Start loading the data: Load the (1) sample size, (2) number of features, and (3) number of ciphertexts per feature.
	// Open this file.
	ifstream ifs_enc_data_matrix(encrypted_data_matrix_path, ios::binary);

	// Set the sample size information.
	int n_samples;
	int n_feats;
	int n_cts_per_feat;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_samples), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_feats), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_cts_per_feat), sizeof(int));

	fprintf(stderr, "Sample size per feature: %d\nn_feats: %d\nn_cts_per_feat: %d\n", n_samples, n_feats, n_cts_per_feat);

	fprintf(stderr, "Loading ciphertexts.\n");
	vector<Ciphertext>** per_feat_ct_vectors = new vector<Ciphertext>*[n_feats];
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		// This is the list of ct's for the current feature.
		per_feat_ct_vectors[i_feat] = new vector<Ciphertext>();

		for (int i_ct = 0; i_ct < n_cts_per_feat; i_ct++)
		{
			Ciphertext cur_ct;
			cur_ct.load(context, ifs_enc_data_matrix);

			// Add this ct to the current feat's list.
			per_feat_ct_vectors[i_feat]->push_back(cur_ct);
		} // i_ct loop.
	} // i_feat.

	fprintf(stderr, "Writing plaintext mask matrix.\n");

	//int n_slots_per_ct = encoder.slot_count();

	FILE* f_plaintext_mask_matrix = open_f(plain_mask_matrix_op_fp, "w");

	fprintf(f_plaintext_mask_matrix, "#Feat_ID");
	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(f_plaintext_mask_matrix, "\tMask_%d", i_s);
	} // i_s loop.
	fprintf(f_plaintext_mask_matrix, "\n");

	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		//vector<Ciphertext>* cur_masking_cts = new vector<Ciphertext>();

		fprintf(f_plaintext_mask_matrix, "Feat_%d", i_feat);

		// Fill the masking cts.
		vector<double> cur_feat_mask_array;
		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			double cur_mask_val = rng->random_double_ran3() * mask_variance;
			cur_feat_mask_array.push_back(cur_mask_val);
			fprintf(f_plaintext_mask_matrix, "\t%.4f", cur_mask_val);
		} // i_s loop.

		fprintf(f_plaintext_mask_matrix, "\n");
	} // i_feat loop.
	close_f(f_plaintext_mask_matrix, plain_mask_matrix_op_fp);
} // generate_mask_matrix_per_data_matrix option.

void write_ciphertext_vital_stats(char* encrypted_data_matrix_path,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* vital_stats_fp)
{
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Setting up the context..\n");
	}

	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	//double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_public_key_path) ||
		!check_file(pooled_relin_key_path) ||
		!check_file(pooled_galois_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s, %s, %s\n", pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path);
		exit(0);
	}

	// Instantiate and load the public key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	// Load galois key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	}

	ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	GaloisKeys pooled_galois_key;
	pooled_galois_key.load(context, ifs_galois_key);
	ifs_galois_key.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded Galois key.\n");
	}

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	//int n_samples_per_ct = encoder.slot_count();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	fprintf(stderr, "Loading encrypted data matrix from %s\n", encrypted_data_matrix_path);
	ifstream ifs_enc_data_matrix(encrypted_data_matrix_path, ios::binary);

	// Set the sample size information.
	int n_samples;
	int n_feats;
	int n_cts_per_feat;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_samples), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_feats), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_cts_per_feat), sizeof(int));

	fprintf(stderr, "Sample size per feature: %d\nn_feats: %d\nn_cts_per_feat: %d\n", n_samples, n_feats, n_cts_per_feat);

	fprintf(stderr, "Loading ciphertexts.\n");
	vector<Ciphertext>** per_feat_ct_vectors = new vector<Ciphertext>*[n_feats];
	FILE* f_vital_stats_op = open_f(vital_stats_fp, "w");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		// This is the list of ct's for the current feature.
		per_feat_ct_vectors[i_feat] = new vector<Ciphertext>();

		for (int i_ct = 0; i_ct < n_cts_per_feat; i_ct++)
		{
			Ciphertext cur_ct;
			cur_ct.load(context, ifs_enc_data_matrix);

			// Add this ct to the current feat's list.
			per_feat_ct_vectors[i_feat]->push_back(cur_ct);

			fprintf(stderr, "Feat: %d / Ciphertext: %d: Chain_Index: %d\tmodulus_degree: %d\n",
				i_feat,
				i_ct,
				(int)(context.get_context_data(cur_ct.parms_id())->chain_index()),
				(int)(cur_ct.coeff_modulus_size()));

			fprintf(f_vital_stats_op, "Feat: %d / Ciphertext: %d: Chain_Index: %d\tmodulus_degree: %d\n",
				i_feat,
				i_ct,
				(int)(context.get_context_data(cur_ct.parms_id())->chain_index()),
				(int)(cur_ct.coeff_modulus_size()));
		} // i_ct loop.
	} // i_feat.

	// Close file.
	close_f(f_vital_stats_op, vital_stats_fp);
} // write_ciphertext_vital_stats function.

void add_remove_encrypted_mask_2_data_matrix(char* encrypted_data_matrix_path, char* encrypted_mask_matrix_fp,
	bool add_mask_flag,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* masked_encrypted_data_matrix_fp)
{
	if (add_mask_flag)
	{
		fprintf(stderr, "MASKING the plaintext mask matrix:\n\
Encrypted matrix file : %s\n\
Text params file : %s\n\
pooled public key path : %s\n\
pooled relin. key path : %s\n\
pooled Galois key path : %s\n\
Output file : %s\n", encrypted_data_matrix_path, text_params_path, pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path, masked_encrypted_data_matrix_fp);
	}
	else
	{
		fprintf(stderr, "UNMASKING the plaintext mask matrix:\n\
Encrypted matrix file : %s\n\
Text params file : %s\n\
pooled public key path : %s\n\
pooled relin. key path : %s\n\
pooled Galois key path : %s\n\
Output file : %s\n", encrypted_data_matrix_path, text_params_path, pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path, masked_encrypted_data_matrix_fp);

	}

	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
	}

	fprintf(stderr, "Setting up the context..\n");
	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	//double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_public_key_path) ||
		!check_file(pooled_relin_key_path) ||
		!check_file(pooled_galois_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s, %s, %s\n", pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path);
		exit(0);
	}

	// Instantiate and load the public key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	// Load galois key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	}

	ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	GaloisKeys pooled_galois_key;
	pooled_galois_key.load(context, ifs_galois_key);
	ifs_galois_key.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded Galois key.\n");
	}

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	//int n_samples_per_ct = encoder.slot_count();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	fprintf(stderr, "Loading encrypted data matrix from %s\n", encrypted_data_matrix_path);
	ifstream ifs_enc_data_matrix(encrypted_data_matrix_path, ios::binary);

	// Set the sample size information.
	int n_samples;
	int n_feats;
	int n_cts_per_feat;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_samples), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_feats), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_cts_per_feat), sizeof(int));

	fprintf(stderr, "Sample size per feature: %d\nn_feats: %d\nn_cts_per_feat: %d\n", n_samples, n_feats, n_cts_per_feat);

	fprintf(stderr, "Loading ciphertexts.\n");
	vector<Ciphertext>** per_feat_ct_vectors = new vector<Ciphertext>*[n_feats];
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		// This is the list of ct's for the current feature.
		per_feat_ct_vectors[i_feat] = new vector<Ciphertext>();

		for (int i_ct = 0; i_ct < n_cts_per_feat; i_ct++)
		{
			Ciphertext cur_ct;
			cur_ct.load(context, ifs_enc_data_matrix);

			// Add this ct to the current feat's list.
			per_feat_ct_vectors[i_feat]->push_back(cur_ct);
		} // i_ct loop.
	} // i_feat.

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// LOAD THE ENCRYPTED MASK CTs:
	fprintf(stderr, "Loading encrypted mask matrix from %s\n", encrypted_mask_matrix_fp);
	ifstream ifs_enc_mask_matrix(encrypted_mask_matrix_fp, ios::binary);

	// Set the sample size information.
	int n_mask_samples;
	int n_mask_feats;
	int n_mask_cts_per_feat;

	ifs_enc_mask_matrix.read(reinterpret_cast<char *>(&n_mask_samples), sizeof(int));
	ifs_enc_mask_matrix.read(reinterpret_cast<char *>(&n_mask_feats), sizeof(int));
	ifs_enc_mask_matrix.read(reinterpret_cast<char *>(&n_mask_cts_per_feat), sizeof(int));

	fprintf(stderr, "Sample size per feature from mask matrix: %d\nn_feats: %d\nn_cts_per_feat: %d\n", n_mask_samples, n_mask_feats, n_mask_cts_per_feat);

	if (n_samples != n_mask_samples ||
		n_feats != n_mask_feats ||
		n_cts_per_feat != n_mask_cts_per_feat)
	{
		fprintf(stderr, "The data matrix size parameters are not matching among data matrix and mask matrix: %d/%d ;; %d/%d ;; %d/%d\n", n_samples, n_mask_samples,
			n_feats, n_mask_feats,
			n_cts_per_feat, n_mask_cts_per_feat);

		exit(0);
	}

	fprintf(stderr, "Loading encrypted mask data matrix:\n");
	vector<Ciphertext>** per_feat_mask_ct_vectors = new vector<Ciphertext>*[n_feats];
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		// This is the list of ct's for the current feature.
		per_feat_mask_ct_vectors[i_feat] = new vector<Ciphertext>();

		for (int i_ct = 0; i_ct < n_cts_per_feat; i_ct++)
		{
			Ciphertext cur_mask_ct;
			cur_mask_ct.load(context, ifs_enc_mask_matrix);

			// Add this ct to the current feat's list.
			per_feat_mask_ct_vectors[i_feat]->push_back(cur_mask_ct);
		} // i_ct loop.
	} // i_feat.

	/////////////////////////////////////////////////////////////////////////////
	// Add the masks and save them.
	ofstream ofs_masked_enc_data_matrix(masked_encrypted_data_matrix_fp, ios::binary);

	// This describes the matrix that is being saved to output.
	int l_result_vec = n_samples; // This parameter has no bearing yet. It should be copied to partially decrypted plaintext file and final pooling file to write the results in the final final text output file. This sets how many values are meaningful (or necessary) in the resulting ciphertext.
	int n_result_vec = n_feats;
	int n_result_cts_per_vec = n_cts_per_feat;
	ofs_masked_enc_data_matrix.write(reinterpret_cast<const char *>(&l_result_vec), sizeof(int));
	ofs_masked_enc_data_matrix.write(reinterpret_cast<const char *>(&(n_result_vec)), sizeof(int));
	ofs_masked_enc_data_matrix.write(reinterpret_cast<const char *>(&n_result_cts_per_vec), sizeof(int));

	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		for (int i_ct = 0; i_ct < n_cts_per_feat; i_ct++)
		{
			Ciphertext cur_data_ct = per_feat_ct_vectors[i_feat]->at(i_ct);
			Ciphertext cur_mask_ct = per_feat_mask_ct_vectors[i_feat]->at(i_ct);

			fprintf(stderr, "Fresh mask ct modulus size: %d, chain index: %d\n", (int)(cur_mask_ct.coeff_modulus_size()),
				(int)(context.get_context_data(cur_mask_ct.parms_id())->chain_index()));

			// Match the modulus of mask and data ct.
			evaluator.mod_switch_to_inplace(cur_mask_ct, cur_data_ct.parms_id());
			cur_mask_ct.scale() = cur_data_ct.scale();

			fprintf(stderr, "Current data ct modulus size: %d, chain index: %d\n", (int)(cur_data_ct.coeff_modulus_size()),
				(int)(context.get_context_data(cur_data_ct.parms_id())->chain_index()));

			// Mask or unmask.
			if (add_mask_flag)
			{
				evaluator.add_inplace(cur_data_ct, cur_mask_ct);
			}
			else
			{
				evaluator.sub_inplace(cur_data_ct, cur_mask_ct);
			}

			cur_data_ct.save(ofs_masked_enc_data_matrix);
		} // i_ct loop.
	} // i_feat loop.

	ofs_masked_enc_data_matrix.close();
} // add_remove_encrypted_mask_2_data_matrix function.

vector<double*>* load_data_processing_vectors_for_encryption(char* data_file_path, int& l_processing_vectors)
{
	bool no_header = false;
	bool no_row_names = false;
	vector<char*>* row_ids = new vector<char*>();
	vector<char*>* col_ids = new vector<char*>();
	int n_cols = 0;
	vector<double*>* data_vectors = load_data_matrix(data_file_path, row_ids, col_ids,
		no_header, no_row_names,
		n_cols);

	fprintf(stderr, "Loaded %d columns id's, %d row id's.\n", (int)col_ids->size(), (int)row_ids->size());

	l_processing_vectors = n_cols;

	return(data_vectors);
}

/***********************************
Example ckks.params file:
####################################
16384
60 40 40 40 60
40
20
***********************************/
t_text_params* load_text_params(char* text_params_path)
{
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading text parameters from %s\n", text_params_path);
	}

	t_text_params* text_params = new t_text_params();

	FILE* f_text_params = open_f(text_params_path, "r");

	// Parse the polynomial degree modulus.
	char* cur_line = getline(f_text_params);
	if (cur_line == NULL)
	{
		fprintf(stderr, "Could not read polynomial degree modulus from %s\n", text_params_path);
		exit(0);
	}
	//sscanf(cur_line, "%d", &(text_params->poly_modulus_degree));
	sscanf(cur_line, "%lu", &(text_params->poly_modulus_degree));

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded requested polynomial degree %d\n", (int)(text_params->poly_modulus_degree));
	}

	// Parse the decomp. bit sizes.
	cur_line = getline(f_text_params);
	if (cur_line == NULL)
	{
		fprintf(stderr, "Could not read decomposition bit sizes from %s\n", text_params_path);
		exit(0);
	}
	t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, " \t");
	text_params->decomp_bit_sizes_ptr = new vector<int>();
	for (int i_dec = 0; i_dec < (int)(toks->size()); i_dec++)
	{
		text_params->decomp_bit_sizes_ptr->push_back(atoi(toks->at(i_dec)->str()));

		if (__DUMP_KEY_SHARING_MSGS__)
		{
			fprintf(stderr, "%d. decomp. : %d bits\n", i_dec, text_params->decomp_bit_sizes_ptr->at(i_dec));
		}

	} // i_dec loop.

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Added %d coefficient decomp bit sizes.\n", (int)(text_params->decomp_bit_sizes_ptr->size()));
	}

	// Parse the CKKS scale.
	cur_line = getline(f_text_params);
	if (cur_line == NULL)
	{
		fprintf(stderr, "Could not read scale from %s\n", text_params_path);
		exit(0);
	}
	sscanf(cur_line, "%d", &(text_params->scale_n_bits));

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded scale: %d bits\n", text_params->scale_n_bits);
	}

	// Parse the per site key noise variance.
	cur_line = getline(f_text_params);
	if (cur_line == NULL)
	{
		fprintf(stderr, "Could not read per site key noise variance bit size from %s\n", text_params_path);
		exit(0);
	}
	sscanf(cur_line, "%d", &(text_params->per_site_key_noise_var_bit_size));

	close_f(f_text_params, text_params_path);

	return(text_params);
} // test param loading function.

vector<double*>* load_data_matrix(char* data_matrix_fp, vector<char*>* row_ids, vector<char*>* col_ids,
	bool no_header, bool no_row_names,
	int& n_cols)
{
	bool processed_header = no_header;

	n_cols = -1;

	vector<double*>* data_matrix_rows = new vector<double*>();

	FILE* f_data_matrix = open_f(data_matrix_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_data_matrix);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");
		if (!processed_header)
		{
			processed_header = true;

			int start_col_i = (int)(!no_row_names);
			for (int col_i = start_col_i; col_i < (int)(toks->size()); col_i++)
			{
				col_ids->push_back(t_string::copy_me_str(toks->at(col_i)->str()));
			} // col_i loop.
		}
		else
		{
			// Add the data.
			int start_col_i = (int)(!no_row_names);
			double* cur_data_vector = new double[(int)(toks->size()) + 2];
			for (int col_i = start_col_i; col_i < (int)(toks->size()); col_i++)
			{
				cur_data_vector[col_i - start_col_i] = atof(toks->at(col_i)->str());
			} // col_i loop.

			// Set the row id if it is processed.
			if (!no_row_names)
			{
				row_ids->push_back(t_string::copy_me_str(toks->at(0)->str()));
			}

			// Add this row.
			data_matrix_rows->push_back(cur_data_vector);

			if (n_cols == -1)
			{
				n_cols = (int)(toks->size()) - start_col_i;
			}
			else
			{
				if (n_cols != ((int)(toks->size()) - start_col_i))
				{
					fprintf(stderr, "Row %d does not have consistent number of columns: %d, %d:\n%s\n", (int)(data_matrix_rows->size()), n_cols, (int)(toks->size() - start_col_i), cur_line);
					exit(0);
				}
			} // column number check.
		} // header check.

		delete[] cur_line;
	} // file reading loop.
	close_f(f_data_matrix, data_matrix_fp);

	return(data_matrix_rows);
}

// Following creates one vector of doubles for testing purposes, used in testing the full protocol. 
// Dont move this function. It is internally used for DSK protocol testing below.
vector<double> generate_random_data_for_testing(size_t slot_count)
{
	vector<double> input;
	input.reserve(slot_count);
	double curr_point = 0;
	double step_size = 1.0 / (static_cast<double>(slot_count) - 1);
	for (size_t i = 0; i < slot_count; i++)
	{
		input.push_back(curr_point);
		curr_point += step_size;
	}
	cout << "Input vector: " << endl;
	print_vector(input, 3, 7);

	return(input);
}

SecretKey get_secret_key_per_NTT_coeff_array(uint64_t* NTT_trans_secret_key_coeff_array, size_t key_coeff_modulus_size, size_t key_coeff_count)
{
	SecretKey sk_per_coeffs = SecretKey();
	sk_per_coeffs.data().resize(mul_safe(key_coeff_count, key_coeff_modulus_size));
	memcpy(sk_per_coeffs.data().data(), NTT_trans_secret_key_coeff_array, sizeof(uint64_t) * (key_coeff_count * key_coeff_modulus_size));

	return(sk_per_coeffs);
}


// This function is used only for testing purposes, not used in the actual protocol.
uint64_t* pool_nonNTT_secret_keys_from_sites(EncryptionParameters parms, int n_sites,
	vector<uint64_t*>* per_site_noisy_keys,
	size_t key_coeff_modulus_size, size_t key_coeff_count,
	uint64_t* sk_data_inv_ntt)
{
	// Do a sanity check by pooling the keys and making sure they add up to the pooled secret key.
	fprintf(stderr, "Validating that the per-site keys sum up to the pooled secret key.\n");
	uint64_t* pooled_noisy_secret_key = new uint64_t[key_coeff_count * key_coeff_modulus_size + 2];
	memset(pooled_noisy_secret_key, 0, sizeof(uint64_t) * key_coeff_count * key_coeff_modulus_size);
	for (int i_mod = 0; i_mod < (int)key_coeff_modulus_size; i_mod++)
	{
		uint64_t* cur_mod_pooled_key = pooled_noisy_secret_key + i_mod * key_coeff_count;

		for (int i_site = 0; i_site < n_sites; i_site++)
		{
			uint64_t* cur_site_cur_mod_key_coeffs = per_site_noisy_keys->at(i_site) + key_coeff_count * i_mod;
			add_poly_coeffmod(cur_mod_pooled_key, cur_site_cur_mod_key_coeffs, key_coeff_count, parms.coeff_modulus().at(i_mod), cur_mod_pooled_key);
		} // i_site loop.

		// Pooling finished.
		if (sk_data_inv_ntt != NULL)
		{
			fprintf(stderr, "Comparing @ mod %d:\n", i_mod);
			uint64_t* cur_mod_sk_coeffs = sk_data_inv_ntt + i_mod * key_coeff_count;
			for (int i_coeff = 0; i_coeff < (int)key_coeff_count; i_coeff++)
			{
				if (cur_mod_pooled_key[i_coeff] != cur_mod_sk_coeffs[i_coeff])
				{
					fprintf(stderr, "Mismatch @ %d: %lu -- %lu\n", i_coeff,
						cur_mod_pooled_key[i_coeff],
						cur_mod_sk_coeffs[i_coeff]);
					exit(0);
				}
			} // i_coeff loop.

			fprintf(stderr, "Looks good..\n");
		}
	} // i_mod loop.

	return(pooled_noisy_secret_key);
}

// TODO::Rejection sampling-based generation of bounded noise::Infinite norm Gaussian sampling may be problematic as it may leak some information about the master key. 
uint64_t** SG_generate_per_site_sk_coeff_selector(int n_sites, int poly_modulus_degree, int coeff_modulus_size, const vector<Modulus> coeff_modulus)
{
	if (n_sites < 2)
	{
		fprintf(stderr, "Need at least 2 sites to generate the noise levels.\n");
		exit(0);
	}

	int n_coeffs = coeff_modulus_size * poly_modulus_degree;
	uint64_t** per_site_coeff_selector = new uint64_t*[n_sites];
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		per_site_coeff_selector[i_site] = new uint64_t[n_coeffs + 2];
	}

	long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
	fprintf(stderr, "Coeff selector seed: %lu\n", rand_seed);
	t_rng* rng = new t_rng(rand_seed);

	for (int i_coeff = 0; i_coeff < poly_modulus_degree; i_coeff++)
	{
		// For each coefficient we select a different site to store the secret key coefficient.
		vector<int>* rand_sites = rng->fast_permute_indices(0, n_sites);
		int top_site_i = rand_sites->at(0);
		int* cur_coeff_per_site_selector = new int[n_sites + 2];
		memset(cur_coeff_per_site_selector, 0, sizeof(int) * n_sites);
		cur_coeff_per_site_selector[top_site_i] = 1;

		// The selectors are shared among different RNS moduli.
		for (int i_mod_degree = 0; i_mod_degree < coeff_modulus_size; i_mod_degree++)
		{
			// Assign the random value for each site to the appropriate coefficient position.
			int cur_coeff_i = (i_mod_degree * poly_modulus_degree) + i_coeff;

			// Assign the coefficient selector to each site for this coefficient.
			for (int i_site = 0; i_site < n_sites; i_site++)
			{
				per_site_coeff_selector[i_site][cur_coeff_i] = cur_coeff_per_site_selector[i_site];
			} // i_site loop.
		} // i_mod_degree loop.

		delete rand_sites;
	} // i_coeff loop.

	return(per_site_coeff_selector);
}

// Following uses SEAL's RNG; we do not need this function, it just calls custom normal generator.
uint64_t* SG_generate_smudging_noise_SEAL_RNG(EncryptionParameters context_parms, double smdging_noise_variance)
{
	fprintf(stderr, "Generating gaussian smudging noise with variance %.5f\n", smdging_noise_variance);

	fprintf(stderr, "Generating smdging noise random values.\n");
	double smdging_noise_std_dev = pow(smdging_noise_variance, .5);
	uint64_t* smdging_noise_unit64 = SEAL_Genomics_Get_Normal_Random_Custom(context_parms, smdging_noise_std_dev, 1024 * 1024);

	return(smdging_noise_unit64);
} // SG_generate_smudging_noise_SEAL_RNG function.

uint64_t* SG_generate_smudging_noise(int poly_modulus_degree, int coeff_modulus_size, double smdging_noise_variance, const vector<Modulus> coeff_modulus)
{
	fprintf(stderr, "Generating gaussian smudging noise with variance %.5f for poly_mod_degree=%d, coeff_mod_size=%d\n", smdging_noise_variance, poly_modulus_degree, coeff_modulus_size);

	// Allocate the noise coefficients vector.
	// poly_modulus_degree :: This is the largest polynomial degree, i.e., number of coefficients, this is usually a large number like 8192.
	// coeff_modulus_size :: This is the size of RNS representation usually a small number representing the decomposition of q.
	int n_coeffs = coeff_modulus_size * poly_modulus_degree;
	uint64_t* smdging_noise_levels = new uint64_t[n_coeffs + 2];

	long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
	fprintf(stderr, "Smdgng noise seed: %lu\n", rand_seed);
	t_rng* rng = new t_rng(rand_seed);

	// Sample each coefficient.
	for (int i_coeff = 0; i_coeff < poly_modulus_degree; i_coeff++)
	{
		double cur_rand = (rng->random_gaussian_double_ran3() * pow(smdging_noise_variance, .5));

		int cur_rand_int = (int)(floor(cur_rand));

		// Go over all coefficients of the secret key and set the value at the correct position in the RNS representation.
		// This sets the noise level for each modulus degree.
		for (int i_mod_degree = 0; i_mod_degree < coeff_modulus_size; i_mod_degree++)
		{
			// Assign the random value for each site to the appropriate coefficient position.
			int cur_coeff_i = (i_mod_degree * poly_modulus_degree) + i_coeff;

			int cur_noise_val_int = cur_rand_int;

			// Negatives should be set with respect to the modulus of the current RNS representation's base, q_k.
			uint64_t cur_mod_noise_val = cur_noise_val_int;
			if (cur_noise_val_int < 0)
			{
				if (__DUMP_KEY_SHARING_MSGS__)
				{
					fprintf(stderr, "Noise is negative (%d), modding with respect to %lu\n", cur_noise_val_int, coeff_modulus.at(i_mod_degree).value());
				}

				cur_mod_noise_val = coeff_modulus.at(i_mod_degree).value() + cur_noise_val_int;

				if (__DUMP_KEY_SHARING_MSGS__)
				{
					fprintf(stderr, "Negative noise is now %lu\n", cur_mod_noise_val);
				}
			}

			smdging_noise_levels[cur_coeff_i] = cur_mod_noise_val;
		} // i_mod_degree loop.
	} // i_coeff loop.

	return(smdging_noise_levels);
} // SG_generate_smudging_noise function.

// This takes the key_context_parms directly and generates the site specific noise levels using SEAL's RNGs.
uint64_t** SG_generate_per_site_complementary_noise_random_complementation_SEAL_RNG(int n_sites, EncryptionParameters key_context_parms, int key_noise_var_bit_size, const vector<Modulus> coeff_modulus)
{
	if (n_sites < 2)
	{
		fprintf(stderr, "Need at least 2 sites to generate the noise levels.\n");
		exit(0);
	}

	// Max the key noise variance bit size at 30-bits.
	if (key_noise_var_bit_size > 30)
	{
		key_noise_var_bit_size = 30;
	}

	// Set the coeff numbers, etc.
	int poly_modulus_degree = key_context_parms.poly_modulus_degree();
	int coeff_modulus_size = key_context_parms.coeff_modulus().size();

	// Sample the modulus noise levels for each site.
	fprintf(stderr, "Sampling noise for each site..\n");
	uint64_t** per_site_base_normal_rand_SEAL_RNG_noise = new uint64_t*[n_sites];
	int** per_site_expanded_normal_rand_SEAL_RNG_noise = new int*[n_sites];
	double key_noise_std_dev = pow(2, key_noise_var_bit_size / 2);
	fprintf(stderr, "Sampling noise for each site using key std-dev of %lu..\n", (uint64_t)key_noise_std_dev);
	double max_dev = 10 * key_noise_std_dev;
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		fprintf(stderr, "Sampling noise for site-%d..\n", i_site);
		per_site_base_normal_rand_SEAL_RNG_noise[i_site] = SEAL_Genomics_Get_Normal_Random_Custom(key_context_parms, key_noise_std_dev, max_dev);

		int n_coeff_count = mul_safe(poly_modulus_degree, coeff_modulus_size);

		per_site_expanded_normal_rand_SEAL_RNG_noise[i_site] = new int[n_coeff_count];

		// Go over the noise and expand it by noise variance.
		for (int i_mod = 0; i_mod < (int)key_context_parms.coeff_modulus().size(); i_mod++)
		{
			size_t cur_modulus = key_context_parms.coeff_modulus().at(i_mod).value();

			for (size_t coeff_i = 0; coeff_i < key_context_parms.poly_modulus_degree(); coeff_i++)
			{
				//fprintf(f_op, "%d\t%d\t%lu\t%lu\n", i_mod, (int)coeff_i, buff[i_mod * key_context_parms.poly_modulus_degree() + coeff_i], key_context_parms.coeff_modulus().at(i_mod).value());
				int cur_coeff_arr_i = i_mod * poly_modulus_degree + coeff_i;

				int cur_signed_expanded_noise_level = per_site_base_normal_rand_SEAL_RNG_noise[i_site][cur_coeff_arr_i];

				if (per_site_base_normal_rand_SEAL_RNG_noise[i_site][cur_coeff_arr_i] > (cur_modulus / 2))
				{
					cur_signed_expanded_noise_level = -1 * (int)(cur_modulus - cur_signed_expanded_noise_level);

					if (__DUMP_KEY_SHARING_MSGS__)
					{
						fprintf(stderr, "Noise is negative: %lu / %lu --> %d\n",
							per_site_base_normal_rand_SEAL_RNG_noise[i_site][cur_coeff_arr_i], cur_modulus,
							cur_signed_expanded_noise_level);
					}
				}
				else
				{
					if (__DUMP_KEY_SHARING_MSGS__)
					{
						fprintf(stderr, "Noise is positive: %lu / %lu --> %d\n",
							per_site_base_normal_rand_SEAL_RNG_noise[i_site][cur_coeff_arr_i], cur_modulus,
							cur_signed_expanded_noise_level);
					}
				}

				//cur_signed_expanded_noise_level *= noise_std_dev;
				//fprintf(stderr, "Expanded noise level ")

				// Store the expanded signed noise.
				per_site_expanded_normal_rand_SEAL_RNG_noise[i_site][cur_coeff_arr_i] = cur_signed_expanded_noise_level;
			} // i loop.
		} // i_mod loop.
	} // i_sitel oop.

	fprintf(stderr, "Generating complementary gaussian noise for %d sites with poly_mod_degree=%d, coeff_mod_size=%d\n", n_sites, poly_modulus_degree, coeff_modulus_size);

	// Allocate the noise coefficients vector.
	// poly_modulus_degree :: This is the largest polynomial degree, i.e., number of coefficients, this is usually a large number like 8192.
	// coeff_modulus_size :: This is the size of RNS representation usually a small number representing the decomposition of q.
	int n_coeffs = coeff_modulus_size * poly_modulus_degree;
	uint64_t** per_site_complementary_noise_levels = new uint64_t*[n_sites];
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		per_site_complementary_noise_levels[i_site] = new uint64_t[n_coeffs + 2];
	} // i_site loop.

	long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
	fprintf(stderr, "Comp. noise site selector seed: %lu\n", rand_seed);
	t_rng* rng = new t_rng(rand_seed);

	// Sample each coefficient.
	for (int i_coeff = 0; i_coeff < poly_modulus_degree; i_coeff++)
	{
		vector<int>* per_site_rand_vals = new vector<int>();
		int total_rand = 0;
		for (int i_site = 0; i_site < n_sites - 1; i_site++)
		{
			//double cur_rand = (rng->random_gaussian_double_ran3() * pow(noise_variance, .5));
			//int cur_rand_int = (int)(floor(cur_rand));
			int cur_rand_int = per_site_expanded_normal_rand_SEAL_RNG_noise[i_site][i_coeff];
			per_site_rand_vals->push_back(cur_rand_int);

			total_rand += cur_rand_int;
		} // i_site loop.

		// Add the total negative to the last site.
		per_site_rand_vals->push_back(-1 * total_rand);

		if ((int)(per_site_rand_vals->size()) != n_sites)
		{
			fprintf(stderr, "Sanity check failed: Number of random values do not match the sites: %d / %d\n",
				(int)(per_site_rand_vals->size()), n_sites);
			exit(0);
		}

		// Following randomizes which site gets assigned the total negative of the random value.
		vector<int>* rand_site_indices = rng->fast_permute_indices(0, n_sites);

		// Go over all coefficients of the secret key and set the value at the correct position in the RNS representation.
		for (int i_mod_degree = 0; i_mod_degree < coeff_modulus_size; i_mod_degree++)
		{
			// Assign the random value for each site to the appropriate coefficient position.
			int cur_coeff_i = (i_mod_degree * poly_modulus_degree) + i_coeff;

			for (int i_site = 0; i_site < n_sites; i_site++)
			{
				int cur_site_rand_site_i = rand_site_indices->at(i_site);
				int cur_noise_val_int = per_site_rand_vals->at(cur_site_rand_site_i);

				// Negatives should be set with respect to the modulus of the current RNS representation's base, q_k.
				uint64_t cur_mod_noise_val = cur_noise_val_int;
				if (cur_noise_val_int < 0)
				{
					if (__DUMP_KEY_SHARING_MSGS__)
					{
						fprintf(stderr, "Noise is negative (%d), modding with respect to %lu\n", cur_noise_val_int, coeff_modulus.at(i_mod_degree).value());
					}

					cur_mod_noise_val = coeff_modulus.at(i_mod_degree).value() + cur_noise_val_int;

					if (__DUMP_KEY_SHARING_MSGS__)
					{
						fprintf(stderr, "Negative noise is now %lu\n", cur_mod_noise_val);
					}
				}

				per_site_complementary_noise_levels[i_site][cur_coeff_i] = cur_mod_noise_val;
			} // i_site loop.
		} // i_mod_degree loop.
	} // i_coeff loop.

	return(per_site_complementary_noise_levels);
} // SG_generate_per_site_complementary_noise_random_complementation_SEAL_RNG function.

uint64_t** SG_generate_per_site_complementary_noise_random_complementation(int n_sites, int poly_modulus_degree, int coeff_modulus_size, double noise_variance, const vector<Modulus> coeff_modulus)
{
	if (n_sites < 2)
	{
		fprintf(stderr, "Need at least 2 sites to generate the noise levels.\n");
		exit(0);
	}

	fprintf(stderr, "Generating complementary gaussian noise for %d sites with poly_mod_degree=%d, coeff_mod_size=%d\n", n_sites, poly_modulus_degree, coeff_modulus_size);

	// Allocate the noise coefficients vector.
	// poly_modulus_degree :: This is the largest polynomial degree, i.e., number of coefficients, this is usually a large number like 8192.
	// coeff_modulus_size :: This is the size of RNS representation usually a small number representing the decomposition of q.
	int n_coeffs = coeff_modulus_size * poly_modulus_degree;
	uint64_t** per_site_complementary_noise_levels = new uint64_t*[n_sites];
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		per_site_complementary_noise_levels[i_site] = new uint64_t[n_coeffs + 2];
	} // i_site loop.

	long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
	fprintf(stderr, "Comp. noise seed: %lu\n", rand_seed);
	t_rng* rng = new t_rng(rand_seed);

	// Sample each coefficient.
	for (int i_coeff = 0; i_coeff < poly_modulus_degree; i_coeff++)
	{
		vector<int>* per_site_rand_vals = new vector<int>();
		int total_rand = 0;
		for (int i_site = 0; i_site < n_sites - 1; i_site++)
		{
			double cur_rand = (rng->random_gaussian_double_ran3() * pow(noise_variance, .5));

			int cur_rand_int = (int)(floor(cur_rand));
			per_site_rand_vals->push_back(cur_rand_int);

			total_rand += cur_rand_int;
		} // i_site loop.

		// Add the total negative to the last site.
		per_site_rand_vals->push_back(-1 * total_rand);

		if ((int)(per_site_rand_vals->size()) != n_sites)
		{
			fprintf(stderr, "Sanity check failed: Number of random values do not match the sites: %d / %d\n",
				(int)(per_site_rand_vals->size()), n_sites);
			exit(0);
		}

		// Following randomizes which site gets assigned the total negative of the random value.
		vector<int>* rand_site_indices = rng->fast_permute_indices(0, n_sites);

		// Go over all coefficients of the secret key and set the value at the correct position in the RNS representation.
		for (int i_mod_degree = 0; i_mod_degree < coeff_modulus_size; i_mod_degree++)
		{
			// Assign the random value for each site to the appropriate coefficient position.
			int cur_coeff_i = (i_mod_degree * poly_modulus_degree) + i_coeff;

			for (int i_site = 0; i_site < n_sites; i_site++)
			{
				int cur_site_rand_site_i = rand_site_indices->at(i_site);
				int cur_noise_val_int = per_site_rand_vals->at(cur_site_rand_site_i);

				// Negatives should be set with respect to the modulus of the current RNS representation's base, q_k.
				uint64_t cur_mod_noise_val = cur_noise_val_int;
				if (cur_noise_val_int < 0)
				{
					if (__DUMP_KEY_SHARING_MSGS__)
					{
						fprintf(stderr, "Noise is negative (%d), modding with respect to %lu\n", cur_noise_val_int, coeff_modulus.at(i_mod_degree).value());
					}

					cur_mod_noise_val = coeff_modulus.at(i_mod_degree).value() + cur_noise_val_int;

					if (__DUMP_KEY_SHARING_MSGS__)
					{
						fprintf(stderr, "Negative noise is now %lu\n", cur_mod_noise_val);
					}
				}

				per_site_complementary_noise_levels[i_site][cur_coeff_i] = cur_mod_noise_val;
			} // i_site loop.
		} // i_mod_degree loop.
	} // i_coeff loop.

	return(per_site_complementary_noise_levels);
} // SG_generate_per_site_complementary_noise_random_complementation function.

// This function encapsulates NTT transformation for an array of coefficients.
uint64_t* NTT_transform_coefficients_array(uint64_t* non_NTT_coeffs_array,
	size_t coeff_modulus_size, size_t coeff_count,
	const NTTTables* small_ntt_tables)
{
	uint64_t* NTT_coeffs_array = new uint64_t[coeff_modulus_size * coeff_count + 2];
	memcpy(NTT_coeffs_array, non_NTT_coeffs_array, sizeof(uint64_t) * coeff_modulus_size * coeff_count);

	RNSIter temp_iter(NTT_coeffs_array, coeff_count);
	ntt_negacyclic_harvey(temp_iter, coeff_modulus_size, small_ntt_tables);

	return(NTT_coeffs_array);
}

void NTT_transform_coefficients_array_in_place(uint64_t* non_NTT_coeffs_array,
	size_t coeff_modulus_size, size_t coeff_count,
	const NTTTables* small_ntt_tables)
{
	RNSIter temp_iter(non_NTT_coeffs_array, coeff_count);
	ntt_negacyclic_harvey(temp_iter, coeff_modulus_size, small_ntt_tables);
}

void inverse_NTT_transform_coefficients_array_in_place(uint64_t* NTT_coeffs_array,
	size_t coeff_modulus_size, size_t coeff_count,
	const NTTTables* small_ntt_tables)
{
	// Invert the ntt representation of the pooled key.
	RNSIter temp_iter(NTT_coeffs_array, coeff_count);
	inverse_ntt_negacyclic_harvey(temp_iter, coeff_modulus_size, small_ntt_tables);
}

uint64_t* inverse_NTT_transform_coefficients_array(uint64_t* NTT_coeffs_array,
	size_t coeff_modulus_size, size_t coeff_count,
	const NTTTables* small_ntt_tables)
{
	uint64_t* inv_NTT_coeffs_array = new uint64_t[coeff_modulus_size * coeff_count + 2];
	memcpy(inv_NTT_coeffs_array, NTT_coeffs_array, sizeof(uint64_t) * coeff_modulus_size * coeff_count);

	// Invert the ntt representation of the pooled key.
	RNSIter temp_iter(inv_NTT_coeffs_array, coeff_count);
	inverse_ntt_negacyclic_harvey(temp_iter, coeff_modulus_size, small_ntt_tables);

	return(inv_NTT_coeffs_array);
}

uint64_t* SEAL_Genomics_Get_Unif_Random_Custom(EncryptionParameters parms, int n_bits_per_max_val)
{
	shared_ptr<UniformRandomGenerator> prng = parms.random_generator()->create();

	// Fill the random buffer.
	size_t n_coeffs = mul_safe(parms.coeff_modulus().size(), parms.poly_modulus_degree());
	uint64_t* destination_orig = new uint64_t[n_coeffs];
	uint64_t* destination = destination_orig;

	// Extract encryption parameters
	auto coeff_modulus = parms.coeff_modulus();
	size_t coeff_modulus_size = coeff_modulus.size();
	size_t coeff_count = parms.poly_modulus_degree();
	size_t dest_byte_count = mul_safe(coeff_modulus_size, coeff_count, sizeof(uint64_t));

	constexpr uint64_t max_random = static_cast<uint64_t>(0xFFFFFFFFFFFFFFFFULL);
	uint64_t abs_max_random = static_cast<uint64_t>((((uint64_t)1) << n_bits_per_max_val) - 1);
	fprintf(stderr, "Max random is %.2f bits\n", log(max_random) / log(2.0));

	// Fill the destination buffer with fresh randomness
	prng->generate(dest_byte_count, reinterpret_cast<seal_byte *>(destination));

	for (size_t j = 0; j < coeff_modulus_size; j++)
	{
		auto &modulus = coeff_modulus[j];
		uint64_t max_multiple = max_random - barrett_reduce_64(max_random, modulus) - 1;

		transform(destination, destination + coeff_count, destination, [&](uint64_t rand) {
			// This ensures uniform distribution
			while (rand >= max_multiple)
			{
				prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand));
			}
			//return barrett_reduce_64(rand, modulus);
			return barrett_reduce_64(rand, abs_max_random);
		});
		destination += coeff_count;
	}

	return(destination_orig);
} // SEAL_Genomics_Get_Unif_Random_Custom option.

// This function wraps a call to sample_poly_uniform.
uint64_t* SEAL_Genomics_Get_Unif_Random(EncryptionParameters parms)
{
	// Do sampling here.
	size_t n_coeffs = mul_safe(parms.coeff_modulus().size(), parms.poly_modulus_degree());
	uint64_t* buff = new uint64_t[n_coeffs];

	fprintf(stderr, "SEAL Sampling..\n");
	RNSIter buff_iter(buff, parms.poly_modulus_degree());
	sample_poly_uniform(parms.random_generator()->create(), parms, buff_iter);
	fprintf(stderr, "Finished SEAL Sampling..\n");

	return(buff);
}

uint64_t* SEAL_Genomics_Get_Normal_Random_Custom(EncryptionParameters parms, double std_dev, double max_dev)
{
	shared_ptr<UniformRandomGenerator> prng = parms.random_generator()->create();
	RandomToStandardAdapter engine(prng);
	ClippedNormalDistribution dist(0, std_dev, max_dev);

	// Fill the random buffer.
	size_t n_coeffs = mul_safe(parms.coeff_modulus().size(), parms.poly_modulus_degree());
	uint64_t* buff = new uint64_t[n_coeffs];

	auto coeff_modulus = parms.coeff_modulus();
	size_t coeff_count = parms.poly_modulus_degree();
	size_t coeff_modulus_size = parms.coeff_modulus().size();
	SEAL_ITERATE(iter(buff), coeff_count, [&](auto &I) {
		int64_t noise = static_cast<int64_t>(dist(engine));
		uint64_t flag = static_cast<uint64_t>(-static_cast<int64_t>(noise < 0));
		SEAL_ITERATE(
			iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
			[&](auto J) { *get<0>(J) = static_cast<uint64_t>(noise) + (flag & get<1>(J).value()); });
	});

	return(buff);
}

uint64_t* SEAL_Genomics_Get_Normal_Random(EncryptionParameters parms)
{
	// Do sampling here.
	size_t n_coeffs = mul_safe(parms.coeff_modulus().size(), parms.poly_modulus_degree());
	uint64_t* buff = new uint64_t[n_coeffs];

	fprintf(stderr, "SEAL Sampling..\n");
	RNSIter buff_iter(buff, parms.poly_modulus_degree());
	sample_poly_normal(parms.random_generator()->create(), parms, buff_iter);
	fprintf(stderr, "Finished SEAL Sampling..\n");

	return(buff);
}

// The main function that generates private key shares, pooled public key for a number of sites.
// 0. Understand the parameter selection and scales: https://blog.openmined.org/ckks-explained-part-5-rescaling/
// What are the conditions for correct pooling of data such that we can pool and type of data? (size-2 ciphertexts is the only requirement.)
// 0.25. Perform some operations on the encrypted data and test decryption of it using pooled keys.
// 1. Save and load parameters
// 2. Save and load keys
// 3. Save and load encrypted data 
// 4. Save and load plaintext data
// 5. Put it all together.
void generate_share_private_keys(int n_sites, char* text_params_path,
	bool perform_key_testing,
	char* op_dir)
{
	// We load the parameters from text file.
	//double secret_key_noise_variance = 20;
	//size_t requested_poly_modulus_degree = 8192 * 2;
	//vector<int> coeff_modulus_bit_sizes = { 60, 40, 40, 40, 60 };

	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find text params file @ %s\n", text_params_path);
		exit(0);
	}

	t_text_params* text_params = load_text_params(text_params_path);

	double secret_key_noise_variance = text_params->per_site_key_noise_var_bit_size;
	double scale = pow(2, text_params->scale_n_bits);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Set the key and ciphertext parameters.
	size_t key_coeff_count = parms.poly_modulus_degree();
	size_t key_coeff_modulus_size = parms.coeff_modulus().size();
	auto &key_context_data = *context.key_context_data();

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Setup the plaintext params.
	//auto &pl_context_data = *context.get_context_data(x1_encrypted.parms_id());
	//auto &pl_context_data = *context.get_context_data(x_plain.parms_id());
	auto &pl_context_data = *context.first_context_data();

	// These are the base parameters, in principle, we should not use these for plaintext as they may have different parameters after processing.
	auto &pl_parms = pl_context_data.parms();
	//auto &pl_coeff_modulus = pl_parms.coeff_modulus();
	size_t pl_coeff_count = pl_parms.poly_modulus_degree();
	size_t pl_coeff_modulus_size = pl_parms.coeff_modulus().size();
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	// This is the master keygenerator that generates the pooled secret/public/relin/galois keys.
	KeyGenerator pooled_keygen(context);

	// This is the main pooled secret key that should be protected, noone will have access to this key after running key generation.
	SecretKey pooled_sk = pooled_keygen.secret_key();
	ofstream ofs_original_sk("pooled.secret_key", ios::binary);
	pooled_sk.save(ofs_original_sk);
	ofs_original_sk.close();

	// Generate the pooled public, relin, and galois keys and save them.
	fprintf(stderr, "Generating public, relin, and galois keys and saving.\n");
	char op_fp[1000];
	PublicKey pooled_pk;
	pooled_keygen.create_public_key(pooled_pk);
	sprintf(op_fp, "%s/pooled.public_key", op_dir);
	ofstream pooled_pk_ofst(op_fp, ios::binary);
	pooled_pk.save(pooled_pk_ofst);
	pooled_pk_ofst.close();

	RelinKeys relin_keys;
	pooled_keygen.create_relin_keys(relin_keys);
	sprintf(op_fp, "%s/pooled.relin_keys", op_dir);
	ofstream pooled_relin_ofst(op_fp, ios::binary);
	relin_keys.save(pooled_relin_ofst);
	pooled_relin_ofst.close();

	GaloisKeys galois_keys;
	pooled_keygen.create_galois_keys(galois_keys);
	sprintf(op_fp, "%s/pooled.galois_keys", op_dir);
	ofstream pooled_galois_ofst(op_fp, ios::binary);
	galois_keys.save(pooled_galois_ofst);
	pooled_galois_ofst.close();

	// Now, share the keys among sites.
	//double noise_variance = secret_key_noise_variance;
	//uint64_t** per_site_complementary_noise = SG_generate_per_site_complementary_noise_random_complementation(n_sites, poly_modulus_degree, key_coeff_modulus_size, secret_key_noise_variance, parms.coeff_modulus());

	uint64_t** per_site_complementary_noise = SG_generate_per_site_complementary_noise_random_complementation_SEAL_RNG(n_sites, key_context_data.parms(), secret_key_noise_variance, parms.coeff_modulus());
	uint64_t** per_site_coeff_selector = SG_generate_per_site_sk_coeff_selector(n_sites, key_coeff_count, key_coeff_modulus_size, parms.coeff_modulus());

	// Generate random values for the current site.
	uint64_t* sk_data_inv_ntt = inverse_NTT_transform_coefficients_array(pooled_sk.data().data(),
																		key_coeff_modulus_size, key_coeff_count,
																		key_context_data.small_ntt_tables());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Save the noise and the masking matrix.
	FILE* f_noise_matrix = NULL;
	FILE* f_noise_plus_key_matrix = NULL;
	FILE* f_coeff_selector_matrix = NULL;

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		f_noise_matrix = open_f("per_site_noise_matrix.txt", "w");
		f_noise_plus_key_matrix = open_f("per_site_noise_plus_key_matrix.txt", "w");
		f_coeff_selector_matrix = open_f("per_site_coeff_selector.txt", "w");
	}

	for (int i_coeff = 0; i_coeff < (int)poly_modulus_degree; i_coeff++)
	{
		for (int i_site = 0; i_site < n_sites; i_site++)
		{
			if (__DUMP_KEY_SHARING_MSGS__)
			{
				if (i_site > 0)
				{
					fprintf(f_noise_matrix, "\t");
					fprintf(f_noise_plus_key_matrix, "\t");
				}
			}

			uint64_t cur_sk_data = sk_data_inv_ntt[i_coeff];
			uint64_t cur_noise = per_site_complementary_noise[i_site][i_coeff];

			double non_mod_sk_data = sk_data_inv_ntt[i_coeff];
			if (cur_sk_data >= parms.coeff_modulus().at(0).value() / 2)
			{
				uint64_t neg_magnitude = (parms.coeff_modulus().at(0).value() - cur_sk_data);
				non_mod_sk_data = (-1 * (double)neg_magnitude);
			}

			double non_mod_noise = cur_noise;
			if (cur_noise >= (parms.coeff_modulus().at(0).value() / 2))
			{
				uint64_t neg_magnitude = (parms.coeff_modulus().at(0).value() - cur_noise);
				non_mod_noise = (-1 * (double)neg_magnitude);

				if (__DUMP_KEY_SHARING_MSGS__)
				{
					fprintf(stderr, "Noise is negative, set the magnitude to %lu (%.0f)\n", neg_magnitude, non_mod_noise);
				}
			}

			if (__DUMP_KEY_SHARING_MSGS__)
			{
				fprintf(stderr, "@ coeff %d (modulus: %lu): sk: %lu, noise: %lu; non_mod-sk: %.2f; non_mod-noise: %.2f\n",
					i_coeff,
					parms.coeff_modulus().at(0).value(),
					cur_sk_data, cur_noise,
					non_mod_sk_data, non_mod_noise);
			}

			if (__DUMP_KEY_SHARING_MSGS__)
			{
				fprintf(f_noise_plus_key_matrix, "%.0f", non_mod_noise + non_mod_sk_data);
				fprintf(f_noise_matrix, "%.0f", non_mod_noise);
				fprintf(f_coeff_selector_matrix, "%lu", per_site_coeff_selector[i_site][i_coeff]);
			}
		} // i_site loop.

		if (__DUMP_KEY_SHARING_MSGS__)
		{
			fprintf(f_noise_matrix, "\n");
			fprintf(f_noise_plus_key_matrix, "\n");
			fprintf(f_coeff_selector_matrix, "\n");
		}
	} // i_coeff loop.

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		close_f(f_noise_plus_key_matrix, "per_site_noise_plus_key_matrix.txt");
		close_f(f_noise_matrix, "per_site_noise_matrix.txt");
		close_f(f_coeff_selector_matrix, "per_site_coeff_selector.txt");
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector<uint64_t*>* per_site_noisy_keys = new vector<uint64_t*>();
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		// For this site, we have the key masking pattern and the noise levels.
		uint64_t* cur_site_noise_levels = per_site_complementary_noise[i_site];
		uint64_t* cur_site_coeff_selector = per_site_coeff_selector[i_site];

		// Allocate the site's noisy key.
		uint64_t* cur_site_noisy_key = new uint64_t[key_coeff_modulus_size * poly_modulus_degree + 2];
		memset(cur_site_noisy_key, 0, sizeof(uint64_t) * key_coeff_modulus_size * key_coeff_count);

		// Allocate and generate the masked key for the current site.
		uint64_t* cur_site_masked_key = new uint64_t[key_coeff_modulus_size * key_coeff_count + 2];
		memcpy(cur_site_masked_key, sk_data_inv_ntt, sizeof(uint64_t) * key_coeff_modulus_size * key_coeff_count);
		for (int i_coeff = 0; i_coeff < (int)(key_coeff_modulus_size * key_coeff_count); i_coeff++)
		{
			if (cur_site_coeff_selector[i_coeff] != 0 &&
				cur_site_coeff_selector[i_coeff] != 1)
			{
				fprintf(stderr, "Sanity check failed, masking array is not as expected: %d. coeff is %lu\n", i_coeff, cur_site_coeff_selector[i_coeff]);
				exit(0);
			}

			cur_site_masked_key[i_coeff] = cur_site_masked_key[i_coeff] * cur_site_coeff_selector[i_coeff];
		} // i_coeff loop.

		// Start processing each modulus.
		for (int i_mod = 0; i_mod < (int)key_coeff_modulus_size; i_mod++)
		{
			uint64_t* cur_mod_sk_coeffs = cur_site_masked_key + i_mod * poly_modulus_degree;
			uint64_t* cur_mod_site_noise_coeffs = cur_site_noise_levels + i_mod * poly_modulus_degree;
			uint64_t* cur_mod_site_noisy_key_coeffs = cur_site_noisy_key + i_mod * poly_modulus_degree;

			// Write the first 5 coefficients.
			fprintf(stderr, "i_mod %d; Before noise addition:\n", i_mod);
			for (int i_coeff = 0; i_coeff < 5; i_coeff++)
			{
				fprintf(stderr, "@ %d: %lu + %lu mod(%lu)\n", i_coeff,
					cur_mod_sk_coeffs[i_coeff],
					cur_mod_site_noise_coeffs[i_coeff],
					parms.coeff_modulus().at(i_mod).value());
			} // i_coeff loop.

			// Write the key values to make sure they are correctly converted from NTT representation: This operation adds the key to all sites.
			add_poly_coeffmod(cur_mod_sk_coeffs, cur_mod_site_noise_coeffs, key_coeff_count, parms.coeff_modulus().at(i_mod), cur_mod_site_noisy_key_coeffs);

			// Write the first 5 coefficients after addition.
			fprintf(stderr, "i_mod %d; After noise addition:\n", i_mod);
			for (int i_coeff = 0; i_coeff < 5; i_coeff++)
			{
				fprintf(stderr, "@ %d: %lu + %lu = %lu mod(%lu)\n", i_coeff,
					cur_mod_sk_coeffs[i_coeff],
					cur_mod_site_noise_coeffs[i_coeff],
					cur_mod_site_noisy_key_coeffs[i_coeff],
					parms.coeff_modulus().at(i_mod).value());
			} // i_coeff loop.
		} // i_mod loop.

		// Add the current noisy key to the list of keys for each site.
		per_site_noisy_keys->push_back(cur_site_noisy_key);

		// Allocate and save the private keys at this point. 
		uint64_t* ntt_transformed_key_coeffs = NTT_transform_coefficients_array(cur_site_noisy_key, key_coeff_modulus_size, key_coeff_count, key_context_data.small_ntt_tables());
		SecretKey cur_site_key = get_secret_key_per_NTT_coeff_array(ntt_transformed_key_coeffs, key_coeff_modulus_size, key_coeff_count);
		cur_site_key.parms_id() = key_context_data.parms_id();

		// Save the file.
		char op_fp[1000];
		sprintf(op_fp, "%s/site_%d.secret_key", op_dir, i_site);
		ofstream key_ofst(op_fp, ios::binary);
		cur_site_key.save(key_ofst);
		key_ofst.close();
	} // i_site loop.

	// We have generated all the keys. After this point, it is about validating the keys.
	if (!perform_key_testing)
	{
		return;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start testing the keys by pooling keys and pooling decryptions.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Do a sanity check by pooling the keys and making sure they add up to the pooled secret key.
	fprintf(stderr, "Validating that the per-site keys sum up to the pooled secret key.\n");
	uint64_t* pooled_noisy_secret_key = pool_nonNTT_secret_keys_from_sites(parms, n_sites,
		per_site_noisy_keys,
		key_coeff_modulus_size, key_coeff_count,
		sk_data_inv_ntt);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Generated %d noisy keys for the sites, testing partial decryption among sites.\n", (int)(per_site_noisy_keys->size()));

	CKKSEncoder encoder(context);
	Encryptor encryptor(context, pooled_pk);
	Evaluator evaluator(context);

	size_t slot_count = encoder.slot_count();

	// Generate some data to encrypt and test.
	vector<double> input = generate_random_data_for_testing(slot_count);

	////////////////////////////////////////////////////////////////////////////////////
	// Following is the data processing portion that analyzes the pooled encrypted data.
	////////////////////////////////////////////////////////////////////////////////////
	Plaintext x_plain;
	print_line(__LINE__);
	cout << "Encode input vectors." << endl;
	encoder.encode(input, scale, x_plain); // This encodes the plain and sets its parameters -- parms.
	Ciphertext x1_encrypted;
	encryptor.encrypt(x_plain, x1_encrypted);

	evaluator.multiply(x1_encrypted, x1_encrypted, x1_encrypted);

	evaluator.relinearize_inplace(x1_encrypted, relin_keys);
	evaluator.rescale_to_next_inplace(x1_encrypted);

	evaluator.multiply(x1_encrypted, x1_encrypted, x1_encrypted);

	evaluator.relinearize_inplace(x1_encrypted, relin_keys);
	evaluator.rescale_to_next_inplace(x1_encrypted);

	evaluator.add_inplace(x1_encrypted, x1_encrypted);

	//decryptor.invariant_noise_budget(x1_encrypted)

	////////////////////////////////////////////////////////////////////////////////////
	// x1_encrypted contains the ciphertext that we will test.
	////////////////////////////////////////////////////////////////////////////////////

	// This is an important value since the ciphertext may not be on the same scale as a fresh PlainText value concordant with initial parms.
	size_t ct_coeff_modulus_size = x1_encrypted.coeff_modulus_size();
	if (context.get_context_data(x1_encrypted.parms_id())->parms().coeff_modulus().size() != ct_coeff_modulus_size)
	{
		fprintf(stderr, "Sanity check failed @ %s(%d); %d/%d\n",
			__FILE__, __LINE__,
			(int)(context.get_context_data(x1_encrypted.parms_id())->parms().coeff_modulus().size()),
			(int)ct_coeff_modulus_size);

		exit(0);
	}

	// Create a key using one of the site's noisy keys.	
	uint64_t* pooled_result_ptr = new uint64_t[ct_coeff_modulus_size * pl_coeff_count];
	memset(pooled_result_ptr, 0, ct_coeff_modulus_size * pl_coeff_count * sizeof(uint64_t));

	// NTT transform, ct and pt are at the same level so we just use the NTTTables of the first level, which is the data level.
	NTT_transform_coefficients_array_in_place(pooled_result_ptr, ct_coeff_modulus_size, pl_coeff_count, pl_context_data.small_ntt_tables());

	// Go over all sites and pool the partial decryptions.
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		// At this point, we have the noisy key for this site in cur_site_noisy_key array. Now, invert the noist key to NTT representation for this site.
		fprintf(stderr, "Converting noisy key for site %d to NTT representation.\n", i_site);

		// NTT transform the current key.
		NTT_transform_coefficients_array_in_place(per_site_noisy_keys->at(i_site), key_coeff_modulus_size, key_coeff_count, key_context_data.small_ntt_tables());

		SecretKey cur_site_sk = get_secret_key_per_NTT_coeff_array(per_site_noisy_keys->at(i_site), key_coeff_modulus_size, key_coeff_count);
		cur_site_sk.parms_id() = key_context_data.parms_id();

		// Allocate the decryptor for the current site.
		Decryptor decryptor(context, cur_site_sk);
		fprintf(stderr, "%s(%d)\n", __FILE__, __LINE__);

		// Perform partial decryption using this key.
		Plaintext cur_site_plain_result;
		decryptor.decrypt(x1_encrypted, cur_site_plain_result);

		auto &pl_parms = context.get_context_data(cur_site_plain_result.parms_id())->parms();
		//auto &pl_coeff_modulus = pl_parms.coeff_modulus();
		size_t pl_coeff_count = pl_parms.poly_modulus_degree();
		size_t pl_coeff_modulus_size = pl_parms.coeff_modulus().size();

		// Add the results from the current partial decryption.
		fprintf(stderr, "%s(%d): pl_coeff_cnt: %d; pl_coeff_mod_size: %d, ciphertext_mod_size: %d;; plaintext coeff length: %d ;; plaintext is_ntt: %d; ciphertext size: %d\n",
			__FILE__, __LINE__,
			(int)(pl_coeff_count),
			(int)(pl_parms.coeff_modulus().size()),
			(int)(x1_encrypted.coeff_modulus_size()),
			(int)(cur_site_plain_result.coeff_count()),
			(int)(cur_site_plain_result.is_ntt_form()),
			(int)(x1_encrypted.size()));

		// Following is the plaintext pooling loop.
		for (int i_mod = 0; i_mod < (int)pl_coeff_modulus_size; i_mod++)
		{
			fprintf(stderr, "Pooling %d. modulus: %lu\n", i_mod, pl_parms.coeff_modulus().at(i_mod).value());
			uint64_t* cur_mod_pooled_result_ptr = pooled_result_ptr + i_mod * pl_coeff_count;
			uint64_t* cur_mod_cur_site_plain_result_ptr = cur_site_plain_result.data(i_mod * pl_coeff_count);
			add_poly_coeffmod(cur_mod_pooled_result_ptr, cur_mod_cur_site_plain_result_ptr, pl_coeff_count, pl_parms.coeff_modulus().at(i_mod), cur_mod_pooled_result_ptr);

			// c0 component must be added only once. Following removes the c0 component from the pooled decryption.
			if (i_site != 0)
			{
				sub_poly_coeffmod(cur_mod_pooled_result_ptr, x1_encrypted.data(0) + i_mod * pl_coeff_count, pl_coeff_count, pl_parms.coeff_modulus().at(i_mod), cur_mod_pooled_result_ptr);
			}
		} // i_mod loop.
	} // i_site loop.

	// Invert the plaintext values for the pooled plaintext values for comparison.
	inverse_NTT_transform_coefficients_array_in_place(pooled_result_ptr, ct_coeff_modulus_size, pl_coeff_count, pl_context_data.small_ntt_tables());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Decrypt using the pooled key: First perform NTT transformation on the key.
	NTT_transform_coefficients_array_in_place(pooled_noisy_secret_key, key_coeff_modulus_size, key_coeff_count, key_context_data.small_ntt_tables());

	// Get the key from the pooled noisy key shares array.
	SecretKey pooled_noisy_sk_obj = get_secret_key_per_NTT_coeff_array(pooled_noisy_secret_key, key_coeff_modulus_size, key_coeff_count);
	pooled_noisy_sk_obj.parms_id() = key_context_data.parms_id();

	// Save the pooled secret key.
	ofstream ofs_pooled_sk("pooled_noisy.secret_key", ios::binary);
	pooled_noisy_sk_obj.save(ofs_pooled_sk);
	ofs_pooled_sk.close();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Instantiate the decryptor using the pooled sk.
	Decryptor pooled_key_decryptor(context, pooled_noisy_sk_obj);

	// Decrypt the data.
	Plaintext pooled_key_decrypted_result;
	pooled_key_decryptor.decrypt(x1_encrypted, pooled_key_decrypted_result);

	vector<double> pooled_key_decreypted_result_array;
	encoder.decode(pooled_key_decrypted_result, pooled_key_decreypted_result_array);
	fprintf(stderr, "Following is from pooled s.k. decryption, should be correct:\n");
	print_vector(pooled_key_decreypted_result_array, 3, 7);

	// Now invert the pooled-key-decrypted plaintext for comparison.
	auto& dec_context_data = *(context.get_context_data(pooled_key_decrypted_result.parms_id()));
	//inverse_NTT_transform_coefficients_array_in_place(pooled_key_decrypted_result.data(), pl_coeff_modulus_size, pl_coeff_count, pl_context_data.small_ntt_tables());
	inverse_NTT_transform_coefficients_array_in_place(pooled_key_decrypted_result.data(),
		dec_context_data.parms().coeff_modulus().size(), // This have to come from the actual decrypted plaintext because it is processed above and it is at a different level.
		pl_coeff_count,
		//dec_context_data.parms().poly_modulus_degree(),
		pl_context_data.small_ntt_tables());
	//dec_context_data.small_ntt_tables());

	fprintf(stderr, "Plaintext parameters: Original pt stats: coeff_mod_size: %d; poly_mod_degree: %d\ncoeff_mod_size: %d; poly_mod_degree: %d; chain_index: %d\n",
		(int)(pl_coeff_modulus_size), (int)(pl_coeff_count),
		(int)(dec_context_data.parms().coeff_modulus().size()), (int)(dec_context_data.parms().poly_modulus_degree()),
		(int)(dec_context_data.chain_index()));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Comparing results from pooled key-based decryption and from the pooled partial decryptions:\n");
	for (int i_mod = 0; i_mod < (int)ct_coeff_modulus_size; i_mod++)
	{
		uint64_t* pooled_key_decrypted_res_ptr = pooled_key_decrypted_result.data() + i_mod * pl_coeff_count;
		uint64_t* pooled_partial_dec_res_ptr = pooled_result_ptr + i_mod * pl_coeff_count;
		for (int i_coeff = 0; i_coeff < (int)pl_coeff_count; i_coeff++)
		{
			//fprintf(stderr, "coeff_i: %d: %lu\t%lu\n", i_coeff, pooled_key_decrypted_res_ptr[i_coeff], pooled_partial_dec_res_ptr[i_coeff]);
			if (pooled_key_decrypted_res_ptr[i_coeff] != pooled_partial_dec_res_ptr[i_coeff])
			{
				fprintf(stderr, "Coeff counts differ @ i_mod %d; i_coeff %d", i_mod, i_coeff);
				exit(0);
			}
		} // i_coeff loop.
	} // i_mod loop.

	// If we are here, the coefficients from the pooling matches the decrpytion using pooled secret key.
	fprintf(stderr, "Plaintext coeffs are matching from partial decryptions.\n");
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Checking final decryption and decoding..\n");

	// Do a final ntt transformation on the pooled partial decryptions since we just inverted it for comparison.
	NTT_transform_coefficients_array_in_place(pooled_result_ptr, ct_coeff_modulus_size, pl_coeff_count, pl_context_data.small_ntt_tables());

	size_t rns_poly_uint64_count = mul_safe(pl_coeff_count, ct_coeff_modulus_size);

	fprintf(stderr, "Plaintext params: %d, %d; %d\n", (int)pl_coeff_count, (int)ct_coeff_modulus_size, (int)rns_poly_uint64_count);

	fprintf(stderr, "Plaintext modulus values: %lu\n%lu\n%lu\n",
		pl_parms.coeff_modulus().at(0).value(),
		pl_parms.coeff_modulus().at(1).value(),
		pl_parms.coeff_modulus().at(2).value());

	fprintf(stderr, "Ciphertext modulus values: %lu\n%lu\n%lu\n%lu\n",
		parms.coeff_modulus().at(0).value(),
		parms.coeff_modulus().at(1).value(),
		parms.coeff_modulus().at(2).value(),
		parms.coeff_modulus().at(3).value());

	fprintf(stderr, "%s(%d)\n", __FILE__, __LINE__);

	Plaintext final_pooled_result;
	final_pooled_result.resize(rns_poly_uint64_count);

	fprintf(stderr, "%s(%d)\n", __FILE__, __LINE__);
	memcpy(final_pooled_result.data(), pooled_result_ptr, sizeof(uint64_t) * (rns_poly_uint64_count));

	fprintf(stderr, "%s(%d): Scale: %.3f\n", __FILE__, __LINE__, log(x1_encrypted.scale()) / log(2.0));
	final_pooled_result.parms_id() = x1_encrypted.parms_id();
	final_pooled_result.scale() = x1_encrypted.scale();

	fprintf(stderr, "%d\n", final_pooled_result.is_ntt_form());

	vector<double> result;
	encoder.decode(final_pooled_result, result);
	cout << "    + Computed result ...... Correct." << endl;
	print_vector(result, 3, 7);
}

/*
There are are two types of data encoding and encryption approaches :
1. Sample - wide CTs(Imputation)
2. Feature - wide CTs(Kinship)

. Sample - wide ciphertexts(CTs) should contain sample size information, and the number of ciphertexts per feature so that evaluation server can easily process this info.each sample - wide CT holds an integral number of individuals(slots).Number of CTs per feature can be computed from sample size and CKKS parameter settings but it is good to include this in encrypted data redundantly to make sure it is concordant among encryptor and evaluator.

. In principle, these can be encrypted using the same function by transposing (using for example, multi_columns_processing_tools) the corresponding data matrix.

. Evaluator functions, however, differs on a case-by - case basis.This function uses the information in encrypted data and generates resulting CTs, which potentially have different scales and decomp - levels.This is already saved with the CT in the output so it is not necessary to store them explicitly anywhere, as this is not trivial.

. Decryptor does not care about sample sizes since results are summary statistics at the end.The CTs that hold the results have different formats.

. Pooling step, similar to decryptor, does not care about sample sizes etc since it is basically pooling of the partial decryptions from sites.This step just takes result PTs, and pools them up.
*/
void encrypt_data_per_pooled_public_key(char* data_file_path, char* text_params_path, char* pooled_pk_path, char* op_fp)
{
	// This is the length of the processing vectors that are loaded from the data that will be used to 
	// divide into slot_count arrays and encrypted into ciphertexts.
	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_pk_path))
	{
		fprintf(stderr, "Could not find public key @ %s\n", pooled_pk_path);
		exit(0);
	}

	// Instantiate and load the public key.
	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_pk_path);
	}

	ifstream pk_ifst(pooled_pk_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, pk_ifst);
	pk_ifst.close();

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	fprintf(stderr, "Instantiating encoder and encryptor.\n");
	CKKSEncoder encoder(context);
	Encryptor encryptor(context, pooled_pk);

	int n_values_per_cipher = encoder.slot_count();
	fprintf(stderr, "%d values per cipher.\n", n_values_per_cipher);

	int l_processing_vectors = 0;
	vector<double*>* data_vectors = load_data_processing_vectors_for_encryption(data_file_path, l_processing_vectors);

	fprintf(stderr, "Loaded %d data vectors of length %d\n", (int)(data_vectors->size()), l_processing_vectors);

	// Load the data file and start encrypting.
	fprintf(stderr, "Encrypting by %d data elements per ciphertext and saving to %s.\n", n_values_per_cipher, op_fp);

	// Open the encrypted data output.
	ofstream ofs_encrypted(op_fp, ios::binary);

	// Write the (1) vector lengths, (2) the number of vectors, (3) number of ciphertexts per vector for this site to the encrypted file first.
	int n_ciphertexts_per_vector = int(ceil((double)l_processing_vectors / n_values_per_cipher));

	// Write the counts, first.
	int n_data_vectors = (int)(data_vectors->size());
	ofs_encrypted.write(reinterpret_cast<const char *>(&l_processing_vectors), sizeof(int));
	ofs_encrypted.write(reinterpret_cast<const char *>(&(n_data_vectors)), sizeof(int));
	ofs_encrypted.write(reinterpret_cast<const char *>(&n_ciphertexts_per_vector), sizeof(int));

	fprintf(stderr, "%d ciphertexts per data vector.\n", n_ciphertexts_per_vector);

	// Now start saving the encrypted data.
	for (int i_data_vector = 0; i_data_vector < (int)(data_vectors->size()); i_data_vector++)
	{
		double* cur_data_vec = data_vectors->at(i_data_vector);

		int data_pt_i = 0;
		int n_cts_per_cur_vec = 0;
		while (data_pt_i < l_processing_vectors)
		{
			fprintf(stderr, "Encrypting data vector %d[%d-%d]\n", i_data_vector, data_pt_i, data_pt_i + n_values_per_cipher - 1);

			vector<double> cur_data_pts;
			for (int cur_cipher_data_pt_i = 0; cur_cipher_data_pt_i < n_values_per_cipher; cur_cipher_data_pt_i++)
			{
				if (data_pt_i < l_processing_vectors)
				{
					cur_data_pts.push_back(cur_data_vec[data_pt_i]);
				}
				else
				{
					cur_data_pts.push_back(0);
				}
				data_pt_i++;
			} // cur_cipher_data_pt_i loop.

			// Encrypt the data points and save to file.
			Plaintext cur_plain;

			// Encode the vector.
			encoder.encode(cur_data_pts, scale, cur_plain);

			// Encrypt the plaintext.
			Ciphertext ct;
			encryptor.encrypt(cur_plain, ct);

			// Save the encrypted ciphertext.
			ct.save(ofs_encrypted);

			n_cts_per_cur_vec++;
		} // data_pt_i loop.

		if (n_ciphertexts_per_vector != n_cts_per_cur_vec)
		{
			fprintf(stderr, "Sanity check failed: # of ciphertexts per vector is not as expected for %d. data vector: %d / %d\n",
				i_data_vector,
				n_ciphertexts_per_vector,
				n_cts_per_cur_vec);
		}
	} // i_data_vector loop.
	ofs_encrypted.close();

	fprintf(stderr, "Finished encryption..");
}

// This function is used for testing the full decryption and rebuilding of the data matrix;; Note that this serves to test the extra parameters that are added at the beginning of the data.
void full_decrypt_encrypted_data_per_pooled_key(char* enc_data_matrix_path, char* text_params_path, char* pooled_secret_key_path, char* op_fp)
{
	// Setup context.
	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	//double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_secret_key_path))
	{
		fprintf(stderr, "Could not find the pooled secret key @ %s\n", pooled_secret_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_secret_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_secret_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
	CKKSEncoder encoder(context);

	// Start loading the ciphertexts, decrypting and saving / formatting them.
	ifstream ifs_enc_data_matrix(enc_data_matrix_path, ios::binary);

	int l_data_vector;
	int n_data_vectors;
	int n_ciphertexts_per_vector;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&l_data_vector), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_data_vectors), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_ciphertexts_per_vector), sizeof(int));

	fprintf(stderr, "Decrypting %d data vectors of length %d by reading %d ciphertexts per vector and saving to %s.\n", n_data_vectors, l_data_vector, n_ciphertexts_per_vector, op_fp);

	FILE* f_plain_op = open_f(op_fp, "w");
	if (f_plain_op == NULL)
	{
		fprintf(stderr, "Could not open %s for writing..\n", op_fp);
		exit(0);
	}

	// Write the header.
	fprintf(f_plain_op, "#COL_ID");
	for (int i_data_col = 0; i_data_col < l_data_vector; i_data_col++)
	{
		fprintf(f_plain_op, "\tCOL_%d", i_data_col);
	} // i_data_col loop.
	fprintf(f_plain_op, "\n");

	// Start decrypting.
	for (int data_vec_i = 0; data_vec_i < n_data_vectors; data_vec_i++)
	{
		fprintf(stderr, "Decrypting %d/%d. vector.\n", data_vec_i, n_data_vectors);
		fprintf(f_plain_op, "Feat_%d", data_vec_i);

		int vector_pt_i = 0;
		for (int cipher_i = 0; cipher_i < n_ciphertexts_per_vector; cipher_i++)
		{
			Ciphertext cur_vec_ct;
			cur_vec_ct.load(context, ifs_enc_data_matrix);

			Plaintext cur_vec_pl;
			decryptor.decrypt(cur_vec_ct, cur_vec_pl);

			vector<double> cur_vec_data;
			encoder.decode(cur_vec_pl, cur_vec_data);

			for (int i_data_pt = 0; i_data_pt < (int)(cur_vec_data.size()); i_data_pt++)
			{
				if (vector_pt_i < l_data_vector)
				{
					fprintf(f_plain_op, "\t%.8f", cur_vec_data.at(i_data_pt));
				}

				vector_pt_i++;
			} // i_data_pt loop.
		} // cipher_i loop.

		fprintf(f_plain_op, "\n");
	} // data_vec_i loop.
}

// site0 flag is necessary for keeping track of c0 term in pooling step.
void partial_decrypt_data_per_noisy_key(char* encrypted_data_path, bool is_site0, char* text_params_path, char* noisy_key_path, double smdg_noise_variance, char* partial_decrypt_data_path)
{
	fprintf(stderr, "Partial decrypting with noisy secret key at %s for site (is_site0=%d), using smdg-ing noise bit size %d and saving to %s\n",
		noisy_key_path, is_site0, (int)smdg_noise_variance, partial_decrypt_data_path);

	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	//double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Load the secret (noisy) key.
	fprintf(stderr, "Loading noisy secret key from %s\n", noisy_key_path);
	ifstream ifs_noisy_key(noisy_key_path, ios::binary);
	SecretKey dec_sk;
	dec_sk.load(context, ifs_noisy_key);
	ifs_noisy_key.close();
	fprintf(stderr, "Loaded noisy secret key.\n");

	// Instantiate the decryptor.
	Decryptor decryptor(context, dec_sk);

	// Just load the ciphertexts from the results and decrypt them with this site's noisy key.
	ifstream ifs_enc_data_matrix(encrypted_data_path, ios::binary);
	ofstream ofs_dec_data_matrix(partial_decrypt_data_path, ios::binary);

	// Read the encrypted data matrix parameters.
	int l_data_vector;
	int n_data_vectors;
	int n_ciphertexts_per_vector;

	// Read the matrix size parameters.
	// l_data_vector is the number of elements in each vector
	// n_data_vectors is the number of vectors
	// n_ciphertexts_per_vector is the number of ct's that we saved for each vector.
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&l_data_vector), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_data_vectors), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_ciphertexts_per_vector), sizeof(int));

	// Write the parameters exactly to the partially decrypted file.
	ofs_dec_data_matrix.write(reinterpret_cast<char *>(&l_data_vector), sizeof(int));
	ofs_dec_data_matrix.write(reinterpret_cast<char *>(&n_data_vectors), sizeof(int));
	ofs_dec_data_matrix.write(reinterpret_cast<char *>(&n_ciphertexts_per_vector), sizeof(int));

	fprintf(stderr, "Decrypting %d data vectors of length %d by reading %d ciphertexts per vector and saving to %s.\n", n_data_vectors, l_data_vector, n_ciphertexts_per_vector, partial_decrypt_data_path);

	// Start loading ciphertexts and decrypt them, then save to the partial decrypted results. For site0, include the c0 term in plaintexts, for others, subtract them from decrypted data.
	for (int i_vector = 0; i_vector < n_data_vectors; i_vector++)
	{
		for (int ct_i = 0; ct_i < n_ciphertexts_per_vector; ct_i++)
		{
			Ciphertext cur_ciphertext;
			cur_ciphertext.load(context, ifs_enc_data_matrix);

			// Decrypt the current ciphertext.
			Plaintext cur_dec_plaintext;
			decryptor.decrypt(cur_ciphertext, cur_dec_plaintext);

			// If this is not the site0, remove the c0 value.
			if (!is_site0)
			{
				size_t ct_coeff_modulus_size = cur_ciphertext.coeff_modulus_size();
				size_t ct_coeff_count = cur_ciphertext.poly_modulus_degree();
				auto &ct_context_data = (*context.get_context_data(cur_ciphertext.parms_id()));
				auto &ct_parms = ct_context_data.parms();

				// 
				// Following is for sanity check; however, note that this may not hold in case we would like to compute custom functions as the modulus may be switched.
				// 
				//auto &pl_context_data = *context.first_context_data(); // Plain parameter context does not have to be at this modulus all the time!
				//auto &pl_parms = pl_context_data.parms();

				//if (pl_parms.coeff_modulus().size() != ct_parms.coeff_modulus().size() ||
				//	pl_parms.poly_modulus_degree() != ct_parms.poly_modulus_degree())
				//{
				//	fprintf(stderr, "Plain vs ciphertext parameters mismatch: %d, %d vs %d, %d\n",
				//			pl_parms.coeff_modulus().size(), ct_parms.coeff_modulus().size(),
				//			pl_parms.poly_modulus_degree(), ct_parms.poly_modulus_degree());

				//	exit(0);
				//}

				fprintf(stderr, "Ciphertext params for vector=%d, ct_i=%d:\ncoeff_modulus_size: %d\n", i_vector, ct_i, (int)(ct_parms.coeff_modulus().size()));

				// We have to remove c0 by going over each decomposition level; the levels are coming from each ciphertext individually.
				for (int i_mod = 0; i_mod < (int)ct_coeff_modulus_size; i_mod++)
				{
					uint64_t* cur_mod_decr_res_ptr = cur_dec_plaintext.data() + i_mod * ct_coeff_count;

					// c0 component must be added only once. Following removes the c0 component from the pooled decryption.
					sub_poly_coeffmod(cur_mod_decr_res_ptr,
						cur_ciphertext.data(0) + i_mod * ct_coeff_count,
						ct_coeff_count, ct_parms.coeff_modulus().at(i_mod),
						cur_mod_decr_res_ptr);
				} // i_mod loop.
			}

			//////////////////////////////////////////////////////////////////////
			// TODO: Do we need a smudging noise term here?
			//////////////////////////////////////////////////////////////////////

			// Save the current partial decrypted data to the output stream.
			cur_dec_plaintext.save(ofs_dec_data_matrix);
		} // ct_i loop.
	} // i_vector loop.

	// Close the files.
	ofs_dec_data_matrix.close();
	ifs_enc_data_matrix.close();
} //  decrypt_data_per_noisy_key

void collaborative_pool_partial_decrypted_plaintext_results(char* partial_decrypted_paths_list_path, char* text_params_path, char* op_fp)
{
	fprintf(stderr, "Pooling partially decrypted data (%s) and saving to %s\n", partial_decrypted_paths_list_path, op_fp);

	t_text_params* text_params = load_text_params(text_params_path);

	vector<int> coeff_modulus_bit_sizes;
	for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
	{
		coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
	} // i_dec loop.

	EncryptionParameters parms(scheme_type::ckks);

	size_t poly_modulus_degree = text_params->poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);
	parms.set_coeff_modulus(CoeffModulus::Create(
		poly_modulus_degree, coeff_modulus_bit_sizes));

	//double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_KEY_SHARING_MSGS__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);
	CKKSEncoder encoder(context);

	// Load the results from sites, pool them and write the final output.
	vector<char*>* per_site_partial_decrypted_paths = buffer_file(partial_decrypted_paths_list_path);
	for (int i_site = 0; i_site < (int)(per_site_partial_decrypted_paths->size()); i_site++)
	{
		if (!check_file(per_site_partial_decrypted_paths->at(i_site)))
		{
			fprintf(stderr, "Could not find site %d's partial decrypted data @ %s\n", i_site, per_site_partial_decrypted_paths->at(i_site));
			exit(0);
		}
	} // i_site loop.

	int n_sites = (int)(per_site_partial_decrypted_paths->size());

	// We have to read the following value from somewhere, possibly the partial decrypted data.
	vector<int>* per_site_l_data_vector = new vector<int>();
	vector<int>* per_site_n_vectors = new vector<int>();
	vector<int>* per_site_n_ct_per_vector = new vector<int>();

	// These three parameters have to match exactly while pooling.
	size_t overall_l_data_vector = 0;
	size_t overall_n_vectors = 0;
	size_t overall_n_ct_per_vector = 0;

	// Load the matrix size parameters for each site.
	vector<ifstream*>* per_site_part_dec_data_ifs = new vector<ifstream*>();
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		ifstream* cur_file_ptr = new ifstream();
		cur_file_ptr->open(per_site_partial_decrypted_paths->at(i_site), ios::binary);
		per_site_part_dec_data_ifs->push_back(cur_file_ptr);

		int l_data_vector;
		int n_vectors;
		int n_ciphertexts_per_vector;

		cur_file_ptr->read(reinterpret_cast<char *>(&l_data_vector), sizeof(int));
		cur_file_ptr->read(reinterpret_cast<char *>(&n_vectors), sizeof(int));
		cur_file_ptr->read(reinterpret_cast<char *>(&n_ciphertexts_per_vector), sizeof(int));

		fprintf(stderr, "i_site: %d: l_data_vector: %d; n_vectors: %d; n_cts_per_vector: %d \n", i_site,
			l_data_vector, n_vectors, n_ciphertexts_per_vector);

		per_site_l_data_vector->push_back(l_data_vector);
		per_site_n_vectors->push_back(n_vectors);
		per_site_n_ct_per_vector->push_back(n_ciphertexts_per_vector);

		if ((int)overall_l_data_vector == 0)
		{
			overall_l_data_vector = l_data_vector;
			overall_n_vectors = n_vectors;
			overall_n_ct_per_vector = n_ciphertexts_per_vector;
		}
		else if (l_data_vector != (int)overall_l_data_vector ||
			n_vectors != (int)overall_n_vectors ||
			n_ciphertexts_per_vector != (int)overall_n_ct_per_vector)
		{
			fprintf(stderr, "Site %d has a different matrix format for pooling: (%d, %d, %d) vs (%d, %d, %d)\n",
				i_site,
				(int)l_data_vector, n_vectors, n_ciphertexts_per_vector,
				(int)overall_l_data_vector, (int)overall_n_vectors, (int)overall_n_ct_per_vector);

			exit(0);
		}
	} // i_site loop.

	// Start reading each plaintext from the sites and pool them, decode and write to final output file.
	FILE* f_op = open_f(op_fp, "w");

	// Write the header.
	fprintf(f_op, "#COL_ID");
	for (int i_data_pt = 0; i_data_pt < (int)overall_l_data_vector; i_data_pt++)
	{
		fprintf(f_op, "\tCOL_%d", i_data_pt);
	} // i_vector loop.
	fprintf(f_op, "\n");

	// Read all vectors and write the final output.
	for (int i_vector = 0; i_vector < (int)overall_n_vectors; i_vector++)
	{
		fprintf(stderr, "Pooling %d. vector.\n", i_vector);

		int cur_vector_pt_i = 0;

		fprintf(f_op, "V%d", i_vector);

		// I am not sure if this is initialized to 0.
		for (int cur_vec_ct_i = 0; cur_vec_ct_i < (int)overall_n_ct_per_vector; cur_vec_ct_i++)
		{
			Plaintext cur_pooled_plaintext;

			// Go over all the sites, load the current plaintext and add it to the final pooled plaintext.
			for (int i_site = 0; i_site < n_sites; i_site++)
			{
				Plaintext cur_site_plaintext;

				cur_site_plaintext.load(context, *(per_site_part_dec_data_ifs->at(i_site)));

				auto &pl_context_data = *context.get_context_data(cur_site_plaintext.parms_id());
				size_t ct_coeff_modulus_size = pl_context_data.parms().coeff_modulus().size();
				auto &pl_parms = pl_context_data.parms();
				size_t pl_coeff_count = pl_parms.poly_modulus_degree();

				// If this is the first site, set the pooled plaintext value.
				if (i_site == 0)
				{
					cur_pooled_plaintext = cur_site_plaintext;
				}
				else
				{
					for (int i_mod = 0; i_mod < (int)ct_coeff_modulus_size; i_mod++)
					{
						add_poly_coeffmod(cur_pooled_plaintext.data() + i_mod * pl_coeff_count,
							cur_site_plaintext.data() + i_mod * pl_coeff_count,
							pl_coeff_count,
							pl_parms.coeff_modulus().at(i_mod),
							cur_pooled_plaintext.data() + i_mod * pl_coeff_count);
					} // i_dim loop.
				}
			} // i_site

			vector<double> pooled_res_array;
			pooled_res_array.resize(encoder.slot_count());
			encoder.decode(cur_pooled_plaintext, pooled_res_array);

			// For the current ciphertext, write the output.
			for (int i_slot = 0; i_slot < (int)(encoder.slot_count()); i_slot++)
			{
				if (cur_vector_pt_i < (int)overall_l_data_vector)
				{
					fprintf(f_op, "\t%.10f", pooled_res_array.at(i_slot));
				}

				// Update the current vector's data point index in the current vector; this counts over all the ciphertexts for this vector.
				cur_vector_pt_i++;
			} // i_slot loop.
		} // cur_vec_ct_i loop.

		fprintf(f_op, "\n");
	} // i_vector loop.
	close_f(f_op, op_fp);
} // collaborative_pool_partial_decrypted_plaintext_results function.