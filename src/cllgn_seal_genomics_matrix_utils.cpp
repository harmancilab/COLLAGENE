#include <stdio.h>
#include <stdlib.h>
#include <seal/seal.h>
#include "cllgn_ansi_string.h"
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_matrix_linalg_utils.h"
#include "cllgn_file_utils.h"
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
#include "cllgn_seal_genomics_key_sharing_utils.h"
#include "cllgn_seal_genomics_matrix_utils.h"

//#include <chrono>
//#include <sys/resource.h>   // check the memory
//#include <random>
//#include <iomanip>          // std::setprecision

#include <vector>

using namespace std;
using namespace seal;
using namespace seal::util;

bool __DUMP_MATRIX_UTILS_MESSAGES__ = false;

#undef __DECRYPT_DEBUG__
//#define __DECRYPT_DEBUG__

void transpose_continuous_encrypted_vector(char* enc_mat_fp, char* text_params_path, char* op_fp)
{
	int loaded_nrow, loaded_ncol;
	vector<Ciphertext>* mat_cts = load_continuous_enc_matrix(enc_mat_fp, loaded_nrow, loaded_ncol, text_params_path);

	if (loaded_nrow != 1 && loaded_ncol != 1)
	{
		fprintf(stderr, "Cannot transpose a non-vector %dx%d matrix..\n", loaded_nrow, loaded_ncol);
		exit(0);
	}

	int nrow = loaded_ncol;
	int ncol = loaded_nrow;

	save_continuous_enc_matrix(mat_cts, nrow, ncol, text_params_path, op_fp);
}

void write_enc_matrix_dimensions(char* enc_mat_fp, char* op_fp)
{
	ifstream ifs_enc_data_matrix(enc_mat_fp, ios::binary);

	// Set the sample size information.
	int n_col;
	int n_row;
	//int n_cts_per_feat;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_col), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_row), sizeof(int));
	ifs_enc_data_matrix.close();

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%d\t%d", n_row, n_col);
	close_f(f_op, op_fp);
}

void write_plain_matrix_dimensions(char* plain_mat_fp, char* op_fp)
{
	int n_col;
	int n_row;
	load_matrix_binary(plain_mat_fp, n_row, n_col);

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%d\t%d", n_row, n_col);
	close_f(f_op, op_fp);
}


void plain_invert_matrix(char* A_mat_fp, char* op_fp)
{
	int loaded_nrow, loaded_ncol;
	double** A_mat = load_matrix_binary(A_mat_fp, loaded_nrow, loaded_ncol, NULL);

	double** inv_A_mat = invert_matrix_GJ(A_mat, loaded_nrow, loaded_ncol, NULL);;

	save_matrix_binary(inv_A_mat, loaded_nrow, loaded_ncol, op_fp);
}

void plain_unpad_matrix_to_size(char* matrix_fp, int new_n_row, int new_n_col, char* op_fp)
{
	int nrows, ncols;
	double** matrix = load_matrix_binary(matrix_fp, nrows, ncols);

	double** new_matrix = allocate_matrix(new_n_row, new_n_col);

	for (int i_row = 0; i_row < new_n_row; i_row++)
	{
		for (int i_col = 0; i_col < new_n_col; i_col++)
		{
			if (i_row < nrows && i_col < ncols)
			{
				new_matrix[i_row][i_col] = matrix[i_row][i_col];
			}
		} // i_col loop.
	} // i_row loop.

	save_matrix_binary(new_matrix, new_n_row, new_n_col, op_fp);
}

void pad_matrix_cols_to_next_power_of_2(char* mat_fp, char* op_fp)
{
	int nRow, nCol;
	double** A_mat = load_matrix_binary(mat_fp, nRow, nCol);
	int n_padded_row = nRow;

	int n_padded_col = 1;
	while (n_padded_col < nCol)
	{
		n_padded_col *= 2;
	}

	fprintf(stderr, "Expanding matrix: [%dx%d]=>[%dx%d]\n", nRow, nCol, n_padded_row, n_padded_col);

	double** padded_mat = allocate_matrix(n_padded_row, n_padded_col);
	copy_matrix(A_mat, nRow, nCol, padded_mat, n_padded_row, n_padded_col);

	save_matrix_binary(padded_mat, n_padded_row, n_padded_col, op_fp);
}

void pad_matrix_rows_to_next_power_of_2(char* mat_fp, char* op_fp)
{
	int nRow, nCol;
	double** A_mat = load_matrix_binary(mat_fp, nRow, nCol);
	int n_padded_row = 1;
	while (n_padded_row < nRow)
	{
		n_padded_row *= 2;
	}

	int n_padded_col = nCol;

	fprintf(stderr, "Expanding matrix: [%dx%d]=>[%dx%d]\n", nRow, nCol, n_padded_row, n_padded_col);

	double** padded_mat = allocate_matrix(n_padded_row, n_padded_col);
	copy_matrix(A_mat, nRow, nCol, padded_mat, n_padded_row, n_padded_col);

	save_matrix_binary(padded_mat, n_padded_row, n_padded_col, op_fp);
}


void pad_matrix_sizee_to_next_power_of_2(char* mat_fp, char* op_fp)
{
	int nRow, nCol;
	double** A_mat = load_matrix_binary(mat_fp, nRow, nCol);
	int n_padded_row = 1;
	while (n_padded_row < nRow)
	{
		n_padded_row *= 2;
	}

	int n_padded_col = 1;
	while (n_padded_col < nCol)
	{
		n_padded_col *= 2;
	}

	fprintf(stderr, "Expanding matrix: [%dx%d]=>[%dx%d]\n", nRow, nCol, n_padded_row, n_padded_col);

	double** padded_mat = allocate_matrix(n_padded_row, n_padded_col);
	copy_matrix(A_mat, nRow, nCol, padded_mat, n_padded_row, n_padded_col);

	save_matrix_binary(padded_mat, n_padded_row, n_padded_col, op_fp);
}

//void plain_add_matrices_per_list(char* matrix_path_list_fp, char* op_fp)
//{
//	vector<char*>* matrix_path_list = buffer_file(matrix_path_list_fp);
//	fprintf(stderr, "Adding %d matrices and writing results to %s\n", matrix_path_list->size(), op_fp);
//
//	double** res_matrix = NULL;
//	int nrows = 0;
//	int ncols = 0;
//
//	for (int i_mat = 0; i_mat < matrix_path_list->size(); i_mat++)
//	{
//		int cur_n_rows, cur_n_cols;
//		double** cur_mat = load_matrix_binary(matrix_path_list->at(i_mat), cur_n_rows, cur_n_cols);
//		if (i_mat == 0)
//		{
//			res_matrix = cur_mat;
//			nrows = cur_n_rows;
//			ncols = cur_n_cols;
//		}
//		else
//		{
//			matrix_add(res_matrix, nrows, ncols, cur_mat, cur_n_rows, cur_n_cols, res_matrix);
//		}
//	} // i_mat loop.
//
//	save_matrix_binary(res_matrix, nrows, ncols, op_fp);
//}


void write_vital_stats_per_continuous_encrypted_matrix(char* enc_mat_fp,
	char* text_params_path,
	char* vital_stats_fp)
{
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading encrypted data matrix from %s\n", enc_mat_fp);
	}

	ifstream ifs_enc_data_matrix(enc_mat_fp, ios::binary);

	// Set the sample size information.
	int n_col;
	int n_row;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_col), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_row), sizeof(int));

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading %dx%d matrix..\n", n_row, n_col);
	}

	double n_data_pts_per_matrix = n_row * n_col;
	int n_cts = (int)(ceil(n_data_pts_per_matrix / n_values_per_ct));

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading ciphertexts.\n");
	}

	FILE* f_vital_stats_op = open_f(vital_stats_fp, "w");
	for (int i_ct = 0; i_ct < n_cts; i_ct++)
	{
		Ciphertext cur_ct;
		cur_ct.load(context, ifs_enc_data_matrix);

		double cur_ct_scale_n_bits = log(cur_ct.scale()) / log(2);

		fprintf(stderr, "Ciphertext: %d: Size: %d\tChain_Index: %d\tmodulus_degree: %d\tpoly_mod_degree: %d\tscale: %.4f\n",
			i_ct,
			(int)(cur_ct.size()),
			(int)(context.get_context_data(cur_ct.parms_id())->chain_index()),
			(int)(cur_ct.coeff_modulus_size()),
			(int)(cur_ct.poly_modulus_degree()),
			(double)cur_ct_scale_n_bits);

		fprintf(f_vital_stats_op, "%d\t%d\t%d\t%d\t%d\t%.4f\n",
			i_ct,
			(int)(cur_ct.size()),
			(int)(context.get_context_data(cur_ct.parms_id())->chain_index()),
			(int)(cur_ct.coeff_modulus_size()),
			(int)(cur_ct.poly_modulus_degree()),
			(double)cur_ct_scale_n_bits);
	} // i_feat.

	// Close file.
	close_f(f_vital_stats_op, vital_stats_fp);
} // write_vital_stats_per_continuous_encrypted_matrix function.

void encrypt_plaintext_matrix(double** matrix, int nrow, int ncol,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_public_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_public_key_path);
		exit(0);
	}

	// Instantiate and load the public key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();
	fprintf(stderr, "Loaded public key.\n");

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	// Load each row to a ciphertext and encrypt.
	ofstream ofs_enc_matrix(op_fp, ios::binary);

	int n_cts_per_row = int(ceil((double)ncol / n_values_per_ct));

	// This describes the matrix that is being saved to output.
	int f_nrow = nrow; // This is the number of rows of the matrix that is being encrypted.
	int f_ncol = ncol; // This is the number of columns for the matrix that is being encrypted.
	//int n_result_cts_per_vec = n_cts_per_row;
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&f_ncol), sizeof(int));
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&(f_nrow)), sizeof(int));
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&n_cts_per_row), sizeof(int));

	// Encrypt each row and write.
	for (int row = 0; row < nrow; row++)
	{
		int i_col = 0;
		int n_cts_per_cur_row = 0;
		while (i_col < ncol)
		{
			fprintf(stderr, "Encrypting data vector %d[%d-%d]\n", row, i_col, i_col + n_values_per_ct - 1);

			vector<double> cur_data_pts;
			for (int cur_cipher_data_pt_i = 0; cur_cipher_data_pt_i < n_values_per_ct; cur_cipher_data_pt_i++)
			{
				if (i_col < ncol)
				{
					cur_data_pts.push_back(matrix[row][i_col]);
				}
				else
				{
					cur_data_pts.push_back(0);
				}

				// Move to the next column for the current row.
				i_col++;
			} // cur_cipher_data_pt_i loop.

			// Encrypt the data points and save to file.
			Plaintext cur_plain;

			// Encode the vector.
			encoder.encode(cur_data_pts, scale, cur_plain);

			// Encrypt the plaintext.
			Ciphertext ct;
			encryptor.encrypt(cur_plain, ct);

			// Save the encrypted ciphertext.
			ct.save(ofs_enc_matrix);

			// Count the number of cts in the current row.
			n_cts_per_cur_row++;
		} // data_pt_i loop.

		if (n_cts_per_row != n_cts_per_cur_row)
		{
			fprintf(stderr, "Sanity check failed: # of ciphertexts per vector is not as expected for %d. row: %d / %d\n",
				row,
				n_cts_per_row,
				n_cts_per_cur_row);
			exit(0);
		}
	} // row loop.

	// Close the file.
	ofs_enc_matrix.close();
} // encrypt_plaintext_matrix

void save_continuous_enc_matrix(vector<Ciphertext>* matrix_cts, int nrow, int ncol,
	char* text_params_path,
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	// Load each row to a ciphertext and encrypt.
	ofstream ofs_enc_matrix(op_fp, ios::binary);

	// This describes the matrix that is being saved to output.
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&ncol), sizeof(int));
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&(nrow)), sizeof(int));

	double n_total_entries_in_matrix = (double)(nrow * ncol);
	int n_cts_per_matrix = (int)ceil(n_total_entries_in_matrix / n_values_per_ct);
	if (n_cts_per_matrix != (int)(matrix_cts->size()))
	{
		fprintf(stderr, "Matrix size is not conformant with the # of cts: %dx%d; %d cts\n", nrow, ncol, n_cts_per_matrix);
	}

	// Encrypt each row and write.
	for (int i_ct = 0; i_ct < (int)(matrix_cts->size()); i_ct++)
	{
		// Encrypt the plaintext.
		Ciphertext ct = matrix_cts->at(i_ct);

		// Save the encrypted ciphertext.
		ct.save(ofs_enc_matrix);
	} // row loop.

	// Close the file.
	ofs_enc_matrix.close();
} // save_continuous_enc_matrix

void save_enc_matrix(vector<Ciphertext>* matrix_cts, int nrow, int ncol,
	char* text_params_path,
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	// Load each row to a ciphertext and encrypt.
	ofstream ofs_enc_matrix(op_fp, ios::binary);

	int n_cts_per_row = int(ceil((double)ncol / n_values_per_ct));

	// This describes the matrix that is being saved to output.
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&ncol), sizeof(int));
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&(nrow)), sizeof(int));
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&n_cts_per_row), sizeof(int));

	// Encrypt each row and write.
	for (int i_ct = 0; i_ct < (int)(matrix_cts->size()); i_ct++)
	{
		// Encrypt the plaintext.
		Ciphertext ct = matrix_cts->at(i_ct);

		// Save the encrypted ciphertext.
		ct.save(ofs_enc_matrix);
	} // row loop.

	// Close the file.
	ofs_enc_matrix.close();
} // save_enc_matrix

vector<Ciphertext>* load_enc_matrix(char* enc_matrix_fp, int& loaded_nrow, int& loaded_ncol, int& loaded_n_cts_per_row, char* text_params_path)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading encrypted matrix from %s..\n", enc_matrix_fp);
	}

	// Load each row to a ciphertext and encrypt.
	ifstream ifs_enc_matrix(enc_matrix_fp, ios::binary);

	// This describes the matrix that is being saved to output.
	int f_nrow = 0; // This parameter has no bearing yet. It should be copied to partially decrypted plaintext file and final pooling file to write the results in the final final text output file. This sets how many values are meaningful (or necessary) in the resulting ciphertext.
	int f_ncol = 0;
	int f_n_cts_per_row = 0;
	ifs_enc_matrix.read(reinterpret_cast<char *>(&f_ncol), sizeof(int));
	ifs_enc_matrix.read(reinterpret_cast<char *>(&(f_nrow)), sizeof(int));
	ifs_enc_matrix.read(reinterpret_cast<char *>(&f_n_cts_per_row), sizeof(int));

	loaded_nrow = f_nrow;
	loaded_ncol = f_ncol;
	loaded_n_cts_per_row = f_n_cts_per_row;

	fprintf(stderr, "Read dimensions %dx%d, %d cts per row..\n", loaded_nrow, loaded_ncol, loaded_n_cts_per_row);

	int n_cts_per_row = int(ceil((double)loaded_ncol / n_values_per_ct));
	if (n_cts_per_row != f_n_cts_per_row)
	{
		fprintf(stderr, "Sanity check failed: $ of cts per row are not conformant while loading %s: %d, %d\n", enc_matrix_fp, n_cts_per_row, f_n_cts_per_row);
		exit(0);
	}

	// Start loading the ct's.
	vector<Ciphertext>* matrix_cts = new vector<Ciphertext>();
	for (int row = 0; row < loaded_nrow; row++)
	{
		int i_ct = 0;
		int n_cts_per_cur_row = 0;
		int i_col = 0;
		while (i_ct < loaded_n_cts_per_row)
		{
			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Loading data vector %d[%d-%d]\n", row, i_col, i_col + n_values_per_ct - 1);
			}

			// Encrypt the plaintext.
			Ciphertext ct;

			// Save the encrypted ciphertext.
			ct.load(context, ifs_enc_matrix);

			// Read the ct.
			matrix_cts->push_back(ct);

			// Count the number of cts in the current row.
			i_col += n_values_per_ct;
			i_ct++;
			n_cts_per_cur_row++;
		} // data_pt_i loop.

		if (n_cts_per_row != n_cts_per_cur_row)
		{
			fprintf(stderr, "Sanity check failed: # of ciphertexts per vector is not as expected for %d. row: %d / %d\n",
				row,
				n_cts_per_row,
				n_cts_per_cur_row);
			exit(0);
		}
	} // row loop.

	// Close the file.
	ifs_enc_matrix.close();

	return(matrix_cts);
} // load_enc_matrix option.

// This is used for partially decrypting other user's data matrices.
void partial_decrypt_matrix(char* encA_mat_fp, bool is_site0, char* text_params_path, char* noisy_sk_path, double smdg_noise_variance, char* partial_decrypt_data_path)
{
	// Use key sharing library.
	partial_decrypt_data_per_noisy_key(encA_mat_fp, is_site0, text_params_path, noisy_sk_path, smdg_noise_variance, partial_decrypt_data_path);
} // partial_decrypt_matrix

// This is used for pooling other site's decryptions into a final result.
void pool_partially_decrypted_matrix(char* partial_decrypted_matrix_path_list_fp, char* text_params_path, char* full_decrypted_matrix_fp)
{
	// Use key sharing library.
	collaborative_pool_partial_decrypted_plaintext_results(partial_decrypted_matrix_path_list_fp, text_params_path, full_decrypted_matrix_fp);
}

void fully_decrypt_matrix(char* enc_matrix_path, char* text_params_path, char* pooled_secret_key_path, char* full_decrypted_matrix_fp)
{
	// Use key sharing library.
	full_decrypt_encrypted_data_per_pooled_key(enc_matrix_path, text_params_path, pooled_secret_key_path, full_decrypted_matrix_fp);
}

void fully_decrypt_continuous_encrypted_matrix(char* enc_matrix_fp, char* text_params_path,
	char* pooled_private_key_path,
	char* full_decrypted_matrix_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_private_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	}

	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	//int n_values_per_ct = encoder.slot_count();

	int loaded_nrows, loaded_ncols;
	vector<Ciphertext>* matrix_cts = load_continuous_enc_matrix(enc_matrix_fp, loaded_nrows, loaded_ncols, text_params_path);

	fprintf(stderr, "Loaded %d cts from %s, processing data points.\n", (int)(matrix_cts->size()), enc_matrix_fp);

	int n_entries_in_matrix = loaded_ncols * loaded_nrows;

	vector<double> all_data_pts;
	for (int i_ct = 0; i_ct < (int)(matrix_cts->size()); i_ct++)
	{
		Plaintext cur_dec_pt;
		decryptor.decrypt(matrix_cts->at(i_ct), cur_dec_pt);

		vector<double> cur_ct_dec_vals;
		encoder.decode(cur_dec_pt, cur_ct_dec_vals);
		all_data_pts.insert(all_data_pts.end(), cur_ct_dec_vals.begin(), cur_ct_dec_vals.end());
	} // i_ct loop.

	fprintf(stderr, "Loaded %d data points from %d continuous encrypted cts (%d entries)\n", (int)(all_data_pts.size()), (int)(matrix_cts->size()), n_entries_in_matrix);

	int i_data = 0;
	FILE* f_op = open_f(full_decrypted_matrix_fp, "w");
	double** dec_matrix = allocate_matrix(loaded_nrows, loaded_ncols);
	for (int i_row = 0; i_row < loaded_nrows; i_row++)
	{
		fprintf(f_op, "ROW_%d", i_row);
		for (int i_col = 0; i_col < loaded_ncols; i_col++)
		{
			fprintf(f_op, "\t%.15f", all_data_pts.at(i_data));
			dec_matrix[i_row][i_col] = all_data_pts.at(i_data);
			i_data++;
		} // i_col loop.

		fprintf(f_op, "\n");
	} // i_row loop.
	close_f(f_op, full_decrypted_matrix_fp);

	char bin_matrix_fp[1000];
	sprintf(bin_matrix_fp, "%s.bin", full_decrypted_matrix_fp);
	save_matrix_binary(dec_matrix, loaded_nrows, loaded_ncols, bin_matrix_fp);
} // fully_decrypt_continuous_encrypted_matrix function.


// We use the algorithm in https://eprint.iacr.org/2018/1041.pdf
// This function multiplies the matrices AxB -- Note that B is given as Bt here.
// The problem with this approach is that each inner product provides one entry and we need to do a lot of shifts to move them into the result matrix.
// Therefore this function will not be efficient for large matrices.
void secure_multiply_matrices(char* encA_mat_fp, char* encBt_mat_fp, char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
		// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}
	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	// Load galois key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	}
	ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	GaloisKeys pooled_galois_key;
	pooled_galois_key.load(context, ifs_galois_key);
	ifs_galois_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded Galois key.\n");
	}

	// Load secret key -- this should be the actual pooled key.
#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not load private key from %s.\n", pooled_private_key_path);
	}

	SecretKey pooled_sk;
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);

	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor.
	Decryptor decryptor(context, pooled_sk);
#endif

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	// Load the matrices and do all pairwise inner products of the rows between A and Bt.
	int loaded_nrow_A, loaded_ncol_A, loaded_n_cts_per_rowA;
	vector<Ciphertext>* A_cts = load_enc_matrix(encA_mat_fp, loaded_nrow_A, loaded_ncol_A, loaded_n_cts_per_rowA, text_params_path);

	int loaded_nrow_Bt, loaded_ncol_Bt, loaded_n_cts_per_row_Bt;
	vector<Ciphertext>* Bt_cts = load_enc_matrix(encBt_mat_fp, loaded_nrow_Bt, loaded_ncol_Bt, loaded_n_cts_per_row_Bt, text_params_path);

	if (loaded_ncol_A != loaded_ncol_Bt)
	{
		fprintf(stderr, "Sanity check failed while secure multiplying matrices: %d, %d; %d, %d\n", loaded_nrow_A, loaded_ncol_A, loaded_ncol_Bt, loaded_nrow_Bt);
		exit(0);
	}

	if (loaded_n_cts_per_rowA != loaded_n_cts_per_row_Bt)
	{
		fprintf(stderr, "Sanity check failed for # cts/per row while secure multiplying matrices (Matrices encrypting with different parameters?): %d, %d\n",
			loaded_n_cts_per_rowA,
			loaded_n_cts_per_row_Bt);
		exit(0);
	}

	fprintf(stderr, "Loaded matrices:\n\
A[%d,%d] @ %d cts/row\n\
Bt[%d,%d] @ %d cts/row\n", loaded_nrow_A, loaded_ncol_A, loaded_n_cts_per_rowA, loaded_nrow_Bt, loaded_ncol_Bt, loaded_n_cts_per_row_Bt);

	//vector<Ciphertext>* res_matrix_cts = new vector<Ciphertext>();
	int res_nrow = loaded_nrow_A;
	int res_ncol = loaded_nrow_Bt;

	// This is a masker array that is used to put the enctries at the correct place.
	vector<double> cur_entry_masker_array;
	for (int cur_cipher_data_pt_i = 0; cur_cipher_data_pt_i < n_values_per_ct; cur_cipher_data_pt_i++)
	{
		cur_entry_masker_array.push_back(0);
	} // cur_cipher_data_pt_i loop.

	// This defines the output matrix's # cts/row, it depends on the column number of result, which is the number of rows of Bt.
	int res_n_cts_per_row = (int)(ceil((double)res_ncol / n_values_per_ct));
	fprintf(stderr, "# cts/row for the result: %d\n", res_n_cts_per_row);

	vector<Ciphertext>* res_cts = new vector<Ciphertext>();
	for (int row_i = 0; row_i < res_nrow; row_i++)
	{
		fprintf(stderr, "Calculating row=%d\n", row_i);

		vector<Ciphertext>* cur_row_cts = new vector<Ciphertext>();
		for (int ct_i = 0; ct_i < res_n_cts_per_row; ct_i++)
		{
			Ciphertext cur_entry_ct;
			encryptor.encrypt_zero(cur_entry_ct);
			cur_entry_ct.scale() = scale;

			cur_row_cts->push_back(cur_entry_ct);
		} // ct_i loop.

		for (int col_i = 0; col_i < res_ncol; col_i++)
		{
			Ciphertext cur_entry_ct;
			encryptor.encrypt_zero(cur_entry_ct);
			cur_entry_ct.scale() = scale;

			for (int ct_i = 0; ct_i < loaded_n_cts_per_rowA; ct_i++)
			{
				// Take the inner product of the ciphertexts: Following assumes there is exactly 1 ct per row.
				// B_ct_i: The index of the ciphertext in Bt
				int A_ct_i = ct_i + row_i * loaded_n_cts_per_rowA;
				int Bt_ct_i = ct_i + col_i * loaded_n_cts_per_row_Bt;

				Ciphertext row_A = A_cts->at(A_ct_i);
				Ciphertext row_B = Bt_cts->at(Bt_ct_i);
				Ciphertext cur_mult;
				evaluator.multiply(row_A, row_B, cur_mult);
				evaluator.relinearize_inplace(cur_mult, pooled_relin_key);
				evaluator.rescale_to_next_inplace(cur_mult);

				int cur_l_rotation = 1;

				// 0th rotation: The final results is stored in cur_mult.
				Ciphertext cur_rotated = cur_mult; // CT is initialized here.
				while (cur_l_rotation < n_values_per_ct)
				{
					fprintf(stderr, "Processing rotation %d\n", cur_l_rotation);

					// Rotate the current by 1; this can be performed much faster.
					evaluator.rotate_vector(cur_mult, cur_l_rotation, pooled_galois_key, cur_rotated);

					// Add the next rotation to the current cumulative sum.
					evaluator.add_inplace(cur_mult, cur_rotated);

					cur_l_rotation *= 2;
				} // i_s loop.

#ifdef __DECRYPT_DEBUG__
				/////////////////////////////////////////////////////////////////////////////
				// Decrypt and check the entries in the summated inner product:
				Plaintext cur_vec_pl;
				decryptor.decrypt(cur_mult, cur_vec_pl);

				vector<double> cur_vec_data;
				encoder.decode(cur_vec_pl, cur_vec_data);
				fprintf(stderr, "Dec[inner_prod(%d,%d)::ct_i[%d]] %.3f\t%.3f\t%.3f\t...\n", row_i, col_i, ct_i,
					cur_vec_data.at(0), cur_vec_data.at(1), cur_vec_data.at(2));
				/////////////////////////////////////////////////////////////////////////////
#endif

				if (__DUMP_MATRIX_UTILS_MESSAGES__)
				{
					fprintf(stderr, "cur entry ct info: scale: %.4f; coeff_modulus_size: %d; poly_modulus_degree: %d; chain index: %d;\n",
						log(cur_entry_ct.scale()) / log(2),
						(int)(cur_entry_ct.coeff_modulus_size()),
						(int)(cur_entry_ct.poly_modulus_degree()),
						(int)((context.get_context_data(cur_entry_ct.parms_id()))->chain_index()));

					fprintf(stderr, "cur mult ct info: scale: %.4f; coeff_modulus_size: %d; poly_modulus_degree: %d; chain index: %d;\n",
						log(cur_mult.scale()) / log(2),
						(int)(cur_entry_ct.coeff_modulus_size()),
						(int)(cur_entry_ct.poly_modulus_degree()),
						(int)((context.get_context_data(cur_entry_ct.parms_id()))->chain_index()));
				}

				// Add the current multiplication to cur entry ct.
				evaluator.mod_switch_to_inplace(cur_entry_ct, cur_mult.parms_id());
				cur_entry_ct.scale() = cur_mult.scale();
				evaluator.add_inplace(cur_entry_ct, cur_mult);
			} // ct_i loop.

			// The important observation is that the above inner product stores the inner product value in every slot of the plaintext that is stored in the resulting ciphertext.
#ifdef __DECRYPT_DEBUG__
			/////////////////////////////////////////////////////////////////////////////
			// Decrypt and check the entries in the summated inner product:
			Plaintext cur_vec_pl;
			decryptor.decrypt(cur_entry_ct, cur_vec_pl);

			vector<double> cur_vec_data;
			encoder.decode(cur_vec_pl, cur_vec_data);
			fprintf(stderr, "Dec[inner_prod(%d,%d)::ct_i[*]] %.3f\t%.3f\t%.3f\t...\n", row_i, col_i,
				cur_vec_data.at(0), cur_vec_data.at(1), cur_vec_data.at(2));
			/////////////////////////////////////////////////////////////////////////////
#endif  // __DECRYPT_DEBUG__

			// Reset the mask.
			std::fill(cur_entry_masker_array.begin(), cur_entry_masker_array.end(), 0);

			// If there are multipl ciphertexts per row, we need to make sure to correct for the columns that are occupied by the columns that are stored in previous ciphertexts.
			cur_entry_masker_array.at(col_i % n_values_per_ct) = 1;

			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "row_i: %d, col_i: %d; ct slot index: %d is set to %.2f (%.2f), n_values_per_ct=%d\n",
					row_i, col_i, col_i % n_values_per_ct, cur_entry_masker_array.at(col_i % n_values_per_ct),
					cur_entry_masker_array.at((col_i % n_values_per_ct) + 1), n_values_per_ct);
			}

			// Encrypt the data points and save to file.
			Plaintext cur_entry_masker_array_pt;

			// Encode the vector.
			encoder.encode(cur_entry_masker_array, scale, cur_entry_masker_array_pt);
			evaluator.mod_switch_to_inplace(cur_entry_masker_array_pt, cur_entry_ct.parms_id());
			cur_entry_masker_array_pt.scale() = cur_entry_ct.scale();

			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "ct: scale: %.4f; coeff_modulus_size: %d; poly_modulus_degree: %d; chain index: %d;\n", log(cur_entry_ct.scale()) / log(2),
					(int)(cur_entry_ct.coeff_modulus_size()),
					(int)(cur_entry_ct.poly_modulus_degree()),
					(int)((context.get_context_data(cur_entry_ct.parms_id()))->chain_index()));
			}

			// Multiple with the mask.
			evaluator.multiply_plain_inplace(cur_entry_ct, cur_entry_masker_array_pt);
			evaluator.relinearize_inplace(cur_entry_ct, pooled_relin_key);
			evaluator.rescale_to_next_inplace(cur_entry_ct);

			// Add this value to the ct for the current row.
			int cur_col_i_ct = (int)(floor((double)col_i / n_values_per_ct));
			evaluator.mod_switch_to_inplace(cur_row_cts->at(cur_col_i_ct), cur_entry_ct.parms_id());
			cur_row_cts->at(cur_col_i_ct).scale() = cur_entry_ct.scale();
			evaluator.add_inplace(cur_row_cts->at(cur_col_i_ct), cur_entry_ct);

#ifdef __DECRYPT_DEBUG__
			/////////////////////////////////////////////////////////////////////////////
			// Decrypt and check the entries in the summated inner product:
			//Plaintext cur_vec_pl;
			decryptor.decrypt(cur_row_cts->at(cur_col_i_ct), cur_vec_pl);

			//vector<double> cur_vec_data;
			encoder.decode(cur_vec_pl, cur_vec_data);
			fprintf(stderr, "Row_i:%d::Dec[row_vec_ct[%d], col_i:%d]]: %.3f\t%.3f\t%.3f\t...\n", row_i, cur_col_i_ct,
				col_i,
				cur_vec_data.at(0), cur_vec_data.at(1), cur_vec_data.at(2));
			/////////////////////////////////////////////////////////////////////////////
#endif
		} // col_i loop.

		// Add the current ct after calculating them. These are objects, not pointers.
		res_cts->insert(res_cts->end(), cur_row_cts->begin(), cur_row_cts->end());
	} // row_i loop.

	fprintf(stderr, "Saving %d cts for a [%d,%d matrix.\n", (int)(res_cts->size()), res_nrow, res_ncol);
	save_enc_matrix(res_cts, res_nrow, res_ncol, text_params_path, op_fp);
} // secure_multiply_matrices

void encrypt_encoded_pt_matrix(char* encoded_matrix_pt_path,
	char* text_params_path,
	char* pooled_public_key_path,
	char* encrypted_data_path)
{
	fprintf(stderr, "Encrypting encoded pt matrix in %s and saving to %s\n",
		encoded_matrix_pt_path, encrypted_data_path);

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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);
	CKKSEncoder encoder(context);

	// Instantiate and load the public key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	int n_entries_per_ct = encoder.slot_count();

	//// Load the secret (noisy) key.
	//fprintf(stderr, "Loading noisy secret key from %s\n", noisy_key_path);
	//ifstream ifs_noisy_key(noisy_key_path, ios::binary);
	//SecretKey dec_sk;
	//dec_sk.load(context, ifs_noisy_key);
	//ifs_noisy_key.close();
	//fprintf(stderr, "Loaded noisy secret key.\n", noisy_key_path);

	//// Instantiate the decryptor.
	//Decryptor decryptor(context, dec_sk);

	// Just load the ciphertexts from the results and decrypt them with this site's noisy key.
	ifstream ifs_dec_data_matrix(encoded_matrix_pt_path, ios::binary);
	ofstream ofs_enc_data_matrix(encrypted_data_path, ios::binary);

	// Read the encrypted data matrix parameters.
	int n_rows;
	int n_cols;

	// Read the matrix size parameters.
	ifs_dec_data_matrix.read(reinterpret_cast<char *>(&n_cols), sizeof(int));
	ifs_dec_data_matrix.read(reinterpret_cast<char *>(&n_rows), sizeof(int));

	ofs_enc_data_matrix.write(reinterpret_cast<const char *>(&n_cols), sizeof(int));
	ofs_enc_data_matrix.write(reinterpret_cast<const char *>(&n_rows), sizeof(int));

	double n_entries_in_matrix = (double)(n_rows * n_cols);
	int n_pts_2_read = (int)ceil(n_entries_in_matrix / n_entries_per_ct);

	fprintf(stderr, "Read %dx%d dimensions, reading %d pt's.\n", n_rows, n_cols, n_pts_2_read);

	// Encrypt each row and write.
	vector<double> cur_ct_data_pts;

	for (int pt_i = 0; pt_i < n_pts_2_read; pt_i++)
	{
		Plaintext cur_pt;

		// Save the encrypted ciphertext.
		cur_pt.load(context, ifs_dec_data_matrix);

		// Encrypt the plaintext.
		Ciphertext cur_ct;
		encryptor.encrypt(cur_pt, cur_ct);

		cur_ct.save(ofs_enc_data_matrix);
	} // data_pt_i loop.

	// Close the file.
	ifs_dec_data_matrix.close();
	ofs_enc_data_matrix.close();
} // encrypt_partial_decrypted_matrix

void encrypt_plaintext_matrix_continuous_ct(double** matrix, int nrow, int ncol,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_fp)
{
	if (nrow == 0 || ncol == 0)
	{
		fprintf(stderr, "Cannot encrypt 0 size matrix.\n");
		exit(1);
	}

	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_public_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_public_key_path);
		exit(0);
	}

	// Instantiate and load the public key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	// Load each row to a ciphertext and encrypt.
	ofstream ofs_enc_matrix(op_fp, ios::binary);

	//int n_cts_per_row = int(ceil((double)ncol / n_values_per_ct));

	// This describes the matrix that is being saved to output.
	int f_nrow = nrow; // This is the number of rows of the matrix that is being encrypted.
	int f_ncol = ncol; // This is the number of columns for the matrix that is being encrypted.
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&f_ncol), sizeof(int));
	ofs_enc_matrix.write(reinterpret_cast<const char *>(&(f_nrow)), sizeof(int));

	// Encrypt each row and write.
	vector<double> cur_ct_data_pts;
	int total_n_saved_entries = 0;
	for (int row = 0; row < nrow; row++)
	{
		for (int col = 0; col < ncol; col++)
		{
			// Push the current point.
			cur_ct_data_pts.push_back(matrix[row][col]);

			if ((int)(cur_ct_data_pts.size()) == n_values_per_ct)
			{
				if (__DUMP_MATRIX_UTILS_MESSAGES__)
				{
					fprintf(stderr, "Writing the ct with data array of length %d\n", (int)(cur_ct_data_pts.size()));
				}

				total_n_saved_entries += cur_ct_data_pts.size();

				// Encrypt the data points and save to file.
				Plaintext cur_plain;

				// Encode the vector.
				encoder.encode(cur_ct_data_pts, scale, cur_plain);

				// Encrypt the plaintext.
				Ciphertext ct;
				encryptor.encrypt(cur_plain, ct);

				// Save the encrypted ciphertext.
				ct.save(ofs_enc_matrix);

				// Empty the data points.
				cur_ct_data_pts.clear();
			}
		} // data_pt_i loop.
	} // row loop.

	// If there are entries, save them to a new ciphertext.
	if (cur_ct_data_pts.size() > 0)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Writing the ct with data array of length %d\n", (int)(cur_ct_data_pts.size()));
		}

		total_n_saved_entries += cur_ct_data_pts.size();

		// Encrypt the data points and save to file.
		Plaintext cur_plain;

		// Encode the vector.
		encoder.encode(cur_ct_data_pts, scale, cur_plain);

		// Encrypt the plaintext.
		Ciphertext ct;
		encryptor.encrypt(cur_plain, ct);

		// Save the encrypted ciphertext.
		ct.save(ofs_enc_matrix);
	}

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Encrypted %d values.\n", total_n_saved_entries);
	}

	// Close the file.
	ofs_enc_matrix.close();
} // encrypt_plaintext_matrix_continuous_ct

void col_expand_dense_encrypt_matrix(char* A_fp, int n_cols_per_expanded_matrix,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_dir)
{
	int nrows, ncols;
	double** matrix = load_matrix_binary(A_fp, nrows, ncols, NULL);
	fprintf(stderr, "Loaded %dx%d matrix from %s\n", nrows, ncols, A_fp);

	// Replicate the column n_cols_per_expanded_matrix times and write the file.
	for (int i_repcol = 0; i_repcol < ncols; i_repcol++)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Generating %d. repcol matrix.\n", i_repcol);
		}

		double** cur_colrep_matrix = allocate_matrix(nrows, n_cols_per_expanded_matrix);

		for (int i_row = 0; i_row < nrows; i_row++)
		{
			for (int i_col = 0; i_col < n_cols_per_expanded_matrix; i_col++)
			{
				cur_colrep_matrix[i_row][i_col] = matrix[i_row][i_repcol];
			} // i_col_loop.
		} // i_row loop.

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			char repcol_mat_fp[1000];
			sprintf(repcol_mat_fp, "%s/repcol_%d.bin", op_dir, i_repcol);
			save_matrix_binary(cur_colrep_matrix, nrows, n_cols_per_expanded_matrix, repcol_mat_fp);
		}

		// Save the encrypted matrix.
		char enc_mat_fp[1000];
		sprintf(enc_mat_fp, "%s/repcol_%d.bin.enc", op_dir, i_repcol);
		encrypt_plaintext_matrix_continuous_ct(cur_colrep_matrix, nrows, n_cols_per_expanded_matrix,
			text_params_path, pooled_public_key_path, enc_mat_fp);
	} // i_repcol loop.
}

void row_expand_continuous_encrypted_matrix(char* enc_A_fp, int n_rows_per_expanded_matrix, char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_dir)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}
	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	// Load galois key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	}
	ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	GaloisKeys pooled_galois_key;
	pooled_galois_key.load(context, ifs_galois_key);
	ifs_galois_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded Galois key.\n");
	}

#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not load private key from %s\n", pooled_private_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor.
	Decryptor decryptor(context, pooled_sk);
#endif 

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	int loaded_nArows, loaded_nAcols;
	vector<Ciphertext>* A_cts = load_continuous_enc_matrix(enc_A_fp, loaded_nArows, loaded_nAcols, text_params_path);

	double n_entries_per_matrix = (double)(loaded_nAcols * loaded_nArows);
	int n_exp_cts = (int)(ceil(n_entries_per_matrix / n_values_per_ct));

	if (n_exp_cts != (int)(A_cts->size()))
	{
		fprintf(stderr, "Expected ct size is not conformant for %dx%d matrices: %d vs %d\n", loaded_nArows, loaded_nAcols, n_exp_cts, (int)A_cts->size());
		exit(0);
	}

	double log2_n_cols = log(loaded_nAcols) / log(2.0);
	if (pow(2, log2_n_cols) != loaded_nAcols)
	{
		fprintf(stderr, "Need a matrix with number of columns that is an exact power of 2.\n");
		exit(0);
	}

	// We currently cannot process matrices whose rows expand over ciphertext boundaries. The reason is that this would require keeping merging ciphertext, etc.
	if (loaded_nAcols > n_values_per_ct)
	{
		fprintf(stderr, "Cannot process %d-column matrices with larger than %d long roaw.\n", loaded_nAcols, n_values_per_ct);
		exit(0);
	}

	if (loaded_nArows > n_values_per_ct)
	{
		fprintf(stderr, "Cannot process matrices with larger than %d rows.\n", loaded_nArows);
		exit(0);
	}

	int n_rows_per_ct = n_values_per_ct / loaded_nAcols;

	vector<double> cur_entry_masker_array;
	for (int cur_cipher_data_pt_i = 0; cur_cipher_data_pt_i < n_values_per_ct; cur_cipher_data_pt_i++)
	{
		if (cur_cipher_data_pt_i < loaded_nAcols)
		{
			cur_entry_masker_array.push_back(1);
		}
		else
		{
			cur_entry_masker_array.push_back(0);
		}
	} // cur_cipher_data_pt_i loop.

	// Encrypt the data points and save to file.
	Plaintext cur_entry_masker_array_pt;

	// Encode the vector.
	encoder.encode(cur_entry_masker_array, scale, cur_entry_masker_array_pt);

	for (int i_rep_row = 0; i_rep_row < loaded_nArows; i_rep_row++)
	{
		Ciphertext shifted_ct;

		// This is the ct index for the current row: This expands to multiple ct's.
		int ct_i_per_cur_row = (i_rep_row * loaded_nAcols) / n_values_per_ct;

		// This is the index of row within ct.
		int row_i_in_ct = i_rep_row % n_rows_per_ct;

		int shift_number = row_i_in_ct * loaded_nAcols;

		// The row is pushed to left.
		evaluator.rotate_vector(A_cts->at(ct_i_per_cur_row), shift_number, pooled_galois_key, shifted_ct);

#ifdef __DECRYPT_DEBUG__
		////////////////
		// Decrypt the current ciphertext.
		Plaintext cur_dec_plaintext;
		decryptor.decrypt(shifted_ct, cur_dec_plaintext);
		vector<double> cur_ct_dec_vals;
		encoder.decode(cur_dec_plaintext, cur_ct_dec_vals);
		fprintf(stderr, "Init shifted for reprow=%d:", i_rep_row);
		for (int i = 0; i < loaded_nAcols * 2; i++)
		{
			fprintf(stderr, " %.3f", cur_ct_dec_vals[i]);
		} // i loop.
		fprintf(stderr, "\n");
		////////////////
#endif

		// Mask this.
		evaluator.mod_switch_to_inplace(cur_entry_masker_array_pt, shifted_ct.parms_id());
		cur_entry_masker_array_pt.scale() = shifted_ct.scale();

		// Multiply with the mask.
		evaluator.multiply_plain_inplace(shifted_ct, cur_entry_masker_array_pt);
		evaluator.relinearize_inplace(shifted_ct, pooled_relin_key);
		evaluator.rescale_to_next_inplace(shifted_ct);

#ifdef __DECRYPT_DEBUG__
		////////////////
		// Decrypt the current ciphertext.
		decryptor.decrypt(shifted_ct, cur_dec_plaintext);
		cur_ct_dec_vals;
		encoder.decode(cur_dec_plaintext, cur_ct_dec_vals);
		fprintf(stderr, "Mask shifted for reprow=%d:", i_rep_row);
		for (int i = 0; i < loaded_nAcols * 2; i++)
		{
			fprintf(stderr, " %.3f", cur_ct_dec_vals[i]);
		} // i loop.
		fprintf(stderr, "\n");
		////////////////
#endif

		// Now we have to shift and add; following loop copies the row to all rows in this ciphertext.
		Ciphertext cur_rotated = shifted_ct; // Rotation 0.
		int l_cur_rotation = loaded_nAcols;
		while (l_cur_rotation < (loaded_nAcols * n_rows_per_ct))
		{
			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Processing %d. rotation.\n", l_cur_rotation);
			}

			// We need to shift this right.
			evaluator.rotate_vector(shifted_ct, -l_cur_rotation, pooled_galois_key, cur_rotated);
			evaluator.add_inplace(shifted_ct, cur_rotated);
			l_cur_rotation *= 2;

#ifdef __DECRYPT_DEBUG__
			////////////////
			decryptor.decrypt(shifted_ct, cur_dec_plaintext);
			cur_ct_dec_vals;
			encoder.decode(cur_dec_plaintext, cur_ct_dec_vals);
			fprintf(stderr, "Log-Copying ct for reprow=%d @ l_rot=%d:", i_rep_row, l_cur_rotation);
			for (int i = 0; i < loaded_nAcols * 2; i++)
			{
				fprintf(stderr, " %.3f", cur_ct_dec_vals[i]);
			} // i loop.
			fprintf(stderr, "\n");
			////////////////
#endif
		} // rotate.

		// Copy the shifted ct to all ct's: Note that row expanded data contains different number of cts than the original data.
		int n_expanded_matrix_cts = ceil(((double)(n_rows_per_expanded_matrix * loaded_nAcols)) / n_values_per_ct);
		vector<Ciphertext>* res_cts = new vector<Ciphertext>();
		for (int i_ct = 0; i_ct < n_expanded_matrix_cts; i_ct++)
		{
			Ciphertext new_ct = shifted_ct;
			res_cts->push_back(new_ct);
		} // i_ct.

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//  Should we fixed the last ct to make sure remaining entries are 0? We should not assume anything further down the entries.  //
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Save the matrix.
		char op_fp[1000];
		sprintf(op_fp, "%s/reprow_%d.bin.enc", op_dir, i_rep_row);
		save_continuous_enc_matrix(res_cts, n_rows_per_expanded_matrix, loaded_nAcols, text_params_path, op_fp);
	} // i_rep_row loop.
} // row_expand_continuous_encrypted_matrix function.

void row_expand_dense_encrypt_matrix(char* A_fp, int n_rows_per_expanded_matrix,
	char* text_params_path,
	char* pooled_public_key_path,
	char* op_dir)
{
	int nrows, ncols;
	double** matrix = load_matrix_binary(A_fp, nrows, ncols, NULL);
	fprintf(stderr, "Loaded %dx%d matrix from %s\n", nrows, ncols, A_fp);

	// Replicate the column n_cols_per_expanded_matrix times and write the file.
	for (int i_reprow = 0; i_reprow < nrows; i_reprow++)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Generating %d. reprow matrix.\n", i_reprow);
		}

		double** cur_rowrep_matrix = allocate_matrix(n_rows_per_expanded_matrix, ncols);

		for (int i_row = 0; i_row < n_rows_per_expanded_matrix; i_row++)
		{
			for (int i_col = 0; i_col < ncols; i_col++)
			{
				cur_rowrep_matrix[i_row][i_col] = matrix[i_reprow][i_col];
			} // i_col_loop.
		} // i_row loop.

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			char reprow_mat_fp[1000];
			sprintf(reprow_mat_fp, "%s/reprow_%d.bin", op_dir, i_reprow);
			save_matrix_binary(cur_rowrep_matrix, n_rows_per_expanded_matrix, ncols, reprow_mat_fp);
		}

		// Save the encrypted matrix.
		char enc_mat_fp[1000];
		sprintf(enc_mat_fp, "%s/reprow_%d.bin.enc", op_dir, i_reprow);
		encrypt_plaintext_matrix_continuous_ct(cur_rowrep_matrix, n_rows_per_expanded_matrix, ncols,
			text_params_path, pooled_public_key_path, enc_mat_fp);
	} // i_repcol loop.
}

vector<Ciphertext>* load_continuous_enc_matrix(char* enc_matrix_fp, int& loaded_nrow, int& loaded_ncol, char* text_params_path)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading encrypted matrix from %s..\n", enc_matrix_fp);
	}

	// Load each row to a ciphertext and encrypt.
	ifstream ifs_enc_matrix(enc_matrix_fp, ios::binary);

	// This describes the matrix that is being saved to output.
	int f_nrow = 0;
	int f_ncol = 0;
	ifs_enc_matrix.read(reinterpret_cast<char *>(&f_ncol), sizeof(int));
	ifs_enc_matrix.read(reinterpret_cast<char *>(&(f_nrow)), sizeof(int));

	loaded_nrow = f_nrow;
	loaded_ncol = f_ncol;

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Read dimensions %dx%d from %s..\n", loaded_nrow, loaded_ncol, enc_matrix_fp);
	}

	double n_entries_in_matrix = (double)(loaded_nrow * loaded_ncol);
	int n_cts_2_read = (int)(ceil(n_entries_in_matrix / n_values_per_ct));

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Reading %d cts.\n", n_cts_2_read);
	}


	// Start loading the ct's.
	vector<Ciphertext>* matrix_cts = new vector<Ciphertext>();
	for (int ct_i = 0; ct_i < n_cts_2_read; ct_i++)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Reading %d. ct..\n", ct_i);
		}

		// Encrypt the plaintext.
		Ciphertext ct;

		// Save the encrypted ciphertext.
		ct.load(context, ifs_enc_matrix);

		// Read the ct.
		matrix_cts->push_back(ct);
	} // ct_i loop.

	// Close the file.
	ifs_enc_matrix.close();

	return(matrix_cts);
}


void secure_add_continuous_encrypted_matrices_in_memory_per_list(char* enc_mat_list_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(1);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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

	//// Instantiate and load the public key.
	//fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	//ifstream ifs_pk(pooled_public_key_path, ios::binary);
	//PublicKey pooled_pk;
	//pooled_pk.load(context, ifs_pk);
	//ifs_pk.close();
	//fprintf(stderr, "Loaded public key.\n");

	// Load relin key.
	fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();
	fprintf(stderr, "Loaded relinearization key.\n");

	// We do not need the galois key here.
	//// Load galois key.
	//fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	//ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	//GaloisKeys pooled_galois_key;
	//pooled_galois_key.load(context, ifs_galois_key);
	//ifs_galois_key.close();
	//fprintf(stderr, "Loaded Galois key.\n");

#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find private key file path @ %s\n", pooled_private_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif

	// Instantiate evaluator.
	Evaluator evaluator(context);

	//Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	//CKKSEncoder encoder(context);
	//int n_values_per_ct = encoder.slot_count();

	vector<char*>* mat_fp_list = buffer_file(enc_mat_list_fp);
	if (mat_fp_list->size() == 0)
	{
		fprintf(stderr, "There are no files in %s, make sure to have at least one file.\n", enc_mat_list_fp);
		exit(1);
	}

	fprintf(stderr, "Pooling %d matrices from %s..\n", (int)(mat_fp_list->size()), enc_mat_list_fp);

	//vector<Ciphertext>* temp_cts = load_continuous_enc_matrix(mat_fp_list->at(0), pooled_nrows, pooled_ncols, text_params_path);
	//delete temp_cts;

	// Load each row to a ciphertext and encrypt.
	ifstream ifs_enc_matrix(mat_fp_list->at(0), ios::binary);

	// This describes the matrix that is being saved to output.
	int f_nrow = 0;
	int f_ncol = 0;
	ifs_enc_matrix.read(reinterpret_cast<char*>(&f_ncol), sizeof(int));
	ifs_enc_matrix.read(reinterpret_cast<char*>(&(f_nrow)), sizeof(int));
	ifs_enc_matrix.close();

	int pooled_nrows = f_nrow;
	int pooled_ncols = f_ncol;

	// Save the pooled matrix.
	//double** pooled_matrix = allocate_matrix(pooled_nrows, pooled_ncols);
	//encrypt_plaintext_matrix_continuous_ct(pooled_matrix, pooled_nrows, pooled_ncols, text_params_path, pooled_public_key_path, op_fp);

	vector<Ciphertext>* res_cts = new vector<Ciphertext>();
	//Ciphertext res_ct;

	for (int i_mat = 0; i_mat < (int)(mat_fp_list->size()); i_mat++)
	{
		int loaded_nArows, loaded_nAcols;
		vector<Ciphertext>* A_cts = load_continuous_enc_matrix(mat_fp_list->at(i_mat), loaded_nArows, loaded_nAcols, text_params_path);

		if (loaded_nAcols != pooled_ncols ||
			loaded_nArows != pooled_nrows)
		{
			fprintf(stderr, "Sanity check failed: Matrix dimensions are not conformant for aggregation on list %s (%s): Agg[%dx%d] vs A[%dx%d]\n",
				enc_mat_list_fp, mat_fp_list->at(i_mat),
				pooled_nrows, pooled_ncols,
				loaded_nArows, loaded_nAcols);

			exit(1);
		}

		fprintf(stderr, "Pooling %s\n", mat_fp_list->at(i_mat));

		for (int i_ct = 0; i_ct < (int)(A_cts->size()); i_ct++)
		{
			// Initialize ct's with first matrix.
			if (i_mat == 0)
			{
				Ciphertext cur_res_ct = A_cts->at(i_ct);
				res_cts->push_back(cur_res_ct);
			}
			else
			{
				if (A_cts->at(i_ct).coeff_modulus_size() > res_cts->at(i_ct).coeff_modulus_size())
				{
					evaluator.mod_switch_to_inplace(A_cts->at(i_ct), res_cts->at(i_ct).parms_id());
					A_cts->at(i_ct).scale() = res_cts->at(i_ct).scale();
				}
				else if (A_cts->at(i_ct).coeff_modulus_size() < res_cts->at(i_ct).coeff_modulus_size())
				{
					evaluator.mod_switch_to_inplace(res_cts->at(i_ct), A_cts->at(i_ct).parms_id());
					res_cts->at(i_ct).scale() = A_cts->at(i_ct).scale();
				}

				// Add the ct.
				evaluator.add_inplace(res_cts->at(i_ct), A_cts->at(i_ct));
			}
		} // i_ct loop.

		delete A_cts;
	} // i_mat loop.

	save_continuous_enc_matrix(res_cts, pooled_nrows, pooled_ncols, text_params_path, op_fp);
} // secure_add_continuous_encrypted_matrices_in_memory_per_list option.

void secure_add_continuous_encrypted_matrices_per_list(char* enc_mat_list_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(1);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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

	//// Instantiate and load the public key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	//}

	//ifstream ifs_pk(pooled_public_key_path, ios::binary);
	//PublicKey pooled_pk;
	//pooled_pk.load(context, ifs_pk);
	//ifs_pk.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded public key.\n");
	//}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	//// Load galois key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	//}
	//ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	//GaloisKeys pooled_galois_key;
	//pooled_galois_key.load(context, ifs_galois_key);
	//ifs_galois_key.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded Galois key.\n");
	//}

#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find private key file path @ %s\n", pooled_private_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif

	// Instantiate evaluator.
	Evaluator evaluator(context);

	//Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	//CKKSEncoder encoder(context);
	//int n_values_per_ct = encoder.slot_count();

	vector<char*>* mat_fp_list = buffer_file(enc_mat_list_fp);
	if (mat_fp_list->size() == 0)
	{
		fprintf(stderr, "There are no files in %s, make sure to have at least one file.\n", enc_mat_list_fp);
		exit(1);
	}

	fprintf(stderr, "Pooling %d matrices from %s..\n", (int)(mat_fp_list->size()), enc_mat_list_fp);

	int pooled_nrows, pooled_ncols;
	vector<Ciphertext>* temp_cts = load_continuous_enc_matrix(mat_fp_list->at(0), pooled_nrows, pooled_ncols, text_params_path);
	delete temp_cts;

	// Save the pooled matrix.
	double** pooled_matrix = allocate_matrix(pooled_nrows, pooled_ncols);
	encrypt_plaintext_matrix_continuous_ct(pooled_matrix, pooled_nrows, pooled_ncols, text_params_path, pooled_public_key_path, op_fp);

	for (int i_mat = 0; i_mat < (int)(mat_fp_list->size()); i_mat++)
	{
		fprintf(stderr, "Pooling %s\n", mat_fp_list->at(i_mat));
		secure_add_continuous_encrypted_matrices(op_fp, mat_fp_list->at(i_mat), text_params_path, pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path, pooled_private_key_path, op_fp);
	} // i_mat loop.	
} // secure_add_continuous_encrypted_matrices_per_list option.

void secure_add_continuous_encrypted_matrices(char* enc_A_fp, char* enc_B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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

	//// Instantiate and load the public key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	//}
	//ifstream ifs_pk(pooled_public_key_path, ios::binary);
	//PublicKey pooled_pk;
	//pooled_pk.load(context, ifs_pk);
	//ifs_pk.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded public key.\n");
	//}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}
	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	//// Load galois key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	//}

	//ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	//GaloisKeys pooled_galois_key;
	//pooled_galois_key.load(context, ifs_galois_key);
	//ifs_galois_key.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded Galois key.\n");
	//}

#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find private key path @ %s\n", pooled_private_key_path);
		exit(0);
	}
	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif

	// Instantiate evaluator.
	Evaluator evaluator(context);

	//Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading A and B matrices.\n");
	}

	int loaded_nArows, loaded_nAcols;
	vector<Ciphertext>* A_cts = load_continuous_enc_matrix(enc_A_fp, loaded_nArows, loaded_nAcols, text_params_path);

	int loaded_nBrows, loaded_nBcols;
	vector<Ciphertext>* B_cts = load_continuous_enc_matrix(enc_B_fp, loaded_nBrows, loaded_nBcols, text_params_path);

	if (A_cts->size() != B_cts->size())
	{
		fprintf(stderr, "Sanity check error: ct size is not conformant among matrices: %d, %d\n", (int)(A_cts->size()), (int)(B_cts->size()));
		exit(1);
	}

	if (loaded_nArows != loaded_nBrows ||
		loaded_nAcols != loaded_nBcols)
	{
		fprintf(stderr, "Sanity check error: Matrix size is not conformant among matrices: %d, %d; %d, %d\n", loaded_nArows, loaded_nAcols, loaded_nBrows, loaded_nBcols);
		exit(1);
	}

	double n_entries_per_matrix = (double)(loaded_nAcols * loaded_nArows);
	int n_exp_cts = (int)(ceil(n_entries_per_matrix / n_values_per_ct));

	if (n_exp_cts != (int)(A_cts->size()))
	{
		fprintf(stderr, "Sanity check error: Expected ct size is not conformant for %dx%d matrices: %d vs %d\n", loaded_nArows, loaded_nAcols, n_exp_cts, (int)A_cts->size());
		exit(1);
	}

	fprintf(stderr, "Adding %dx%d matrices over %d cts..\n", loaded_nArows, loaded_nAcols, (int)(A_cts->size()));

	vector<Ciphertext>* res_cts = new vector<Ciphertext>();
	for (int i_ct = 0; i_ct < (int)(A_cts->size()); i_ct++)
	{
		Ciphertext res_ct;

		if (A_cts->at(i_ct).coeff_modulus_size() > B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(A_cts->at(i_ct), B_cts->at(i_ct).parms_id());
			A_cts->at(i_ct).scale() = B_cts->at(i_ct).scale();
		}
		else if (A_cts->at(i_ct).coeff_modulus_size() < B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(B_cts->at(i_ct), A_cts->at(i_ct).parms_id());
			B_cts->at(i_ct).scale() = A_cts->at(i_ct).scale();
		}

		evaluator.add(A_cts->at(i_ct), B_cts->at(i_ct), res_ct);

		res_cts->push_back(res_ct);
	} // i_ct loop.

	save_continuous_enc_matrix(res_cts, loaded_nArows, loaded_nAcols, text_params_path, op_fp);
} // secure_add_continuous_encrypted_matrices option.

void secure_subtract_continuous_encrypted_matrices(char* A_fp, char* B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	//	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	if (!check_file(pooled_public_key_path) ||
		!check_file(pooled_relin_key_path) ||
		!check_file(pooled_galois_key_path) ||
		!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s, %s, %s\n", pooled_public_key_path, pooled_relin_key_path, pooled_galois_key_path);
		exit(0);
	}

	//// Instantiate and load the public key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	//}

	//ifstream ifs_pk(pooled_public_key_path, ios::binary);
	//PublicKey pooled_pk;
	//pooled_pk.load(context, ifs_pk);
	//ifs_pk.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded public key.\n");
	//}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	//// Load galois key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	//}

	//ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	//GaloisKeys pooled_galois_key;
	//pooled_galois_key.load(context, ifs_galois_key);
	//ifs_galois_key.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded Galois key.\n");
	//}

	// Load secret key -- this should be the actual pooled key.
#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_private_key_path);
		exit(0);
	}

	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif

	// Instantiate evaluator.
	Evaluator evaluator(context);

	//Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading A and B matrices.\n");
	}

	int loaded_nArows, loaded_nAcols;
	vector<Ciphertext>* A_cts = load_continuous_enc_matrix(A_fp, loaded_nArows, loaded_nAcols, text_params_path);

	int loaded_nBrows, loaded_nBcols;
	vector<Ciphertext>* B_cts = load_continuous_enc_matrix(B_fp, loaded_nBrows, loaded_nBcols, text_params_path);

	if (A_cts->size() != B_cts->size())
	{
		fprintf(stderr, "Sanity check error: ct size is not conformant among matrices: %d, %d\n", (int)(A_cts->size()), (int)(B_cts->size()));
		exit(0);
	}

	if (loaded_nArows != loaded_nBrows ||
		loaded_nAcols != loaded_nBcols)
	{
		fprintf(stderr, "Sanity check error: matrix size is not conformant among matrices: %d, %d; %d, %d\n", loaded_nArows, loaded_nAcols, loaded_nBrows, loaded_nBcols);
		exit(0);
	}

	double n_entries_per_matrix = (double)(loaded_nAcols * loaded_nArows);
	int n_exp_cts = (int)(ceil(n_entries_per_matrix / n_values_per_ct));

	if (n_exp_cts != (int)(A_cts->size()))
	{
		fprintf(stderr, "Sanity check error: Expected ct size is not conformant for %dx%d matrices: %d vs %d\n", loaded_nArows, loaded_nAcols, n_exp_cts, (int)A_cts->size());
		exit(0);
	}

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Adding %dx%d matrices over %d cts..\n", loaded_nArows, loaded_nAcols, (int)(A_cts->size()));
	}

	vector<Ciphertext>* res_cts = new vector<Ciphertext>();
	for (int i_ct = 0; i_ct < (int)(A_cts->size()); i_ct++)
	{
		Ciphertext res_ct;

		if (A_cts->at(i_ct).coeff_modulus_size() > B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(A_cts->at(i_ct), B_cts->at(i_ct).parms_id());
			A_cts->at(i_ct).scale() = B_cts->at(i_ct).scale();
		}
		else if (A_cts->at(i_ct).coeff_modulus_size() < B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(B_cts->at(i_ct), A_cts->at(i_ct).parms_id());
			B_cts->at(i_ct).scale() = A_cts->at(i_ct).scale();
		}

		evaluator.sub(A_cts->at(i_ct), B_cts->at(i_ct), res_ct);

		res_cts->push_back(res_ct);
	} // i_ct loop.

	save_continuous_enc_matrix(res_cts, loaded_nArows, loaded_nAcols, text_params_path, op_fp);
} // secure_subtract_continuous_encrypted_matrices option.

void secure_multiply_elementwise_continuous_encrypted_matrices(char* A_fp, char* B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
	// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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

	//// Instantiate and load the public key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	//}

	//ifstream ifs_pk(pooled_public_key_path, ios::binary);
	//PublicKey pooled_pk;
	//pooled_pk.load(context, ifs_pk);
	//ifs_pk.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded public key.\n");
	//}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	//// Load galois key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	//}

	//ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	//GaloisKeys pooled_galois_key;
	//pooled_galois_key.load(context, ifs_galois_key);
	//ifs_galois_key.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded Galois key.\n");
	//}

	// Load secret key -- this should be the actual pooled key.
#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_private_key_path);
		exit(0);
	}

	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif 

	// Instantiate evaluator.
	Evaluator evaluator(context);

	//Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading A and B matrices.\n");
	}

	int loaded_nArows, loaded_nAcols;
	vector<Ciphertext>* A_cts = load_continuous_enc_matrix(A_fp, loaded_nArows, loaded_nAcols, text_params_path);

	int loaded_nBrows, loaded_nBcols;
	vector<Ciphertext>* B_cts = load_continuous_enc_matrix(B_fp, loaded_nBrows, loaded_nBcols, text_params_path);

	if (A_cts->size() != B_cts->size())
	{
		fprintf(stderr, "Sanity check failed: ct size is not conformant among matrices: %d, %d\n", (int)(A_cts->size()), (int)(B_cts->size()));
		exit(0);
	}

	if (loaded_nArows != loaded_nBrows ||
		loaded_nAcols != loaded_nBcols)
	{
		fprintf(stderr, "matrix size is not conformant among matrices: %d, %d; %d, %d\n", loaded_nArows, loaded_nAcols, loaded_nBrows, loaded_nBcols);
		exit(0);
	}

	double n_entries_per_matrix = (double)(loaded_nAcols * loaded_nArows);
	int n_exp_cts = (int)(ceil(n_entries_per_matrix / n_values_per_ct));

	if (n_exp_cts != (int)(A_cts->size()))
	{
		fprintf(stderr, "Sanity check error: Expected ct size is not conformant for %dx%d matrices: %d vs %d\n", loaded_nArows, loaded_nAcols, n_exp_cts, (int)A_cts->size());
		exit(0);
	}

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Multiplying %dx%d matrices over %d cts..\n", loaded_nArows, loaded_nAcols, (int)(A_cts->size()));
	}

	vector<Ciphertext>* res_cts = new vector<Ciphertext>();
	for (int i_ct = 0; i_ct < (int)(A_cts->size()); i_ct++)
	{
		// First, make sure the ct's match in scale and mod chain.
		if (A_cts->at(i_ct).coeff_modulus_size() > B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(A_cts->at(i_ct), B_cts->at(i_ct).parms_id());
			A_cts->at(i_ct).scale() = B_cts->at(i_ct).scale();
		}
		else if (A_cts->at(i_ct).coeff_modulus_size() < B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(B_cts->at(i_ct), A_cts->at(i_ct).parms_id());
			B_cts->at(i_ct).scale() = A_cts->at(i_ct).scale();
		}

		Ciphertext res_ct;
		evaluator.multiply(A_cts->at(i_ct), B_cts->at(i_ct), res_ct);

		// We must relinearize the result.
		evaluator.relinearize_inplace(res_ct, pooled_relin_key);
		evaluator.rescale_to_next_inplace(res_ct);

		res_cts->push_back(res_ct);
	} // i_ct loop.

	save_continuous_enc_matrix(res_cts, loaded_nArows, loaded_nAcols, text_params_path, op_fp);
} // secure_multiply_elementwise_continuous_encrypted_matrices option.

void secure_row2row_inner_prod_continuous_encrypted_matrices(char* A_fp, char* B_fp,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Sanity check error: Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	// Load galois key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	}

	ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	GaloisKeys pooled_galois_key;
	pooled_galois_key.load(context, ifs_galois_key);
	ifs_galois_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded Galois key.\n");
	}

#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_private_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif 

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading A and B matrices.\n");
	}

	int loaded_nArows, loaded_nAcols;
	vector<Ciphertext>* A_cts = load_continuous_enc_matrix(A_fp, loaded_nArows, loaded_nAcols, text_params_path);

	int loaded_nBrows, loaded_nBcols;
	vector<Ciphertext>* B_cts = load_continuous_enc_matrix(B_fp, loaded_nBrows, loaded_nBcols, text_params_path);

	if (A_cts->size() != B_cts->size())
	{
		fprintf(stderr, "Sanity check error: ct size is not conformant among matrices: %d, %d\n", (int)(A_cts->size()), (int)(B_cts->size()));
		exit(0);
	}

	if (loaded_nArows != loaded_nBrows ||
		loaded_nAcols != loaded_nBcols)
	{
		fprintf(stderr, "Sanity check error: matrix size is not conformant among matrices: %d, %d; %d, %d\n", loaded_nArows, loaded_nAcols, loaded_nBrows, loaded_nBcols);
		exit(0);
	}

	double n_entries_per_matrix = (double)(loaded_nAcols * loaded_nArows);
	int n_exp_cts = (int)(ceil(n_entries_per_matrix / n_values_per_ct));

	if (n_exp_cts != (int)(A_cts->size()))
	{
		fprintf(stderr, "Sanity check error: Expected ct size is not conformant for %dx%d matrices: %d vs %d\n", loaded_nArows, loaded_nAcols, n_exp_cts, (int)A_cts->size());
		exit(0);
	}

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Adding %dx%d matrices over %d cts..\n", loaded_nArows, loaded_nAcols, (int)(A_cts->size()));
	}

	// Make sure that the number of columns is a power of two.
	double log2_nAcol = log(loaded_nAcols) / log(2);

	if (pow(2.0, log2_nAcol) != loaded_nAcols)
	{
		fprintf(stderr, "Sanity check error: # of columns must be a power of two, read %d.. (Pad with 0's?)\n", loaded_nAcols);
		exit(0);
	}

	// We currently cannot process matrices whose rows expand over ciphertext boundaries. The reason is that this would require summing upto mid ciphertexts, etc.
	if (loaded_nAcols > n_values_per_ct)
	{
		fprintf(stderr, "Sanity check error: Cannot process %d-column matrices with larger than %d long roaw.\n", loaded_nAcols, n_values_per_ct);
		exit(0);
	}

	if (loaded_nArows > n_values_per_ct)
	{
		fprintf(stderr, "Sanity check error: Cannot process matrices with larger than %d rows.\n", loaded_nArows);
		exit(0);
	}

	// This is a masker array that is used to put the enctries at the correct place.
	vector<double> masker_array;
	for (int cur_cipher_data_pt_i = 0; cur_cipher_data_pt_i < n_values_per_ct; cur_cipher_data_pt_i++)
	{
		masker_array.push_back(0);
	} // cur_cipher_data_pt_i loop.

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Performing elementwise row2row product for %d rows..\n", loaded_nArows);
	}

	Ciphertext final_row2row_multsum_res_ct; // We assume there is only one ct that is holding the result, i.e., nrows < n_slots.
	encryptor.encrypt_zero(final_row2row_multsum_res_ct);
	final_row2row_multsum_res_ct.scale() = scale;

	int row_res_ct_data_pt_i = 0; // This holds the row index we are setting in the result ct.
	for (int i_ct = 0; 
		(row_res_ct_data_pt_i < loaded_nArows) && (i_ct < (int)(A_cts->size()));
		i_ct++)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "At %d. ct..\n", i_ct);
		}

		if (A_cts->at(i_ct).coeff_modulus_size() > B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(A_cts->at(i_ct), B_cts->at(i_ct).parms_id());
			A_cts->at(i_ct).scale() = B_cts->at(i_ct).scale();
		}
		else if (A_cts->at(i_ct).coeff_modulus_size() < B_cts->at(i_ct).coeff_modulus_size())
		{
			evaluator.mod_switch_to_inplace(B_cts->at(i_ct), A_cts->at(i_ct).parms_id());
			B_cts->at(i_ct).scale() = A_cts->at(i_ct).scale();
		}

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Elementwise multiplying and relinearizing ct_i=%d.\n", i_ct);
		}

		Ciphertext elem_mult_ct;
		evaluator.multiply(A_cts->at(i_ct), B_cts->at(i_ct), elem_mult_ct);
		evaluator.relinearize_inplace(elem_mult_ct, pooled_relin_key);
		evaluator.rescale_to_next_inplace(elem_mult_ct);

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Rotating and summating over %d columns..\n", loaded_nAcols);
		}

		Ciphertext cur_rotated = elem_mult_ct; // 2^n rotations.
		int l_rotation = 1; // 2^0 rotations.
		while (l_rotation < loaded_nAcols)
		{
			evaluator.rotate_vector(elem_mult_ct, l_rotation, pooled_galois_key, cur_rotated);

			evaluator.add_inplace(elem_mult_ct, cur_rotated);

			l_rotation *= 2;
		} // rotation loop.

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Processing summations from rotations..\n");
		}

		// Go over all the rows in the res_ct and mask/assign them.
		for (int mask_i = 0; 
			(row_res_ct_data_pt_i < loaded_nArows) && (mask_i < n_values_per_ct);
			mask_i += loaded_nAcols)
		{
			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "@ %d. masking index..\n", mask_i);
			}

			// This check does an early stop if we are past the number of rows in the matrix, i.e., we are processing unused values.
			if (row_res_ct_data_pt_i > loaded_nArows)
			{
				break;
			}

			// Reset the mask.
			std::fill(masker_array.begin(), masker_array.end(), 0);
			masker_array.at(mask_i) = 1;

			// Encrypt the data points and save to file.
			Plaintext cur_entry_masker_array_pt;

			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Encoding mask array @ mask_i=%d, res_row_i: %d\n", mask_i, row_res_ct_data_pt_i);
			}

			// Encode the vector.
			encoder.encode(masker_array, scale, cur_entry_masker_array_pt);
			evaluator.mod_switch_to_inplace(cur_entry_masker_array_pt, elem_mult_ct.parms_id());
			cur_entry_masker_array_pt.scale() = elem_mult_ct.scale();

			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Multiplying with the mask array @ mask_i=%d, res_row_i: %d\n", mask_i, row_res_ct_data_pt_i);
			}

			// Multiple with the mask.
			Ciphertext masked_res_ct;
			evaluator.multiply_plain(elem_mult_ct, cur_entry_masker_array_pt, masked_res_ct);
			evaluator.relinearize_inplace(masked_res_ct, pooled_relin_key);
			evaluator.rescale_to_next_inplace(masked_res_ct);

			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Rotating result to %d to access %d..\n", mask_i - row_res_ct_data_pt_i, row_res_ct_data_pt_i);
			}

			// Add the current masked value to the current result: First, move it to the correct position.
			evaluator.rotate_vector_inplace(masked_res_ct, mask_i - row_res_ct_data_pt_i, pooled_galois_key);
			evaluator.mod_switch_to_inplace(final_row2row_multsum_res_ct, masked_res_ct.parms_id());
			final_row2row_multsum_res_ct.scale() = masked_res_ct.scale();

			//fprintf(stderr, "Adding result.\n");
			evaluator.add_inplace(final_row2row_multsum_res_ct, masked_res_ct);

			row_res_ct_data_pt_i++;
		} // mask_i loop.
	} // i_ct loop.

	vector<Ciphertext>* final_res_cts = new vector<Ciphertext>();
	final_res_cts->push_back(final_row2row_multsum_res_ct);
	save_continuous_enc_matrix(final_res_cts, loaded_nArows, 1, text_params_path, op_fp);
} // secure_row2row_inner_prod_continuous_encrypted_matrices

/*
A[axb], B[bxc]
Take each column of A, [ax1] colwise-copy them c times to generate axc -- this generates b matrices
Take each row of B, [1xc] rowwise-copy them a times to generate axc -- this generates b matrices

Encoding:
Each ct contains the concatenation of the rows of expanded matrices.
In computation, multiply the cts, then shift and copy the va

. We do not need to store the number of ct's or number of rows per ct.
*/
void secure_multiply_matrices_Acol_Brow_expansions(char* A_dir, char* B_dir,
	char* text_params_path,
	char* pooled_public_key_path,
	char* pooled_relin_key_path,
	char* pooled_galois_key_path,
	char* pooled_private_key_path, // To be removed.
	char* op_fp)
{
	// Load the per site matrices, compute the 4th power of each matrix element, save the results.
		// Setup context.
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find the text parameters @ %s\n", text_params_path);
		exit(0);
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

	double scale = pow(2, text_params->scale_n_bits);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
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
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading public key from %s.\n", pooled_public_key_path);
	}

	ifstream ifs_pk(pooled_public_key_path, ios::binary);
	PublicKey pooled_pk;
	pooled_pk.load(context, ifs_pk);
	ifs_pk.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded public key.\n");
	}

	// Load relin key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading relinearization key from %s.\n", pooled_relin_key_path);
	}

	ifstream ifs_relin_key(pooled_relin_key_path, ios::binary);
	RelinKeys pooled_relin_key;
	pooled_relin_key.load(context, ifs_relin_key);
	ifs_relin_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded relinearization key.\n");
	}

	//// Load galois key.
	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loading galois key from %s.\n", pooled_galois_key_path);
	//}

	//ifstream ifs_galois_key(pooled_galois_key_path, ios::binary);
	//GaloisKeys pooled_galois_key;
	//pooled_galois_key.load(context, ifs_galois_key);
	//ifs_galois_key.close();

	//if (__DUMP_MATRIX_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Loaded Galois key.\n");
	//}

#ifdef __DECRYPT_DEBUG__
	if (!check_file(pooled_private_key_path))
	{
		fprintf(stderr, "Could not find one of the keys @ %s\n", pooled_private_key_path);
		exit(0);
	}

	// Load secret key -- this should be the actual pooled key.
	fprintf(stderr, "Loading pooled secret key from %s.\n", pooled_private_key_path);
	SecretKey pooled_sk;
	ifstream ifs_pooled_sk(pooled_private_key_path, ios::binary);
	pooled_sk.load(context, ifs_pooled_sk);
	ifs_pooled_sk.close();
	fprintf(stderr, "Loaded pooled secret key.\n");

	// Instantiate the decryptor and encoder.
	Decryptor decryptor(context, pooled_sk);
#endif 

	// Instantiate evaluator.
	Evaluator evaluator(context);

	Encryptor encryptor(context, pooled_pk);

	//// We need the encoder to get the number of samples per ct.
	//CKKSEncoder encoder(context);
	////int n_values_per_ct = encoder.slot_count();

	// Load first matrix and get dimensions of the multiplication.
	char first_mat_fp[1000];
	sprintf(first_mat_fp, "%s/repcol_0.bin.enc", A_dir);
	if (!check_file(first_mat_fp))
	{
		fprintf(stderr, "Could not find first matrix file @ %s\n", first_mat_fp);
		exit(0);
	}

	int loaded_nrow, loaded_ncol;
	vector<Ciphertext>* cts = load_continuous_enc_matrix(first_mat_fp, loaded_nrow, loaded_ncol, text_params_path);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded %d cts per %dx%d matrix.\n", (int)(cts->size()), loaded_nrow, loaded_ncol);
	}

	int n_cts_per_exp_matrix = (int)(cts->size());

	vector<Ciphertext>* res_cts = new vector<Ciphertext>();
	for (int i_ct = 0; i_ct < (int)(cts->size()); i_ct++)
	{
		Ciphertext cur_entry_ct;
		encryptor.encrypt_zero(cur_entry_ct);
		cur_entry_ct.scale() = scale;

		res_cts->push_back(cur_entry_ct);
	}

	for (int i_exp = 0; ; i_exp++)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Processing %d. expansion for A and B.\n", i_exp);
		}

		char cur_Acol_exp_fp[1000];
		sprintf(cur_Acol_exp_fp, "%s/repcol_%d.bin.enc", A_dir, i_exp);

		if (!check_file(cur_Acol_exp_fp))
		{
			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Processed %d expanded matrices.\n", i_exp);
			}
			break;
		}

		int loaded_nrow_A, loaded_ncol_A;
		vector<Ciphertext>* A_cts = load_continuous_enc_matrix(cur_Acol_exp_fp, loaded_nrow_A, loaded_ncol_A, text_params_path);

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Loaded %d cts from %s for matrix [%dx%d]\n", (int)(A_cts->size()), cur_Acol_exp_fp, loaded_nrow_A, loaded_ncol_A);
		}

		char cur_Brow_exp_fp[1000];
		sprintf(cur_Brow_exp_fp, "%s/reprow_%d.bin.enc", B_dir, i_exp);

		int loaded_nrow_B, loaded_ncol_B;
		vector<Ciphertext>* B_cts = load_continuous_enc_matrix(cur_Brow_exp_fp, loaded_nrow_B, loaded_ncol_B, text_params_path);

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Loaded %d cts from %s for matrix [%dx%d]\n", (int)(B_cts->size()), cur_Brow_exp_fp, loaded_nrow_B, loaded_ncol_B);
		}

		if (loaded_nrow_A != loaded_nrow_B ||
			loaded_ncol_A != loaded_ncol_B)
		{
			fprintf(stderr, "Sanity check failed while processing %s, %s: %d, %d; %d, %d\n", cur_Acol_exp_fp, cur_Brow_exp_fp, loaded_nrow_A, loaded_ncol_A, loaded_nrow_B, loaded_ncol_B);
			exit(0);
		}

		if ((int)(A_cts->size()) != n_cts_per_exp_matrix ||
			(int)(B_cts->size()) != n_cts_per_exp_matrix)
		{
			fprintf(stderr, "The ct count is not conformant: %d, %d, %d\n", n_cts_per_exp_matrix, (int)(A_cts->size()), (int)(B_cts->size()));
			exit(0);
		}

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Loaded %s, %s; %d cts\n", cur_Acol_exp_fp, cur_Brow_exp_fp, (int)(B_cts->size()));
		}

		for (int i_ct = 0; i_ct < (int)(A_cts->size()); i_ct++)
		{
			Ciphertext cur_ct_mult;

			// If the ct's are not at the same level, we need to relevel them.
			if (A_cts->at(i_ct).coeff_modulus_size() > B_cts->at(i_ct).coeff_modulus_size())
			{
				evaluator.mod_switch_to_inplace(A_cts->at(i_ct), B_cts->at(i_ct).parms_id());
			}
			else if (A_cts->at(i_ct).coeff_modulus_size() < B_cts->at(i_ct).coeff_modulus_size())
			{
				evaluator.mod_switch_to_inplace(B_cts->at(i_ct), A_cts->at(i_ct).parms_id());
			}

			evaluator.multiply(A_cts->at(i_ct), B_cts->at(i_ct), cur_ct_mult);
			evaluator.relinearize_inplace(cur_ct_mult, pooled_relin_key);
			evaluator.rescale_to_next_inplace(cur_ct_mult);

			evaluator.mod_switch_to_inplace(res_cts->at(i_ct), cur_ct_mult.parms_id());
			res_cts->at(i_ct).scale() = cur_ct_mult.scale();

			evaluator.add_inplace(res_cts->at(i_ct), cur_ct_mult);
		} // i_ct loop.
	} // i_mats loop.

	// Write the result.
	save_continuous_enc_matrix(res_cts, loaded_nrow, loaded_ncol, text_params_path, op_fp);
} // secure_multiply_matrices_Acol_Brow_expansions option.

// site0 flag is necessary for keeping track of c0 term in pooling step.
void partial_decrypt_continuous_enc_per_noisy_secretkey_w_smdgng_noise(char* encrypted_data_path, bool is_site0, char* text_params_path, char* noisy_key_path, int n_bits_per_max_smdging_noise, char* partial_decrypt_data_path)
{
	if (!check_file(text_params_path))
	{
		fprintf(stderr, "Could not find text params file @ %s\n", text_params_path);
		exit(0);
	}

	t_text_params* text_params = load_text_params(text_params_path);

	// If the smdging noise is set to 0, reset it to default value of key noise variance.
	const int DEFAULT_SMDGING_NOISE_BIT_SIZE = 40; // This is not the size of uniform noise, it defines the variance of the discrete normal distribution.
	if (n_bits_per_max_smdging_noise == 0)
	{
		fprintf(stderr, "Setting smdging noise level to default of %d-bits.\n", DEFAULT_SMDGING_NOISE_BIT_SIZE);
		n_bits_per_max_smdging_noise = DEFAULT_SMDGING_NOISE_BIT_SIZE;
	}

	// Max the smdging noise size at 60-bits.
	if (n_bits_per_max_smdging_noise > 60)
	{
		n_bits_per_max_smdging_noise = 60;
	}

	fprintf(stderr, "Partial decrypting with noisy secret key at %s for site (is_site0=%d), using uniform smdg-ing noise bit size of %d and saving to %s\n",
		noisy_key_path, is_site0, (int)n_bits_per_max_smdging_noise, partial_decrypt_data_path);

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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);
	CKKSEncoder encoder(context);

	// We store both of these but use only one.
	auto &first_context_data = *context.first_context_data();
	auto &key_context_data = *context.key_context_data();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "first_ctx: %d; key_ctx: %d\n", (int)first_context_data.chain_index(), (int)key_context_data.chain_index());
	}

	fprintf(stderr, "Generating the smudging noise @ %d-bits.\n", n_bits_per_max_smdging_noise);
	double smdg_noise_variance = pow(2, n_bits_per_max_smdging_noise);
	double MAX_SMDGNG_ERROR = 10 * smdg_noise_variance; // 10 times the noise variance.
	//uint64_t* SG_generate_smudging_noise(int poly_modulus_degree, int coeff_modulus_size, double smdging_noise_variance, const vector<Modulus> coeff_modulus);
	//uint64_t* smdging_noise_levels = SG_generate_smudging_noise(poly_modulus_degree, parms.coeff_modulus().size(), smdg_noise_variance, parms.coeff_modulus());
	uint64_t* smdging_noise_levels = SEAL_Genomics_Get_Normal_Random_Custom(first_context_data.parms(), pow(smdg_noise_variance, 0.5), MAX_SMDGNG_ERROR);
	//uint64_t* smdging_noise_levels = SEAL_Genomics_Get_Unif_Random_Custom(first_context_data.parms(), n_bits_per_max_smdging_noise);

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Smudging noise levels: ");
		//for (size_t coeff_i = 0; coeff_i < poly_modulus_degree; coeff_i++)
		for (size_t coeff_i = 0; coeff_i < 20; coeff_i++)
		{
			for (size_t mod_i = 0; mod_i < first_context_data.parms().coeff_modulus().size(); mod_i++)
			{
				fprintf(stderr, "coeff: %lu / mod: %lu: %lu\n", coeff_i, mod_i, smdging_noise_levels[mod_i * poly_modulus_degree + coeff_i]);
			} // mod_i loop.
		} // coeff_i
	}

	// Convert the smdging noise into NTT representation.
	// It is very important to match the modulus sizes to the NTT tables.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Converting the smudging noise to NTT representation: %d; %d; %lu\n", (int)parms.coeff_modulus().size(), (int)poly_modulus_degree, (long unsigned int)(void*)((key_context_data.small_ntt_tables())));
	}
	//NTT_transform_coefficients_array_in_place(smdging_noise_levels, key_context_data.parms().coeff_modulus().size(), poly_modulus_degree, key_context_data.small_ntt_tables());
	NTT_transform_coefficients_array_in_place(smdging_noise_levels, first_context_data.parms().coeff_modulus().size(), poly_modulus_degree, first_context_data.small_ntt_tables());

	int n_entries_per_ct = encoder.slot_count();

	// Load the secret (noisy) key.
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading noisy secret key from %s\n", noisy_key_path);
	}

	ifstream ifs_noisy_key(noisy_key_path, ios::binary);
	SecretKey dec_sk;
	dec_sk.load(context, ifs_noisy_key);
	ifs_noisy_key.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded noisy secret key.\n");
	}

	// Instantiate the decryptor.
	Decryptor decryptor(context, dec_sk);

	// Just load the ciphertexts from the results and decrypt them with this site's noisy key.
	ifstream ifs_enc_data_matrix(encrypted_data_path, ios::binary);
	ofstream ofs_dec_data_matrix(partial_decrypt_data_path, ios::binary);

	// Read the encrypted data matrix parameters.
	int n_rows;
	int n_cols;

	// Read the matrix size parameters.
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_cols), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_rows), sizeof(int));

	// Write the parameters exactly to the partially decrypted file.
	ofs_dec_data_matrix.write(reinterpret_cast<char *>(&n_cols), sizeof(int));
	ofs_dec_data_matrix.write(reinterpret_cast<char *>(&n_rows), sizeof(int));

	double total_entries_per_site = (double)(n_rows * n_cols);
	int n_cts_per_matrix = (int)(ceil(total_entries_per_site / n_entries_per_ct));

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Decrypting %d data points for a %dx%d matrix by reading every ciphertext and saving to %s.\n", (int)total_entries_per_site, n_rows, n_cols, partial_decrypt_data_path);
	}

	// Start loading ciphertexts and decrypt them, then save to the partial decrypted results. For site0, include the c0 term in plaintexts, for others, subtract them from decrypted data.
	for (int ct_i = 0; ct_i < n_cts_per_matrix; ct_i++)
	{
		Ciphertext cur_ciphertext;
		cur_ciphertext.load(context, ifs_enc_data_matrix);

		if (cur_ciphertext.size() != 2)
		{
			fprintf(stderr, "Must make sure the ct's are of size 2, it seems like previous processing did not relinearize these ciphertexts. Currently the CT has size of %d.\n", (int)(cur_ciphertext.size()));
			exit(0);
		}

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

			if (__DUMP_MATRIX_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Ciphertext params for ct_i=%d:\ncoeff_modulus_size: %d\n", ct_i, (int)(ct_parms.coeff_modulus().size()));
			}

			// We have to remove c0 by going over each decomposition level; the levels are coming from each ciphertext individually.
			for (int i_mod = 0; i_mod < (int)ct_coeff_modulus_size; i_mod++)
			{
				uint64_t* cur_mod_decr_res_ptr = cur_dec_plaintext.data() + i_mod * ct_coeff_count;

				// c0 component must be added only once. Following removes the c0 component from the pooled decryption.
				sub_poly_coeffmod(cur_mod_decr_res_ptr,
					cur_ciphertext.data(0) + i_mod * ct_coeff_count,
					ct_coeff_count, ct_parms.coeff_modulus().at(i_mod),
					cur_mod_decr_res_ptr);

				// Add the smdging noise here for this decryption instance:
				// This step aims at flooding or at least disrupting the current decryption to a certain extent independent of the secret key of this site that is performing decryption.
				add_poly_coeffmod(cur_mod_decr_res_ptr,
					smdging_noise_levels + i_mod * ct_coeff_count,
					ct_coeff_count, ct_parms.coeff_modulus().at(i_mod),
					cur_mod_decr_res_ptr);
			} // i_mod loop.
		} // is_site0 check.

		// Save the current partial decrypted data to the output stream.
		cur_dec_plaintext.save(ofs_dec_data_matrix);
	} // ct_i loop.

	// Close the files.
	ofs_dec_data_matrix.close();
	ifs_enc_data_matrix.close();
} //  partial_decrypt_continuous_enc_per_noisy_secretkey

void collaborative_pool_partial_decrypted_continuous_enc_plaintext_results(char* partial_decrypted_paths_list_path, char* text_params_path, char* op_fp)
{
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Pooling partially decrypted data (%s) and saving to %s\n", partial_decrypted_paths_list_path, op_fp);
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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);
	CKKSEncoder encoder(context);

	int n_entries_per_ct = encoder.slot_count();

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
	vector<int>* per_site_n_row = new vector<int>();
	vector<int>* per_site_n_col = new vector<int>();

	// These three parameters have to match exactly while pooling.
	size_t overall_n_row = 0;
	size_t overall_n_col = 0;

	// Load the matrix size parameters for each site.
	vector<ifstream*>* per_site_part_dec_data_ifs = new vector<ifstream*>();
	for (int i_site = 0; i_site < n_sites; i_site++)
	{
		ifstream* cur_file_ptr = new ifstream();
		cur_file_ptr->open(per_site_partial_decrypted_paths->at(i_site), ios::binary);
		per_site_part_dec_data_ifs->push_back(cur_file_ptr);

		int cur_site_n_row;
		int cur_site_n_col;

		cur_file_ptr->read(reinterpret_cast<char *>(&cur_site_n_col), sizeof(int));
		cur_file_ptr->read(reinterpret_cast<char *>(&cur_site_n_row), sizeof(int));

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "i_site: %d: n_row: %d; n_col: %d;\n", i_site, cur_site_n_row, cur_site_n_col);
		}

		per_site_n_col->push_back(cur_site_n_col);
		per_site_n_row->push_back(cur_site_n_row);

		if (overall_n_row == 0)
		{
			overall_n_row = cur_site_n_row;
			overall_n_col = cur_site_n_col;
		}
		else if ((int)overall_n_row != cur_site_n_row ||
			(int)overall_n_col != cur_site_n_col)
		{
			fprintf(stderr, "Sanity check error: Site %d has a different matrix format for pooling: (%d, %d) vs (%d, %d)\n",
				i_site,
				cur_site_n_row, cur_site_n_col,
				(int)overall_n_row, (int)overall_n_col);

			exit(0);
		}
	} // i_site loop.

	double total_entries_per_site = (double)(overall_n_col * overall_n_row);
	int n_pts_2_read_per_site = (int)(ceil(total_entries_per_site / n_entries_per_ct));

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Open the encoded plaintext matrix file for saving.
	char encoded_pt_matrix_op_fp[1000];
	sprintf(encoded_pt_matrix_op_fp, "%s_encoded.bin", op_fp);
	ofstream ofs_enc_data_matrix(encoded_pt_matrix_op_fp, ios::binary);

	ofs_enc_data_matrix.write(reinterpret_cast<const char *>(&overall_n_col), sizeof(int));
	ofs_enc_data_matrix.write(reinterpret_cast<const char *>(&overall_n_row), sizeof(int));
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read all vectors and write the final output.
	//for (int i_vector = 0; i_vector < overall_n_row; i_vector++)
	//vector<Plaintext>* pooled_pts = new vector<Plaintext>();
	vector<double>* pooled_data_values = new vector<double>();
	for (int i_pt = 0; i_pt < n_pts_2_read_per_site; i_pt++)
	{
		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Pooling %d. pt.\n", i_pt);
		}

		//int cur_vector_pt_i = 0;

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

			// If this is the first site, set the pooled plaintext value; this re-allocates the pooled plaintext.
			if (i_site == 0)
			{
				cur_pooled_plaintext = cur_site_plaintext;
			}
			else
			{
				for (int i_mod = 0; i_mod < (int)ct_coeff_modulus_size; i_mod++)
				{
					//fprintf(stderr, "Pooling %d. modulus: coeff_count: %d\n", i_mod, pl_coeff_count);

					add_poly_coeffmod(cur_pooled_plaintext.data() + i_mod * pl_coeff_count,
						cur_site_plaintext.data() + i_mod * pl_coeff_count,
						pl_coeff_count,
						pl_parms.coeff_modulus().at(i_mod),
						cur_pooled_plaintext.data() + i_mod * pl_coeff_count);
				} // i_dim loop.
			}
		} // i_site

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cur_pooled_plaintext.save(ofs_enc_data_matrix);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		vector<double> pooled_res_array;
		pooled_res_array.resize(encoder.slot_count());
		encoder.decode(cur_pooled_plaintext, pooled_res_array);

		// For the current ciphertext, write the output.
		pooled_data_values->insert(pooled_data_values->end(), pooled_res_array.begin(), pooled_res_array.end());
	} // i_pt loop.

	ofs_enc_data_matrix.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Pooled %d plaintext results, saving them into %dx%d matrix in %s\n", (int)(pooled_data_values->size()), (int)overall_n_row, (int)overall_n_col, op_fp);
	}

	if ((int)(overall_n_row * overall_n_col) > (int)(pooled_data_values->size()))
	{
		fprintf(stderr, "Sanity check error: Matrix has a size larger than the pooled points; %dx%d=%d vs %d data points.\n",
			(int)overall_n_row, (int)overall_n_col,
			(int)(overall_n_row * overall_n_col), (int)(pooled_data_values->size()));
		exit(0);
	}

	// Start reading each plaintext from the sites and pool them, decode and write to final output file.
	//FILE* f_op = open_f(op_fp, "w");

	// Write the header.
	double** dec_matrix = allocate_matrix(overall_n_row, overall_n_col);
	//fprintf(f_op, "#COL_ID");
	//for (int i_col = 0; i_col < overall_n_col; i_col++)
	//{
	//	fprintf(f_op, "\tCOL_%d", i_col);
	//} // i_vector loop.
	//fprintf(f_op, "\n");

	int data_pt_i = 0;
	for (int i_row = 0; i_row < (int)overall_n_row; i_row++)
	{
		//fprintf(f_op, "V%d", i_row);

		for (int i_col = 0; i_col < (int)overall_n_col; i_col++)
		{
			//fprintf(f_op, "\t%.10f", pooled_data_values->at(data_pt_i));

			dec_matrix[i_row][i_col] = pooled_data_values->at(data_pt_i);

			data_pt_i++;
		} // i_col loop.

		//fprintf(f_op, "\n");
	} // i_row loop.

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Saving %dx%d matrix to %s\n", (int)overall_n_row, (int)overall_n_col, op_fp);
	}

	save_matrix_binary(dec_matrix, overall_n_row, overall_n_col, op_fp);

	//close_f(f_op, op_fp);
} // collaborative_pool_partial_decrypted_continuous_enc_plaintext_results function.

void add_remove_encrypted_mask_2_continuous_encrypted_data_matrix(char* encrypted_data_matrix_path, char* encrypted_mask_matrix_fp,
	bool add_mask_flag,
	char* text_params_path,
	char* masked_encrypted_data_matrix_fp)
{
	if (add_mask_flag)
	{
		fprintf(stderr, "MASKING the plaintext mask matrix:\n\
Encrypted matrix file : %s\n\
Text params file : %s\n\
Output file : %s\n", encrypted_data_matrix_path, text_params_path, masked_encrypted_data_matrix_fp);
	}
	else
	{
		fprintf(stderr, "UNMASKING the plaintext mask matrix:\n\
Encrypted matrix file : %s\n\
Text params file : %s\n\
Output file : %s\n", encrypted_data_matrix_path, text_params_path, masked_encrypted_data_matrix_fp);

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

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Max coefficient modulus bit count: %d\n", CoeffModulus::MaxBitCount(poly_modulus_degree));
	}

	SEALContext context(parms);

	// Instantiate evaluator.
	Evaluator evaluator(context);

	// We need the encoder to get the number of samples per ct.
	CKKSEncoder encoder(context);
	int n_values_per_ct = encoder.slot_count();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading encrypted data matrix from %s\n", encrypted_data_matrix_path);
	}

	ifstream ifs_enc_data_matrix(encrypted_data_matrix_path, ios::binary);

	// Set the sample size information.
	int n_data_row;
	int n_data_col;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_data_col), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_data_row), sizeof(int));

	double n_entries_in_matrix = (double)(n_data_col * n_data_row);
	int n_cts_in_matrix = (int)(ceil(n_entries_in_matrix / n_values_per_ct));

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded %dx%d data matrix.\n", n_data_row, n_data_col);
	}

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading ciphertexts.\n");
	}

	vector<Ciphertext>* data_cts = new vector<Ciphertext>();
	for (int i_ct = 0; i_ct < n_cts_in_matrix; i_ct++)
	{
		Ciphertext cur_ct;
		cur_ct.load(context, ifs_enc_data_matrix);

		// Add this ct to the current feat's list.
		data_cts->push_back(cur_ct);
	} // i_ct loop.

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// LOAD THE ENCRYPTED MASK CTs:
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading encrypted mask matrix from %s\n", encrypted_mask_matrix_fp);
	}

	ifstream ifs_enc_mask_matrix(encrypted_mask_matrix_fp, ios::binary);

	// Set the sample size information.
	int n_mask_row;
	int n_mask_col;

	ifs_enc_mask_matrix.read(reinterpret_cast<char *>(&n_mask_col), sizeof(int));
	ifs_enc_mask_matrix.read(reinterpret_cast<char *>(&n_mask_row), sizeof(int));

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loaded %dx%d mask matrix\n", n_mask_row, n_mask_col);
	}

	if (n_data_row != n_mask_row ||
		n_data_col != n_mask_col)
	{
		fprintf(stderr, "The data matrix and mask matrix are not conformant: %dx%d ;; %dx%d\n",
			n_data_row, n_data_col,
			n_mask_row, n_mask_col);

		exit(0);
	}

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading encrypted mask data matrix:\n");
	}

	vector<Ciphertext>* mask_cts = new vector<Ciphertext>();
	for (int i_ct = 0; i_ct < n_cts_in_matrix; i_ct++)
	{
		Ciphertext cur_ct;
		cur_ct.load(context, ifs_enc_mask_matrix);

		// Add this ct to the current feat's list.
		mask_cts->push_back(cur_ct);
	} // i_ct loop.

	/////////////////////////////////////////////////////////////////////////////
	// Add the masks and save them.
	ofstream ofs_masked_enc_data_matrix(masked_encrypted_data_matrix_fp, ios::binary);

	// This describes the matrix that is being saved to output.
	int n_res_row = n_data_row; // This parameter has no bearing yet. It should be copied to partially decrypted plaintext file and final pooling file to write the results in the final final text output file. This sets how many values are meaningful (or necessary) in the resulting ciphertext.
	int n_res_col = n_data_col;
	ofs_masked_enc_data_matrix.write(reinterpret_cast<const char *>(&n_res_col), sizeof(int));
	ofs_masked_enc_data_matrix.write(reinterpret_cast<const char *>(&(n_res_row)), sizeof(int));

	for (int i_ct = 0; i_ct < n_cts_in_matrix; i_ct++)
	{
		Ciphertext cur_data_ct = data_cts->at(i_ct);
		Ciphertext cur_mask_ct = mask_cts->at(i_ct);

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Fresh mask ct modulus size: %d, chain index: %d\n", (int)(cur_mask_ct.coeff_modulus_size()),
				(int)(context.get_context_data(cur_mask_ct.parms_id())->chain_index()));
		}

		// Match the modulus of mask and data ct.
		evaluator.mod_switch_to_inplace(cur_mask_ct, cur_data_ct.parms_id());
		cur_mask_ct.scale() = cur_data_ct.scale();

		if (__DUMP_MATRIX_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Current data ct modulus size: %d, chain index: %d\n", (int)(cur_data_ct.coeff_modulus_size()),
				(int)(context.get_context_data(cur_data_ct.parms_id())->chain_index()));
		}

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

	ofs_masked_enc_data_matrix.close();
} // add_remove_encrypted_mask_2_data_matrix function.

// Get the data sizes and generate masking matrix: For each subject's value, generate a masking value and add it.
void generate_mask_matrix_per_data_matrix(char* encrypted_data_matrix_path,
	double mask_variance,
	char* plain_mask_matrix_op_fp)
{
	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Generating the plaintext mask matrix:\n\
Encrypted matrix file : %s\n\
Output file : %s\n", encrypted_data_matrix_path, plain_mask_matrix_op_fp);
	}

	long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Mask matrix random seed: %lu\n", rand_seed);
	}

	t_rng* rng = new t_rng(rand_seed);

	// Start loading the data: Load the (1) sample size, (2) number of features, and (3) number of ciphertexts per feature.
	// Open this file.
	ifstream ifs_enc_data_matrix(encrypted_data_matrix_path, ios::binary);

	// Read the dimensions.
	int n_row;
	int n_col;

	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_col), sizeof(int));
	ifs_enc_data_matrix.read(reinterpret_cast<char *>(&n_row), sizeof(int));

	ifs_enc_data_matrix.close();

	if (__DUMP_MATRIX_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Loading %dx%d matrix..\n", n_row, n_col);
	}

	double** mask_matrix = allocate_matrix(n_row, n_col);

	for (int i_row = 0; i_row < n_row; i_row++)
	{
		for (int i_col = 0; i_col < n_col; i_col++)
		{
			double cur_mask_val = rng->random_double_ran3() * mask_variance;
			mask_matrix[i_row][i_col] = cur_mask_val;
		} // i_s loop.
	} // i_feat loop.

	// Save mask matrix.
	save_matrix_binary(mask_matrix, n_row, n_col, plain_mask_matrix_op_fp);
} // generate_mask_matrix_per_data_matrix option.

