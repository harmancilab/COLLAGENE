#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_file_utils.h"
#include "cllgn_histogram.h"
#include "cllgn_ansi_string.h"
#include <math.h>
#include <string.h>
#include "cllgn_features_weight_utils.h"
#include "cllgn_x_sigmoid.h"
#include "cllgn_x_chisqr.h"
#include "cllgn_LR_model_stats.h"
#include "cllgn_LR_utils.h"
#include "cllgn_fedLR_utils.h"
#include "cllgn_fedLR_secure_convertible_protocol_utils.h"
#include "cllgn_matrix_linalg_utils.h"

bool __DUMP_CRYPTABLE_FEDLR_MESSAGES__ = false;

/*
HE functions: 
-------------
- Transpose.
- Extract diagonal vector from a matrix
- Get matrix from vector.
- Matrix multiplication.
*/

using namespace std;

void plain_add_matrices_per_list(char* matrix_path_list_fp, char* op_fp)
{
	vector<char*>* matrix_path_list = buffer_file(matrix_path_list_fp);
	fprintf(stderr, "Adding %d matrices and writing results to %s\n", (int)(matrix_path_list->size()), op_fp);

	double** res_matrix = NULL;
	int nrows = 0;
	int ncols = 0;

	for (int i_mat = 0; i_mat < (int)(matrix_path_list->size()); i_mat++)
	{
		int cur_n_rows, cur_n_cols;
		double** cur_mat = load_matrix_binary(matrix_path_list->at(i_mat), cur_n_rows, cur_n_cols);
		if (i_mat == 0)
		{
			res_matrix = cur_mat;
			nrows = cur_n_rows;
			ncols = cur_n_cols;
		}
		else
		{
			matrix_add(res_matrix, nrows, ncols, cur_mat, cur_n_rows, cur_n_cols, res_matrix);
		}
	} // i_mat loop.

	save_matrix_binary(res_matrix, nrows, ncols, op_fp);
}

// This is a placeholder for real partial decryption function; it implements a simple partializing of the matrix given decrypting client index and total clients.
double** partial_decrypt_matrix(double** matrix, int nrow, int ncol, int dec_i_client, int n_clients)
{
	double num = (double)(dec_i_client + 1);
	double den = (double)(n_clients * (n_clients + 1) / 2);

	double part_dec_factor = num / den;
	double** part_dec_epoc_weights = matrix_multiply_by_scalar(matrix, nrow, ncol, part_dec_factor, NULL);
	return(part_dec_epoc_weights);
}

// Ensures that the mult noise is non-singular.
double** generate_multiplicative_full_noise_matrix(int nrow, int ncol)
{
	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	double** mult_noise_matrix = allocate_matrix(nrow, ncol);

	// This is the full random matrix.
	for(int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			mult_noise_matrix[i_row][i_col] = rng->random_gaussian_double_ran3();
		} // i_col loop.
	} // i_row loop.

	return(mult_noise_matrix);
}

double** generate_constant_diagonal_matrix(int nrow, int ncol, double diagonal_val)
{
	double** mult_noise_matrix = allocate_matrix(nrow, ncol);
	int min_dim = (ncol < nrow) ? (ncol) : (nrow);
	for (int i_diag = 0; i_diag < min_dim; i_diag++)
	{
		mult_noise_matrix[i_diag][i_diag] = diagonal_val;
	} // i_diag loop.

	return(mult_noise_matrix);
}

double** generate_multiplicative_diagonal_noise_matrix(int nrow, int ncol)
{
	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	double** mult_noise_matrix = allocate_matrix(nrow, ncol);
	for (int i_diag = 0; i_diag < nrow; i_diag++)
	{
		mult_noise_matrix[i_diag][i_diag] = rng->random_gaussian_double_ran3();
	} // i_diag loop.

	return(mult_noise_matrix);
}

void get_fulldec_pval_stats_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/full_dec_pval_matrix_iter_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

void get_partdec_pval_stats_matrix_fp(char* shared_working_dir, int iter_i, int var_block_start_i, int var_block_end_i, int client_i, int dec_client_i, char* buffer)
{
	sprintf(buffer, "%s/part_dec_pval_matrix_iter_%d_block_%d_%d_client_%d_by_client_%d.bin", shared_working_dir, iter_i, var_block_start_i, var_block_end_i, client_i, dec_client_i);
}

void get_pval_stat_matrix_op_fp(char* shared_working_dir, int iter_i, int var_block_start_i, int var_block_end_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/PVAL_STATS_ITER_%d_BLOCK_%d_%d_client_%d.bin", shared_working_dir, iter_i, var_block_start_i, var_block_end_i, client_i);
}

void get_inv_XtWX_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/inv_XtWX_iter_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

void get_site_specific_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/site_specific_all_site_noise_XtWX_iter_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

void get_sitewise_pooled_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/sitewise_pooled_all_site_noise_XtWX_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

void get_partdec_sitewise_pooled_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, int dec_client_i, char* buffer)
{
	sprintf(buffer, "%s/part_dec_sitewise_pooled_all_site_noise_XtWX_iter_%d_client_%d_by_client_%d.bin", shared_working_dir, iter_i, client_i, dec_client_i);
}

void get_partdec_beta_fp(char* shared_working_dir, int iter_i, int client_i, int dec_client_i, char* buffer)
{
	sprintf(buffer, "%s/part_dec_beta_%d_client_%d_by_client_%d.bin", shared_working_dir, iter_i, client_i, dec_client_i);
}

void get_fulldec_beta_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/full_dec_beta_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

void get_fulldec_sitewise_pooled_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/full_dec_sitewise_pooled_all_site_noise_XtWX_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

void get_mult_noise_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/mult_noise_matrix_iter_%d_client_%d.bin", shared_working_dir, iter_i, client_i);
}

/* 
This function loads and processes the current null model stats using IRLS fitting model.
1- Load the current beta matrix
2- Update mu, nu, z, W, ...
3- Save the intermediate matrix.
*/
void cryptable_client_calculate_save_XtWX_XtWz(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "@ iteration %d, Client %d/%d -- Calculating XtWX and XtWz and saving to shared working space %s\n", i_iter, client_i, n_clients, shared_working_dir);

	// We need the two sigmoid functions here.
	double(*GET_SIGMOID_PER_FEAT_COMB)(double);

	if (sigmoid_approx_type == SIGMOID_APPROX_NATIVE)
	{
		fprintf(stderr, "Using native sigmoid..\n");
		GET_SIGMOID_PER_FEAT_COMB = get_sigmoid_val_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_KIM_ETAL)
	{
		fprintf(stderr, "Using Kim. etal sigmoid..\n");
		GET_SIGMOID_PER_FEAT_COMB = get_Kim_etal_poly_approx_sigmoid_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_TENSEAL)
	{
		fprintf(stderr, "Using TenSEAL sigmoid..\n");
		GET_SIGMOID_PER_FEAT_COMB = get_TenSeal_poly_approx_sigmoid_per_feat_comb;
	}
	else
	{
		fprintf(stderr, "Unknown sigmoid approximation: %d\n", sigmoid_approx_type);
		exit(0);
	}

	// The training data/model.
	int n_feats = 0;
	int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loading/generation.
	// Load the features if they exist.
	fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
	int n_covars_per_feats_file = 0;
	vector<double*>* per_ind_feat_matrix = new vector<double*>();
	vector<double>* per_ind_pheno_vector = new vector<double>();
	load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

	n_covars = n_covars_per_feats_file;
	n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

	// Copy the phenotype vector to an array.
	double* per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
	for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
	{
		per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
	}

	if (__DUMP_CRYPTABLE_FEDLR_MESSAGES__)
	{
		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno_vector->at(i_s));
		} // i_s loop.
	}

	// Sample size must be fixed after loading.
	sample_size = (int)(per_ind_feat_matrix->size());

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Calculating XtWX and XtWz with local sample size of %d:\n\
n_covars: %d\n\
n_epoch: %d\n\
shared working directory: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, n_epoch, shared_working_dir, LL_EPSILON);

	if (__DUMP_CRYPTABLE_FEDLR_MESSAGES__)
	{
		// Write the feats and pheno.
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			fprintf(stderr, "Individual %d features:\n", i_s);
			for (int feat_i = 0; feat_i < n_feats; feat_i++)
			{
				fprintf(stderr, "Feat %d: %.4f;", feat_i, per_ind_feat_matrix->at(i_s)[feat_i]);
			} // feat_i loop.

			fprintf(stderr, "Pheno: %.2f\n", per_ind_pheno_vector->at(i_s));
		} // i_s loop.
	}

	double** cur_epoch_weights = NULL;

	// This is the main point of loading of the beta values.
	// Load the current epoch weights, if they exist.
	int prev_i_iter = (i_iter > 0) ? (i_iter - 1) : (n_epoch);
	char beta_matrix_fp[1000];

	// Check for fully decrypted beta.
	get_fulldec_beta_fp(shared_working_dir, prev_i_iter, client_i, beta_matrix_fp);
	if (!check_file(beta_matrix_fp))
	{
		// This automatically sets all weights to 0, i.e., mu=0.5.
		cur_epoch_weights = allocate_matrix(n_feats, 1);
	}
	else
	{
		fprintf(stderr, "Loading current weights from %s\n", beta_matrix_fp);

		int n_loaded_rows;
		int n_loaded_cols;
		cur_epoch_weights = load_matrix_binary(beta_matrix_fp, n_loaded_rows, n_loaded_cols);

		if (n_loaded_rows != n_feats || n_loaded_cols != 1)
		{
			fprintf(stderr, "Loaded weight matrix is not conformant with the sizes; %d, %d ; %d, %d.\n", n_loaded_rows, n_loaded_cols, n_feats, 1);
		}
	}
	
	// We don't need all of these for this step.
	double** nu = allocate_matrix(sample_size, 1);
	double** mu = allocate_matrix(sample_size, 1);
	double** one_min_mu = allocate_matrix(sample_size, 1);
	double** mu_one_min_mu = allocate_matrix(sample_size, 1);
	double** z = allocate_matrix(sample_size, 1);
	double** Yt = allocate_copy_matrix(&(per_ind_pheno), 1, sample_size);
	double** Y = transpose_matrix(Yt, 1, sample_size, NULL);

	double** W_diag = allocate_matrix(sample_size, 1);
	double** W = allocate_matrix(sample_size, sample_size);

	double** Y_min_mu = allocate_matrix(sample_size, 1);
	double** Y_min_mu_by_mu_one_min_mu = allocate_matrix(sample_size, 1);

	double** one_vec_sample = allocate_matrix(sample_size, 1, 1.0);

	// Copy the features matrix and allocat the related matrices.
	double** X = get_matrix_per_row_vector_list(per_ind_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);
	double** XtW = allocate_matrix(n_feats, sample_size);
	double** XtWX = allocate_matrix(n_feats, n_feats);
	double** XtWz = allocate_matrix(n_feats, 1);

	double max_abs_nu = 0;

	// Calculate nu and mu, and z=nu + (y[i]-mu[i])/(mu[i]*(1-mu[i])
	matrix_multiply(X, sample_size, n_feats, cur_epoch_weights, n_feats, 1, nu);

	// Calculate mu from nu
	process_matrix_elementwise_by_callback(nu, sample_size, 1, GET_SIGMOID_PER_FEAT_COMB, mu);

	// Calculate z. This is required for updating beta.
	matrix_subtract(one_vec_sample, sample_size, 1, mu, sample_size, 1, one_min_mu);

	// Calculate mu(1-mu) elementwise.
	matrix_multiply_elementwise(one_min_mu, sample_size, 1, mu, sample_size, 1, mu_one_min_mu);

	// Calculate Y-mu
	matrix_subtract(Y, sample_size, 1, mu, sample_size, 1, Y_min_mu);

	// Divide Y-mu by mu(1-mu).
	matrix_divide_elementwise(Y_min_mu, sample_size, 1, mu_one_min_mu, sample_size, 1, Y_min_mu_by_mu_one_min_mu);

	// Add nu to (Y-mu)/(mu(1-mu)).
	matrix_add(nu, sample_size, 1, Y_min_mu_by_mu_one_min_mu, sample_size, 1, z);

	// Generate W_diag; this is already calcualted above.
	//matrix_multiply_elementwise(mu, sample_size, 1, one_min_mu, sample_size, 1, W_diag);
	copy_matrix(mu_one_min_mu, sample_size, 1, W_diag, sample_size, 1);

	// Convert W_diag to W.
	get_diag_matrix(W_diag, sample_size, 1, W);

	// Calculate X'WX and X'Wz
	matrix_multiply(Xt, n_feats, sample_size, W, sample_size, sample_size, XtW);

	matrix_multiply(XtW, n_feats, sample_size, X, sample_size, n_feats, XtWX);

	// We compute this XtWX and save it. This needs to be pooled, inverted and used.
	char XtWX_matrix_fp[1000];
	get_XtWX_matrix_fp(private_working_dir, i_iter, client_i, XtWX_matrix_fp);
	save_matrix_binary(XtWX, n_feats, n_feats, XtWX_matrix_fp);

	// Calculate XtWz and store it in the shared directory.
	matrix_multiply(XtW, n_feats, sample_size, z, sample_size, 1, XtWz);

	// We compute this XtWX and save it. This needs to be pooled, inverted and used.
	char XtWz_matrix_fp[1000];
	get_XtWz_matrix_fp(shared_working_dir, i_iter, client_i, XtWz_matrix_fp);
	save_matrix_binary(XtWz, n_feats, 1, XtWz_matrix_fp);

	// Compute the LL and save it for the current client. This calculate the LL for the previous beta values, not with the new values.
	double cur_epoch_updated_LL = get_log_likelihood_per_mu_matrix(mu, per_ind_pheno, sample_size);

	double** LL_matrix = allocate_matrix(1, 1);
	LL_matrix[0][0] = cur_epoch_updated_LL;

	char LL_matrix_fp[1000];
	get_LL_matrix_fp(shared_working_dir, i_iter, client_i, LL_matrix_fp);
	save_matrix_binary(LL_matrix, 1, 1, LL_matrix_fp);

	// Get the largest nu, this is necessary for sigmoid approximations.
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_subject_nu = nu[i_s][0];
		if (max_abs_nu < fabs(cur_subject_nu))
		{
			max_abs_nu = fabs(cur_subject_nu);
		}
	} // i_s loop.

	// Generate the multiplicative noise that will be used to calculate inv(XtWX).
	char mult_noise_fp[1000];
	get_mult_noise_fp(shared_working_dir, i_iter, client_i, mult_noise_fp);
	double** mult_noise = generate_multiplicative_diagonal_noise_matrix(n_feats, n_feats);
	save_matrix_binary(mult_noise, n_feats, n_feats, mult_noise_fp);

	if (__DUMP_CRYPTABLE_FEDLR_MESSAGES__)
	{
		char op_fp[1000];
		sprintf(op_fp, "mult_noise_iter_%d_client_%d.txt", i_iter, client_i);
		save_matrix_plain(mult_noise, n_feats, n_feats, op_fp);
	}
} // client_calculate_save_XtWX_XtWz option.

void cryptable_client_add_mult_noise_2_XtWX(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	//int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		//sample_size = per_ind_feat_matrix->size();
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	char XtWX_matrix_fp[1000];
	get_XtWX_matrix_fp(private_working_dir, i_iter, client_i, XtWX_matrix_fp);

	int loaded_nrows, loaded_ncols;
	double** XtWX = load_matrix_binary(XtWX_matrix_fp, loaded_nrows, loaded_ncols);
	
	// Load the noise matrices from each site from shared directory and HE-multiply.
	double** noisy_XtWX = allocate_matrix(n_feats, n_feats);
	copy_matrix(XtWX, n_feats, n_feats, noisy_XtWX, n_feats, n_feats);

	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		char mult_noise_fp[1000];
		get_mult_noise_fp(shared_working_dir, i_iter, cur_client_i, mult_noise_fp);
		int loaded_nrows, loaded_ncols;
		double** cur_client_noise_matrix = load_matrix_binary(mult_noise_fp, loaded_nrows, loaded_ncols, NULL);
		if (loaded_nrows != n_feats ||
			loaded_ncols != n_feats)
		{
			fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", mult_noise_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
			exit(0);
		}

		double** cur_mat = matrix_multiply(noisy_XtWX, n_feats, n_feats, cur_client_noise_matrix, n_feats, n_feats, NULL);
		noisy_XtWX = cur_mat;
	} // cur_client_i loop.

	// Save the noisy encrypted matrix to shared directory.
	// This matrix is site-specific XtWX with noise levels from all sites.
	char site_spec_pooled_noise_XtWX_fp[1000];
	get_site_specific_all_site_noise_XtWX_fp(shared_working_dir, i_iter, client_i, site_spec_pooled_noise_XtWX_fp);
	save_matrix_binary(noisy_XtWX, n_feats, n_feats, site_spec_pooled_noise_XtWX_fp);
}

void cryptable_client_pool_site_specific_all_site_noise_XtWX(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	//int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		//sample_size = per_ind_feat_matrix->size();
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double** pooled_noisy_XtWX = allocate_matrix(n_feats, n_feats);

	// Pool the all-site-noisy-XtWX from all sites.
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		char site_spec_pooled_noise_XtWX_fp[1000];
		get_site_specific_all_site_noise_XtWX_fp(shared_working_dir, i_iter, cur_client_i, site_spec_pooled_noise_XtWX_fp);
		int loaded_nrows, loaded_ncols;
		double** site_specific_all_site_noise_XtWX = load_matrix_binary(site_spec_pooled_noise_XtWX_fp, loaded_nrows, loaded_ncols, NULL);

		if (loaded_nrows != n_feats ||
			loaded_ncols != n_feats)
		{
			fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", site_spec_pooled_noise_XtWX_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
			exit(0);
		}

		matrix_add(pooled_noisy_XtWX, n_feats, n_feats, site_specific_all_site_noise_XtWX, n_feats, n_feats, pooled_noisy_XtWX);
	} // cur_client_i loop.

	// This is the matrix that contains noise added and pooled among sites.
	char sitewise_pooled_all_site_noise_XtWX_fp[1000];
	get_sitewise_pooled_all_site_noise_XtWX_fp(shared_working_dir, i_iter, client_i, sitewise_pooled_all_site_noise_XtWX_fp);
	save_matrix_binary(pooled_noisy_XtWX, n_feats, n_feats, sitewise_pooled_all_site_noise_XtWX_fp);
}

void cryptable_collaborative_decrypt_pooled_noisy_XtWX(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	//int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		//sample_size = per_ind_feat_matrix->size();
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Pool the all-site-noisy-XtWX from all sites.
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{	
		// Load the current client's matrix.
		int loaded_nrows, loaded_ncols;
		char sitewise_pooled_all_site_noise_XtWX_fp[1000];
		get_sitewise_pooled_all_site_noise_XtWX_fp(shared_working_dir, i_iter, cur_client_i, sitewise_pooled_all_site_noise_XtWX_fp);
		double** sitewise_pooled_all_site_noise_XtWX = load_matrix_binary(sitewise_pooled_all_site_noise_XtWX_fp, loaded_nrows, loaded_ncols, NULL);

		if (loaded_nrows != n_feats ||
			loaded_ncols != n_feats)
		{
			fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", sitewise_pooled_all_site_noise_XtWX_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
			exit(0);
		}

		// Just divide by the # of clients.
		//matrix_multiply_by_scalar(sitewise_pooled_all_site_noise_XtWX, n_feats, n_feats, (double)1.0/(double)n_clients, sitewise_pooled_all_site_noise_XtWX);
		double** part_dec_sitewise_pooled_all_site_noise_XtWX = partial_decrypt_matrix(sitewise_pooled_all_site_noise_XtWX, n_feats, n_feats, client_i, n_clients);

		// Save the partdec by client_i.
		char partdec_sitewise_pooled_all_site_noise_XtWX_fp[1000];
		get_partdec_sitewise_pooled_all_site_noise_XtWX_fp(shared_working_dir, i_iter, cur_client_i, client_i, partdec_sitewise_pooled_all_site_noise_XtWX_fp);
		save_matrix_binary(part_dec_sitewise_pooled_all_site_noise_XtWX, n_feats, n_feats, partdec_sitewise_pooled_all_site_noise_XtWX_fp);
	} // cur_client_i loop.
}

// This is the last step: Pool the partial decryptions then remove noise, which can only be done in encrypted domain.
void cryptable_pool_partially_decrypted_pooled_noisy_XtWx_remove_noise(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	//int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		//sample_size = per_ind_feat_matrix->size();
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double** dec_pooled_noisy_XtWX = allocate_matrix(n_feats, n_feats);

	// Pool all the partially decrypted results from decrypting sites.
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		int loaded_nrows, loaded_ncols;
		char partdec_sitewise_pooled_all_site_noise_XtWX_fp[1000];
		get_partdec_sitewise_pooled_all_site_noise_XtWX_fp(shared_working_dir, i_iter, client_i, cur_client_i, partdec_sitewise_pooled_all_site_noise_XtWX_fp);
		double** part_dec_sitewise_pooled_all_site_noise_XtWX = load_matrix_binary(partdec_sitewise_pooled_all_site_noise_XtWX_fp, loaded_nrows, loaded_ncols, NULL);

		if (loaded_nrows != n_feats ||
			loaded_ncols != n_feats)
		{
			fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", partdec_sitewise_pooled_all_site_noise_XtWX_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
			exit(0);
		}

		// This is simple addition.
		matrix_add(dec_pooled_noisy_XtWX, n_feats, n_feats, part_dec_sitewise_pooled_all_site_noise_XtWX, n_feats, n_feats, dec_pooled_noisy_XtWX);
	} // cur_client_i loop.

	// Remove noise.
	fprintf(stderr, "Inverting all-site noisy XtWX matrix.\n");
	double** inv_dec_pooled_noisy_XtWX = invert_matrix_GJ(dec_pooled_noisy_XtWX, n_feats, n_feats, NULL);

	// We can remove the noise at this point.
	for (int cur_client_i = n_clients-1; cur_client_i >= 0 ; cur_client_i--)
	{
		char mult_noise_fp[1000];
		get_mult_noise_fp(shared_working_dir, i_iter, cur_client_i, mult_noise_fp);
		int loaded_nrows, loaded_ncols;
		double** cur_client_noise_matrix = load_matrix_binary(mult_noise_fp, loaded_nrows, loaded_ncols, NULL);
		if (loaded_nrows != n_feats ||
			loaded_ncols != n_feats)
		{
			fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", mult_noise_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
			exit(0);
		}

		// In practice, this is an encrypted product and we can only get an encrypted matrix at this point, so noise removal can be done here.
		matrix_multiply(cur_client_noise_matrix, n_feats, n_feats, inv_dec_pooled_noisy_XtWX, n_feats, n_feats, inv_dec_pooled_noisy_XtWX);
	} // cur_client_i loop.

	// Save the final XtWX.
	char inv_XtWX_fp[1000];
	get_inv_XtWX_matrix_fp(private_working_dir, i_iter, client_i, inv_XtWX_fp);
	save_matrix_binary(inv_dec_pooled_noisy_XtWX, n_feats, n_feats, inv_XtWX_fp);
}

// This function takes encrypted XtWz and multiplies it with XtWX
void cryptable_client_pool_XtWX_XtWz_update_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "@ iteration %d, Client %d/%d -- Pooling XtWX and XtWz and updating beta and saving to shared working space %s\n", i_iter, client_i, n_clients, shared_working_dir);

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loeading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		sample_size = (int)(per_ind_feat_matrix->size());
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Training plaintext LR with IRLS with sample size of %d:\n\
n_covars: %d\n\
n_epoch: %d\n\
shared working directory: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, n_epoch, shared_working_dir, LL_EPSILON);

	if (__DUMP_CRYPTABLE_FEDLR_MESSAGES__)
	{
		// Write the feats and pheno.
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			fprintf(stderr, "Individual %d features:\n", i_s);
			for (int feat_i = 0; feat_i < n_feats; feat_i++)
			{
				fprintf(stderr, "Feat %d: %.4f;", feat_i, per_ind_feat_matrix->at(i_s)[feat_i]);
			} // feat_i loop.

			fprintf(stderr, "Pheno: %.2f\n", per_ind_pheno[i_s]);
		} // i_s loop.
	}

	double** XtWX_inv = allocate_matrix(n_feats, n_feats);
	double** XtWz = allocate_matrix(n_feats, 1);
	double** cur_epoch_weights = allocate_matrix(n_feats, 1);

	// Load the inverted XtWX matrix for this site -- this is coming from the above functions that perform collaborative inversion.

	// Load the XtWX matrices from all clients and pool them.
	fprintf(stderr, "Pooling XtWX and XtWz matrices..\n");
	char inv_XtWX_fp[1000];
	get_inv_XtWX_matrix_fp(private_working_dir, i_iter, client_i, inv_XtWX_fp);
	int loaded_nrows, loaded_ncols;
	load_matrix_binary(inv_XtWX_fp, loaded_nrows, loaded_ncols, XtWX_inv);
	if (loaded_nrows != n_feats ||
		loaded_ncols != n_feats)
	{
		fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", inv_XtWX_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
		exit(0);
	}

	// Pool XtWz matrix from all sites; this is encrypted.
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		int n_loaded_row, n_loaded_col;

		// Load and pool XtWz matrix.
		char cur_client_XtWz_matrix_fp[1000];
		get_XtWz_matrix_fp(shared_working_dir, i_iter, cur_client_i, cur_client_XtWz_matrix_fp);
		fprintf(stderr, "Loading XtWz matrix from %s\n", cur_client_XtWz_matrix_fp);

		double** cur_XtWz = load_matrix_binary(cur_client_XtWz_matrix_fp, n_loaded_row, n_loaded_col);
		if (n_loaded_col != 1 || n_loaded_row != n_feats)
		{
			fprintf(stderr, "Loaded XtWz is non-conformant for client %d @ iteration %d: %d, %d ; %d, %d\n", client_i, i_iter, n_feats, 1, n_loaded_row, n_loaded_col);
			exit(0);
		}

		// Pool this matrix, this can be done in place.
		fprintf(stderr, "Pooling XtWz.\n");
		matrix_add(XtWz, n_feats, 1, cur_XtWz, n_feats, 1, XtWz);
	} // client_i loop.

	// This is the update step: We have the pooled XtWz matrix from all sites. Just do multiplication, this is done in secure domain.
	matrix_multiply(XtWX_inv, n_feats, n_feats, XtWz, n_feats, 1, cur_epoch_weights);

	// Save the contribution of this site to the epoch weights.
	char beta_matrix_fp[1000];
	get_beta_matrix_fp(shared_working_dir, i_iter, client_i, beta_matrix_fp);
	save_matrix_binary(cur_epoch_weights, n_feats, 1, beta_matrix_fp);

	char params_txt_fp[1000];
	sprintf(params_txt_fp, "current_params_iter_%d_client_%d.txt", i_iter, client_i);
	save_matrix_plain(cur_epoch_weights, n_feats, 1, params_txt_fp);
} // client_pool_XtWX_XtWz_update_beta

void cryptable_client_collaborative_decrypt_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "@ iteration %d, Client %d/%d -- Collaboratively decryption beta and saving to shared working space %s\n", i_iter, client_i, n_clients, shared_working_dir);

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loeading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		sample_size = (int)(per_ind_feat_matrix->size());
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Training plaintext LR with IRLS with sample size of %d:\n\
n_covars: %d\n\
n_epoch: %d\n\
shared working directory: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, n_epoch, shared_working_dir, LL_EPSILON);

	// Decrypt other sites' current encrypted beta.
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		char beta_matrix_fp[1000];
		get_beta_matrix_fp(shared_working_dir, i_iter, cur_client_i, beta_matrix_fp);
		int loaded_nrow, loaded_ncol;
		double** cur_epoch_weights = load_matrix_binary(beta_matrix_fp, loaded_nrow, loaded_ncol, NULL);
		if (loaded_nrow != n_feats || loaded_ncol != 1)
		{
			fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", beta_matrix_fp, loaded_nrow, loaded_ncol, n_feats, 1);
			exit(0);
		}

		double** part_dec_epoc_weights = partial_decrypt_matrix(cur_epoch_weights, n_feats, 1, client_i, n_clients);

		char partdec_beta_fp[1000];
		get_partdec_beta_fp(shared_working_dir, i_iter, cur_client_i, client_i, partdec_beta_fp);
		save_matrix_binary(part_dec_epoc_weights, n_feats, 1, partdec_beta_fp);
	} // cur_client_i loop.
}

void cryptable_client_pool_partially_decrypted_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "@ iteration %d, Client %d/%d -- Pooling XtWX and XtWz and updating beta and saving to shared working space %s\n", i_iter, client_i, n_clients, shared_working_dir);

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	int n_feats = 0;
	int sample_size = 0;
	int n_covars = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start feature loeading/generation.
	// Load the features if they exist.
	if (check_file(subject_per_row_feats_fp))
	{
		fprintf(stderr, "Loading the training features data from %s.\n", subject_per_row_feats_fp);
		int n_covars_per_feats_file = 0;
		per_ind_feat_matrix = new vector<double*>();
		vector<double>* per_ind_pheno_vector = new vector<double>();
		load_LR_data_row_samples(subject_per_row_feats_fp, subject_per_row_pheno_fp, n_covars_per_feats_file, per_ind_feat_matrix, per_ind_pheno_vector);

		n_covars = n_covars_per_feats_file;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[(int)(per_ind_feat_matrix->size()) + 2];
		for (int i_s = 0; i_s < (int)per_ind_feat_matrix->size(); i_s++)
		{
			per_ind_pheno[i_s] = per_ind_pheno_vector->at(i_s);
		}

		fprintf(stderr, "First 10 subjects' %d features/phenotypes:\n", n_feats);
		for (int i_s = 0; i_s < 10; i_s++)
		{
			fprintf(stderr, "Subject_%d: ", i_s);
			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(stderr, "%d: %.4f;", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		sample_size = (int)(per_ind_feat_matrix->size());
	} // feats file check.
	else
	{
		fprintf(stderr, "Could not find the features file @ %s\n", subject_per_row_feats_fp);
		exit(0);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Training plaintext LR with IRLS with sample size of %d:\n\
n_covars: %d\n\
n_epoch: %d\n\
shared working directory: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, n_epoch, shared_working_dir, LL_EPSILON);

	// Decrypt other sites' current encrypted beta.
	double** pooled_beta_matrix = allocate_matrix(n_feats, 1);
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		char partdec_beta_fp[1000];
		get_partdec_beta_fp(shared_working_dir, i_iter, client_i, cur_client_i, partdec_beta_fp);
		int loaded_nrows, loaded_ncols;
		double** cur_partdec_beta_matrix = load_matrix_binary(partdec_beta_fp, loaded_nrows, loaded_ncols, NULL);
		if (loaded_nrows != n_feats || loaded_ncols != 1)
		{
			fprintf(stderr, "Sanity check failed loading %s: %d, %d; %d, %d\n", partdec_beta_fp, loaded_nrows, loaded_ncols, n_feats, 1);
			exit(0);
		}

		// Pooling is simple addition for now.
		matrix_add(pooled_beta_matrix, n_feats, 1, cur_partdec_beta_matrix, n_feats, 1, pooled_beta_matrix);
	} // cur_client_i loop.

	// Save the fully decrypted beta to file; beta should be read from this file next time.
	char full_dec_beta_fp[1000];
	get_fulldec_beta_fp(shared_working_dir, i_iter, client_i, full_dec_beta_fp);
	save_matrix_binary(pooled_beta_matrix, n_feats, 1, full_dec_beta_fp);
}

// This function stays as is, it just needs to be collaboratively decrypted in case there are a lot of iterations.
void cryptable_client_check_convergence_per_updated_beta(int i_iter, int client_i, int n_clients, char* shared_working_dir, double LL_EPSILON)
{
	if (i_iter > 0)
	{
		fprintf(stderr, "Pooling previous iteration's LL stats from all sites.\n");
		int prev_i_iter = i_iter - 1;
		double prev_iter_pooled_LL = 0;
		for (int cur_i_client = 0; cur_i_client < n_clients; cur_i_client++)
		{
			char LL_matrix_fp[1000];
			get_LL_matrix_fp(shared_working_dir, prev_i_iter, client_i, LL_matrix_fp);

			int nrow, ncol;
			double** cur_client_LL = load_matrix_binary(LL_matrix_fp, nrow, ncol);

			prev_iter_pooled_LL += cur_client_LL[0][0];
		} // cur_i_client loop.

		fprintf(stderr, "Pooling current iteration's LL stats from all sites.\n");
		double cur_iter_pooled_LL = 0;
		for (int cur_i_client = 0; cur_i_client < n_clients; cur_i_client++)
		{
			char LL_matrix_fp[1000];
			get_LL_matrix_fp(shared_working_dir, i_iter, client_i, LL_matrix_fp);

			int nrow, ncol;
			double** cur_client_LL = load_matrix_binary(LL_matrix_fp, nrow, ncol);

			cur_iter_pooled_LL += cur_client_LL[0][0];
		} // cur_i_client loop.

		if (fabs(cur_iter_pooled_LL - prev_iter_pooled_LL) < LL_EPSILON)
		{
			fprintf(stderr, "iter: %d -- CONVERGENCE IS REACHED!!!\n", i_iter);
			fprintf(stderr, "iter: %d -- CONVERGENCE IS REACHED!!!\n", i_iter);
			fprintf(stderr, "iter: %d -- CONVERGENCE IS REACHED!!!\n", i_iter);
			fprintf(stderr, "iter: %d -- CONVERGENCE IS REACHED!!!\n", i_iter);
		}
	}
	else
	{
		fprintf(stderr, "Skipping convergence check.\n");
	}
} // check convergence.

// This function computes GtWG (vx1) and GtWX (vxp)
void cryptable_client_calculate_save_pvalue_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, 
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "Client %d: Calculating p-value stats..\n", client_i);

	fprintf(stderr, "Loading the features data from %s.\n", subject_per_row_feats_fp);
	int n_covars_per_feats_file = 0;
	vector<double*>* per_subject_feat_matrix = new vector<double*>();
	vector<double>* per_subject_obs_pheno = new vector<double>();
	load_LR_data_row_samples(subject_per_row_feats_fp, obs_pheno_fp, n_covars_per_feats_file, per_subject_feat_matrix, per_subject_obs_pheno);

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars_per_feats_file);
	fprintf(stderr, "%d features loaded from covariates file.\n", n_feats);

	int sample_size = (int)(per_subject_obs_pheno->size());

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	double** X = get_matrix_per_row_vector_list(per_subject_feat_matrix, n_feats, NULL);
	double** nu = allocate_matrix(sample_size, 1);
	double** Y0 = allocate_matrix(sample_size, 1);

	char beta_matrix_fp[1000];
	//get_beta_matrix_fp(shared_working_dir, i_iter, client_i, beta_matrix_fp);
	get_fulldec_beta_fp(shared_working_dir, i_iter, client_i, beta_matrix_fp);
	int loaded_nrow, loaded_ncol;
	double** cur_epoch_weights = load_matrix_binary(beta_matrix_fp, loaded_nrow, loaded_ncol);
	if (loaded_nrow != n_feats || loaded_ncol != 1)
	{
		fprintf(stderr, "Sanity check failed while loading %s: %d, %d ; %d, %d\n", beta_matrix_fp, loaded_nrow, loaded_ncol, n_feats, 1);
		exit(0);
	}

	fprintf(stderr, "Loaded beta values:\n");
	for (int i_row = 0; i_row < n_feats; i_row++)
	{
		fprintf(stderr, "%.10f\n", cur_epoch_weights[i_row][0]);
	} // i_row loop.

	// Calculate nu and null phenotype matrix.
	matrix_multiply(X, sample_size, n_feats, cur_epoch_weights, n_feats, 1, nu);

	// Calculate mu from nu
	process_matrix_elementwise_by_callback(nu, sample_size, 1, get_sigmoid_val_per_feat_comb, Y0);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)per_subject_obs_pheno->size(), (int)per_subject_feat_matrix->size());
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\t%.4f\t%.4f\n", i_s, Y0[i_s][0] - per_subject_obs_pheno->at(i_s),
				Y0[i_s][0], 
				per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	// Compute the projection matrix: 
	fprintf(stderr, "Setting Sigma inv..\n");
	double* W_vec = new double[sample_size];
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		W_vec[i_s] = Y0[i_s][0] * (1 - Y0[i_s][0]);
	} // i_s loop.

	// Get the sigma inverse matrix; this is basically W matrix of size N by N. It is in fact a vector.
	double** W = get_diag_matrix(W_vec, sample_size, NULL);
	
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);
	double** XtW = matrix_multiply(Xt, n_feats, sample_size, W, sample_size, sample_size, NULL);

	// This will be needed to compute the G'WX matrices.
	double** WX = transpose_matrix(XtW, n_feats, sample_size, NULL);

	// This is the matrix that we will compute and store for the scale parameter pooling.
	double** GtWX = allocate_matrix(var_block_size, n_feats);

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);

	// Obs and Null Pheno; these are needed to save the chisqr stat.
	double** Y = get_matrix_per_col_vector(per_subject_obs_pheno, NULL);
	double** Y_min_Y0 = matrix_subtract(Y, sample_size, 1, Y0, sample_size, 1, NULL);

	// These are needed for chisqr stat in the pooling stage.
	double** Gt_Y_min_Y0 = allocate_matrix(var_block_size, 1);

	// These are needed for the first part of the P matrix in pooling.
	double** GtW = allocate_matrix(var_block_size, sample_size);
	double** GtWG = allocate_matrix(var_block_size, 1);

	fprintf(stderr, "Processing genotypes by %d-variant blocks..\n", var_block_size);
	FILE* f_geno = open_f(GMMAT_text_genotype_matrix_fp, "r");
	int cur_var_i = 0;
	bool file_ended = false;
	while (!file_ended)
	{
		// read and load the next block.
		int cur_block_loaded_n_vars = 0;
		int var_block_start_i = cur_var_i;
		vector<char*>* cur_block_var_ids = new vector<char*>();
		for (int i_var = cur_var_i;
			i_var < cur_var_i + var_block_size;
			i_var++)
		{
			char* cur_line = getline(f_geno);
			if (cur_line == NULL)
			{
				file_ended = true;
				break;
			}

			n_processed_vars++;
			if (n_processed_vars % 100 == 0)
			{
				fprintf(stderr, "@ %d. variant.           \r", n_processed_vars);
			}

			char buff[100];
			int i_cur_char = 0;
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			char* var_id = t_string::copy_me_str(buff);
			cur_block_var_ids->push_back(var_id);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* ref_all = t_string::copy_me_str(buff);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* alt_all = t_string::copy_me_str(buff);

			for (int i_s = 0; i_s < sample_size; i_s++)
			{
				t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
				double cur_geno = atof(buff);
				if (cur_geno != 0 &&
					cur_geno != 1 &&
					cur_geno != 2)
				{
					fprintf(stderr, "Found invalid geno: %.0f for %s\n", cur_geno, var_id);
					exit(0);
				}

				G[i_s][cur_block_loaded_n_vars] = cur_geno;
				Gt[cur_block_loaded_n_vars][i_s] = cur_geno;
			} // i_t loop.

			// Update the # of loaded vars here.
			cur_block_loaded_n_vars++;
			delete[] cur_line;
		} // block loading option.

		fprintf(stderr, "Processed %d variants.\n", cur_block_loaded_n_vars);

		// Make sure that we have block size many variants processed.
		if (cur_block_loaded_n_vars != var_block_size)
		{
			break;
		}

		int var_block_end_i = cur_var_i + cur_block_loaded_n_vars - 1;
		cur_var_i += cur_block_loaded_n_vars;

		fprintf(stderr, "Saving Block[%d-%d]\n", var_block_start_i, var_block_end_i);

		// Compute chisqr stat for the current block of variants.
		matrix_multiply(Gt, var_block_size, sample_size, Y_min_Y0, sample_size, 1, Gt_Y_min_Y0);

		// Save the absolute chi-square stats matrix.
		char cur_block_Gt_Ymin_Y0_fp[1000];
		char var_block_ids_fp[1000];
		get_chisq_stat_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_Gt_Ymin_Y0_fp, var_block_ids_fp);
		save_matrix_binary(Gt_Y_min_Y0, var_block_size, 1, cur_block_Gt_Ymin_Y0_fp);

		// Write the variant identifiers.
		FILE* f_var_block_ids = open_f(var_block_ids_fp, "w");
		for (int i_var = 0; i_var < (int)cur_block_var_ids->size(); i_var++)
		{
			fprintf(f_var_block_ids, "%s\n", cur_block_var_ids->at(i_var));
		} // i_var loop.
		close_f(f_var_block_ids, var_block_ids_fp);

		// Save plaintext for debugging.
		char cur_block_Gt_Ymin_Y0_fp_txt_fp[1000];
		sprintf(cur_block_Gt_Ymin_Y0_fp_txt_fp, "chisq_stat_site_%d.txt", client_i);
		save_matrix_plain(Gt_Y_min_Y0, var_block_size, 1, cur_block_Gt_Ymin_Y0_fp_txt_fp);

		// Compute GtWX matrix.
		matrix_multiply(Gt, var_block_size, sample_size, WX, sample_size, n_feats, GtWX);

		// Save the GtWX.
		char cur_block_GtWX_fp[1000];
		get_GtWX_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_GtWX_fp);
		save_matrix_binary(GtWX, var_block_size, n_feats, cur_block_GtWX_fp);

		// We also need G'WG term: This is computed using more optimal multiplications using diagonal structure of W.
		// This is not done in encrypted domain, it is done at each site locally.
		matrix_right_multiply_with_diag_matrix(Gt, var_block_size, sample_size, W, sample_size, sample_size, GtW);
		matrix_row_by_row_inner_product(GtW, var_block_size, sample_size, Gt, var_block_size, sample_size, GtWG);

		// Save GtWG.
		char cur_block_GtWG_fp[1000];
		get_GtWG_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_GtWG_fp);
		save_matrix_binary(GtWG, var_block_size, 1, cur_block_GtWG_fp);
	} // genotype file reading loop.
	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
} // client_calculate_save_pvalue_stats option.

// This function does update and pooling. Note that we have all the encrypted data for pooling and p-value estimation at this stage.
void cryptable_client_pool_pvalue_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, 
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "Client %d: Updating chisqr-scale stats..\n", client_i);

	fprintf(stderr, "Loading the features data from %s.\n", subject_per_row_feats_fp);
	int n_covars_per_feats_file = 0;
	vector<double*>* per_subject_feat_matrix = new vector<double*>();
	vector<double>* per_subject_obs_pheno = new vector<double>();
	load_LR_data_row_samples(subject_per_row_feats_fp, obs_pheno_fp, n_covars_per_feats_file, per_subject_feat_matrix, per_subject_obs_pheno);

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars_per_feats_file);
	fprintf(stderr, "%d features loaded from covariates file.\n", n_feats);

	double** X = get_matrix_per_row_vector_list(per_subject_feat_matrix, n_feats, NULL);

	int sample_size = (int)(per_subject_feat_matrix->size());

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculate null phenotypes.
	double** nu = allocate_matrix(sample_size, 1);
	double** Y0 = allocate_matrix(sample_size, 1);

	// Load the fully decrypted beta, not the original beta, which is supposed to be encrypted.
	char beta_matrix_fp[1000];
	get_fulldec_beta_fp(shared_working_dir, i_iter, client_i, beta_matrix_fp);
	int loaded_nrow, loaded_ncol;
	double** cur_epoch_weights = load_matrix_binary(beta_matrix_fp, loaded_nrow, loaded_ncol);
	if (loaded_nrow != n_feats || loaded_ncol != 1)
	{
		fprintf(stderr, "Sanity check failed while loading %s: %d, %d ; %d, %d\n", beta_matrix_fp, loaded_nrow, loaded_ncol, n_feats, 1);
		exit(0);
	}

	// Calculate nu and null phenotype matrix.
	matrix_multiply(X, sample_size, n_feats, cur_epoch_weights, n_feats, 1, nu);

	// Calculate mu from nu
	process_matrix_elementwise_by_callback(nu, sample_size, 1, get_sigmoid_val_per_feat_comb, Y0);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)per_subject_obs_pheno->size(), (int)per_subject_feat_matrix->size());
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\n", i_s, Y0[i_s][0] - per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	// Pool XtWX matrix and invert it.
	double** XtWX_inv = allocate_matrix(n_feats, n_feats);

	fprintf(stderr, "Loading inv(XtWX)..\n");
	char inv_XtWX_fp[1000];
	get_inv_XtWX_matrix_fp(private_working_dir, i_iter, client_i, inv_XtWX_fp);
	int loaded_nrows, loaded_ncols;
	load_matrix_binary(inv_XtWX_fp, loaded_nrows, loaded_ncols, XtWX_inv);
	if (loaded_nrows != n_feats ||
		loaded_ncols != n_feats)
	{
		fprintf(stderr, "Sanity check failed while loading %s: %d, %d; %d, %d\n", inv_XtWX_fp, loaded_nrows, loaded_ncols, n_feats, n_feats);
		exit(0);
	}

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);

	double** cur_block_cur_client_GtWX_matrix = allocate_matrix(var_block_size, n_feats);
	double** cur_block_pooled_GtWX_matrix = allocate_matrix(var_block_size, n_feats);

	double** cur_block_cur_client_GtWG_matrix = allocate_matrix(var_block_size, 1);
	double** cur_block_pooled_GtWG_matrix = allocate_matrix(var_block_size, 1);

	double** cur_block_cur_client_Gt_Y_min_Y0_matrix = allocate_matrix(var_block_size, 1);
	double** cur_block_pooled_Gt_Y_min_Y0_matrix = allocate_matrix(var_block_size, 1);

	double** cur_block_pooled_GtWX_Xi = allocate_matrix(var_block_size, n_feats);

	fprintf(stderr, "Processing genotypes by %d-variant blocks..\n", var_block_size);
	FILE* f_geno = open_f(GMMAT_text_genotype_matrix_fp, "r");
	int cur_var_i = 0;
	bool file_ended = false;
	while (!file_ended)
	{
		// read and load the next block.
		int cur_block_loaded_n_vars = 0;
		int var_block_start_i = cur_var_i;
		vector<char*>* cur_block_var_ids = new vector<char*>();
		for (int i_var = cur_var_i;
			i_var < cur_var_i + var_block_size;
			i_var++)
		{
			char* cur_line = getline(f_geno);
			if (cur_line == NULL)
			{
				file_ended = true;
				break;
			}

			n_processed_vars++;
			if (n_processed_vars % 100 == 0)
			{
				fprintf(stderr, "@ %d. variant.           \r", n_processed_vars);
			}

			char buff[100];
			int i_cur_char = 0;
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			char* var_id = t_string::copy_me_str(buff);
			cur_block_var_ids->push_back(var_id);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* ref_all = t_string::copy_me_str(buff);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* alt_all = t_string::copy_me_str(buff);

			for (int i_s = 0; i_s < sample_size; i_s++)
			{
				t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
				double cur_geno = atof(buff);
				if (cur_geno != 0 &&
					cur_geno != 1 &&
					cur_geno != 2)
				{
					fprintf(stderr, "Found invalid geno: %.0f for %s\n", cur_geno, var_id);
					exit(0);
				}

				G[i_s][cur_block_loaded_n_vars] = cur_geno;
				Gt[cur_block_loaded_n_vars][i_s] = cur_geno;
			} // i_t loop.

			// Update the # of loaded vars here.
			cur_block_loaded_n_vars++;
			delete[] cur_line;
		} // block loading option.

		fprintf(stderr, "Processed %d variants.\n", cur_block_loaded_n_vars);

		// Make sure that we have block size many variants processed.
		if (cur_block_loaded_n_vars != var_block_size)
		{
			break;
		}

		int var_block_end_i = cur_var_i + cur_block_loaded_n_vars - 1;
		cur_var_i += cur_block_loaded_n_vars;

		fprintf(stderr, "Saving Block[%d-%d]\n", var_block_start_i, var_block_end_i);

		// Reset the pooled matrix for this block.
		set_matrix_val(cur_block_pooled_GtWX_matrix, var_block_size, n_feats, 0);
		set_matrix_val(cur_block_pooled_GtWG_matrix, var_block_size, 1, 0);
		set_matrix_val(cur_block_pooled_Gt_Y_min_Y0_matrix, var_block_size, 1, 0);
		vector<char*>* var_ids = NULL;

		// This pools data from other sites, which means it will be performed in encrypted domain.
		for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
		{
			// Pool GWX matrices.
			char cur_client_GtWX_fp[1000];
			get_GtWX_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, cur_client_i, cur_client_GtWX_fp);
			load_matrix_binary(cur_client_GtWX_fp, loaded_nrow, loaded_ncol, cur_block_cur_client_GtWX_matrix);
			if (loaded_nrow != var_block_size || loaded_ncol != n_feats)
			{
				fprintf(stderr, "Sanity check failed loading %d. client's GtWX matrix: %d, %d ; %d, %d\n", cur_client_i, loaded_nrow, loaded_ncol, var_block_size, n_feats);
				exit(0);
			}

			matrix_add(cur_block_pooled_GtWX_matrix, var_block_size, n_feats, cur_block_cur_client_GtWX_matrix, var_block_size, n_feats, cur_block_pooled_GtWX_matrix);

			// Pool the GtWG matrix.
			char cur_client_GtWG_fp[1000];
			get_GtWG_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, cur_client_i, cur_client_GtWG_fp);
			load_matrix_binary(cur_client_GtWG_fp, loaded_nrow, loaded_ncol, cur_block_cur_client_GtWG_matrix);
			if (loaded_nrow != var_block_size || loaded_ncol != 1)
			{
				fprintf(stderr, "Sanity check failed loading %d. client's GtWX matrix: %d, %d ; %d, %d\n", cur_client_i, loaded_nrow, loaded_ncol, var_block_size, 1);
				exit(0);
			}

			matrix_add(cur_block_pooled_GtWG_matrix, var_block_size, 1, cur_block_cur_client_GtWG_matrix, var_block_size, 1, cur_block_pooled_GtWG_matrix);

			// Pool the Gt(y-y0) matrix.
			char cur_client_Gt_YminY0_fp[1000];
			char var_ids_fp[1000];
			get_chisq_stat_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, cur_client_i, cur_client_Gt_YminY0_fp, var_ids_fp);
			load_matrix_binary(cur_client_Gt_YminY0_fp, loaded_nrow, loaded_ncol, cur_block_cur_client_Gt_Y_min_Y0_matrix);
			if (loaded_nrow != var_block_size || loaded_ncol != 1)
			{
				fprintf(stderr, "Sanity check failed loading %d. client's Gt_YminY0 matrix: %d, %d ; %d, %d\n", cur_client_i, loaded_nrow, loaded_ncol, var_block_size, 1);
				exit(0);
			}

			// Load the var ids.
			if (var_ids == NULL)
			{
				var_ids = buffer_file(var_ids_fp);
			}

			matrix_add(cur_block_pooled_Gt_Y_min_Y0_matrix, var_block_size, 1, cur_block_cur_client_Gt_Y_min_Y0_matrix, var_block_size, 1, cur_block_pooled_Gt_Y_min_Y0_matrix);
		} // cur_client_i loop.

		// Save the pooled GtWX matrix.
		char cur_block_pooled_GtWX_matrix_fp[1000];
		sprintf(cur_block_pooled_GtWX_matrix_fp, "%s/pooled_GtWX_block_%d_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i);
		save_matrix_binary(cur_block_pooled_GtWX_matrix, var_block_size, n_feats, cur_block_pooled_GtWX_matrix_fp);

		// Multiply on the right by Z=XtWX_inv; use row expansion of Z, which is calculated before.
		matrix_multiply(cur_block_pooled_GtWX_matrix, var_block_size, n_feats, XtWX_inv, n_feats, n_feats, cur_block_pooled_GtWX_Xi);

		// Save the pooled GtWX matrix.
		char cur_block_pooled_GtWXZ_matrix_fp[1000];
		sprintf(cur_block_pooled_GtWXZ_matrix_fp, "%s/pooled_GtWXZ_block_%d_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i);
		save_matrix_binary(cur_block_pooled_GtWX_Xi, var_block_size, n_feats, cur_block_pooled_GtWXZ_matrix_fp);

		// We have everything we need to calculate p-values: pooled(Gt(y-y0)) ; pooled(GtWG) - pooled(GtWXZXWG)
		// Following return a vx1 matrix.
		double** cur_block_pooled_GtWx_Xi_XtWG = matrix_row_by_row_inner_product(cur_block_pooled_GtWX_Xi, var_block_size, n_feats, cur_block_pooled_GtWX_matrix, var_block_size, n_feats, NULL);

		char cur_block_pooled_GtWXZXtWG_matrix_fp[1000];
		sprintf(cur_block_pooled_GtWXZXtWG_matrix_fp, "%s/pooled_GtWXZXtWG_block_%d_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i);
		save_matrix_binary(cur_block_pooled_GtWx_Xi_XtWG, var_block_size, n_feats, cur_block_pooled_GtWXZXtWG_matrix_fp);

		char cur_block_pooled_GtWG_matrix_fp[1000];
		sprintf(cur_block_pooled_GtWG_matrix_fp, "%s/pooled_GtWG_block_%d_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i);
		save_matrix_binary(cur_block_pooled_GtWG_matrix, var_block_size, n_feats, cur_block_pooled_GtWG_matrix_fp);

		double** cur_block_pooled_chi_sq_scale_stats = matrix_subtract(cur_block_pooled_GtWG_matrix, var_block_size, 1, cur_block_pooled_GtWx_Xi_XtWG, var_block_size, 1, NULL);
		
		// Square the stat.
		double** cur_block_pooled_chi_sq_stat = matrix_multiply_elementwise(cur_block_pooled_Gt_Y_min_Y0_matrix, var_block_size, 1, cur_block_pooled_Gt_Y_min_Y0_matrix, var_block_size, 1, NULL);

		// Calculate normalized stats.
		double** chi_sq_norm_stats = matrix_divide_elementwise(cur_block_pooled_chi_sq_stat, var_block_size, 1, cur_block_pooled_chi_sq_scale_stats, var_block_size, 1, NULL);

		// Get p-values.
		double** p_values = process_matrix_elementwise_by_callback(chi_sq_norm_stats, var_block_size, 1, get_1DOF_chisqr_pval, NULL);

		double** p_val_stats_matrix = allocate_matrix(var_block_size, 5);
		for (int var_i = 0; var_i < var_block_size; var_i++)
		{
			p_val_stats_matrix[var_i][0] = p_values[var_i][0];
			p_val_stats_matrix[var_i][1] = chi_sq_norm_stats[var_i][0];
			p_val_stats_matrix[var_i][2] = cur_block_pooled_chi_sq_stat[var_i][0];
			p_val_stats_matrix[var_i][3] = cur_block_pooled_chi_sq_scale_stats[var_i][0];
			p_val_stats_matrix[var_i][4] = p_values[var_i][0];
		} // var_i loop.

		// Save the p-values stats for this client.
		char p_val_stats_fp[1000];
		get_pval_stat_matrix_op_fp(shared_working_dir, i_iter, var_block_start_i, var_block_end_i, client_i, p_val_stats_fp);
		save_matrix_binary(p_val_stats_matrix, var_block_size, 5, p_val_stats_fp);
	} // genotype file reading loop.

	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
} // client_update_pvalue_scale_stats_per_block

// This function does update and pooling. Note that we have all the encrypted data for pooling and p-value estimation at this stage.
void cryptable_client_collaborative_decrypt_pval_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp,
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "Client %d: Updating chisqr-scale stats..\n", client_i);

	fprintf(stderr, "Loading the features data from %s.\n", subject_per_row_feats_fp);
	int n_covars_per_feats_file = 0;
	vector<double*>* per_subject_feat_matrix = new vector<double*>();
	vector<double>* per_subject_obs_pheno = new vector<double>();
	load_LR_data_row_samples(subject_per_row_feats_fp, obs_pheno_fp, n_covars_per_feats_file, per_subject_feat_matrix, per_subject_obs_pheno);

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars_per_feats_file);
	fprintf(stderr, "%d features loaded from covariates file.\n", n_feats);

	double** X = get_matrix_per_row_vector_list(per_subject_feat_matrix, n_feats, NULL);

	int sample_size = (int)(per_subject_feat_matrix->size());

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculate null phenotypes.
	double** nu = allocate_matrix(sample_size, 1);
	double** Y0 = allocate_matrix(sample_size, 1);

	// Load the fully decrypted beta, not the original beta, which is supposed to be encrypted.
	char beta_matrix_fp[1000];
	get_fulldec_beta_fp(shared_working_dir, i_iter, client_i, beta_matrix_fp);
	int loaded_nrow, loaded_ncol;
	double** cur_epoch_weights = load_matrix_binary(beta_matrix_fp, loaded_nrow, loaded_ncol);
	if (loaded_nrow != n_feats || loaded_ncol != 1)
	{
		fprintf(stderr, "Sanity check failed while loading %s: %d, %d ; %d, %d\n", beta_matrix_fp, loaded_nrow, loaded_ncol, n_feats, 1);
		exit(0);
	}

	// Calculate nu and null phenotype matrix.
	matrix_multiply(X, sample_size, n_feats, cur_epoch_weights, n_feats, 1, nu);

	// Calculate mu from nu
	process_matrix_elementwise_by_callback(nu, sample_size, 1, get_sigmoid_val_per_feat_comb, Y0);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)per_subject_obs_pheno->size(), (int)per_subject_feat_matrix->size());
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\n", i_s, Y0[i_s][0] - per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	//int sample_size = per_subject_feat_matrix->size();

	// Compute the projection matrix: 
	fprintf(stderr, "Setting Sigma inv..\n");
	double* W_vec = new double[sample_size];
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		W_vec[i_s] = Y0[i_s][0] * (1 - Y0[i_s][0]);
	} // i_s loop.

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);

	fprintf(stderr, "Processing genotypes by %d-variant blocks..\n", var_block_size);
	FILE* f_geno = open_f(GMMAT_text_genotype_matrix_fp, "r");
	int cur_var_i = 0;
	bool file_ended = false;
	while (!file_ended)
	{
		// read and load the next block.
		int cur_block_loaded_n_vars = 0;
		int var_block_start_i = cur_var_i;
		vector<char*>* cur_block_var_ids = new vector<char*>();
		for (int i_var = cur_var_i;
			i_var < cur_var_i + var_block_size;
			i_var++)
		{
			char* cur_line = getline(f_geno);
			if (cur_line == NULL)
			{
				file_ended = true;
				break;
			}

			n_processed_vars++;
			if (n_processed_vars % 100 == 0)
			{
				fprintf(stderr, "@ %d. variant.           \r", n_processed_vars);
			}

			char buff[100];
			int i_cur_char = 0;
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			char* var_id = t_string::copy_me_str(buff);
			cur_block_var_ids->push_back(var_id);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* ref_all = t_string::copy_me_str(buff);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* alt_all = t_string::copy_me_str(buff);

			for (int i_s = 0; i_s < sample_size; i_s++)
			{
				t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
				double cur_geno = atof(buff);
				if (cur_geno != 0 &&
					cur_geno != 1 &&
					cur_geno != 2)
				{
					fprintf(stderr, "Found invalid geno: %.0f for %s\n", cur_geno, var_id);
					exit(0);
				}

				G[i_s][cur_block_loaded_n_vars] = cur_geno;
				Gt[cur_block_loaded_n_vars][i_s] = cur_geno;
			} // i_t loop.

			// Update the # of loaded vars here.
			cur_block_loaded_n_vars++;
			delete[] cur_line;
		} // block loading option.

		fprintf(stderr, "Processed %d variants.\n", cur_block_loaded_n_vars);

		// Make sure that we have block size many variants processed.
		if (cur_block_loaded_n_vars != var_block_size)
		{
			break;
		}

		int var_block_end_i = cur_var_i + cur_block_loaded_n_vars - 1;
		cur_var_i += cur_block_loaded_n_vars;

		for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
		{
			fprintf(stderr, "Partially decrypting p-value stats matrix for block[%d-%d]\n", var_block_start_i, var_block_end_i);

			// Save the p-values stats for this client.
			int loaded_nrows, loaded_ncols;
			char p_val_stats_fp[1000];
			get_pval_stat_matrix_op_fp(shared_working_dir, i_iter, var_block_start_i, var_block_end_i, cur_client_i, p_val_stats_fp);
			double** p_val_stats_matrix = load_matrix_binary(p_val_stats_fp, loaded_nrows, loaded_ncols, NULL);
			if (loaded_nrows != var_block_size || loaded_ncols != 5)
			{
				fprintf(stderr, "Sanity check failed while loading %s: %d, %d ; %d, %d\n", p_val_stats_fp, loaded_nrows, loaded_ncols, var_block_size, 5);
				exit(0);
			}

			double** partdec_pval_stats_matrix = partial_decrypt_matrix(p_val_stats_matrix, var_block_size, 5, client_i, n_clients);
			char partdec_pval_matrix_fp[1000];
			get_partdec_pval_stats_matrix_fp(shared_working_dir, i_iter, var_block_start_i, var_block_end_i, cur_client_i, client_i, partdec_pval_matrix_fp);
			save_matrix_binary(partdec_pval_stats_matrix, var_block_size, 5, partdec_pval_matrix_fp);
		} // cur_client_i loop.
	} // genotype file reading loop.

	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
}

void cryptable_client_pool_partially_decrypted_pval_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp,
	char* private_working_dir,
	char* shared_working_dir)
{
	fprintf(stderr, "Client %d: Updating chisqr-scale stats..\n", client_i);

	fprintf(stderr, "Loading the features data from %s.\n", subject_per_row_feats_fp);
	int n_covars_per_feats_file = 0;
	vector<double*>* per_subject_feat_matrix = new vector<double*>();
	vector<double>* per_subject_obs_pheno = new vector<double>();
	load_LR_data_row_samples(subject_per_row_feats_fp, obs_pheno_fp, n_covars_per_feats_file, per_subject_feat_matrix, per_subject_obs_pheno);

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars_per_feats_file);
	fprintf(stderr, "%d features loaded from covariates file.\n", n_feats);

	int sample_size = (int)(per_subject_feat_matrix->size());

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)per_subject_obs_pheno->size(), (int)per_subject_feat_matrix->size());
		exit(0);
	}

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);

	// Read partial decryptions from all sites, aggregate them and save.
	char full_dec_pval_stats_matrix_fp[1000];
	get_fulldec_pval_stats_matrix_fp(shared_working_dir, i_iter, client_i, full_dec_pval_stats_matrix_fp);
	FILE* f_fulldec_pval_stats = open_f(full_dec_pval_stats_matrix_fp, "w");

	fprintf(stderr, "Processing genotypes by %d-variant blocks..\n", var_block_size);
	FILE* f_geno = open_f(GMMAT_text_genotype_matrix_fp, "r");
	int cur_var_i = 0;
	bool file_ended = false;
	while (!file_ended)
	{
		// read and load the next block.
		int cur_block_loaded_n_vars = 0;
		int var_block_start_i = cur_var_i;
		vector<char*>* cur_block_var_ids = new vector<char*>();
		for (int i_var = cur_var_i;
			i_var < cur_var_i + var_block_size;
			i_var++)
		{
			char* cur_line = getline(f_geno);
			if (cur_line == NULL)
			{
				file_ended = true;
				break;
			}

			n_processed_vars++;
			if (n_processed_vars % 100 == 0)
			{
				fprintf(stderr, "@ %d. variant.           \r", n_processed_vars);
			}

			char buff[100];
			int i_cur_char = 0;
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			char* var_id = t_string::copy_me_str(buff);
			cur_block_var_ids->push_back(var_id);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* ref_all = t_string::copy_me_str(buff);
			t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
			//char* alt_all = t_string::copy_me_str(buff);

			for (int i_s = 0; i_s < sample_size; i_s++)
			{
				t_string::get_next_token(cur_line, buff, 100, "\t", i_cur_char);
				double cur_geno = atof(buff);
				if (cur_geno != 0 &&
					cur_geno != 1 &&
					cur_geno != 2)
				{
					fprintf(stderr, "Found invalid geno: %.0f for %s\n", cur_geno, var_id);
					exit(0);
				}

				G[i_s][cur_block_loaded_n_vars] = cur_geno;
				Gt[cur_block_loaded_n_vars][i_s] = cur_geno;
			} // i_t loop.

			// Update the # of loaded vars here.
			cur_block_loaded_n_vars++;
			delete[] cur_line;
		} // block loading option.

		fprintf(stderr, "Processed %d variants.\n", cur_block_loaded_n_vars);

		// Make sure that we have block size many variants processed.
		if (cur_block_loaded_n_vars != var_block_size)
		{
			break;
		}

		int var_block_end_i = cur_var_i + cur_block_loaded_n_vars - 1;
		cur_var_i += cur_block_loaded_n_vars;

		// Aggregate partial decryptions of the p-value matrix from other sites.
		double** fulldec_cur_block_p_val_stats_matrix = allocate_matrix(var_block_size, 5);
		for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
		{
			fprintf(stderr, "Partially decrypting p-value stats matrix for block[%d-%d]\n", var_block_start_i, var_block_end_i);

			// Save the p-values stats for this client.
			int loaded_nrows, loaded_ncols;
			char partdec_pval_matrix_fp[1000];
			get_partdec_pval_stats_matrix_fp(shared_working_dir, i_iter, var_block_start_i, var_block_end_i, client_i, cur_client_i, partdec_pval_matrix_fp);
			double** cur_block_partdec_p_val_stats_matrix = load_matrix_binary(partdec_pval_matrix_fp, loaded_nrows, loaded_ncols, NULL);
			if (loaded_nrows != var_block_size || loaded_ncols != 5)
			{
				fprintf(stderr, "Sanity check failed while loading %s: %d, %d ; %d, %d\n", partdec_pval_matrix_fp, loaded_nrows, loaded_ncols, var_block_size, 5);
				exit(0);
			}
			
			matrix_add(fulldec_cur_block_p_val_stats_matrix, var_block_size, 5, cur_block_partdec_p_val_stats_matrix, var_block_size, 5, fulldec_cur_block_p_val_stats_matrix);
		} // cur_client_i loop.

		for (int var_i = 0; var_i < var_block_size; var_i++)
		{
			fprintf(f_fulldec_pval_stats, "%s\t%.5f\t%5f\t%5f\t%5f\n", cur_block_var_ids->at(var_i), 
					fulldec_cur_block_p_val_stats_matrix[var_i][0],
					fulldec_cur_block_p_val_stats_matrix[var_i][1],
					fulldec_cur_block_p_val_stats_matrix[var_i][2], 
					fulldec_cur_block_p_val_stats_matrix[var_i][3]);
		} // var_i loop.
	} // genotype file reading loop.

	close_f(f_fulldec_pval_stats, full_dec_pval_stats_matrix_fp);
	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
} // cryptable_client_pool_partially_decrypted_pval_stats
