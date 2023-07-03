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
#include "cllgn_matrix_linalg_utils.h"

bool __DUMP_FEDLR_MESSAGES__ = false;

/*
	- This is the SGD implementation of LR for testing with genotype data.
	. https://github.com/OpenMined/TenSEAL/blob/main/tutorials/Tutorial%201%20-%20Training%20and%20Evaluation%20of%20Logistic%20Regression%20on%20Encrypted%20Data.ipynb
	. https://cseweb.ucsd.edu/~elkan/250B/logreg.pdf

	- TODOs::
	. 1-bit SGD good for parallelization?
	.
*/

using namespace std;

double get_1DOF_chisqr_pval(double norm_chisqr_scale)
{
	double pval = chisqr(1, norm_chisqr_scale);
	return(pval);
}

void get_pval_stat_op_fp(char* shared_working_dir, int client_i, char* p_val_stats_fp)
{
	sprintf(p_val_stats_fp, "%s/p_val_stats_%d.txt", shared_working_dir, client_i);
}

void get_GtWX_Xi_XtWG_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/GtWX_Xi_XtWG_VAR_BLOCK_%d_%d_CLIENT_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i, client_i);
}

void get_GtWX_matrix_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/GtWX_VAR_BLOCK_%d_%d_CLIENT_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i, client_i);
}

void get_GtWG_matrix_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/GtWG_VAR_BLOCK_%d_%d_CLIENT_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i, client_i);
}

void get_chisq_stat_matrix_fp(char* shared_working_dir, 
	int var_block_start_i, int var_block_end_i, 
	int client_i, 
	char* chisq_buffer, char* var_block_ids_buffer)
{
	sprintf(chisq_buffer, "%s/P_CHISQ_STAT_VAR_BLOCK_%d_%d_CLIENT_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i, client_i);
	sprintf(var_block_ids_buffer, "%s/VAR_BLOCK_%d_%d_IDS_CLIENT_%d.list", shared_working_dir, var_block_start_i, var_block_end_i, client_i);
}

void get_scale_matrix_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/P_CHISQ_SCALE_VAR_BLOCK_%d_%d_CLIENT_%d.bin", shared_working_dir, var_block_start_i, var_block_end_i, client_i);
}

void get_LL_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/LL_MATRIX_%d_ITER_%d.bin", shared_working_dir, client_i, iter_i);
}

void get_mu_matrix_fp(char* shared_working_dir, int i_iter, int client_i, char* buffer)
{
	sprintf(buffer, "%s/MU_MATRIX_%d_ITER_%d.bin", shared_working_dir, client_i, i_iter);
}

void get_beta_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/BETA_MATRIX_CLIENT_%d_ITER_%d.bin", shared_working_dir, client_i, iter_i);
}

void get_XtWX_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/XtWX_CLIENT_%d_ITER_%d.bin", shared_working_dir, client_i, iter_i);
}

void get_XtWz_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer)
{
	sprintf(buffer, "%s/XtWz_CLIENT_%d_ITER_%d.bin", shared_working_dir, client_i, iter_i);
}

/* 
This function loads and processes the current null model stats using IRLS fitting model.
1- Load the current beta matrix
2- Update mu, nu, z, W, ...
3- Save the intermediate matrix.
*/
void client_calculate_save_XtWX_XtWz(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* shared_working_dir)
{
	fprintf(stderr, "@ iteration %d, Client %d/%d -- Calculating XtWX and XtWz and saving to shared working space %s\n", i_iter, client_i, n_clients, shared_working_dir);

	// We need the two sigmoid functions here.
	//double(*GET_SIGMOID_PER_FEATS_WEIGHTS)(double*, int, double*);
	double(*GET_SIGMOID_PER_FEAT_COMB)(double);

	if (sigmoid_approx_type == SIGMOID_APPROX_NATIVE)
	{
		fprintf(stderr, "Using native sigmoid..\n");
		//GET_SIGMOID_PER_FEATS_WEIGHTS = get_sigmoid_val_per_feats_weights;
		GET_SIGMOID_PER_FEAT_COMB = get_sigmoid_val_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_KIM_ETAL)
	{
		fprintf(stderr, "Using Kim. etal sigmoid..\n");
		//GET_SIGMOID_PER_FEATS_WEIGHTS = get_Kim_etal_poly_approx_sigmoid_per_feats_weights;
		GET_SIGMOID_PER_FEAT_COMB = get_Kim_etal_poly_approx_sigmoid_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_TENSEAL)
	{
		fprintf(stderr, "Using TenSEAL sigmoid..\n");
		//GET_SIGMOID_PER_FEATS_WEIGHTS = get_TenSeal_poly_approx_sigmoid_per_feats_weights;
		GET_SIGMOID_PER_FEAT_COMB = get_TenSeal_poly_approx_sigmoid_per_feat_comb;
	}
	else
	{
		fprintf(stderr, "Unknown sigmoid approximation: %d\n", sigmoid_approx_type);
		exit(0);
	}

	// We need the RNG for different tasks.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	//double* per_feat_weights = NULL;
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
		for (int i_s = 0; i_s < (int)(per_ind_feat_matrix->size()); i_s++)
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

	if (__DUMP_FEDLR_MESSAGES__)
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

	double** cur_epoch_weights = NULL;

	// This is the main point of loading of the beta values.
	// Load the current epoch weights, if they exist.
	int prev_i_iter = (i_iter > 0) ? (i_iter - 1) : (n_epoch);
	char beta_matrix_fp[1000];

	// Load the beta from previous iteration.
	get_beta_matrix_fp(shared_working_dir, prev_i_iter, client_i, beta_matrix_fp);
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

	//double cur_LL = 0;

	// Copy the features matrix and allocat the related matrices.
	double** X = get_matrix_per_row_vector_list(per_ind_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);
	double** XtW = allocate_matrix(n_feats, sample_size);
	double** XtWX = allocate_matrix(n_feats, n_feats);
	//double** XtWX_inv = allocate_matrix(n_feats, n_feats);
	double** XtWz = allocate_matrix(n_feats, 1);

	double max_abs_nu = 0;

	// Calculate nu and mu, and z=nu + (y[i]-mu[i])/(mu[i]*(1-mu[i])
	matrix_multiply(X, sample_size, n_feats, cur_epoch_weights, n_feats, 1, nu);

	// Calculate mu from nu
	process_matrix_elementwise_by_callback(nu, sample_size, 1, GET_SIGMOID_PER_FEAT_COMB, mu);

	// Save the mu matrix.
	char mu_matrix_fp[1000];
	get_mu_matrix_fp(shared_working_dir, i_iter, client_i, mu_matrix_fp);
	save_matrix_binary(mu, sample_size, 1, mu_matrix_fp);

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
	get_XtWX_matrix_fp(shared_working_dir, i_iter, client_i, XtWX_matrix_fp);
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
} // client_calculate_save_XtWX_XtWz option.

void client_pool_XtWX_XtWz_update_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* shared_working_dir)
{
	fprintf(stderr, "@ iteration %d, Client %d/%d -- Pooling XtWX and XtWz and updating beta and saving to shared working space %s\n", i_iter, client_i, n_clients, shared_working_dir);

	//// We need the two sigmoid functions here.
	////double(*GET_SIGMOID_PER_FEATS_WEIGHTS)(double*, int, double*);
	//double(*GET_SIGMOID_PER_FEAT_COMB)(double);

	//if (sigmoid_approx_type == SIGMOID_APPROX_NATIVE)
	//{
	//	fprintf(stderr, "Using native sigmoid..\n");
	//	//GET_SIGMOID_PER_FEATS_WEIGHTS = get_sigmoid_val_per_feats_weights;
	//	GET_SIGMOID_PER_FEAT_COMB = get_sigmoid_val_per_feat_comb;
	//}
	//else if (sigmoid_approx_type == SIGMOID_APPROX_KIM_ETAL)
	//{
	//	fprintf(stderr, "Using Kim. etal sigmoid..\n");
	//	//GET_SIGMOID_PER_FEATS_WEIGHTS = get_Kim_etal_poly_approx_sigmoid_per_feats_weights;
	//	GET_SIGMOID_PER_FEAT_COMB = get_Kim_etal_poly_approx_sigmoid_per_feat_comb;
	//}
	//else if (sigmoid_approx_type == SIGMOID_APPROX_TENSEAL)
	//{
	//	fprintf(stderr, "Using TenSEAL sigmoid..\n");
	//	//GET_SIGMOID_PER_FEATS_WEIGHTS = get_TenSeal_poly_approx_sigmoid_per_feats_weights;
	//	GET_SIGMOID_PER_FEAT_COMB = get_TenSeal_poly_approx_sigmoid_per_feat_comb;
	//}
	//else
	//{
	//	fprintf(stderr, "Unknown sigmoid approximation: %d\n", sigmoid_approx_type);
	//	exit(0);
	//}

	// We need the RNG for different tasks.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	//double* per_feat_weights = NULL;
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
		for (int i_s = 0; i_s < (int)(per_ind_feat_matrix->size()); i_s++)
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

	if (__DUMP_FEDLR_MESSAGES__)
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

	// These will be updated now.
	double** cur_epoch_weights = allocate_matrix(n_feats, 1);

	//double** nu = allocate_matrix(sample_size, 1);
	//double** mu = allocate_matrix(sample_size, 1);
	//double** one_min_mu = allocate_matrix(sample_size, 1);
	//double** mu_one_min_mu = allocate_matrix(sample_size, 1);
	//double** z = allocate_matrix(sample_size, 1);
	//double** Yt = allocate_copy_matrix(&(per_ind_pheno), 1, sample_size);
	//double** Y = transpose_matrix(Yt, 1, sample_size, NULL);

	//double** W_diag = allocate_matrix(sample_size, 1);
	//double** W = allocate_matrix(sample_size, sample_size);

	//double** Y_min_mu = allocate_matrix(sample_size, 1);
	//double** Y_min_mu_by_mu_one_min_mu = allocate_matrix(sample_size, 1);

	//double** one_vec_sample = allocate_matrix(sample_size, 1, 1.0);

	//double cur_LL = 0;

	// Copy the features matrix and allocat the related matrices.
	//double** X = get_matrix_per_row_vector_list(per_ind_feat_matrix, n_feats, NULL);
	//double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);
	//double** XtW = allocate_matrix(n_feats, sample_size);
	double** XtWX = allocate_matrix(n_feats, n_feats);
	double** XtWX_inv = allocate_matrix(n_feats, n_feats);
	double** XtWz = allocate_matrix(n_feats, 1);

	// Load the XtWX matrices from all clients and pool them.
	fprintf(stderr, "Pooling XtWX and XtWz matrices..\n");
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		int n_loaded_row, n_loaded_col;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Load and pool XtWX matrix.
		char cur_client_XtWX_matrix_fp[1000];
		get_XtWX_matrix_fp(shared_working_dir, i_iter, cur_client_i, cur_client_XtWX_matrix_fp);
		fprintf(stderr, "Loading XtWX matrix from %s\n", cur_client_XtWX_matrix_fp);

		double** cur_XtWX = load_matrix_binary(cur_client_XtWX_matrix_fp, n_loaded_row, n_loaded_col);
		if (n_loaded_col != n_feats || n_loaded_row != n_feats)
		{
			fprintf(stderr, "Loaded XtWX is non-conformant for client %d @ iteration %d: %d, %d ; %d, %d\n", cur_client_i, i_iter, n_feats, n_feats, n_loaded_row, n_loaded_col);
			exit(0);
		}

		// Pool this matrix, this can be done in place.
		fprintf(stderr, "Pooling XtWX.\n");
		matrix_add(XtWX, n_feats, n_feats, cur_XtWX, n_feats, n_feats, XtWX);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Load and pool XtWz matrix.
		char cur_client_XtWz_matrix_fp[1000];
		get_XtWz_matrix_fp(shared_working_dir, i_iter, cur_client_i, cur_client_XtWz_matrix_fp);
		fprintf(stderr, "Loading XtWz matrix from %s\n", cur_client_XtWz_matrix_fp);

		double** cur_XtWz = load_matrix_binary(cur_client_XtWz_matrix_fp, n_loaded_row, n_loaded_col);
		if (n_loaded_col != 1 || n_loaded_row != n_feats)
		{
			fprintf(stderr, "Loaded XtWz is non-conformant for client %d @ iteration %d: %d, %d ; %d, %d\n", cur_client_i, i_iter, n_feats, 1, n_loaded_row, n_loaded_col);
			exit(0);
		}

		// Pool this matrix, this can be done in place.
		fprintf(stderr, "Pooling XtWz.\n");
		matrix_add(XtWz, n_feats, 1, cur_XtWz, n_feats, 1, XtWz);
	} // client_i loop.

	//double max_abs_nu = 0;

	//FILE* f_per_epoch_info = open_f("per_epoch_info.txt", "w");

	// Process current iteration:
	// Invert the pooled XtWX matrix.
	invert_matrix_GJ(XtWX, n_feats, n_feats, XtWX_inv);

	// Test the inversion step.
	double** inv_tester = matrix_multiply(XtWX_inv, n_feats, n_feats, XtWX, n_feats, n_feats, NULL);
	fprintf(stderr, "Inversion test:\n");
	print_matrix(inv_tester, n_feats, n_feats);

	// This is the update step: We have the pooled XtWz matrix from all sites. Just do multiplication.
	matrix_multiply(XtWX_inv, n_feats, n_feats, XtWz, n_feats, 1, cur_epoch_weights);

	// Save the contribution of this site to the epoch weights.
	char beta_matrix_fp[1000];
	get_beta_matrix_fp(shared_working_dir, i_iter, client_i, beta_matrix_fp);
	save_matrix_binary(cur_epoch_weights, n_feats, 1, beta_matrix_fp);

	char params_txt_fp[1000];
	sprintf(params_txt_fp, "current_params_iter_%d_client_%d.txt", i_iter, client_i);
	save_matrix_plain(cur_epoch_weights, n_feats, 1, params_txt_fp);
} // client_pool_XtWX_XtWz_update_beta

void client_check_convergence_per_updated_beta(int i_iter, int client_i, int n_clients, char* shared_working_dir, double LL_EPSILON)
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

void client_calculate_save_pvalue_stats(int client_i, int n_clients, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_matrix_fp,
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

	// Laod the null model phenotypes.
	//vector<char*>* null_model_pheno_lines = buffer_file(null_model_pheno_fp);
	int loaded_nrow, loaded_ncol;
	double** null_pheno_matrix = load_matrix_binary(null_model_pheno_matrix_fp, loaded_nrow, loaded_ncol);
	if (loaded_nrow != (int)(per_subject_obs_pheno->size()) ||
		loaded_ncol != 1)
	{
		fprintf(stderr, "null phenotype matrix dimensions non-conformant: %d, %d\n", loaded_nrow, loaded_ncol);
		exit(0);
	}

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)(per_subject_obs_pheno->size()), (int)(per_subject_feat_matrix->size()));
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\n", i_s, null_pheno_matrix[i_s][0] - per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	int sample_size = (int)(per_subject_feat_matrix->size());

	// Compute the projection matrix: 
	fprintf(stderr, "Setting Sigma inv..\n");
	double* sigma_inv_vec = new double[sample_size];
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		sigma_inv_vec[i_s] = null_pheno_matrix[i_s][0] * (1 - null_pheno_matrix[i_s][0]);
	} // i_s loop.

	// Get the sigma inverse matrix.
	double** sigma_inv = get_diag_matrix(sigma_inv_vec, sample_size, NULL);
	double** X = get_matrix_per_row_vector_list(per_subject_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);

	double** Xtsigma_inv = matrix_multiply(Xt, n_feats, sample_size, sigma_inv, sample_size, sample_size, NULL);

	// This will be needed to compute the G'WX matrices.
	double** sigma_invX = transpose_matrix(Xtsigma_inv, n_feats, sample_size, NULL);

	// This is the matrix that we will compute and store for the scale parameter pooling.
	double** Gtsigma_invX = allocate_matrix(var_block_size, n_feats);

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);

	// Obs and Null Pheno; these are needed to save the chisqr stat.
	double** Y = get_matrix_per_col_vector(per_subject_obs_pheno, NULL);
	//double** Yt = transpose_matrix(Y, sample_size, 1, NULL);
	double** Y0 = null_pheno_matrix;
	//double** Y0t = transpose_matrix(Y0, sample_size, 1, NULL);
	double** Y_min_Y0 = matrix_subtract(Y, sample_size, 1, Y0, sample_size, 1, NULL);

	// These are needed for chisqr stat in the pooling stage.
	double** Gt_Y_min_Y0 = allocate_matrix(var_block_size, 1);

	// These are needed for the first part of the P matrix in pooling.
	double** Gtsigma_inv = allocate_matrix(var_block_size, sample_size);
	double** Gtsigma_invG = allocate_matrix(var_block_size, var_block_size);

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
		for (int i_var = 0; i_var < (int)(cur_block_var_ids->size()); i_var++)
		{
			fprintf(f_var_block_ids, "%s\n", cur_block_var_ids->at(i_var));
		} // i_var loop.
		close_f(f_var_block_ids, var_block_ids_fp);

		// Save plaintext for debugging.
		char cur_block_Gt_Ymin_Y0_fp_txt_fp[1000];
		sprintf(cur_block_Gt_Ymin_Y0_fp_txt_fp, "chisq_stat_site_%d.txt", client_i);
		save_matrix_plain(Gt_Y_min_Y0, var_block_size, 1, cur_block_Gt_Ymin_Y0_fp_txt_fp);

		// Compute GtWX matrix.
		matrix_multiply(Gt, var_block_size, sample_size, sigma_invX, sample_size, n_feats, Gtsigma_invX);

		// Save the matrix.
		char cur_block_Gtsigma_invX_fp[1000];
		get_GtWX_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_Gtsigma_invX_fp);
		save_matrix_binary(Gtsigma_invX, var_block_size, n_feats, cur_block_Gtsigma_invX_fp);

		// We also need G'WG term.
		matrix_multiply(Gt, var_block_size, sample_size, sigma_inv, sample_size, sample_size, Gtsigma_inv);
		matrix_multiply(Gtsigma_inv, var_block_size, sample_size, G, sample_size, var_block_size, Gtsigma_invG);

		char cur_block_GtWG_fp[1000];
		get_GtWG_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_GtWG_fp);
		save_matrix_binary(Gtsigma_invG, var_block_size, var_block_size, cur_block_GtWG_fp);
	} // genotype file reading loop.
	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
} // client_calculate_save_pvalue_stats option.

void client_update_pvalue_scale_stats_per_block(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_matrix_fp,
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

	// Laod the null model phenotypes.
	//vector<char*>* null_model_pheno_lines = buffer_file(null_model_pheno_fp);
	int loaded_nrow, loaded_ncol;
	double** null_pheno_matrix = load_matrix_binary(null_model_pheno_matrix_fp, loaded_nrow, loaded_ncol);
	if (loaded_nrow != (int)(per_subject_obs_pheno->size()) ||
		loaded_ncol != 1)
	{
		fprintf(stderr, "null phenotype matrix dimensions non-conformant: %d, %d\n", loaded_nrow, loaded_ncol);
		exit(0);
	}

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)(per_subject_obs_pheno->size()), (int)(per_subject_feat_matrix->size()));
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\n", i_s, null_pheno_matrix[i_s][0] - per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	int sample_size = (int)(per_subject_feat_matrix->size());

	// Compute the projection matrix: 
	fprintf(stderr, "Setting Sigma inv..\n");
	double* sigma_inv_vec = new double[sample_size];
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		sigma_inv_vec[i_s] = null_pheno_matrix[i_s][0] * (1 - null_pheno_matrix[i_s][0]);
	} // i_s loop.

	// Get the sigma inverse matrix.
	double** sigma_inv = get_diag_matrix(sigma_inv_vec, sample_size, NULL);
	double** X = get_matrix_per_row_vector_list(per_subject_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);

	double** Xtsigma_inv = matrix_multiply(Xt, n_feats, sample_size, sigma_inv, sample_size, sample_size, NULL);

	// This will be needed to compute the G'WX matrices.
	//double** sigma_invX = transpose_matrix(Xtsigma_inv, n_feats, sample_size, NULL);

	// Pool XtWX matrix and invert it.
	double** XtWX = allocate_matrix(n_feats, n_feats);
	double** XtWX_inv = allocate_matrix(n_feats, n_feats);

	// Load the XtWX matrices from all clients and pool them.
	fprintf(stderr, "Pooling XtWX and XtWz matrices..\n");
	for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
	{
		int n_loaded_row, n_loaded_col;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Load and pool XtWX matrix.
		char cur_client_XtWX_matrix_fp[1000];
		get_XtWX_matrix_fp(shared_working_dir, i_iter, cur_client_i, cur_client_XtWX_matrix_fp);
		fprintf(stderr, "Loading XtWX matrix from %s\n", cur_client_XtWX_matrix_fp);

		double** cur_XtWX = load_matrix_binary(cur_client_XtWX_matrix_fp, n_loaded_row, n_loaded_col);
		if (n_loaded_col != n_feats || n_loaded_row != n_feats)
		{
			fprintf(stderr, "Loaded XtWX is non-conformant for client %d @ iteration %d: %d, %d ; %d, %d\n", cur_client_i, i_iter, n_feats, n_feats, n_loaded_row, n_loaded_col);
			exit(0);
		}

		// Pool this matrix, this can be done in place.
		fprintf(stderr, "Pooling XtWX.\n");
		matrix_add(XtWX, n_feats, n_feats, cur_XtWX, n_feats, n_feats, XtWX);
	} // XtWX pooling loop.

	// Invert the pooled XtWX matrix.
	invert_matrix_GJ(XtWX, n_feats, n_feats, XtWX_inv);

	// Test the inversion step.
	double** inv_tester = matrix_multiply(XtWX_inv, n_feats, n_feats, XtWX, n_feats, n_feats, NULL);
	fprintf(stderr, "Inversion test:\n");
	print_matrix(inv_tester, n_feats, n_feats);

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);

	//// Obs and Null Pheno; these are needed to save the chisqr stat.
	//double** Y = get_matrix_per_col_vector(per_subject_obs_pheno, NULL);
	//double** Yt = transpose_matrix(Y, sample_size, 1, NULL);
	//double** Y0 = null_pheno_matrix;
	//double** Y0t = transpose_matrix(Y0, sample_size, 1, NULL);
	//double** Y_min_Y0 = matrix_subtract(Y, sample_size, 1, Y0, sample_size, 1, NULL);

	// These are needed for the first part of the P matrix in pooling.
	double** cur_block_cur_client_GtWX_matrix = allocate_matrix(var_block_size, n_feats);
	double** cur_block_pooled_GtWX_matrix = allocate_matrix(var_block_size, n_feats);

	double** GtWX_Xi = allocate_matrix(var_block_size, n_feats);
	double** GtWX_Xi_XtW = allocate_matrix(var_block_size, sample_size);
	double** GtWX_Xi_XtWG = allocate_matrix(var_block_size, var_block_size);

	double** cur_block_GtWG = allocate_matrix(var_block_size, var_block_size);

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
		} // cur_client_i loop.

		// Multiply on the right by XtWX_inv.
		matrix_multiply(cur_block_pooled_GtWX_matrix, var_block_size, n_feats, XtWX_inv, n_feats, n_feats, GtWX_Xi);

		// Multiply on the right by Xtsigma_inv
		matrix_multiply(GtWX_Xi, var_block_size, n_feats, Xtsigma_inv, n_feats, sample_size, GtWX_Xi_XtW);

		// Finally, multiply with G on the right to get the final contribution.
		matrix_multiply(GtWX_Xi_XtW, var_block_size, sample_size, G, sample_size, var_block_size, GtWX_Xi_XtWG);

		// Load the GtWG for this site and save it, it is the final P contribution of this site.
		char GtWG_fp[1000];
		get_GtWG_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, GtWG_fp);
		int loaded_nrow, loaded_ncol;
		load_matrix_binary(GtWG_fp, loaded_nrow, loaded_ncol, cur_block_GtWG);
		if (loaded_nrow != var_block_size || loaded_ncol != var_block_size)
		{
			fprintf(stderr, "Sanity check failed loading %d. client's GtWG matrix (%s): %d, %d ; %d, %d\n", client_i, GtWG_fp,
				loaded_nrow, loaded_ncol, var_block_size, var_block_size);
			exit(0);
		}

		// This adds G'WG - G'WX_Xi_X'WG ; which is the scale parameter for this site that pools all information from other sites.
		matrix_subtract(cur_block_GtWG, var_block_size, var_block_size, GtWX_Xi_XtWG, var_block_size, var_block_size, cur_block_GtWG);

		// Save this site's contribution to 2nd term of GP.
		char cur_site_scale_fp[1000];
		get_scale_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_site_scale_fp);
		save_matrix_binary(cur_block_GtWG, var_block_size, var_block_size, cur_site_scale_fp);
	} // genotype file reading loop.
	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
} // client_update_pvalue_scale_stats_per_block

void client_pool_pvalue_stats(int client_i, int n_clients, int var_block_size, char* shared_working_dir)
{
	fprintf(stderr, "Client %d: Pooling p-value stats..\n", client_i);

	double** pooled_Gt_Y_min_Y0 = allocate_matrix(var_block_size, 1);
	double** pooled_chi_sqr_scale_vector = allocate_matrix(var_block_size, 1);
	double** cur_block_chi_sqr_norm_stat_vector = allocate_matrix(var_block_size, 1);

	fprintf(stderr, "Processing genotypes by %d-variant blocks..\n", var_block_size);

	char p_val_stats_fp[1000];
	get_pval_stat_op_fp(shared_working_dir, client_i, p_val_stats_fp);
	FILE* f_p_val_stats = open_f(p_val_stats_fp, "w");

	for (int var_block_start_i = 0; ; var_block_start_i += var_block_size)
	{
		int var_block_end_i = var_block_start_i + var_block_size - 1;
		fprintf(stderr, "Pooling stats for block[%d-%d]\n", var_block_start_i, var_block_end_i);
		char cur_block_Gt_Ymin_Y0_fp[1000];
		char var_block_ids_fp[1000];
		get_chisq_stat_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_Gt_Ymin_Y0_fp, var_block_ids_fp);
		if (!check_file(cur_block_Gt_Ymin_Y0_fp))
		{
			fprintf(stderr, "Could not find block %d-%d\n", var_block_start_i, var_block_end_i);
			break;
		}

		vector<char*>* var_block_ids = buffer_file(var_block_ids_fp);

		// Go over all clients.
		set_matrix_val(pooled_chi_sqr_scale_vector, var_block_size, 1, 0);
		set_matrix_val(pooled_Gt_Y_min_Y0, var_block_size, 1, 0);
		for (int cur_client_i = 0; cur_client_i < n_clients; cur_client_i++)
		{
			char cur_block_Gt_Ymin_Y0_fp[1000];
			char var_block_ids_buffer[1000];
			get_chisq_stat_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, cur_client_i, cur_block_Gt_Ymin_Y0_fp, var_block_ids_buffer);

			int loaded_nrow, loaded_ncol;
			double** cur_block_Gt_Y_min_Y0 = load_matrix_binary(cur_block_Gt_Ymin_Y0_fp, loaded_nrow, loaded_ncol);
			if (var_block_size != loaded_nrow ||
				loaded_ncol != 1)
			{
				fprintf(stderr, "Sanity check failed.\n");
			}

			// Pool the stat matrix.
			matrix_add(cur_block_Gt_Y_min_Y0, var_block_size, 1, pooled_Gt_Y_min_Y0, var_block_size, 1, pooled_Gt_Y_min_Y0);

			// Load the scale values for the current block.
			char cur_block_chisq_scale_fp[1000];
			get_scale_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, cur_client_i, cur_block_chisq_scale_fp);
			double** cur_block_chi_sqr_scale_matrix = load_matrix_binary(cur_block_chisq_scale_fp, loaded_nrow, loaded_ncol);
			if (var_block_size != loaded_nrow ||
				loaded_ncol != var_block_size)
			{
				fprintf(stderr, "Sanity check failed.\n");
			}

			// Get the diagonal entry of this matrix -- must prove that this is what we need here.
			double** cur_block_chi_sqr_scale_vector = get_diagonal_of_matrix(cur_block_chi_sqr_scale_matrix, var_block_size, var_block_size, NULL);

			// Pool the scale matrix.
			matrix_add(cur_block_chi_sqr_scale_vector, var_block_size, 1, pooled_chi_sqr_scale_vector, var_block_size, 1, pooled_chi_sqr_scale_vector);
		} // cur_client_i loop.

		double** cur_block_chi_sqr_stat = matrix_multiply_elementwise(pooled_Gt_Y_min_Y0, var_block_size, 1, pooled_Gt_Y_min_Y0, var_block_size, 1, NULL);

		// Get the p-values.
		matrix_divide_elementwise(cur_block_chi_sqr_stat, var_block_size, 1, pooled_chi_sqr_scale_vector, var_block_size, 1, cur_block_chi_sqr_norm_stat_vector);
		double** cur_var_block_chi_sqr_pvals = process_matrix_elementwise_by_callback(cur_block_chi_sqr_norm_stat_vector, var_block_size, 1, get_1DOF_chisqr_pval, NULL);
		//double pval = chisqr(1, cur_chi_sqr_stat / chisqr_scale);

		for (int i_var = 0; i_var < (int)(var_block_ids->size()); i_var++)
		{
			fprintf(f_p_val_stats, "%s\t%.5f\t%.5f\t%.5f\t%.5f\n", var_block_ids->at(i_var), cur_var_block_chi_sqr_pvals[i_var][0],
				cur_block_chi_sqr_norm_stat_vector[i_var][0], cur_block_chi_sqr_stat[i_var][0], pooled_chi_sqr_scale_vector[i_var][0]);
		} // i_var loop.		
	} // genotype file reading loop.

	close_f(f_p_val_stats, p_val_stats_fp);
} // client_pool_pvalue_stats

