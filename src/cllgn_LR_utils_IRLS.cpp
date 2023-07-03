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
#include "cllgn_LR_model_stats.h"
#include "cllgn_LR_utils.h"
#include "cllgn_matrix_linalg_utils.h"

#include "cllgn_x_chisqr.h"
#include "cllgn_x_sigmoid.h"

/*
	- This is the SGD implementation of LR for testing with genotype data.
	. https://github.com/OpenMined/TenSEAL/blob/main/tutorials/Tutorial%201%20-%20Training%20and%20Evaluation%20of%20Logistic%20Regression%20on%20Encrypted%20Data.ipynb
	. https://cseweb.ucsd.edu/~elkan/250B/logreg.pdf

	- TODOs::
	. 1-bit SGD good for parallelization?
	.
*/

using namespace std;

// Uses the IRLS iterations here:: https://towardsdatascience.com/iterated-reweighted-least-squares-and-glms-explained-9c0cc0063526
void train_modified_LR_per_IRLS(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* op_dir)
{
	// We need the two sigmoid functions here.
	//double(*GET_SIGMOID_PER_FEATS_WEIGHTS)(double*, int, double*);
	double(*GET_SIGMOID_PER_FEAT_COMB)(double);
	//double(*GET_LOGIT_PER_WEIGHTS)(double*, int, double*);
	//double(*GET_LOGIT_PER_P)(double);

	//GET_LOGIT_PER_P = get_logit;

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

	// this is the threshold at which the training is stopped early.
	//double LL_EPSILON = 0.01;

	// We need the RNG for different tasks.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	double* per_feat_weights = NULL;
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
				fprintf(stderr, "\t%d: %.4f", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
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
op_dir: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, n_epoch, op_dir, LL_EPSILON);

	// Save the pheno and model data.
	save_LR_data_row_samples(per_ind_feat_matrix, per_ind_pheno, per_feat_weights, n_covars, op_dir);

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

	//// Allocate the weights for the current epoch.
	//double* cur_epoch_weights = new double[n_feats + 2];
	//memset(cur_epoch_weights, 0, sizeof(double) * (n_feats + 2));

	//// Allocate and initialize the nu and mu vectors.
	//double* mu = new double[sample_size];
	//memset(mu, 0, sizeof(double) * sample_size);
	//double* nu = new double[sample_size];
	//memset(nu, 0, sizeof(double) * sample_size);

	double** cur_epoch_weights = allocate_matrix(n_feats, 1);
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

	double cur_LL = 0;

	// Copy the features matrix and allocat the related matrices.
	double** X = get_matrix_per_row_vector_list(per_ind_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);
	double** XtW = allocate_matrix(n_feats, sample_size);
	double** XtWX = allocate_matrix(n_feats, n_feats);
	double** XtWX_inv = allocate_matrix(n_feats, n_feats);
	double** XtWz = allocate_matrix(n_feats, 1);

	double max_abs_nu = 0;

	FILE* f_per_epoch_info = open_f("per_epoch_info.txt", "w");
	for (int i_iter = 0; i_iter < n_epoch; i_iter++)
	{
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

		invert_matrix_GJ(XtWX, n_feats, n_feats, XtWX_inv);

		// Test the inversion step.
		double** inv_tester = matrix_multiply(XtWX_inv, n_feats, n_feats,
												XtWX, n_feats, n_feats, NULL);
		fprintf(stderr, "Inversion test:\n");
		print_matrix(inv_tester, n_feats, n_feats);

		matrix_multiply(XtW, n_feats, sample_size, z, sample_size, 1, XtWz);

		//// Update current beta: This is n_featsx1 matrix.
		//double** new_beta_matrix = matrix_multiply(XtWX_inv, n_feats, n_feats,
		//	XtWz, n_feats, 1);
		matrix_multiply(XtWX_inv, n_feats, n_feats, XtWz, n_feats, 1, cur_epoch_weights);

		//// Copy the beta to the weights vector.
		//for (int i_feat = 0; i_feat < n_feats; i_feat++)
		//{
		//	cur_epoch_weights[i_feat] = new_beta_matrix[i_feat][0];
		//} // i_feat loop.

		// Get the largest nu, this is necessary for sigmoid approximations.
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			double cur_subject_nu = nu[i_s][0];
			if (max_abs_nu < fabs(cur_subject_nu))
			{
				max_abs_nu = fabs(cur_subject_nu);
			}
		} // i_s loop.
		//	mu[i_s] = GET_SIGMOID_PER_FEATS_WEIGHTS(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);

		//	//double test_mu = GET_SIGMOID_PER_FEAT_COMB(cur_subject_nu);
		//	//fprintf(stderr, "%.5f // %.5f\n", test_mu, mu[i_s]);
		//} // i_S loop.

		//double cur_epoch_updated_LL = get_log_likelihood(per_ind_feat_matrix, per_ind_pheno, cur_epoch_weights, n_feats);
		double cur_epoch_updated_LL = get_log_likelihood_per_mu_matrix(mu, per_ind_pheno, sample_size);

		fprintf(f_per_epoch_info, "%d\t%.10f\n", i_iter, cur_epoch_updated_LL);

		if (fabs(cur_LL - cur_epoch_updated_LL) < LL_EPSILON)
		{
			break;
		}

		cur_LL = cur_epoch_updated_LL;
	} // iteration loop.
	fclose(f_per_epoch_info);

	fprintf(stderr, "Max nu: %.5f\n", max_abs_nu);

	// Calculate the phenotype expectations for each subject.
	FILE* f_per_subject_pheno = open_f("per_subject_phenotypes.txt", "w");
	fprintf(f_per_subject_pheno, "SAMPLE_ID\tPHENOTYPE\n");
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		//double cur_pred_pheno = GET_SIGMOID_PER_FEATS_WEIGHTS(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);
		double cur_pred_pheno = mu[i_s][0];
		fprintf(f_per_subject_pheno, "%d\t%.5f\n", i_s, cur_pred_pheno);
	} // i_s loop.

	FILE* f_trained_model = open_f("trained_model.txt", "w");
	fprintf(f_trained_model, "N_COVARS\t%d\n", n_covars);
	int covar_i = 0;
	fprintf(stderr, "Real/Trained weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		//fprintf(stderr, "%.4f\t%.4f\n", per_feat_weights[i_feat], cur_epoch_weights[i_feat]);
		fprintf(stderr, "Feat %d: %.4f\n", i_feat, cur_epoch_weights[i_feat][0]);

		// Write the current feature weight.
		if (i_feat == GET_INTERCEPT_FEAT_I())
		{
			fprintf(f_trained_model, "INTERCEPT\t%lf\n", cur_epoch_weights[i_feat][0]);
		}
		else
		{
			fprintf(f_trained_model, "COVAR_%d\t%lf\n", covar_i, cur_epoch_weights[i_feat][0]);
			covar_i++;
		}
	} // i_feat loop.
	close_f(f_trained_model, "trained_model.txt");
} // train_modified_LR_per_IRLS function.

// Uses the IRLS iterations here:: https://towardsdatascience.com/iterated-reweighted-least-squares-and-glms-explained-9c0cc0063526
void train_baseline_LR_per_IRLS(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp, 
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* op_dir)
{
	// We need the two sigmoid functions here.
	double(*GET_SIGMOID_PER_FEATS_WEIGHTS)(double*, int, double*);
	//double(*GET_SIGMOID_PER_FEAT_COMB)(double);
	//double(*GET_LOGIT_PER_WEIGHTS)(double*, int, double*);
	double(*GET_LOGIT_PER_P)(double);

	GET_LOGIT_PER_P = get_logit;

	if (sigmoid_approx_type == SIGMOID_APPROX_NATIVE)
	{
		fprintf(stderr, "Using native sigmoid..\n");
		GET_SIGMOID_PER_FEATS_WEIGHTS = get_sigmoid_val_per_feats_weights;
		//GET_SIGMOID_PER_FEAT_COMB = get_sigmoid_val_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_KIM_ETAL)
	{
		fprintf(stderr, "Using Kim. etal sigmoid..\n");
		GET_SIGMOID_PER_FEATS_WEIGHTS = get_Kim_etal_poly_approx_sigmoid_per_feats_weights;
		//GET_SIGMOID_PER_FEAT_COMB = get_Kim_etal_poly_approx_sigmoid_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_TENSEAL)
	{
		fprintf(stderr, "Using TenSEAL sigmoid..\n");
		GET_SIGMOID_PER_FEATS_WEIGHTS = get_TenSeal_poly_approx_sigmoid_per_feats_weights;
		//GET_SIGMOID_PER_FEAT_COMB = get_TenSeal_poly_approx_sigmoid_per_feat_comb;
	}
	else
	{
		fprintf(stderr, "Unknown sigmoid approximation: %d\n", sigmoid_approx_type);
		exit(0);
	}

	// this is the threshold at which the training is stopped early.
	//double LL_EPSILON = 0.01;

	// We need the RNG for different tasks.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	double* per_feat_weights = NULL;
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
				fprintf(stderr, "\t%d: %.4f", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
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
op_dir: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, n_epoch, op_dir, LL_EPSILON);

	// Save the pheno and model data.
	save_LR_data_row_samples(per_ind_feat_matrix, per_ind_pheno, per_feat_weights, n_covars, op_dir);

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

	//// Write the model.
	//fprintf(stderr, "Model weights:\n");
	//for (int i_feat = 0; i_feat < n_feats; i_feat++)
	//{
	//	fprintf(stderr, "Feat %d: %.5f\n", i_feat, per_feat_weights[i_feat]);
	//} // i_feat loop.

	// Allocate the weights for the current epoch.
	double* cur_epoch_weights = new double[n_feats + 2];
	memset(cur_epoch_weights, 0, sizeof(double) * (n_feats + 2));

	// Mu vector.
	double* mu = new double[sample_size];
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		mu[i_s] = 0.5;
	} // i_s loop.

	double cur_LL = 0;

	// Copy the features matrix.
	double** X = get_matrix_per_row_vector_list(per_ind_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);

	double max_abs_nu = 0;

	FILE* f_per_epoch_info = open_f("per_epoch_info.txt", "w");
	for(int i_iter = 0; i_iter < n_epoch; i_iter++)
	{
		// This is necessary for 
		double* z_vec = new double[sample_size];
		double* W_diag = new double[sample_size];
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			// mu[i_s]=sigmoid(nu); thus, we actually do not need the logit computation here. We can simply compute nu.
			z_vec[i_s] = GET_LOGIT_PER_P(mu[i_s]) + (per_ind_pheno[i_s] - mu[i_s]) / (mu[i_s] * (1 - mu[i_s]));

			W_diag[i_s] = (mu[i_s] * (1 - mu[i_s]));
		} // i_s loop.

		double** z = get_matrix_per_col_vector(z_vec, sample_size, NULL);
		//fprintf(stderr, "Matrix: z:\n");
		//print_matrix(z, sample_size, 1);

		double** W = get_diag_matrix(W_diag, sample_size, NULL);

		// Calculate X'WX and X'Wz
		double** XtW = matrix_multiply(Xt, n_feats, sample_size,
										W, sample_size, sample_size, NULL);

		double** XtWX = matrix_multiply(XtW, n_feats, sample_size,
										X, sample_size, n_feats, NULL);

		double** XtWX_inv = invert_matrix_GJ(XtWX, n_feats, n_feats, NULL);

		// Test the inversion step.
		double** inv_tester = matrix_multiply(XtWX_inv, n_feats, n_feats, 
												XtWX, n_feats, n_feats, NULL);
		fprintf(stderr, "Inversion test:\n");
		print_matrix(inv_tester, n_feats, n_feats);

		double** XtWz = matrix_multiply(XtW, n_feats, sample_size,
										z, sample_size, 1, NULL);

		// Update current beta: This is n_featsx1 matrix.
		double** new_beta_matrix = matrix_multiply(XtWX_inv, n_feats, n_feats,
													XtWz, n_feats, 1, NULL);

		// Copy the beta to the weights vector.
		for (int i_feat = 0; i_feat < n_feats; i_feat++)
		{
			cur_epoch_weights[i_feat] = new_beta_matrix[i_feat][0];
		} // i_feat loop.

		// Update mu.
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			double cur_subject_nu = get_vec2vec_inner_product(per_ind_feat_matrix->at(i_s), cur_epoch_weights, n_feats);
			if (max_abs_nu < fabs(cur_subject_nu))
			{
				max_abs_nu = fabs(cur_subject_nu);
			}

			mu[i_s] = GET_SIGMOID_PER_FEATS_WEIGHTS(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);

			//double test_mu = GET_SIGMOID_PER_FEAT_COMB(cur_subject_nu);
			//fprintf(stderr, "%.5f // %.5f\n", test_mu, mu[i_s]);
		} // i_S loop.

		double cur_epoch_updated_LL = get_log_likelihood(per_ind_feat_matrix, per_ind_pheno, cur_epoch_weights, n_feats);

		fprintf(f_per_epoch_info, "%d\t%.10f\n", i_iter, cur_epoch_updated_LL);

		if (fabs(cur_LL - cur_epoch_updated_LL) < LL_EPSILON)
		{
			break;
		}

		cur_LL = cur_epoch_updated_LL;
	} // iteration loop.
	fclose(f_per_epoch_info);

	fprintf(stderr, "Max nu: %.5f\n", max_abs_nu);

	// Calculate the phenotype expectations for each subject.
	FILE* f_per_subject_pheno = open_f("per_subject_phenotypes.txt", "w");
	fprintf(f_per_subject_pheno, "SAMPLE_ID\tPHENOTYPE\n");
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_pred_pheno = GET_SIGMOID_PER_FEATS_WEIGHTS(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);
		fprintf(f_per_subject_pheno, "%d\t%.5f\n", i_s, cur_pred_pheno);
	} // i_s loop.

	FILE* f_trained_model = open_f("trained_model.txt", "w");
	fprintf(f_trained_model, "N_COVARS\t%d\n", n_covars);
	int covar_i = 0;
	fprintf(stderr, "Real/Trained weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		//fprintf(stderr, "%.4f\t%.4f\n", per_feat_weights[i_feat], cur_epoch_weights[i_feat]);
		fprintf(stderr, "Feat %d: %.4f\n", i_feat, cur_epoch_weights[i_feat]);

		// Write the current feature weight.
		if (i_feat == GET_INTERCEPT_FEAT_I())
		{
			fprintf(f_trained_model, "INTERCEPT\t%lf\n", cur_epoch_weights[i_feat]);
		}
		else
		{
			fprintf(f_trained_model, "COVAR_%d\t%lf\n", covar_i, cur_epoch_weights[i_feat]);
			covar_i++;
		}
	} // i_feat loop.
	close_f(f_trained_model, "trained_model.txt");
} // train_modified_LR_per_IRLS function.

void assign_ChiSqr_p_values_per_fit_model_pheno_per_GMMAT_text_geno(char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_fp,
	char* op_fp)
{
	fprintf(stderr, "Computing chi-squared p-values of variants in %s using null model phenotypes in %s and saving to %s\n", GMMAT_text_genotype_matrix_fp, null_model_pheno_fp, op_fp);

	fprintf(stderr, "Loading the features data from %s.\n", subject_per_row_feats_fp);
	int n_covars_per_feats_file = 0;
	vector<double*>* per_subject_feat_matrix = new vector<double*>();
	vector<double>* per_subject_obs_pheno = new vector<double>();
	load_LR_data_row_samples(subject_per_row_feats_fp, obs_pheno_fp, n_covars_per_feats_file, per_subject_feat_matrix, per_subject_obs_pheno);

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars_per_feats_file);
	fprintf(stderr, "%d features loaded from covariates file.\n", n_feats);

	// Laod the null model phenotypes.
	vector<char*>* null_model_pheno_lines = buffer_file(null_model_pheno_fp);
	vector<double>* per_subject_null_pheno = new vector<double>();

	// MAke sre to skip the header.
	for (int i_l = 1; i_l < (int)(null_model_pheno_lines->size()); i_l++)
	{
		double cur_subj_pheno;
		sscanf(null_model_pheno_lines->at(i_l), "%*s %lf", &cur_subj_pheno);
		per_subject_null_pheno->push_back(cur_subj_pheno);
	} // i_l loop.
	fprintf(stderr, "Loaded %d null phenotypes.\n", (int)(per_subject_null_pheno->size()));

	if (per_subject_null_pheno->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between observed and null phenotypes: %d, %d\n", (int)(per_subject_null_pheno->size()), (int)(per_subject_obs_pheno->size()));
		exit(0);
	}

	if (per_subject_feat_matrix->size() != per_subject_null_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)(per_subject_null_pheno->size()), (int)(per_subject_feat_matrix->size()));
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\n", i_s, per_subject_null_pheno->at(i_s) - per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	int sample_size = (int)(per_subject_feat_matrix->size());

	// Compute the projection matrix: 
	fprintf(stderr, "Setting Sigma inv..\n");
	double* sigma_inv_vec = new double[sample_size];
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		sigma_inv_vec[i_s] = per_subject_null_pheno->at(i_s) * (1 - per_subject_null_pheno->at(i_s));
	} // i_s loop.

	// Get the sigma inverse matrix.
	double** sigma_inv = get_diag_matrix(sigma_inv_vec, sample_size, NULL);
	double** X = get_matrix_per_row_vector_list(per_subject_feat_matrix, n_feats, NULL);
	double** Xt = transpose_matrix(X, sample_size, n_feats, NULL);

	double** Xtsigma_inv = matrix_multiply(Xt, n_feats, sample_size, sigma_inv, sample_size, sample_size, NULL);

	// Transpose X'Sigma-1 = Sigma-1X
	double** sigmainv_X = transpose_matrix(Xtsigma_inv, n_feats, sample_size, NULL);

	double** Xtsigma_invX = matrix_multiply(Xtsigma_inv, n_feats, sample_size, X, sample_size, n_feats, NULL);

	double** inv_term = invert_matrix_GJ(Xtsigma_invX, n_feats, n_feats, NULL);

	fprintf(stderr, "Inversion test for XtS-1X:\n");
	print_matrix(matrix_multiply(inv_term, n_feats, n_feats, Xtsigma_invX, n_feats, n_feats, NULL), n_feats, n_feats);

	double** right1 = matrix_multiply(inv_term, n_feats, n_feats, Xtsigma_inv, n_feats, sample_size, NULL);

	double** right = matrix_multiply(sigmainv_X, sample_size, n_feats, right1, n_feats, sample_size, NULL);

	double** P = matrix_subtract(sigma_inv, sample_size, sample_size, right, sample_size, sample_size, NULL);

	fprintf(stderr, "Processing genotypes..\n");
	FILE* f_geno = open_f(GMMAT_text_genotype_matrix_fp, "r");
	FILE* f_op = open_f(op_fp, "w");
	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size,1);
	double** Gt = allocate_matrix(1, sample_size);
	while (1)
	{
		char* cur_line = getline(f_geno);
		if (cur_line == NULL)
		{
			break;
		}

		n_processed_vars++;
		if (n_processed_vars % 100 == 0)
		{
			fprintf(stderr, "@ %d. variant.           \r", n_processed_vars);
		}

		t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, "\t");
		vector<char*>* cur_line_tok_strs = t_string::copy_tokens_2_strs(cur_line_toks);
		char* var_id = cur_line_tok_strs->at(0);
		//char* ref_all = cur_line_tok_strs->at(1);
		//char* alt_all = cur_line_tok_strs->at(2);

		if ((int)(cur_line_tok_strs->size()) != (sample_size + 3))
		{
			fprintf(stderr, "Genotype sample size is not as expected: %d, %d\n",
				(int)(cur_line_tok_strs->size()), (int)(per_subject_null_pheno->size()));
			exit(0);
		}

		double cur_chi_stat = 0;

		for (int i_t = 3; i_t < (int)(cur_line_tok_strs->size()); i_t++)
		{
			int i_s = i_t - 3;

			double cur_geno = atof(cur_line_tok_strs->at(i_t));
			if (cur_geno != 0 &&
				cur_geno != 1 &&
				cur_geno != 2)
			{
				fprintf(stderr, "Found invalid geno: %.0f for %s\n", cur_geno, var_id);
				exit(0);
			}

			G[i_s][0] = cur_geno;
			Gt[0][i_s] = cur_geno;

			double cur_sqi_sqr_cont = (cur_geno * (per_subject_obs_pheno->at(i_s) - per_subject_null_pheno->at(i_s)));
			cur_chi_stat += cur_sqi_sqr_cont;
		} // i_t loop.

		// Take square.
		double cur_chi_sqr_stat = pow(cur_chi_stat, 2);

		double** first_mul = matrix_multiply(Gt, 1, sample_size, P, sample_size, sample_size, NULL);
		double** second_mul = matrix_multiply(first_mul, 1, sample_size, G, sample_size, 1, NULL);
		double chisqr_scale = second_mul[0][0];

		double pval = chisqr(1, cur_chi_sqr_stat / chisqr_scale);

		fprintf(f_op, "%s\t%.5f\t%.5f\t%.5f\t%.5f\n", var_id, pval, cur_chi_sqr_stat / chisqr_scale, cur_chi_sqr_stat, chisqr_scale);

		t_string::clean_string_list(cur_line_tok_strs);
		t_string::clean_tokens(cur_line_toks);
		delete[] cur_line;
	} // genotype file reading loop.
	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
	close_f(f_op, op_fp);
}


