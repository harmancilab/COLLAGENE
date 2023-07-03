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
#include "cllgn_matrix_linalg_utils.h"

bool __DUMP_FEDLR_META_TYPE_MESSAGES__ = false;

using namespace std;

void client_calculate_save_pvalue_stats_for_meta_analysis(int client_i, int n_clients, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
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
	if (loaded_nrow != (int)per_subject_obs_pheno->size() ||
		loaded_ncol != 1)
	{
		fprintf(stderr, "null phenotype matrix dimensions non-conformant: %d, %d\n", loaded_nrow, loaded_ncol);
		exit(0);
	}

	if (per_subject_feat_matrix->size() != per_subject_obs_pheno->size())
	{
		fprintf(stderr, "The sample size does not match between null phenotypes and feature matrix: %d, %d\n", (int)per_subject_obs_pheno->size(), (int)per_subject_feat_matrix->size());
		exit(0);
	}

	for (int i_s = 0; i_s < 10; i_s++)
	{
		fprintf(stderr, "Subject %d: %.4f\n", i_s, null_pheno_matrix[i_s][0] - per_subject_obs_pheno->at(i_s));
	} // i_s loop.

	int sample_size = (int)per_subject_feat_matrix->size();

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

	// Transpose X'Sigma-1 = Sigma-1X
	double** sigmainv_X = transpose_matrix(Xtsigma_inv, n_feats, sample_size, NULL);

	// This matrix should be the pooled matrix from all sites, we are current using the site-specific Xtsigma_invX matrix, which is not the correct one.
	// This can be done by simply saving this matrix while training. This is a pxp matrix.
	double** Xtsigma_invX = matrix_multiply(Xtsigma_inv, n_feats, sample_size, X, sample_size, n_feats, NULL);

	double** inv_term = invert_matrix_GJ(Xtsigma_invX, n_feats, n_feats, NULL);

	fprintf(stderr, "Inversion test for XtS-1X:\n");
	print_matrix(matrix_multiply(inv_term, n_feats, n_feats, Xtsigma_invX, n_feats, n_feats, NULL), n_feats, n_feats);

	double** right1 = matrix_multiply(inv_term, n_feats, n_feats, Xtsigma_inv, n_feats, sample_size, NULL);

	double** right = matrix_multiply(sigmainv_X, sample_size, n_feats, right1, n_feats, sample_size, NULL);

	// This is necessary for estimating chi-squared scale parameter;; Can we not compute it here but below to make multiplications more efficient?
	// We need G'PG = G' (sigma_inv - sigmainvX inv_term X'sigmainv) G
	// sigmainv is a diagonal matrix. Multiplication only multiplies rows / columns on left / right; this does not require n^3.
	// Just implement matrix multiplication routines into custom processing step of the SEAL_Genomics and test with different sizes of inputs.
	// 
	double** P = matrix_subtract(sigma_inv, sample_size, sample_size, right, sample_size, sample_size, NULL);

	int n_processed_vars = 0;
	double** G = allocate_matrix(sample_size, var_block_size);
	double** Gt = allocate_matrix(var_block_size, sample_size);
	double** first_mul = allocate_matrix(var_block_size, sample_size);
	//double** chisqr_scale = allocate_matrix(var_block_size, 1);

	// Obs and Null Pheno.
	double** Y = get_matrix_per_col_vector(per_subject_obs_pheno, NULL);
	//double** Yt = transpose_matrix(Y, sample_size, 1, NULL);

	double** Y0 = null_pheno_matrix;
	//double** Y0t = transpose_matrix(Y0, sample_size, 1, NULL);

	double** Y_min_Y0 = matrix_subtract(Y, sample_size, 1, Y0, sample_size, 1, NULL);

	double** Gt_Y_min_Y0 = allocate_matrix(var_block_size, 1);
	double** cur_block_chi_sqr_scale_matrix = allocate_matrix(var_block_size, var_block_size);
	double** cur_block_chi_sqr_scale_vector = allocate_matrix(var_block_size, 1);
	//double** cur_block_chi_sqr_norm_stat_vector = allocate_matrix(var_block_size, 1);

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

				//double cur_chi_sqr_cont = (cur_geno * (per_subject_obs_pheno->at(i_s) - per_subject_null_pheno->at(i_s)));
				//cur_chi_stat += cur_chi_sqr_cont;
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

		// Compute scale for the current block.
		matrix_multiply(Gt, var_block_size, sample_size, P, sample_size, sample_size, first_mul);
		matrix_multiply(first_mul, var_block_size, sample_size, G, sample_size, var_block_size, cur_block_chi_sqr_scale_matrix);
		get_diagonal_of_matrix(cur_block_chi_sqr_scale_matrix, var_block_size, var_block_size, cur_block_chi_sqr_scale_vector);

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

		char cur_block_chisq_scale_fp[1000];
		get_scale_matrix_fp(shared_working_dir, var_block_start_i, var_block_end_i, client_i, cur_block_chisq_scale_fp);
		save_matrix_binary(cur_block_chi_sqr_scale_vector, var_block_size, 1, cur_block_chisq_scale_fp);

		char chisq_scale_matrix_txt_fp[1000];
		sprintf(chisq_scale_matrix_txt_fp, "chisq_scale_site_%d.txt", client_i);
		save_matrix_plain(cur_block_chi_sqr_scale_vector, var_block_size, 1, chisq_scale_matrix_txt_fp);

		//matrix_divide_elementwise(cur_block_chi_sqr_stat, var_block_size, 1, cur_block_chi_sqr_scale_vector, var_block_size, 1, cur_block_chi_sqr_norm_stat_vector);
		//double** cur_var_block_chi_sqr_pvals = process_matrix_elementwise_by_callback(cur_block_chi_sqr_norm_stat_vector, var_block_size, 1, get_1DOF_chisqr_pval, NULL);
		//double pval = chisqr(1, cur_chi_sqr_stat / chisqr_scale);

		//fprintf(f_op, "%s\t%.5f\t%.5f\t%.5f\t%.5f\n", var_id, pval, cur_chi_sqr_stat / chisqr_scale, cur_chi_sqr_stat, chisqr_scale);		
	} // genotype file reading loop.
	close_f(f_geno, GMMAT_text_genotype_matrix_fp);
} // client_calculate_save_pvalue_stats_for_meta_analysis option.

void client_pool_pvalue_stats_for_meta_analysis(int client_i, int n_clients, int var_block_size, char* shared_working_dir)
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
			double** cur_block_chi_sqr_scale_vector = load_matrix_binary(cur_block_chisq_scale_fp, loaded_nrow, loaded_ncol);
			if (var_block_size != loaded_nrow ||
				loaded_ncol != 1)
			{
				fprintf(stderr, "Sanity check failed.\n");
			}

			// Pool the scale matrix.
			matrix_add(cur_block_chi_sqr_scale_vector, var_block_size, 1, pooled_chi_sqr_scale_vector, var_block_size, 1, pooled_chi_sqr_scale_vector);
		} // cur_client_i loop.

		double** cur_block_chi_sqr_stat = matrix_multiply_elementwise(pooled_Gt_Y_min_Y0, var_block_size, 1, pooled_Gt_Y_min_Y0, var_block_size, 1, NULL);

		// Get the p-values.
		matrix_divide_elementwise(cur_block_chi_sqr_stat, var_block_size, 1, pooled_chi_sqr_scale_vector, var_block_size, 1, cur_block_chi_sqr_norm_stat_vector);
		double** cur_var_block_chi_sqr_pvals = process_matrix_elementwise_by_callback(cur_block_chi_sqr_norm_stat_vector, var_block_size, 1, get_1DOF_chisqr_pval, NULL);
		//double pval = chisqr(1, cur_chi_sqr_stat / chisqr_scale);

		for (int i_var = 0; i_var < (int)var_block_ids->size(); i_var++)
		{
			fprintf(f_p_val_stats, "%s\t%.5f\t%.5f\t%.5f\t%.5f\n", var_block_ids->at(i_var), cur_var_block_chi_sqr_pvals[i_var][0],
				cur_block_chi_sqr_norm_stat_vector[i_var][0], cur_block_chi_sqr_stat[i_var][0], pooled_chi_sqr_scale_vector[i_var][0]);
		} // i_var loop.		
	} // genotype file reading loop.

	close_f(f_p_val_stats, p_val_stats_fp);
} // client_pool_pvalue_stats
