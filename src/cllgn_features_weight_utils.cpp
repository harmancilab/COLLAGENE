#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_ansi_string.h"
#include "cllgn_file_utils.h"
#include "cllgn_histogram.h"
#include <math.h>
#include <string.h>
#include "cllgn_LR_model_stats.h"
#include "cllgn_x_sigmoid.h"
#include "cllgn_x_chisqr.h"
#include "cllgn_features_weight_utils.h"

void generate_save_random_feats_model(int sample_size, int n_covars, double pheno_noise_std_dev, char* op_dir)
{
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	fprintf(stderr, "Generating model.\n");
	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);
	double* per_feat_weights = generate_per_feat_weights(rng, n_feats);

	fprintf(stderr, "Generating the training data for %d sample and %d features.\n", sample_size, n_feats);

	n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

	// We first generate a random feature matrix.
	fprintf(stderr, "Generating feature matrix.\n");
	vector<double*>* per_ind_feat_matrix = new vector<double*>();
	double* per_ind_pheno = new double[sample_size + 2];
	generate_LR_geno_covar_pheno_data(rng, per_ind_feat_matrix, per_ind_pheno, per_feat_weights, pheno_noise_std_dev, sample_size, n_covars);

	fprintf(stderr, "Saving features and model.\n");
	save_LR_data_row_samples(per_ind_feat_matrix, per_ind_pheno, per_feat_weights, n_covars, op_dir);
} // generate_save_random_feats_model option.

double* load_LR_model_weights(char* model_fp, int& model_loaded_n_covars)
{
	fprintf(stderr, "Loading model weights from %s..\n", model_fp);
	vector<char*>* model_file_lines = buffer_file(model_fp);

	int n_covars_from_file;
	char n_covars_col_str[100];
	sscanf(model_file_lines->at(0), "%s %d", n_covars_col_str, &n_covars_from_file);

	if (!t_string::compare_strings(n_covars_col_str, "N_COVARS"))
	{
		fprintf(stderr, "Model file does not start as expected: %s\n", model_file_lines->at(0));
		exit(0);
	}

	model_loaded_n_covars = n_covars_from_file;
	int n_feats = GET_N_FEATS_PER_N_COVARS(model_loaded_n_covars);

	double* per_feat_weights = new double[n_feats + 2];
	for (int i_l = 1; i_l < (int)(model_file_lines->size()); i_l++)
	{
		int feat_i = i_l - 1;

		vector<char*>* model_line_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(model_file_lines->at(i_l), "\t"));
		if ((int)(model_line_toks->size()) != 2)
		{
			fprintf(stderr, "Model file line is not as expected: %s\n", model_file_lines->at(i_l));
			exit(0);
		}

		if (feat_i == GET_INTERCEPT_FEAT_I())
		{
			if (!t_string::compare_strings(model_line_toks->at(0), "INTERCEPT"))
			{
				fprintf(stderr, "Intercept column id is not as expected in %s: %s\n", model_fp, model_file_lines->at(i_l));
				exit(0);
			}
		}

		// Copy the weight for this feature.
		per_feat_weights[feat_i] = atof(model_line_toks->at(1));

		fprintf(stderr, "Feat %d: Weight: %.5f\n", feat_i, per_feat_weights[feat_i]);
	} // i_l loop.

	return(per_feat_weights);
} 

void load_LR_data_row_samples(char* LR_data_fp, char* pheno_fp, int& n_covars, vector<double*>* per_ind_feat_matrix, vector<double>* per_ind_pheno_vector)
{
	fprintf(stderr, "Loading LR data with samples in the rows from %s.\n", LR_data_fp);

	// One of the feature names is expected to be we need to have INTERCEPT TERMS.
	vector<char*>* LR_data_lines = buffer_file(LR_data_fp);

	// Parse the header and validate the two columns.
	vector<char*>* header_tok_strs = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(LR_data_lines->at(0), "\t"));
	for (int tok_i = 1; tok_i < (int)(header_tok_strs->size()); tok_i++)
	{
		int feat_i = tok_i - 1;

		if (feat_i == GET_INTERCEPT_FEAT_I())
		{
			if (!t_string::compare_strings(header_tok_strs->at(tok_i), "INTERCEPT"))
			{
				fprintf(stderr, "Intercept column id in the header is not valid: %s\n", header_tok_strs->at(feat_i));
				exit(0);
			}
		}
	} // col_i loop.

	// # of features excludes the subject identifier.
	int n_feats = (int)(header_tok_strs->size()) - 1;
	
	// Set the number of covariates from the data file.
	n_covars = GET_N_COVARS_PER_N_FEATS(n_feats);

	// Start processing the data.
	for (int i_l = 1; i_l < (int)(LR_data_lines->size()); i_l++)
	{
		vector<char*>* cur_sample_feat_tok_strs = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(LR_data_lines->at(i_l), "\t"));

		double* cur_sample_feats = new double[n_feats];
		for (int feat_i = 0; feat_i < n_feats; feat_i++)
		{
			int tok_i = feat_i + 1;
			cur_sample_feats[feat_i] = atof(cur_sample_feat_tok_strs->at(tok_i));
		} // feat_i loop.

		per_ind_feat_matrix->push_back(cur_sample_feats);
	} // i_l loop.

	fprintf(stderr, "Loaded %d subjects' features.\n", (int)(per_ind_feat_matrix->size()));

	// Load the phenotypes.
	if (!check_file(pheno_fp))
	{
		fprintf(stderr, "Could not find phenotype file @ %s\n", pheno_fp);
		exit(0);
	}

	vector<char*>* pheno_lines = buffer_file(pheno_fp);
	vector<char*>* pheno_header_tok_str = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(pheno_lines->at(0), "\t"));
	if ((int)(pheno_header_tok_str->size()) != 2 ||
		!t_string::compare_strings(pheno_header_tok_str->at(1), "PHENOTYPE"))
	{
		fprintf(stderr, "Phenotype file header is not as expected: %s\n", pheno_lines->at(0));
		exit(0);
	}
	t_string::clean_string_list(pheno_header_tok_str);

	fprintf(stderr, "Loading phenotypes from %s..\n", pheno_fp);
	for (int i_l = 1; i_l < (int)(pheno_lines->size()); i_l++)
	{
		vector<char*>* cur_sample_pheno_tok_str = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(pheno_lines->at(i_l), "\t"));

		per_ind_pheno_vector->push_back(atof(cur_sample_pheno_tok_str->at(1)));

		t_string::clean_string_list(cur_sample_pheno_tok_str);
	} // i_l loop.

	fprintf(stderr, "Loaded phenotypes..\n");
} // load_LR_data_row_samples function.

void copy_weights(double* copy_per_feat_weights, double* per_feat_weights, int n_feats)
{
	for (int feat_i = 0; feat_i < n_feats; feat_i++)
	{
		copy_per_feat_weights[feat_i] = per_feat_weights[feat_i];
	} // feat_i loop
}

void save_LR_data_row_samples(vector<double*>* per_ind_feat_matrix, double* per_ind_pheno, double* per_feat_weights, int n_covars, char* op_dir)
{
	fprintf(stderr, "Saving pheno/covar data for subjects and model data to under \"%s/\"\n", op_dir);

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

	// Save each sample on the row.
	if (per_ind_feat_matrix != NULL)
	{
		fprintf(stderr, "Saving the pheno/covar data for %d subjects.\n", (int)(per_ind_feat_matrix->size()));

		char feat_matrix_fp[1000];
		sprintf(feat_matrix_fp, "%s/feat_matrix.txt", op_dir);
		FILE* f_feat_matrix = open_f(feat_matrix_fp, "w");
		if (f_feat_matrix == NULL)
		{
			fprintf(stderr, "Could not open %s for writing.\n", feat_matrix_fp);
			exit(0);
		}

		// Write the header.
		int covar_i = 0;
		fprintf(f_feat_matrix, "SUBJECT_ID");
		for (int i_feat = 0; i_feat < n_feats; i_feat++)
		{
			if (i_feat == GET_INTERCEPT_FEAT_I())
			{
				fprintf(f_feat_matrix, "\tINTERCEPT");
			}
			else
			{
				fprintf(f_feat_matrix, "\tCOVAR_%d", covar_i);
				covar_i++;
			}
		} // i_feat loop.

		// Finish the header line.
		fprintf(f_feat_matrix, "\n");

		// Write the data per subject.
		for (int i_s = 0; i_s < (int)(per_ind_feat_matrix->size()); i_s++)
		{
			fprintf(f_feat_matrix, "SUBJECT_%d", i_s);

			for (int i_feat = 0; i_feat < n_feats; i_feat++)
			{
				fprintf(f_feat_matrix, "\t%.4f", per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(f_feat_matrix, "\n");
		} // i_s loop.

		close_f(f_feat_matrix, feat_matrix_fp);
	} // data saving check.

	if (per_ind_pheno != NULL)
	{
		char pheno_matrix_fp[1000];
		sprintf(pheno_matrix_fp, "%s/phenotypes.txt", op_dir);
		FILE* f_pheno_matrix = open_f(pheno_matrix_fp, "w");
		if (f_pheno_matrix == NULL)
		{
			fprintf(stderr, "Could not open %s for writing.\n", pheno_matrix_fp);
			exit(0);
		}

		fprintf(f_pheno_matrix, "SUBJECT_ID\tPHENOTYPE\n");

		for (int i_s = 0; i_s < (int)(per_ind_feat_matrix->size()); i_s++)
		{
			fprintf(f_pheno_matrix, "SUBJECT_%d\t%.4f\n", i_s, per_ind_pheno[i_s]);
		} // i_s loop.

		close_f(f_pheno_matrix, pheno_matrix_fp);
	} // phenotype data check.

	// Write the feature weights, if they are provided.
	if (per_feat_weights != NULL)
	{
		fprintf(stderr, "Saving model parameters.\n");

		char model_fp[1000];
		sprintf(model_fp, "%s/model.txt", op_dir);
		FILE* f_model = open_f(model_fp, "w");
		if (f_model == NULL)
		{
			fprintf(stderr, "Could not open %s for writing.\n", model_fp);
			exit(0);
		}

		// Write the number of covariates.
		fprintf(f_model, "N_COVARS\t%d\n", n_covars);

		int covar_i = 0;
		for (int i_feat = 0; i_feat < n_feats; i_feat++)
		{
			if (i_feat == GET_INTERCEPT_FEAT_I())
			{
				fprintf(f_model, "INTERCEPT\t%lf\n", per_feat_weights[i_feat]);
			}
			else
			{
				fprintf(f_model, "COVAR_%d\t%lf\n", covar_i, per_feat_weights[i_feat]);
				covar_i++;
			}
		} // i_feat loop.

		close_f(f_model, model_fp);
	} // model saving check.
} // data/model saving function.

void normalize_feats(vector<double*>* per_ind_feature_matrix, int n_covars)
{
	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);
	//int sample_size = (int)(per_ind_feature_matrix->size());

	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		vector<double>* cur_feat_vals = new vector<double>();
		for (int i_s = 0; i_s < (int)(per_ind_feature_matrix->size()); i_s++)
		{
			cur_feat_vals->push_back(per_ind_feature_matrix->at(i_s)[i_feat]);
		} // i_s loop.

		double feat_mean;
		double feat_std_dev;
		get_stats(cur_feat_vals, feat_mean, feat_std_dev);

		for (int i_s = 0; i_s < (int)(per_ind_feature_matrix->size()); i_s++)
		{
			per_ind_feature_matrix->at(i_s)[i_feat] = (per_ind_feature_matrix->at(i_s)[i_feat] - feat_mean) / feat_std_dev;
		} // i_s loop.
	} // i_feat loop.
}

// Random weights for initializing the algorithm.
double* generate_per_feat_weights(t_rng* rng, int n_covars)
{
	// The total number of features is covariate count plus 1: Intercept.
	int n_feats = n_covars + 2;

	double* per_feat_weights = new double[n_feats + 1];

	// Randomly select the weights.
	for (int feat_i = 0; feat_i < n_feats; feat_i++)
	{
		if (feat_i == GET_INTERCEPT_FEAT_I())
		{
			per_feat_weights[feat_i] = rng->random_gaussian_double_ran3() * 0.1;
		}
		else
		{
			per_feat_weights[feat_i] = rng->random_gaussian_double_ran3() * 3;
		}
	} // feat_i loop.

	return(per_feat_weights);
}

void generate_LR_geno_covar_pheno_data(t_rng* rng,
	vector<double*>* per_ind_feat_matrix,
	double* per_ind_pheno,
	double* per_feat_weights,
	double noise_std_dev,
	int sample_size, int n_covars)
{
	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

	fprintf(stderr, "Generating model weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		fprintf(stderr, "Feat %d: %.6f\n", i_feat, per_feat_weights[i_feat]);
	} // i_feat loop.

	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		fprintf(stderr, "Generating data for sample %d\r", i_s);

		// Generate the covariates.
		double* cur_sample_feats = new double[n_feats + 2];

		// Add the intercept to the first position.
		cur_sample_feats[GET_INTERCEPT_FEAT_I()] = 1; // This is the intercept.

		for (int feat_i = 0; feat_i < n_feats; feat_i++)
		{
			// Make sure we do not overwrite intercept.
			if (feat_i != GET_INTERCEPT_FEAT_I())
			{
				//double cur_covar = (rng->random_double_ran3() - 0.5) * (feat_i + 1);
				double cur_covar = (rng->random_double_ran3() - 0.5);

				cur_sample_feats[feat_i] = cur_covar;
			}
		} // feat_i loop.

		// Get the pheno for this sample: Make sure to have the genotype impact here.
		double cur_residual_noise = rng->random_gaussian_double_ran3() * noise_std_dev;

		// Add the noise.
		double cur_subject_feat_comb = cur_residual_noise;
		for (int feat_i = 0; feat_i < n_feats; feat_i++)
		{
			cur_subject_feat_comb += cur_sample_feats[feat_i] * per_feat_weights[feat_i];
		} // feat_i loop.

		//double cur_sample_pheno = (get_sigmoid_val(cur_sample_feats, n_feats, per_feat_weights) > 0.5) ? (1.0) : (0.0);
		double cur_sample_pheno = (get_sigmoid_val_per_feat_comb(cur_subject_feat_comb) > 0.5) ? (1.0) : (0.0);
		per_ind_pheno[i_s] = cur_sample_pheno;

		per_ind_feat_matrix->push_back(cur_sample_feats);
	} // i_s loop.

	fprintf(stderr, "Generated sample data of size %d.\n", (int)(per_ind_feat_matrix->size()));
}