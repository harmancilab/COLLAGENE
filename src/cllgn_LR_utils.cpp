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
#include "cllgn_x_chisqr.h"
#include "cllgn_x_sigmoid.h"
#include "cllgn_LR_model_stats.h"
#include "cllgn_LR_utils.h"

/*
	- This is the SGD implementation of LR for testing with genotype data.
	. https://github.com/OpenMined/TenSEAL/blob/main/tutorials/Tutorial%201%20-%20Training%20and%20Evaluation%20of%20Logistic%20Regression%20on%20Encrypted%20Data.ipynb
	. https://cseweb.ucsd.edu/~elkan/250B/logreg.pdf

	- TODOs::
	. 1-bit SGD good for parallelization?
	.
*/

using namespace std;

// This is the R compatible LR training function -- the baseline model.
void train_LR_baseline(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp, char* model_fp,
	int sample_size, int n_covars,
	double pheno_noise_std_dev,
	double update_step,
	double regularization_weight,
	int n_epoch,
	char* op_dir)
{
	// this is the threshold at which the training is stopped early.
	double LL_EPSILON = 0.01;

	// We need the RNG for different tasks.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	double* per_feat_weights = NULL;
	int n_feats = 0;

	// Load or generate model.
	if (check_file(model_fp))
	{
		int model_loaded_covars = 0;
		per_feat_weights = load_LR_model_weights(model_fp, model_loaded_covars);

		n_covars = model_loaded_covars;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);
	}
	else
	{
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Generate the random model data.
		per_feat_weights = generate_per_feat_weights(rng, n_feats);
		
	} // model loading check.

	// Model loading is finished, at this stage, we should have the n_feats, n_covars, and model setup.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// At this point, n_feats and n_covars must be set.
	if (n_feats == 0 || 
		n_covars == 0)
	{
		fprintf(stderr, "Could not set n_feats after loading/setting the model..\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Loaded the model with %d covariates (%d features) into memory.\n", n_covars, n_feats);
	}

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

		if (n_covars != n_covars_per_feats_file)
		{
			fprintf(stderr, "Could not match the covariate number between model and features file: %d, %d\n", n_covars, n_covars_per_feats_file);
			exit(0);
		}

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[per_ind_feat_matrix->size() + 2];
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
				fprintf(stderr, "\t%d: %.4f", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		sample_size = per_ind_feat_matrix->size();
	} // feats file check.
	else
	{
		fprintf(stderr, "Generating the training data for %d sample and %d features.\n", sample_size, n_feats);

		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);		

		// We first generate a random feature matrix.
		per_ind_feat_matrix = new vector<double*>();
		per_ind_pheno = new double[sample_size + 2];
		generate_LR_geno_covar_pheno_data(rng, per_ind_feat_matrix, per_ind_pheno, per_feat_weights, pheno_noise_std_dev, sample_size, n_covars);
	} // feats file check.

	// End of feature/genotype loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Training plaintext LR with sample size of %d:\n\
n_covars: %d\n\
pheno_noise_std_dev: %.5f\n\
update_step: %.5f\n\
regularization_weight: %.5f\n\
n_epoch: %d\n\
op_dir: %s\n\
LL_EPSILON: %.6f\n", sample_size, 
n_covars, pheno_noise_std_dev,  update_step, regularization_weight, n_epoch, op_dir, LL_EPSILON);

	// Save the pheno/geno and model data.
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

	// Write the model.
	fprintf(stderr, "Model weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		fprintf(stderr, "Feat %d: %.5f\n", i_feat, per_feat_weights[i_feat]);
	} // i_feat loop.

	// Allocate the weights for the current epoch.
	double* cur_epoch_weights = new double[n_feats + 2];
	memset(cur_epoch_weights, 0, sizeof(double) * (n_feats + 2));

	// TODO::Keep track of the sigmoid parameters; this is necessary to make sure sigmoid approximations work.
	FILE* f_per_epoch_info = open_f("per_epoch_info.txt", "w");
	double* cur_epoc_cur_sample_per_feat_gradients = new double[n_feats + 2];
	vector<double>* per_epoch_LL = new vector<double>();

	// Compute the first LL.
	double init_epoch_LL = get_log_likelihood(per_ind_feat_matrix, per_ind_pheno, cur_epoch_weights, n_feats);
	per_epoch_LL->push_back(init_epoch_LL);

	for (int epoch_i = 0; epoch_i < n_epoch; epoch_i++)
	{
		fprintf(stderr, "@ Epoch %d\n", epoch_i);

		// This is the per sample gradient for SGD; we keep it for each feature because the weights are updated after all gradients are computed.
		memset(cur_epoc_cur_sample_per_feat_gradients, 0, sizeof(double) * n_feats);

		vector<int>* sgd_rand_subj_i = rng->fast_permute_indices(0, sample_size);

		// Each sample updates all parameters once.
		for (int sgd_i_s = 0; sgd_i_s < sample_size; sgd_i_s++)
		{
			int i_s = sgd_rand_subj_i->at(sgd_i_s);

			fprintf(stderr, "@ Sample %d (%d)             \r", sgd_i_s, i_s);

			for (int feat_i = 0; feat_i < n_feats; feat_i++)
			{
				// This computes the inner product of the features, including the intercept and generates the predicted phenotype.
				double cur_pred_pheno = get_sigmoid_val_per_feats_weights(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);

				double inv_sample_size = ((double)(1.0) / sample_size);

				double cur_sample_cur_feat_gradient = 0;

				// Do not regularize the intercept.
				if (feat_i == GET_INTERCEPT_FEAT_I())
				{
					if (per_ind_feat_matrix->at(i_s)[feat_i] != 1)
					{
						fprintf(stderr, "Sanity check failed, intercept feature value is not 1.0 for subject %d: %.4f\n", i_s, per_ind_feat_matrix->at(i_s)[feat_i]);
						exit(0);
					}

					// Update step should be tuned to correct for inverse sample size factor.
					cur_sample_cur_feat_gradient = update_step * inv_sample_size *
						(cur_pred_pheno - per_ind_pheno[i_s]) * per_ind_feat_matrix->at(i_s)[feat_i];
				}
				else
				{
					cur_sample_cur_feat_gradient = update_step * inv_sample_size *
						((cur_pred_pheno - per_ind_pheno[i_s]) * per_ind_feat_matrix->at(i_s)[feat_i] +
						(regularization_weight * cur_epoch_weights[feat_i]));
				}
				
				cur_epoc_cur_sample_per_feat_gradients[feat_i] = cur_sample_cur_feat_gradient;
				//cur_epoch_weights[feat_i] = cur_epoch_weights[feat_i] - cur_sample_cur_feat_gradient;
			} // feat_i loop.

			// Update the weights with the current sample's gradient.
			for (int feat_i = 0; feat_i < n_feats; feat_i++)
			{
				cur_epoch_weights[feat_i] = cur_epoch_weights[feat_i] - cur_epoc_cur_sample_per_feat_gradients[feat_i];
			} // feat_i loop.
		} // i_s loop.

		// Write the current model's likelihood.
		double cur_epoch_updated_LL = get_log_likelihood(per_ind_feat_matrix, per_ind_pheno, cur_epoch_weights, n_feats);
		fprintf(f_per_epoch_info, "%d\t%.6f\n", epoch_i, cur_epoch_updated_LL);	

		// Do an early stop check.
		if (per_epoch_LL->size() > 0)
		{
			double prev_epoch_LL = per_epoch_LL->back();
			per_epoch_LL->push_back(cur_epoch_updated_LL);

			fprintf(stderr, "\nLL: %.5f (%.5f)\n", cur_epoch_updated_LL, prev_epoch_LL);

			if (fabs(cur_epoch_updated_LL - prev_epoch_LL) < LL_EPSILON)
			{
				break;
			}
		}
	} // epoch_i loop

	// Clean gradient memory.
	delete[] cur_epoc_cur_sample_per_feat_gradients;

	close_f(f_per_epoch_info, "per_epoch_info.txt");

	FILE* f_trained_model = open_f("trained_model.txt", "w");
	fprintf(f_trained_model, "N_COVARS\t%d\n", n_covars);
	int covar_i = 0;
	fprintf(stderr, "Real/Trained weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		fprintf(stderr, "%.4f\t%.4f\n", per_feat_weights[i_feat], cur_epoch_weights[i_feat]);

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

	// Generate testing data.
	vector<double*>* per_testing_ind_feat_matrix = new vector<double*>();
	double* per_testing_ind_pheno = new double[sample_size + 2];
	generate_LR_geno_covar_pheno_data(rng, per_testing_ind_feat_matrix, per_testing_ind_pheno, per_feat_weights, pheno_noise_std_dev, sample_size, n_covars);

	int n_correct = 0;
	double total_prob_error = 0;
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_pred_pheno = (get_sigmoid_val_per_feats_weights(per_testing_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights) > 0.5) ? (1) : (0);

		if (cur_pred_pheno == per_testing_ind_pheno[i_s])
		{
			n_correct++;
		}

		total_prob_error += abs(cur_pred_pheno - per_testing_ind_pheno[i_s]);

		fprintf(stderr, "Sample %d: %.0f / %.0f\n", i_s, per_testing_ind_pheno[i_s], cur_pred_pheno);
	} // i_s loop.

	fprintf(stderr, "%d/%d correct, Prob error: %.5f.\n", n_correct, sample_size, total_prob_error);

	double* per_feat_LRT_score = new double[n_feats + 2];
	compute_LRT_p_value(per_feat_LRT_score, per_ind_feat_matrix, per_ind_pheno, per_feat_weights, n_covars);

	//for(int )
	//per_feat_LRT_score
} // train_LR_baseline function.

void train_modified_LR(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp, char* model_fp,
	int sample_size, int n_covars,
	double pheno_noise_std_dev,
	double update_step,
	double regularization_weight,
	int n_epoch,
	int gradient_addition_order_option,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* op_dir)
{
	// We need the two sigmoid functions here.
	double(*GET_SIGMOID_PER_FEATS_WEIGHTS)(double*, int, double*);
	//double(*GET_SIGMOID_PER_FEAT_COMB)(double);

	if (sigmoid_approx_type == SIGMOID_APPROX_NATIVE)
	{
		fprintf(stderr, "Using native sigmoid..\n");
		GET_SIGMOID_PER_FEATS_WEIGHTS = get_sigmoid_val_per_feats_weights;
		//GET_SIGMOID_PER_FEAT_COMB = get_sigmoid_val_per_feat_comb;
	}
	else if (sigmoid_approx_type == SIGMOID_APPROX_KIM_ETAL)
	{
		fprintf(stderr, "Using native Kim. etal sigmoid..\n");
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

	if (gradient_addition_order_option == GRADIENT_ADD_PER_FEAT)
	{
		fprintf(stderr, "Using per feature gradient update.\n");
	}
	else if (gradient_addition_order_option == GRADIENT_ADD_PER_SUBJECT)
	{
		fprintf(stderr, "Using per subject gradient update.\n");
	}
	else if (gradient_addition_order_option == GRADIENT_ADD_PER_SAMPLE)
	{
		fprintf(stderr, "Using per sample gradient update.\n");
	}
	else
	{
		fprintf(stderr, "Unknown gradient order update: %d\n", gradient_addition_order_option);
		exit(0);
	}

	// this is the threshold at which the training is stopped early.
	//double LL_EPSILON = 0.01;

	// We need the RNG for different tasks.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// The training data/model.
	vector<double*>* per_ind_feat_matrix = NULL;
	double* per_ind_pheno = NULL;
	double* per_feat_weights = NULL;
	int n_feats = 0;

	// Load or generate model.
	if (check_file(model_fp))
	{
		int model_loaded_covars = 0;
		per_feat_weights = load_LR_model_weights(model_fp, model_loaded_covars);

		n_covars = model_loaded_covars;
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);
	}
	else
	{
		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// Generate the random model data.
		per_feat_weights = generate_per_feat_weights(rng, n_feats);

	} // model loading check.

	// Model loading is finished, at this stage, we should have the n_feats, n_covars, and model setup.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// At this point, n_feats and n_covars must be set.
	if (n_feats == 0 ||
		n_covars == 0)
	{
		fprintf(stderr, "Could not set n_feats after loading/setting the model..\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Loaded the model with %d covariates (%d features) into memory.\n", n_covars, n_feats);
	}

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

		if (n_covars != n_covars_per_feats_file)
		{
			fprintf(stderr, "Could not match the covariate number between model and features file: %d, %d\n", n_covars, n_covars_per_feats_file);
			exit(0);
		}

		// Copy the phenotype vector to an array.
		per_ind_pheno = new double[per_ind_feat_matrix->size() + 2];
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
				fprintf(stderr, "\t%d: %.4f", i_feat, per_ind_feat_matrix->at(i_s)[i_feat]);
			} // i_feat loop.

			fprintf(stderr, "Pheno: %.3f\n", per_ind_pheno[i_s]);
		} // i_s loop.

		// Sample size must be fixed after loading.
		sample_size = per_ind_feat_matrix->size();
	} // feats file check.
	else
	{
		fprintf(stderr, "Generating the training data for %d sample and %d features.\n", sample_size, n_feats);

		n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

		// We first generate a random feature matrix.
		per_ind_feat_matrix = new vector<double*>();
		per_ind_pheno = new double[sample_size + 2];
		generate_LR_geno_covar_pheno_data(rng, per_ind_feat_matrix, per_ind_pheno, per_feat_weights, pheno_noise_std_dev, sample_size, n_covars);
	} // feats file check.

	// End of feature loading/generation step.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Training plaintext LR with sample size of %d:\n\
n_covars: %d\n\
pheno_noise_std_dev: %.5f\n\
update_step: %.5f\n\
regularization_weight: %.5f\n\
n_epoch: %d\n\
op_dir: %s\n\
LL_EPSILON: %.6f\n", sample_size,
n_covars, pheno_noise_std_dev, update_step, regularization_weight, n_epoch, op_dir, LL_EPSILON);

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

	// Write the model.
	fprintf(stderr, "Model weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		fprintf(stderr, "Feat %d: %.5f\n", i_feat, per_feat_weights[i_feat]);
	} // i_feat loop.

	// Allocate the weights for the current epoch.
	double* cur_epoch_weights = new double[n_feats + 2];
	memset(cur_epoch_weights, 0, sizeof(double) * (n_feats + 2));

	// TODO::Keep track of the sigmoid parameters; this is necessary to make sure sigmoid approximations work.
	FILE* f_per_epoch_info = open_f("per_epoch_info.txt", "w");
	double* cur_epoc_cur_subject_per_feat_gradients = new double[n_feats + 2]; // This is the default SGD gradient update value, where each gradient is updated for each subject at once.
	vector<double>* per_epoch_LL = new vector<double>();

	double* cur_epoc_samplewide_aggregated_per_feat_gradient_update = new double[n_feats + 2];

	// Compute the first LL.
	double init_epoch_LL = get_log_likelihood(per_ind_feat_matrix, per_ind_pheno, cur_epoch_weights, n_feats);
	per_epoch_LL->push_back(init_epoch_LL);

	for (int epoch_i = 0; epoch_i < n_epoch; epoch_i++)
	{
		fprintf(stderr, "@ Epoch %d\n", epoch_i);

		// Set the samplewide aggregated gradients for the current epoch.
		memset(cur_epoc_samplewide_aggregated_per_feat_gradient_update, 0, sizeof(double) * n_feats);

		vector<int>* sgd_rand_subj_i = rng->fast_permute_indices(0, sample_size);

		// Each sample updates all parameters once.
		for (int sgd_i_s = 0; sgd_i_s < sample_size; sgd_i_s++)
		{
			// This is the per subject gradient for SGD; we keep it for each feature because the weights are updated after all gradients are computed.
			memset(cur_epoc_cur_subject_per_feat_gradients, 0, sizeof(double) * n_feats);

			int i_s = sgd_rand_subj_i->at(sgd_i_s);

			fprintf(stderr, "@ Sample %d (%d)             \r", sgd_i_s, i_s);

			for (int feat_i = 0; feat_i < n_feats; feat_i++)
			{
				// This computes the inner product of the features, including the intercept and generates the predicted phenotype.
				// The inference utilizes the sigmoid approximation selected by the argument.
				double cur_pred_pheno = GET_SIGMOID_PER_FEATS_WEIGHTS(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);

				double inv_sample_size = ((double)(1.0) / sample_size);

				double cur_sample_cur_feat_gradient = 0;

				// Do not regularize the intercept.
				if (feat_i == GET_INTERCEPT_FEAT_I())
				{
					if (per_ind_feat_matrix->at(i_s)[feat_i] != 1)
					{
						fprintf(stderr, "Sanity check failed, intercept feature value is not 1.0 for subject %d: %.4f\n", i_s, per_ind_feat_matrix->at(i_s)[feat_i]);
						exit(0);
					}

					// Update step should be tuned to correct for inverse sample size factor.
					cur_sample_cur_feat_gradient = update_step * inv_sample_size *
						(cur_pred_pheno - per_ind_pheno[i_s]) * per_ind_feat_matrix->at(i_s)[feat_i];
				}
				else
				{
					cur_sample_cur_feat_gradient = update_step * inv_sample_size *
						((cur_pred_pheno - per_ind_pheno[i_s]) * per_ind_feat_matrix->at(i_s)[feat_i] +
						(regularization_weight * cur_epoch_weights[feat_i]));
				}

				cur_epoc_cur_subject_per_feat_gradients[feat_i] = cur_sample_cur_feat_gradient;

				// Update the samplewide aggregated gradient for this feature.
				cur_epoc_samplewide_aggregated_per_feat_gradient_update[feat_i] += cur_sample_cur_feat_gradient;

				// If gradient updates are done per feature, update here.
				if (gradient_addition_order_option == GRADIENT_ADD_PER_FEAT)
				{
					cur_epoch_weights[feat_i] = cur_epoch_weights[feat_i] - cur_sample_cur_feat_gradient;
				}
			} // feat_i loop.

			// If gradient updates are done per sample, compute all feat gradients and set them here.
			if (gradient_addition_order_option == GRADIENT_ADD_PER_SUBJECT)
			{
				// Update the weights with the current sample's gradient.
				for (int feat_i = 0; feat_i < n_feats; feat_i++)
				{
					cur_epoch_weights[feat_i] = cur_epoch_weights[feat_i] - cur_epoc_cur_subject_per_feat_gradients[feat_i];
				} // feat_i loop.
			}
		} // sgd_i_s loop.

		// Per sample weight updates are performed after all subjects are processed.
		if (gradient_addition_order_option == GRADIENT_ADD_PER_SAMPLE)
		{
			for (int feat_i = 0; feat_i < n_feats; feat_i++)
			{
				// This makes sure that the amount of update is small.
				cur_epoch_weights[feat_i] = cur_epoch_weights[feat_i] - (cur_epoc_samplewide_aggregated_per_feat_gradient_update[feat_i] / sample_size);
				//cur_epoch_weights[feat_i] = cur_epoch_weights[feat_i] - (cur_epoc_samplewide_aggregated_per_feat_gradient_update[feat_i] / update_step);
			}
		}

		// Write the current model's likelihood.
		double cur_epoch_updated_LL = get_log_likelihood(per_ind_feat_matrix, per_ind_pheno, cur_epoch_weights, n_feats);
		fprintf(f_per_epoch_info, "%d\t%.6f\n", epoch_i, cur_epoch_updated_LL);

		// Do an early stop check.
		if (per_epoch_LL->size() > 0)
		{
			double prev_epoch_LL = per_epoch_LL->back();
			per_epoch_LL->push_back(cur_epoch_updated_LL);

			fprintf(stderr, "\nLL: %.5f (%.5f)\n", cur_epoch_updated_LL, prev_epoch_LL);

			if (fabs(cur_epoch_updated_LL - prev_epoch_LL) < LL_EPSILON)
			{
				break;
			}
		}
	} // epoch_i loop

	// Clean gradient memory.
	delete[] cur_epoc_cur_subject_per_feat_gradients;
	delete[] cur_epoc_samplewide_aggregated_per_feat_gradient_update;

	close_f(f_per_epoch_info, "per_epoch_info.txt");

	// Calculate the phenotype expectations for each subject.
	FILE* f_per_subject_pheno = open_f("per_subject_phenotypes.txt", "w");
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_pred_pheno = GET_SIGMOID_PER_FEATS_WEIGHTS(per_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights);
		fprintf(f_per_subject_pheno, "%d\t%.5f\n", i_s, cur_pred_pheno);
	} // i_s loop.

	close_f(f_per_subject_pheno, "per_subject_phenotypes.txt");

	FILE* f_trained_model = open_f("trained_model.txt", "w");
	fprintf(f_trained_model, "N_COVARS\t%d\n", n_covars);
	int covar_i = 0;
	fprintf(stderr, "Real/Trained weights:\n");
	for (int i_feat = 0; i_feat < n_feats; i_feat++)
	{
		fprintf(stderr, "%.4f\t%.4f\n", per_feat_weights[i_feat], cur_epoch_weights[i_feat]);

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

	// Generate testing data.
	vector<double*>* per_testing_ind_feat_matrix = new vector<double*>();
	double* per_testing_ind_pheno = new double[sample_size + 2];
	generate_LR_geno_covar_pheno_data(rng, per_testing_ind_feat_matrix, per_testing_ind_pheno, per_feat_weights, pheno_noise_std_dev, sample_size, n_covars);

	int n_correct = 0;
	double total_prob_error = 0;
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_pred_pheno = (get_sigmoid_val_per_feats_weights(per_testing_ind_feat_matrix->at(i_s), n_feats, cur_epoch_weights) > 0.5) ? (1) : (0);

		if (cur_pred_pheno == per_testing_ind_pheno[i_s])
		{
			n_correct++;
		}

		total_prob_error += abs(cur_pred_pheno - per_testing_ind_pheno[i_s]);

		fprintf(stderr, "Sample %d: %.0f / %.0f\n", i_s, per_testing_ind_pheno[i_s], cur_pred_pheno);
	} // i_s loop.

	fprintf(stderr, "%d/%d correct, Prob error: %.5f.\n", n_correct, sample_size, total_prob_error);

	double* per_feat_LRT_score = new double[n_feats + 2];
	compute_LRT_p_value(per_feat_LRT_score, per_ind_feat_matrix, per_ind_pheno, per_feat_weights, n_covars);
} // train_modified_LR function.



