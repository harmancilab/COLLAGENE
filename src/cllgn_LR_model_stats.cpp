#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_file_utils.h"
#include "cllgn_histogram.h"
#include <math.h>
#include <string.h>
#include "cllgn_features_weight_utils.h"
#include "cllgn_LR_model_stats.h"
#include "cllgn_x_chisqr.h"
#include "cllgn_x_sigmoid.h"
#include "cllgn_LR_utils.h"

/*
From : https://web.stanford.edu/class/archive/stats/stats200/stats200.1172/Lecture26.pdf
For the last weight, compute
	D = 2log( l(X|\beta) ) - 2log( l(X|\beta_0) )
*/
void compute_LRT_p_value(double* per_feat_LRT_score, vector<double*>* per_ind_feature_matrix, double* per_ind_pheno, double* weights, int n_covars)
{
	//int sample_size = (int)(per_ind_feature_matrix->size());

	int n_feats = GET_N_FEATS_PER_N_COVARS(n_covars);

	double base_model_LL = get_log_likelihood(per_ind_feature_matrix, per_ind_pheno, weights, n_feats);

	for (int feat_i = 0; feat_i < n_feats; feat_i++)
	{
		double* cur_model_weights = new double[n_feats];
		copy_weights(cur_model_weights, weights, n_feats);

		// set the current feature's weight to 0.
		cur_model_weights[feat_i] = 0;

		double cur_feat_excluded_model_model_LL = get_log_likelihood(per_ind_feature_matrix, per_ind_pheno, cur_model_weights, n_feats);

		double LRT_score = -2 * (cur_feat_excluded_model_model_LL - base_model_LL);

		double chsqr_pval = chisqr(1, LRT_score);

		fprintf(stderr, "Feat %d: LRT-score: %.5f (%lf)\n", feat_i, LRT_score, log(chsqr_pval) / log(10.0));

		delete[] cur_model_weights;
	} // feat_i loop.
}

double get_log_likelihood(vector<double*>* per_ind_feat_matrix, double* per_ind_pheno, double* cur_model_weights, int n_feats)
{
	int sample_size = per_ind_feat_matrix->size();

	// Compute the log likelihood at this stage.
	double cur_epoch_log_likelihood = 0;
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_pred_pheno = get_sigmoid_val_per_feats_weights(per_ind_feat_matrix->at(i_s), n_feats, cur_model_weights);

		cur_epoch_log_likelihood += (per_ind_pheno[i_s] * log(cur_pred_pheno) + (1 - per_ind_pheno[i_s]) * log(1 - cur_pred_pheno));
	} // i_s loop.

	return(cur_epoch_log_likelihood);
}

double get_log_likelihood_per_mu_matrix(double** mu, double* per_ind_pheno, int sample_size)
{
	// Compute the log likelihood at this stage.
	double cur_epoch_log_likelihood = 0;
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		double cur_pred_pheno = mu[i_s][0];

		cur_epoch_log_likelihood += (per_ind_pheno[i_s] * log(cur_pred_pheno) + (1 - per_ind_pheno[i_s]) * log(1 - cur_pred_pheno));
	} // i_s loop.

	return(cur_epoch_log_likelihood);
}