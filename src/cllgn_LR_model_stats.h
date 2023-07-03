#ifndef __LR_model_stats__
#define __LR_model_stats__

void compute_LRT_p_value(double* per_feat_LRT_score, vector<double*>* per_ind_feature_matrix, double* per_ind_pheno, double* weights, int n_covars);

double get_log_likelihood(vector<double*>* per_ind_feat_matrix, double* per_ind_pheno, double* cur_model_weights, int n_feats);

double get_log_likelihood_per_mu_matrix(double** mu, double* per_ind_pheno, int sample_size);

#endif // __LR_model_stats__