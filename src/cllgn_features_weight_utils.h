#ifndef __FEATURES_WEIGHTS_UTILS__
#define __FEATURES_WEIGHTS_UTILS__

#include <vector>
using namespace std;

class t_rng;

//#define GET_GENO_FEAT_I(n_feats) ((n_feats)-1)
#define GET_INTERCEPT_FEAT_I() (0)
#define GET_N_FEATS_PER_N_COVARS(n_covars) ((n_covars) + 1)
#define GET_N_COVARS_PER_N_FEATS(n_feats) ((n_feats) - 1)

void generate_save_random_feats_model(int sample_size, int n_covars, double pheno_noise_std_dev, char* op_dir);

void generate_LR_geno_covar_pheno_data(t_rng* rng,
	vector<double*>* per_ind_feat_matrix,
	double* per_ind_pheno,
	double* per_feat_weights,
	double noise_std_dev,
	int sample_size, int n_covars);

void generate_save_random_feats_model(int sample_size, int n_covars, double pheno_noise_std_dev, char* op_dir);

void normalize_feats(vector<double*>* per_ind_feature_matrix, int n_covars);

double* load_LR_model_weights(char* model_fp, int& model_loaded_n_covars);

void load_LR_data_row_samples(char* LR_data_fp, char* pheno_fp, int& n_covars, vector<double*>* per_ind_feat_matrix, vector<double>* per_ind_pheno_vector);

void save_LR_data_row_samples(vector<double*>* per_ind_feat_matrix, double* per_ind_pheno, double* per_feat_weights, int n_covars, char* op_dir);

double* generate_per_feat_weights(t_rng* rng, int n_covars);

void copy_weights(double* copy_per_feat_weights, double* per_feat_weights, int n_feats);

#endif // __FEATURES_WEIGHTS_UTILS__
