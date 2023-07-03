#ifndef __SEAL_PLAIN_LR_UTILS__
#define __SEAL_PLAIN_LR_UTILS__

enum { GRADIENT_ADD_PER_FEAT, GRADIENT_ADD_PER_SUBJECT, GRADIENT_ADD_PER_SAMPLE };

enum {SIGMOID_APPROX_NATIVE, 
SIGMOID_APPROX_TENSEAL,
SIGMOID_APPROX_KIM_ETAL};

// Use this function to evaluate LR training parameters such as sigmoid approximation, parameter initialization, and # of iterations.
void train_LR_baseline(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp, char* model_fp,
	int sample_size, int n_covars,
	double pheno_noise_std_dev,
	double update_step,
	double regularization_weight,
	int n_epoch,
	char* op_dir);

void train_modified_LR(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp, char* model_fp,
	int sample_size, int n_covars,
	double pheno_noise_std_dev,
	double update_step,
	double regularization_weight,
	int n_epoch,
	int gradient_addition_order_option,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* op_dir);

void train_modified_LR_per_IRLS(char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* op_dir);

void save_LR_data_row_samples(vector<double*>* per_ind_feat_matrix, double* per_ind_pheno, double* per_feat_weights, int n_covars, char* op_dir);

void assign_ChiSqr_p_values_per_fit_model_pheno_per_GMMAT_text_geno(char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_fp,
	char* op_fp);
#endif // __SEAL_PLAIN_LR_UTILS__
