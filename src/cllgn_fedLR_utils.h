#ifndef __FEDLR_UTILS__
#define __FEDLR_UTILS__

// This defines the global set of parameters, makes it easier to copy parameters between functions.
struct t_FedLR_environment_params
{
	int i_iter;
	int client_i;
	int n_clients;
	char* subject_per_row_feats_fp;
	char* subject_per_row_pheno_fp;
	int n_epoch;
	int sigmoid_approx_type;
	double LL_EPSILON;
	char* private_working_dir;
	char* shared_working_dir;;
};

double get_1DOF_chisqr_pval(double norm_chisqr_scale);

void get_pval_stat_op_fp(char* shared_working_dir, int client_i, char* p_val_stats_fp);

void get_GtWX_Xi_XtWG_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer);

void get_GtWX_matrix_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer);

void get_GtWG_matrix_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer);

void get_chisq_stat_matrix_fp(char* shared_working_dir,
	int var_block_start_i, int var_block_end_i,
	int client_i,
	char* chisq_buffer, char* var_block_ids_buffer);

void get_scale_matrix_fp(char* shared_working_dir, int var_block_start_i, int var_block_end_i, int client_i, char* buffer);

void get_LL_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_mu_matrix_fp(char* shared_working_dir, int i_iter, int client_i, char* buffer);

void get_beta_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_XtWX_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_XtWz_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void client_pool_XtWX_XtWz_update_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* shared_working_dir);

void client_calculate_save_XtWX_XtWz(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* shared_working_dir);

void client_check_convergence_per_updated_beta(int i_iter, int client_i, int n_clients, char* shared_working_dir, double LL_EPSILON);

void client_calculate_save_pvalue_stats_for_meta_analysis(int client_i, int n_clients, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_fp,
	char* shared_working_dir);

void client_calculate_save_pvalue_stats(int client_i, int n_clients, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_matrix_fp,
	char* shared_working_dir);

void client_update_pvalue_scale_stats_per_block(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, char* null_model_pheno_matrix_fp,
	char* shared_working_dir);

void client_pool_pvalue_stats_for_meta_analysis(int client_i, int n_clients, int var_block_size, char* shared_working_dir);
void client_pool_pvalue_stats(int client_i, int n_clients, int var_block_size, char* shared_working_dir);

#endif // __FEDLR_UTILS__