#ifndef __FEDLR_SECURE_CONVERTIBLE_PROTOCOL__
#define __FEDLR_SECURE_CONVERTIBLE_PROTOCOL__

#include <vector>
using namespace std;

void plain_add_matrices_per_list(char* matrix_path_list_fp, char* op_fp);
double** generate_multiplicative_full_noise_matrix(int nrow, int ncol);
double** generate_multiplicative_diagonal_noise_matrix(int nrow, int ncol);
double** generate_constant_diagonal_matrix(int nrow, int ncol, double diagonal_val);

void get_inv_XtWX_matrix_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_site_specific_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_sitewise_pooled_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_partdec_sitewise_pooled_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, int dec_client_i, char* buffer);

void get_fulldec_sitewise_pooled_all_site_noise_XtWX_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void get_mult_noise_fp(char* shared_working_dir, int iter_i, int client_i, char* buffer);

void cryptable_client_calculate_save_XtWX_XtWz(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_add_mult_noise_2_XtWX(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_pool_site_specific_all_site_noise_XtWX(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_collaborative_decrypt_pooled_noisy_XtWX(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

// This is the last step: Pool the partial decryptions then remove noise, which can only be done in encrypted domain.
void cryptable_pool_partially_decrypted_pooled_noisy_XtWx_remove_noise(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_pool_XtWX_XtWz_update_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_collaborative_decrypt_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_pool_partially_decrypted_beta(int i_iter, int client_i, int n_clients,
	char* subject_per_row_feats_fp, char* subject_per_row_pheno_fp,
	int n_epoch,
	int sigmoid_approx_type,
	double LL_EPSILON,
	char* private_working_dir,
	char* shared_working_dir);

// This function is skipped for now, the sites can just run 10 iterations and assume convergence.
void cryptable_client_check_convergence_per_updated_beta(int i_iter, int client_i, int n_clients, char* shared_working_dir, double LL_EPSILON);

void cryptable_client_calculate_save_pvalue_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_pool_pvalue_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp, 
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_collaborative_decrypt_pval_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp,
	char* private_working_dir,
	char* shared_working_dir);

void cryptable_client_pool_partially_decrypted_pval_stats(int client_i, int n_clients, int i_iter, int var_block_size, char* GMMAT_text_genotype_matrix_fp,
	char* subject_per_row_feats_fp,
	char* obs_pheno_fp,
	char* private_working_dir,
	char* shared_working_dir);



#endif // __FEDLR_SECURE_CONVERTIBLE_PROTOCOL__

