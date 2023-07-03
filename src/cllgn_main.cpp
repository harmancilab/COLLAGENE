#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cllgn_file_utils.h"
#include "cllgn_ansi_string.h"
#include "cllgn_annot_region_tools.h"
#include "cllgn_variation_tools.h"
#include "cllgn_multicolumn_processing.h"

#include "cllgn_nomenclature.h"

#include "cllgn_features_weight_utils.h"
#include "cllgn_LR_model_stats.h"
#include "cllgn_LR_utils.h"

#include "cllgn_x_sigmoid.h"
#include "cllgn_x_chisqr.h"
#include "cllgn_matrix_linalg_utils.h"
#include "cllgn_fedLR_secure_convertible_protocol_utils.h"
#include "cllgn_fedLR_utils.h"

#include "cllgn_multicolumn_processing.h"

#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_seal_genomics_matrix_utils.h"
#include "SEAL-main/native/examples/examples.h"
#include <seal/util/common.h>
#include "seal/keygenerator.h"
#include "seal/randomtostd.h"
#include "seal/util/common.h"
#include "seal/util/galois.h"
#include "seal/util/ntt.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include "seal/util/rlwe.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintcore.h"

// Functions that are essentially used for key sharing protocols: Encrypt, Partial Decrypt, Pooling.
#include "cllgn_seal_genomics_key_sharing_utils.h"

#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Data Processing/Formatting:\n\
	-separate_VCF_2_chroms\n\
	-convert_VCF_2_matbed\n\
	-extract_variant_regions\n\
	-extract_subset_subjects\n\
	-write_GMMAT_formatted_geno\n\
	-impute_missing_GMMAT_formatted_geno\n\
	-extract_subsample_GMMAT_genotypes_from_GMMAT_genotype_matrix\n\
	-extract_rows_per_query_column_preserve_query_order\n\
Utilities:\n\
	-save_matrix_plain_2_bin\n\
	-dump_matrix_plain\n\
	-assign_chisqr_pvals_per_ST_stats\n\
	-generate_mult_diagonal_noise_matrix\n\
	-generate_mult_full_noise_matrix\n\
Matrix Library Options:\n\
	-scalar_multiply_matrix_plain\n\
	-plain_unpad_matrix_to_size\n\
	-plain_pad_matrix_to_next_power_of_2\n\
	-plain_pad_matrix_rows_to_next_power_of_2\n\
	-plain_pad_matrix_cols_to_next_power_of_2\n\
	-plain_invert_matrix\n\
	-plain_add_matrices\n\
	-plain_add_matrices_per_list\n\
	-generate_random_pt_matrix\n\
	-generate_constant_full_matrix\n\
	-generate_constant_diagonal_matrix\n\
	-multiply_matrices_pt\n\
	-multiply_matrices_elementwise_pt\n\
	-row2row_multiply_pt\n\
	-transpose_pt_matrix\n\
	-transform_elementwise_per_callback_pt\n\
Continuous encrypted Secure Matrix Operations:\n\
	-validate_ckks_text_params\n\
	-transpose_continuous_encrypted_vector\n\
	-write_enc_matrix_dimensions\n\
	-write_plain_matrix_dimensions\n\
	-continuous_encrypt_data_matrix\n\
	-encrypt_encoded_pt_matrix\n\
	-col_expand_dense_encrypt_matrix\n\
	-row_expand_dense_encrypt_matrix\n\
	-row_expand_continuous_encrypted_matrix\n\
	-secure_multiply_matrices_Acol_Brow_expansions\n\
	-secure_add_cont_ct_matrices\n\
	-secure_add_cont_ct_matrices_per_list\n\
	-secure_sub_cont_ct_matrices\n\
	-secure_elementwise_mul_cont_ct_matrices\n\
	-secure_row2row_inner_prod_continuous_encrypted_matrices\n\
	-get_uniform_rand\n\
	-get_normal_rand\n\
Encrypt/Decrypt Continuous Encrypted Matrices:\n\
	-fully_decrypt_continuous_encrypted_matrix\n\
	-partial_decrypt_continuous_enc_per_noisy_secretkey\n\
	-pool_partial_decrypted_continuous_enc_data\n\
CT-Refresh:\n\
	-generate_plaintext_mask_per_continuous_encrypted_data\n\
	-additive_mask_continuous_encrypted_data\n\
	-additive_unmask_continuous_encrypted_data\n\
	-write_continuous_encrypted_ciphertext_vital_stats\n\
DSK Options:\n\
	-generate_per_site_noisy_keys\n", argv[0]);
	exit(0);
	}

	// Initialize the global ID's.
	t_g_Column_IDs::init_global_IDs();

	try
	{
		if (t_string::compare_strings(argv[1], "-transform_elementwise_per_callback_pt"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "USAGE: %s %s [Matrix file] [Function (log/exp/sigmoid)] [Output matrix file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			char* callback_name = argv[3];
			char* op_matrix_fp = argv[4];

			int nrows = 0;
			int ncols = 0;
			double** input_matrix = load_matrix_binary(matrix_fp, nrows, ncols);

			double** res_matrix = NULL;
			if (t_string::compare_strings(callback_name, "sigmoid"))
			{
				res_matrix = process_matrix_elementwise_by_callback(input_matrix, nrows, ncols, get_sigmoid_val_per_feat_comb, NULL);
			}
			else if (t_string::compare_strings(callback_name, "exp"))
			{
				res_matrix = process_matrix_elementwise_by_callback(input_matrix, nrows, ncols, exp, NULL);
			}
			else if (t_string::compare_strings(callback_name, "log"))
			{
				res_matrix = process_matrix_elementwise_by_callback(input_matrix, nrows, ncols, get_log_val_protected, NULL);
			}
			
			if (res_matrix != NULL)
			{
				save_matrix_binary(res_matrix, nrows, ncols, op_matrix_fp);
			}
		} // -transform_elementwise_per_callback_pt option.
		else if (t_string::compare_strings(argv[1], "-get_uniform_rand"))
		{
			if (argc != 3)
			{
				fprintf(stderr, "USAGE: %s %s [text params file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* text_params_path = argv[2];

			t_text_params* text_params = load_text_params(text_params_path);

			//double secret_key_noise_variance = text_params->per_site_key_noise_variance;
			//double scale = pow(2, text_params->scale_n_bits);

			vector<int> coeff_modulus_bit_sizes;
			for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
			{
				coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
			} // i_dec loop.

			EncryptionParameters parms(scheme_type::ckks);

			size_t poly_modulus_degree = text_params->poly_modulus_degree;
			parms.set_poly_modulus_degree(poly_modulus_degree);
			parms.set_coeff_modulus(CoeffModulus::Create(
				poly_modulus_degree, coeff_modulus_bit_sizes));

			SEALContext context(parms);
			auto &context_data = *context.key_context_data();
			//auto &context_data = *context.first_context_data();
			auto key_context_parms = context_data.parms();

			uint64_t* buff = SEAL_Genomics_Get_Unif_Random_Custom(key_context_parms, 30);

			// At this point, buff should be filled with random normal values from the clipped normal distribution.
			FILE* f_op = open_f("SEAL_unif_rands.txt", "w");
			for (int i_mod = 0; i_mod < (int)key_context_parms.coeff_modulus().size(); i_mod++)
			{
				for (size_t coeff_i = 0; coeff_i < key_context_parms.poly_modulus_degree(); coeff_i++)
				{
					fprintf(f_op, "%d\t%d\t%lu\t%lu\n", i_mod, (int)coeff_i, buff[i_mod * key_context_parms.poly_modulus_degree() + coeff_i], key_context_parms.coeff_modulus().at(i_mod).value());
				} // i loop.
			}
			close_f(f_op, "SEAL_rands.txt");
		} // get_uniform_rand
		else if (t_string::compare_strings(argv[1], "-get_normal_rand"))
		{
			if (argc != 3)
			{
				fprintf(stderr, "USAGE: %s %s [text params file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* text_params_path = argv[2];

			t_text_params* text_params = load_text_params(text_params_path);

			int secret_key_noise_variance_n_bits = text_params->per_site_key_noise_var_bit_size;
			double noise_std_dev = pow(2, secret_key_noise_variance_n_bits / 2);

			vector<int> coeff_modulus_bit_sizes;
			for (int i_dec = 0; i_dec < (int)(text_params->decomp_bit_sizes_ptr->size()); i_dec++)
			{
				coeff_modulus_bit_sizes.push_back(text_params->decomp_bit_sizes_ptr->at(i_dec));
			} // i_dec loop.

			EncryptionParameters parms(scheme_type::ckks);

			size_t poly_modulus_degree = text_params->poly_modulus_degree;
			parms.set_poly_modulus_degree(poly_modulus_degree);
			parms.set_coeff_modulus(CoeffModulus::Create(
				poly_modulus_degree, coeff_modulus_bit_sizes));

			SEALContext context(parms);
			auto &context_data = *context.key_context_data();
			//auto &context_data = *context.first_context_data();
			auto key_context_parms = context_data.parms();

			// Generate random values at the key context level.
			uint64_t* buff = SEAL_Genomics_Get_Normal_Random_Custom(key_context_parms, noise_std_dev, 1000000);

			// At this point, buff should be filled with random normal values from the clipped normal distribution.
			FILE* f_op = open_f("SEAL_normal_rands.txt", "w");
			for (int i_mod = 0; i_mod < (int)key_context_parms.coeff_modulus().size(); i_mod++)
			{
				for (size_t coeff_i = 0; coeff_i < key_context_parms.poly_modulus_degree(); coeff_i++)
				{
					fprintf(f_op, "%d\t%d\t%lu\t%lu\n", i_mod, (int)coeff_i, buff[i_mod * key_context_parms.poly_modulus_degree() + coeff_i], key_context_parms.coeff_modulus().at(i_mod).value());
				} // i loop.
			}
			close_f(f_op, "SEAL_normal_rands.txt");
		} // get_normal_rand
		else if (t_string::compare_strings(argv[1], "-separate_VCF_2_chroms"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "USAGE: %s %s [VCF file path] [chrom id's list path] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}

			char* vcf_fp = argv[2];
			char* chrom_ids_list_fp = argv[3];
			char* op_dir = argv[4];

			vector<char*>* chr_ids = buffer_file(chrom_ids_list_fp);
			fprintf(stderr, "Loaded %d chromosome..\n", (int)chr_ids->size());

			vector<FILE*>* per_chrom_VCF_f_ptrs = new vector<FILE*>();
			for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
			{
				// Normalize the chromosome identifier first.
				normalize_chr_id(chr_ids->at(i_chr));

				char chr_fp[1000];
				sprintf(chr_fp, "%s/%s.vcf.gz", op_dir, chr_ids->at(i_chr));
				FILE* f_vcf = open_f(chr_fp, "w");
				per_chrom_VCF_f_ptrs->push_back(f_vcf);
			} // i_chr loop.

			int n_vars_processed = 0;
			FILE* f_VCF = open_f(vcf_fp, "r");
			while (1)
			{
				char* cur_line = getline(f_VCF);
				if (cur_line == NULL)
				{
					break;
				}

				// Write each header to all files.
				if (cur_line[0] == '#')
				{
					for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
					{
						fprintf(per_chrom_VCF_f_ptrs->at(i_chr), "%s\n", cur_line);
					} // i_chr loop.
					delete[] cur_line;
					continue;
				}

				if (n_vars_processed % 1000 == 0)
				{
					fprintf(stderr, "Processed %d variants..          \r", n_vars_processed);
				}

				n_vars_processed++;

				char cur_var_chrom[100];
				sscanf(cur_line, "%s", cur_var_chrom);
				normalize_chr_id(cur_var_chrom); // make sure to normalize the chromosome identifier from the vcf file.
				int i_chr = t_string::get_i_str(chr_ids, cur_var_chrom);
				if (i_chr < (int)chr_ids->size())
				{
					fprintf(per_chrom_VCF_f_ptrs->at(i_chr), "%s\n", cur_line);
				}

				delete[] cur_line;
			} // VCF reading loop.
			close_f(f_VCF, vcf_fp);
			fprintf(stderr, "Done, closing files..\n");

			for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
			{
				char chr_fp[1000];
				sprintf(chr_fp, "%s/%s.vcf.gz", op_dir, chr_ids->at(i_chr));
				close_f(per_chrom_VCF_f_ptrs->at(i_chr), chr_fp);
			} // i_chr loop.
		} // separate_VCF_2_chroms option.
		else if (t_string::compare_strings(argv[1], "-convert_VCF_2_matbed"))
		{
			// Load VCF and save it in matbed format.
			if (argc != 11)
			{
				fprintf(stderr, "USAGE: %s %s [VCF file path] [VCF sample ids list file path (Use EpiLeak to extract)] \
[Variant regions BED file path] [chromosome id to process] \
[Binary sequence directory (Necessary for ref matching)] [Match reference allele? (0/1)] \
[Match region names? (0/1)] \
[Haplotype specific encoding (0/1)] \
[Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* vcf_fp = argv[2];
			char* vcf_sample_ids_list_fp = argv[3];
			char* var_regions_BED_fp = argv[4];
			char* chr_id_2_process = argv[5];
			char* bin_seq_dir = argv[6];
			bool match_ref_alleles_flag = (argv[7][0] == '1');
			bool match_region_names_flag = (argv[8][0] == '1');
			bool haplotype_specific_encoding = (argv[9][0] == '1');
			char* op_fp = argv[10];

			extract_genotype_signals_per_VCF_memoptimized(vcf_fp,
				vcf_sample_ids_list_fp,
				var_regions_BED_fp,
				chr_id_2_process,
				bin_seq_dir,
				match_ref_alleles_flag,
				match_region_names_flag,
				haplotype_specific_encoding,
				op_fp);
		} // -import_VCF option.
		else if (t_string::compare_strings(argv[1], "-extract_subset_subjects"))
		{
			// Extract variants from the chromosome separated directory.
			if (argc != 6)
			{
				fprintf(stderr, "USAGE: %s %s [matbed file path] [sample ids list file path] [variant regions BED path] [Output matbed path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* geno_sig_matbed_fp = argv[2];
			char* sample_ids_fp = argv[3];
			char* subset_subject_ids_list_fp = argv[4];
			char* op_geno_sig_fp = argv[5];

			extract_genotype_signals_per_subsample_list(geno_sig_matbed_fp, sample_ids_fp, subset_subject_ids_list_fp, op_geno_sig_fp);
		} // -extract_subset_subjects option.
		else if (t_string::compare_strings(argv[1], "-extract_variant_regions"))
		{
			// Extract variants from the chromosome separated directory.
			if (argc != 6)
			{
				fprintf(stderr, "USAGE: %s %s [matbed file path] [sample ids list file path] [variant regions BED path] [Output matbed path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* geno_sig_matbed_fp = argv[2];
			char* sample_ids_fp = argv[3];
			char* var_regs_BED_fp = argv[4];
			char* op_geno_sig_fp = argv[5];

			vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(geno_sig_matbed_fp, sample_ids_fp);
			vector<char*>* sample_ids = buffer_file(sample_ids_fp);
			fprintf(stderr, "Loaded %d subject ids..\n", (int)(sample_ids->size()));
			vector<t_annot_region*>* var_regs = load_BED(var_regs_BED_fp);

			vector<t_annot_region*>* selected_var_geno_regs = intersect_annot_regions(geno_regs, var_regs, false);
			fprintf(stderr, "Found %d genotyped variants out of %d variants.\n", (int)(selected_var_geno_regs->size()), (int)(var_regs->size()));

			fprintf(stderr, "Saving to %s..\n", op_geno_sig_fp);
			binarize_variant_genotype_signal_regions(selected_var_geno_regs, NULL, sample_ids, op_geno_sig_fp);
		} // -extract_variant_regions option.
		else if (t_string::compare_strings(argv[1], "-write_GMMAT_formatted_geno"))
		{
			// Given the variant filtered per chromosome formatted matbed files, write the GMMAT formatted genotypes for using in FEDGWAS.
			// Load all genotype regions from all chromosomes.
			// Write the GMMAT formatted results.
			if (argc != 5)
			{
				fprintf(stderr, "USAGE: %s %s [Geno signal regs files list path] [sample ids list file path] [Output GMMAT genotypes path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* geno_sig_files_list_fp = argv[2];
			char* genotype_sample_ids_fp = argv[3];
			char* geno_GDS_fp = argv[4];

			vector<char*>* geno_sig_files_list = buffer_file(geno_sig_files_list_fp);

			vector<char*>* genotype_sample_ids = buffer_file(genotype_sample_ids_fp);

			// Write the genotypes in GDS formatted file.
			FILE* f_gds_geno = open_f(geno_GDS_fp, "w");

			for(int i_geno_f = 0; i_geno_f < (int)geno_sig_files_list->size(); i_geno_f++)
			{
				// Genotype sample is the standard sample that we take identifiers from.
				char genotype_regs_fp[1000];
				strcpy(genotype_regs_fp, geno_sig_files_list->at(i_geno_f));
				if (!check_file(genotype_regs_fp))
				{
					fprintf(stderr, "Skipping %s..\n", genotype_regs_fp);
					continue;
				}

				fprintf(stderr, "Processing %s...\n", genotype_regs_fp);
				vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(genotype_regs_fp, genotype_sample_ids_fp);
				if (geno_regs->size() == 0)
				{
					fprintf(stderr, "Could not find any variants on %s\n", genotype_regs_fp);
					continue;
				}
				// Don't write a header to the text file. This is not GDS.
				//fprintf(f_gds_geno, "SNP\tREF\tALT");
				//for (int i_s = 0; i_s < genotype_sample_ids->size(); i_s++)
				//{
				//	fprintf(f_gds_geno, "\t%s", genotype_sample_ids->at(i_s));
				//} // i_s loop.

				//fprintf(f_gds_geno, "\n");

				for (int i_var = 0; i_var < (int)geno_regs->size(); i_var++)
				{
					void** cur_var_info = (void**)(geno_regs->at(i_var)->data);
					char* cur_var_geno = (char*)(cur_var_info[0]);

					char* var_name = geno_regs->at(i_var)->name;
					t_string_tokens* toks = t_string::tokenize_by_chars(var_name, "_");
					char ref_all = (toks->at(1)->str()[0]);
					char alt_all = (toks->at(2)->str()[0]);

					fprintf(f_gds_geno, "%s\t%c\t%c", var_name, ref_all, alt_all);
					for (int i_s = 0; i_s < (int)genotype_sample_ids->size(); i_s++)
					{
						fprintf(f_gds_geno, "\t%d", (int)(cur_var_geno[i_s]));
					} // i_s loop.

					fprintf(f_gds_geno, "\n");

					// Free memory.
					delete[] cur_var_geno;
				} // i_var loop.

				// Free genotype signal memory.
				delete_annot_regions(geno_regs);
			} // chr_ids list.
			close_f(f_gds_geno, geno_GDS_fp);
		} // write_GMMAT_formatted_geno option.
		else if (t_string::compare_strings(argv[1], "-impute_missing_GMMAT_formatted_geno"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "USAGE: %s %s [GMMAT with missing genotypes path] [Output GMMAT genotypes path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* GMMAT_fp = argv[2];
			char* imp_GMMAT_fp = argv[3];

			int n_processed_vars = 0;

			double* cur_var_geno_sig = NULL;;
			int sample_size = 0;
			char col_buff[1000];
			FILE* f_gmmat = open_f(GMMAT_fp, "r");
			FILE* f_imp_gmmat = open_f(imp_GMMAT_fp, "w");
			FILE* f_missing_rates = open_f("missingness.txt", "w");
			while (1)
			{
				char* cur_gmmat_line = getline(f_gmmat);
				if (cur_gmmat_line == NULL)
				{
					break;
				}

				// Initialize the sample size and genotypes.
				if (sample_size == 0)
				{
					t_string_tokens* toks = t_string::tokenize_by_chars(cur_gmmat_line, "\t");

					sample_size = toks->size() - 3;
					cur_var_geno_sig = new double[sample_size];
					fprintf(stderr, "Set sample size to %d\n", sample_size);

					t_string::clean_tokens(toks);
				}

				int i_cur_char = 0;
				t_string::get_next_token(cur_gmmat_line, col_buff, 1000, "\t", i_cur_char);
				char* var_id = t_string::copy_me_str(col_buff);
				t_string::get_next_token(cur_gmmat_line, col_buff, 1000, "\t", i_cur_char);
				char* ref_all = t_string::copy_me_str(col_buff);
				t_string::get_next_token(cur_gmmat_line, col_buff, 1000, "\t", i_cur_char);
				char* alt_all = t_string::copy_me_str(col_buff);

				double tot_sig = 0;
				int n_valid_subjects = 0;
				for (int i_s = 0; i_s < sample_size; i_s++)
				{
					t_string::get_next_token(cur_gmmat_line, col_buff, 1000, "\t", i_cur_char);
					double cur_geno = atof(col_buff);
					cur_var_geno_sig[i_s] = cur_geno;

					if (cur_geno == 0 ||
						cur_geno == 1 ||
						cur_geno == 2)
					{
						tot_sig += cur_geno;
						n_valid_subjects += 1;
					}
				} // i loop.

				int mean_geno = (int)(floor(tot_sig / n_valid_subjects));

				// Impute the missing or invalid genotypes to mean genotype..
				fprintf(f_imp_gmmat, "%s\t%s\t%s", var_id, ref_all, alt_all);
				int n_corrected = 0;
				for (int i_s = 0; i_s < sample_size; i_s++)
				{
					// Write valid genotypes as they are, fix the invalid/missing ones.
					if (cur_var_geno_sig[i_s] == 0 ||
						cur_var_geno_sig[i_s] == 1 ||
						cur_var_geno_sig[i_s] == 2)
					{
						fprintf(f_imp_gmmat, "\t%d", (int)cur_var_geno_sig[i_s]);
					}
					else
					{
						fprintf(f_imp_gmmat, "\t%d", mean_geno);

						//fprintf(stderr, "%s: %d: %d -> %d\n", toks->at(0)->str(), i_s, cur_var_geno_sig[i_s], mean_geno);
						n_corrected++;
					}
				} // i_s loop.

				fprintf(f_imp_gmmat, "\n");

				if (n_processed_vars % 1000 == 0)
				{
					fprintf(stderr, "%s: Mean: %d; # Corrected: %d           \r", var_id, mean_geno, n_corrected);
				}

				fprintf(f_missing_rates, "%s\t%d\t%d\n", var_id, n_corrected, sample_size);

				n_processed_vars++;

				// Free memory.
				delete[] cur_gmmat_line;
				delete[] var_id;
				delete[] ref_all;
				delete[] alt_all;
			} // GMMAT reading loop.

			close_f(f_missing_rates, "missingness.txt");
			close_f(f_gmmat, GMMAT_fp);
			close_f(f_imp_gmmat, imp_GMMAT_fp);

			exit(0);
		} // -impute_missing_GMMAT_formatted_geno option.
		else if (t_string::compare_strings(argv[1], "-extract_rows_per_query_column_preserve_query_order"))
		{
			if (argc != 7)
			{
				fprintf(stderr, "USAGE: %s %s [Query column values path] [Multicolumn file path] [Column index (number) to match] [Write all matching rows? (0/1)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* queries_col_fp = argv[2];
			char* multicol_fp = argv[3];
			int col_i = atoi(argv[4]);
			bool write_all_matches = argv[5][0] == '1';
			char* op_fp = argv[6];
			extract_rows_per_query_column_preserve_query_order(queries_col_fp, col_i, multicol_fp, write_all_matches, op_fp);

			exit(0);
		}
		else if (t_string::compare_strings(argv[1], "-extract_subsample_GMMAT_genotypes_from_GMMAT_genotype_matrix"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [GMMAT genotype matrix] [GMMAT genotype matrix sample ids list] [Subsample ids list] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* gmmat_geno_fp = argv[2];
			char* gmmat_sample_ids_fp = argv[3];
			char* subsample_ids_fp = argv[4];
			char* op_fp = argv[5];

			vector<char*>* gmmat_sample_ids = buffer_file(gmmat_sample_ids_fp);
			vector<char*>* subsample_ids = buffer_file(subsample_ids_fp);
			fprintf(stderr, "Loaded %d GMMAT sample id's and %d subsample id's, mapping them..\n", (int)gmmat_sample_ids->size(), (int)subsample_ids->size());

			vector<int>* per_subsample_gmmat_sample_id = new vector<int>();
			for (int i_sub = 0; i_sub < (int)subsample_ids->size(); i_sub++)
			{
				int subsample_gmmat_sample_i = t_string::get_i_str(gmmat_sample_ids, subsample_ids->at(i_sub));
				if (subsample_gmmat_sample_i == (int)gmmat_sample_ids->size())
				{
					fprintf(stderr, "Could not find %s in GMMAT sample, make sure all subjects are in GMMAT sample ids in %s\n", gmmat_sample_ids_fp, gmmat_sample_ids_fp);
				}
				else
				{
					per_subsample_gmmat_sample_id->push_back(subsample_gmmat_sample_i);
				}
			} // i_sub loop.

			FILE* f_gmmat = open_f(gmmat_geno_fp, "r");
			FILE* f_op = open_f(op_fp, "w");
			while (1)
			{
				char* cur_line = getline(f_gmmat);
				if (cur_line == NULL)
				{
					break;
				}

				t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");
				fprintf(f_op, "%s\t%s\t%s", toks->at(0)->str(), toks->at(1)->str(), toks->at(2)->str());

				for (int i_sub = 0; i_sub < (int)per_subsample_gmmat_sample_id->size(); i_sub++)
				{
					if (per_subsample_gmmat_sample_id->at(i_sub) == (int)gmmat_sample_ids->size())
					{
						fprintf(stderr, "Sanity check failed: %d. subsample id (%s) is not in GMMAT sample ids.\n", i_sub, subsample_ids->at(i_sub));
						exit(0);
					}

					int i_sub_tok_i = per_subsample_gmmat_sample_id->at(i_sub) + 3;
					fprintf(f_op, "\t%s", toks->at(i_sub_tok_i)->str());
				} // i_sub loop.

				fprintf(f_op, "\n");

				// Free memory.
				t_string::clean_tokens(toks);
				delete[] cur_line;
			} // file reading loop.
			close_f(f_op, op_fp);
			close_f(f_gmmat, gmmat_geno_fp);
		} // -extract_subsample_GMMAT_genotypes_from_GMMAT_genotype_matrix option
		else if (t_string::compare_strings(argv[1], "-test_linux_srandom_seeded"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Rand type: Uniform/Gaussian] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* rand_type = argv[2];
			char* op_fp = argv[3];

			long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
			t_rng* rng = new t_rng(rand_seed);

			fprintf(stderr, "Seed is %lu\n", rand_seed);

			FILE* f_op = open_f(op_fp, "w");
			if (t_string::compare_strings(rand_type, "Gaussian"))
			{
				fprintf(stderr, "Generating Gaussian random values..\n");
				for (int i = 0; i < 10000; i++)
				{
					double cur_rand = rng->random_gaussian_double_ran3();
					fprintf(f_op, "%.7f\n", cur_rand);
				} // i loop.
			} // Gaussian check.
			else if (t_string::compare_strings(rand_type, "Uniform"))
			{
				fprintf(stderr, "Generating Uniform random values..\n");
				for (int i = 0; i < 10000; i++)
				{
					double cur_rand = rng->random_double_ran3();
					fprintf(f_op, "%.7f\n", cur_rand);
				} // i loop.
			} // UNiform check.
			close_f(f_op, op_fp);
		} // -test_linux_srandom_seeded option.
		else if (t_string::compare_strings(argv[1], "-encrypt_encoded_pt_matrix"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Encoded plaintext matrix file] [Text params file] [Pooled public key file] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encoded_matrix_pt_path = argv[2];
			char* text_params_path = argv[3];
			char* pooled_public_key_path = argv[4];
			char* encrypted_data_path = argv[5];

			encrypt_encoded_pt_matrix(encoded_matrix_pt_path,
				text_params_path,
				pooled_public_key_path,
				encrypted_data_path);
		} // -encrypt_encoded_pt_matrix option.
		else if (t_string::compare_strings(argv[1], "-validate_ckks_text_params"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Text params file] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* text_params_fp = argv[2];
			char* op_fp = argv[3];

			if (validate_coeff_modulus_bit_size(text_params_fp))
			{
				fprintf(stderr, "Parameters seem to be VALID.\n");
				FILE* f_op = open_f(op_fp, "w");
				fprintf(f_op, "VALID\n");
				close_f(f_op, op_fp);
			}
			else
			{
				fprintf(stderr, "Parameters seem to be INVALID.\n");
				FILE* f_op = open_f(op_fp, "w");
				fprintf(f_op, "INVALID\n");
				close_f(f_op, op_fp);
			}

			exit(0);
		} // -validate_ckks_text_params option.
		else if (t_string::compare_strings(argv[1], "-transpose_continuous_encrypted_vector"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Encrypted Matrix file] [Text params file] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_mat_fp = argv[2];
			char* text_params_fp = argv[3];
			char* op_fp = argv[4];

			transpose_continuous_encrypted_vector(enc_mat_fp, text_params_fp, op_fp);
		} // -transpose_continuous_encrypted_vector option.
		else if (t_string::compare_strings(argv[1], "-scalar_multiply_matrix_plain"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Scalar factor] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* pt_mat_fp = argv[2];
			double scalar_factor = atof(argv[3]);
			char* op_fp = argv[4];

			int loaded_nrows, loaded_ncols;
			double** matrix = load_matrix_binary(pt_mat_fp, loaded_nrows, loaded_ncols);

			for (int i_row = 0; i_row < loaded_nrows; i_row++)
			{
				for (int i_col = 0; i_col < loaded_ncols; i_col++)
				{
					matrix[i_row][i_col] *= scalar_factor;
				} // i_col loop. 
			} // i_row loop.

			save_matrix_binary(matrix, loaded_nrows, loaded_ncols, op_fp);
		} // -write_enc_matrix_dimensions option.
		else if (t_string::compare_strings(argv[1], "-write_enc_matrix_dimensions"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_mat_fp = argv[2];
			char* op_fp = argv[3];

			write_enc_matrix_dimensions(enc_mat_fp, op_fp);
		} // -write_enc_matrix_dimensions option.
		else if (t_string::compare_strings(argv[1], "-write_plain_matrix_dimensions"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* plain_mat_fp = argv[2];
			char* op_fp = argv[3];

			write_plain_matrix_dimensions(plain_mat_fp, op_fp);
		} // -write_plain_matrix_dimensions option.
		else if (t_string::compare_strings(argv[1], "-plain_pad_matrix_to_next_power_of_2"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			char* op_fp = argv[3];

			pad_matrix_sizee_to_next_power_of_2(matrix_fp, op_fp);
		} // -plain_pad_matrix_to_next_power_of_2 option.
		else if (t_string::compare_strings(argv[1], "-plain_pad_matrix_rows_to_next_power_of_2"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			char* op_fp = argv[3];

			pad_matrix_rows_to_next_power_of_2(matrix_fp, op_fp);
		} // -plain_pad_matrix_rows_to_next_power_of_2 option.	
		else if (t_string::compare_strings(argv[1], "-plain_pad_matrix_cols_to_next_power_of_2"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			char* op_fp = argv[3];

			pad_matrix_cols_to_next_power_of_2(matrix_fp, op_fp);
		} // -plain_pad_matrix_cols_to_next_power_of_2 option.	
		else if (t_string::compare_strings(argv[1], "-plain_unpad_matrix_to_size"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Matrix path list] [New # rows] [New # cols] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			int new_n_row = atoi(argv[3]);
			int new_n_col = atoi(argv[4]);
			char* op_fp = argv[5];

			plain_unpad_matrix_to_size(matrix_fp, new_n_row, new_n_col, op_fp);
		} // -plain_unpad_matrix_to_size option.
		else if (t_string::compare_strings(argv[1], "-plain_invert_matrix"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			char* op_fp = argv[3];

			plain_invert_matrix(matrix_fp, op_fp);
		} // -plain_invert_matrix option.
		else if (t_string::compare_strings(argv[1], "-plain_add_matrices"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* A_mat_fp = argv[2];
			char* B_mat_fp = argv[3];
			char* op_fp = argv[4];

			int loaded_nArow, loaded_nAcol;
			int loaded_nBrow, loaded_nBcol;

			double** A_mat = load_matrix_binary(A_mat_fp, loaded_nArow, loaded_nAcol, NULL);
			double** B_mat = load_matrix_binary(B_mat_fp, loaded_nBrow, loaded_nBcol, NULL);

			double** sum_mat = matrix_add(A_mat, loaded_nArow, loaded_nAcol, B_mat, loaded_nBrow, loaded_nBcol, NULL);

			save_matrix_binary(sum_mat, loaded_nArow, loaded_nAcol, op_fp);
		} // -add_matrices_per_list option.
		else if (t_string::compare_strings(argv[1], "-plain_add_matrices_per_list"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [Matrix path list] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_path_list_fp = argv[2];
			char* op_fp = argv[3];

			plain_add_matrices_per_list(matrix_path_list_fp, op_fp);
		} // -add_matrices_per_list option.
		else if (t_string::compare_strings(argv[1], "-secure_add_cont_ct_matrices_per_list"))
		{
			if (argc != 9)
			{
				fprintf(stderr, "%s %s [Continuous encrypted matrix path list file] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_mat_list_fp = argv[2];
			char* text_params_path = argv[3];
			char* pooled_public_key_path = argv[4];
			char* pooled_relin_key_path = argv[5];
			char* pooled_galois_key_path = argv[6];
			char* pooled_private_key_path = argv[7];
			char* op_fp = argv[8];

			secure_add_continuous_encrypted_matrices_per_list(enc_mat_list_fp,
				text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_fp);
		} // -secure_add_cont_ct_matrices_per_list option.
		else if (t_string::compare_strings(argv[1], "-secure_add_cont_ct_matrices"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [Matrix A file] [Matrix B file] [Text parameters file path] \
[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_A_fp = argv[2];
			char* enc_B_fp = argv[3];
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_fp = argv[9];

			secure_add_continuous_encrypted_matrices(enc_A_fp, enc_B_fp,
				text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_fp);
		}
		else if (t_string::compare_strings(argv[1], "-secure_sub_cont_ct_matrices"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [Matrix A col expansion dir.] [Matrix B row expansion dir.] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_A_fp = argv[2];
			char* enc_B_fp = argv[3];
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_fp = argv[9];

			secure_subtract_continuous_encrypted_matrices(enc_A_fp, enc_B_fp,
				text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_fp);
		}
		else if (t_string::compare_strings(argv[1], "-secure_elementwise_mul_cont_ct_matrices"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [Matrix A path] [Matrix B path] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_A_fp = argv[2];
			char* enc_B_fp = argv[3];
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_fp = argv[9];

			secure_multiply_elementwise_continuous_encrypted_matrices(enc_A_fp, enc_B_fp,
				text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_fp);
		}
		else if (t_string::compare_strings(argv[1], "-generate_plaintext_mask_per_continuous_encrypted_data"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Matrix A col expansion dir.] [Mask variance] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_A_fp = argv[2];
			double mask_variance = atof(argv[3]);
			char* plain_mask_matrix_op_fp = argv[4];

			generate_mask_matrix_per_data_matrix(enc_A_fp,
				mask_variance,
				plain_mask_matrix_op_fp);
		} // -generate_plaintext_mask_per_continuous_encrypted_data option.
		else if (t_string::compare_strings(argv[1], "-additive_mask_continuous_encrypted_data"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Cont. encrypted matrix file] [Cont. encrypted mask file] [Text parameters file path] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encrypted_data_matrix_path = argv[2];
			char* encrypted_mask_matrix_fp = argv[3];
			char* text_params_path = argv[4];
			char* masked_encrypted_data_matrix_fp = argv[5];

			bool add_mask_flag = true;

			add_remove_encrypted_mask_2_continuous_encrypted_data_matrix(encrypted_data_matrix_path, encrypted_mask_matrix_fp,
				add_mask_flag,
				text_params_path,
				masked_encrypted_data_matrix_fp);
		} // -additive_mask_continuous_encrypted_data option
		else if (t_string::compare_strings(argv[1], "-additive_unmask_continuous_encrypted_data"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Cont. encrypted matrix file] [Cont. encrypted mask file] [Text parameters file path] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encrypted_data_matrix_path = argv[2];
			char* encrypted_mask_matrix_fp = argv[3];
			char* text_params_path = argv[4];
			char* masked_encrypted_data_matrix_fp = argv[5];

			bool add_mask_flag = false;

			add_remove_encrypted_mask_2_continuous_encrypted_data_matrix(encrypted_data_matrix_path, encrypted_mask_matrix_fp,
				add_mask_flag,
				text_params_path,
				masked_encrypted_data_matrix_fp);
		} // -additive_unmask_continuous_encrypted_data option.
		else if (t_string::compare_strings(argv[1], "-write_continuous_encrypted_ciphertext_vital_stats"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Data matrix path] [Text params path] [Output path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_mat_fp = argv[2];
			char* text_params_path = argv[3];
			char* vital_stats_fp = argv[4];

			write_vital_stats_per_continuous_encrypted_matrix(enc_mat_fp,
				text_params_path,
				vital_stats_fp);
		} // --write_continuous_encrypted_ciphertext_vital_stats option.
		else if (t_string::compare_strings(argv[1], "-continuous_encrypt_data_matrix"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Data matrix path] [Text params path] [Public key path] [Output path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* matrix_fp = argv[2];
			char* text_params_path = argv[3];
			char* public_key_path = argv[4];
			char* op_fp = argv[5];

			int loaded_nrow, loaded_ncol;
			double** data_matrix = load_matrix_binary(matrix_fp, loaded_nrow, loaded_ncol);

			encrypt_plaintext_matrix_continuous_ct(data_matrix, loaded_nrow, loaded_ncol, text_params_path, public_key_path, op_fp);
		} // -continuous_encrypt_data_matrix
		else if (t_string::compare_strings(argv[1], "-partial_decrypt_continuous_enc_per_noisy_secretkey"))
		{
			if (argc != 8)
			{
				fprintf(stderr, "%s %s [Continuous Encrypted data path] [is site 0? (0/1)] [Text params path] [Noisy secret key path] [Smdg-ing noise bit size (Default is 40-bits for variance)] [Partial decrypted data output path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encrypted_data_path = argv[2];
			bool is_site0 = (argv[3][0] == '1');
			char* text_params_path = argv[4];
			char* noisy_key_path = argv[5];
			int n_bits_per_max_smdging_noise = atoi(argv[6]);
			char* partial_decrypt_data_path = argv[7];

			partial_decrypt_continuous_enc_per_noisy_secretkey_w_smdgng_noise(encrypted_data_path, is_site0, text_params_path, noisy_key_path, n_bits_per_max_smdging_noise, partial_decrypt_data_path);
		} // -partial_decrypt_per_noisy_secretkey option.
		else if (t_string::compare_strings(argv[1], "-pool_partial_decrypted_continuous_enc_data"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "USAGE: %s %s [Partial decrypted continuous encrypted matrix paths list path] [Text params path] [Output path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* partial_decrypted_paths_list_path = argv[2];
			char* text_params_path = argv[3];
			char* op_fp = argv[4];

			collaborative_pool_partial_decrypted_continuous_enc_plaintext_results(partial_decrypted_paths_list_path, text_params_path, op_fp);
		} // -pool_partial_decrypted_data option.
		else if (t_string::compare_strings(argv[1], "-secure_multiply_matrices_Acol_Brow_expansions"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [Matrix A col expansion dir.] [Matrix B row expansion dir.] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* A_dir = argv[2];
			char* B_dir = argv[3];
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_fp = argv[9];

			secure_multiply_matrices_Acol_Brow_expansions(A_dir, B_dir,
				text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_fp);
		} // -secure_multiply_matrices_Acol_Brow_expansions option.
		else if (t_string::compare_strings(argv[1], "-fully_decrypt_continuous_encrypted_matrix"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Matrix path] [Text parameters file path] [Pooled Secret Key file path] [Encrypted output path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* enc_matrix_path = argv[2];
			char* text_params_path = argv[3];
			char* pooled_secret_key_path = argv[4];
			char* full_decrypted_matrix_fp = argv[5];

			fully_decrypt_continuous_encrypted_matrix(enc_matrix_path, text_params_path, pooled_secret_key_path, full_decrypted_matrix_fp);
		} // fully_decrypt_continuous_encrypted_matrix option.
		else if (t_string::compare_strings(argv[1], "-col_expand_dense_encrypt_matrix"))
		{
			if (argc != 7)
			{
				fprintf(stderr, "%s %s [Matrix A path] [# of column replications] [Text parameters file path] [Pooled public key path] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}
			char* matA_fp = argv[2];
			int n_cols_per_rep = atoi(argv[3]);
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* op_dir = argv[6];

			col_expand_dense_encrypt_matrix(matA_fp, n_cols_per_rep, text_params_path, pooled_public_key_path, op_dir);
		} // -col_expand_dense_encrypt_matrix option.
		else if (t_string::compare_strings(argv[1], "-row_expand_dense_encrypt_matrix"))
		{
			if (argc != 7)
			{
				fprintf(stderr, "%s %s [Matrix A path] [# of row replications] [Text parameters file path] [Pooled public key path] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}
			char* matA_fp = argv[2];
			int n_rows_per_rep = atoi(argv[3]);
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* op_dir = argv[6];

			row_expand_dense_encrypt_matrix(matA_fp, n_rows_per_rep, text_params_path, pooled_public_key_path, op_dir);
		} // -row_expand_dense_encrypt_matrix option.
		else if (t_string::compare_strings(argv[1], "-row2row_multiply_pt"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Matrix A path] [Matrix B path] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}
			char* matA_fp = argv[2];
			char* matB_fp = argv[3];
			char* op_fp = argv[4];

			int loaded_nrowA, loaded_ncolA;
			double** matA = load_matrix_binary(matA_fp, loaded_nrowA, loaded_ncolA, NULL);
			int loaded_nrowB, loaded_ncolB;
			double** matB = load_matrix_binary(matB_fp, loaded_nrowB, loaded_ncolB, NULL);

			if (loaded_nrowA != loaded_nrowB ||
				loaded_ncolA != loaded_ncolB)
			{
				fprintf(stderr, "Dimensions not conformant for row2row inner multiplication: %d, %d; %d, %d\n", loaded_nrowA, loaded_ncolA, loaded_nrowB, loaded_ncolB);
				exit(0);
			}

			fprintf(stderr, "Row2Row multiplying A[%dx%d]xB[%dx%d] and saving to %s\n", loaded_nrowA, loaded_ncolA, loaded_nrowB, loaded_ncolB, op_fp);

			double** mat_AB = matrix_row_by_row_inner_product(matA, loaded_nrowA, loaded_ncolA, matB, loaded_nrowB, loaded_ncolB, NULL);
			save_matrix_binary(mat_AB, loaded_nrowA, 1, op_fp);
		} // -multiply_matrices_pt option.	
		else if (t_string::compare_strings(argv[1], "-multiply_matrices_pt"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Matrix A path] [Matrix B path] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}
			char* matA_fp = argv[2];
			char* matB_fp = argv[3];
			char* op_fp = argv[4];

			int loaded_nrowA, loaded_ncolA;
			double** matA = load_matrix_binary(matA_fp, loaded_nrowA, loaded_ncolA, NULL);
			int loaded_nrowB, loaded_ncolB;
			double** matB = load_matrix_binary(matB_fp, loaded_nrowB, loaded_ncolB, NULL);

			if (loaded_ncolA != loaded_nrowB)
			{
				fprintf(stderr, "Dimensions not conformant for multiplication: %d, %d; %d, %d\n", loaded_nrowA, loaded_ncolA, loaded_nrowB, loaded_ncolB);
				exit(0);
			}

			fprintf(stderr, "Calculating A[%dx%d]xB[%dx%d] and saving to %s\n", loaded_nrowA, loaded_ncolA, loaded_nrowB, loaded_ncolB, op_fp);

			double** mat_AB = matrix_multiply(matA, loaded_nrowA, loaded_ncolA, matB, loaded_nrowB, loaded_ncolB, NULL);
			save_matrix_binary(mat_AB, loaded_nrowA, loaded_ncolB, op_fp);
		} // -multiply_matrices_pt option.
		else if (t_string::compare_strings(argv[1], "-multiply_matrices_elementwise_pt"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [Matrix A path] [Matrix B path] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}
			char* matA_fp = argv[2];
			char* matB_fp = argv[3];
			char* op_fp = argv[4];

			int loaded_nrowA, loaded_ncolA;
			double** matA = load_matrix_binary(matA_fp, loaded_nrowA, loaded_ncolA, NULL);
			int loaded_nrowB, loaded_ncolB;
			double** matB = load_matrix_binary(matB_fp, loaded_nrowB, loaded_ncolB, NULL);

			if (loaded_ncolA != loaded_ncolB ||
				loaded_nrowA != loaded_nrowB)
			{
				fprintf(stderr, "Dimensions not conformant for multiplication: %d, %d; %d, %d\n", loaded_nrowA, loaded_ncolA, loaded_nrowB, loaded_ncolB);
				exit(0);
			}

			double** res_matrix = matrix_multiply_elementwise(matA, loaded_nrowA, loaded_ncolA, matB, loaded_nrowB, loaded_ncolB, NULL);
			save_matrix_binary(res_matrix, loaded_nrowA, loaded_ncolB, op_fp);
		} // -multiply_matrices_elementwise_pt option.
		else if (t_string::compare_strings(argv[1], "-row_expand_continuous_encrypted_matrix"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [Matrix A path] [# rows in expanded matrix] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encA_mat_fp = argv[2];
			int n_rows_per_expanded_matrix = atoi(argv[3]);
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_dir = argv[9];

			row_expand_continuous_encrypted_matrix(encA_mat_fp, n_rows_per_expanded_matrix, text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_dir);
		} // -row_expand_continuous_encrypted_matrix option.
		else if (t_string::compare_strings(argv[1], "-secure_row2row_inner_prod_continuous_encrypted_matrices"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "USAGE: %s %s [Matrix A path] [Matrix B path] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encA_mat_fp = argv[2];
			char* encB_mat_fp = argv[3];
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_fp = argv[9];

			secure_row2row_inner_prod_continuous_encrypted_matrices(encA_mat_fp, encB_mat_fp,
				text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path, // To be removed.
				op_fp);
		} // -secure_multiply_matrices option.
		else if (t_string::compare_strings(argv[1], "-secure_multiply_matrices"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [Matrix A path] [Matrix Bt path] [Text parameters file path] \
	[Pooled public key path] [Pooled relinearization key path] [Pooled Galois key path] [Pooled SK path (To be removed...)] [Output file path]\n", argv[0], argv[1]);
				exit(0);
			}

			char* encA_mat_fp = argv[2];
			char* encBt_mat_fp = argv[3];
			char* text_params_path = argv[4];
			char* pooled_public_key_path = argv[5];
			char* pooled_relin_key_path = argv[6];
			char* pooled_galois_key_path = argv[7];
			char* pooled_private_key_path = argv[8];
			char* op_fp = argv[9];

			secure_multiply_matrices(encA_mat_fp, encBt_mat_fp, text_params_path,
				pooled_public_key_path,
				pooled_relin_key_path,
				pooled_galois_key_path,
				pooled_private_key_path,
				op_fp);
		} // -secure_multiply_matrices option.
		else if (t_string::compare_strings(argv[1], "-generate_random_pt_matrix"))
		{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [# rows] [# cols] [Output path]\n", argv[0], argv[1]);
			exit(0);
		}

		int nrow = atoi(argv[2]);
		int ncol = atoi(argv[3]);
		char* op_fp = argv[4];

		fprintf(stderr, "Generating random %dx%d matrix and saving to %s.\n", nrow, ncol, op_fp);

		long int rand_seed = (long)(t_seed_manager::seed_me_getrandom());
		t_rng* rng = new t_rng(rand_seed);

		double** matrix = allocate_matrix(nrow, ncol);
		for (int row = 0; row < nrow; row++)
		{
			for (int col = 0; col < ncol; col++)
			{
				matrix[row][col] = rng->random_double_ran3() * 5;
			} // col loop.
		} // row loop.

		save_matrix_binary(matrix, nrow, ncol, op_fp);
		} // -generate_random_pt_matrix option.
		else if (t_string::compare_strings(argv[1], "-transpose_pt_matrix"))
		{
		if (argc != 4)
		{
			fprintf(stderr, "%s %s [Matrix file path] [Output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* matrix_fp = argv[2];
		char* op_fp = argv[3];

		int loaded_nrow, loaded_ncol;
		double** matrix = load_matrix_binary(matrix_fp, loaded_nrow, loaded_ncol);

		fprintf(stderr, "Loaded %dx%d matrix, transposing and saving to %s..\n", loaded_nrow, loaded_ncol, op_fp);

		double** trans_matrix = transpose_matrix(matrix, loaded_nrow, loaded_ncol, NULL);

		save_matrix_binary(trans_matrix, loaded_ncol, loaded_nrow, op_fp);
		} // -transpose_pt_matrix option.
		else if (t_string::compare_strings(argv[1], "-generate_per_site_noisy_keys"))
		{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [# sites] [Text parameters file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		int n_sites = atoi(argv[2]);
		char* text_params_path = argv[3];
		char* op_dir = argv[4];

		// This is now implemented in utils.
		generate_share_private_keys(n_sites,
			text_params_path,
			true, // Do key testing.
			op_dir);
		} // -generate_per_site_noisy_keys option.
		else if (t_string::compare_strings(argv[1], "-generate_mult_diagonal_noise_matrix"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [# rows/cols] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			int l_diag = atoi(argv[2]);
			char* op_fp = argv[3];

			// Generate the multiplicative noise that will be used to calculate inv(XtWX).
			double** mult_noise = generate_multiplicative_diagonal_noise_matrix(l_diag, l_diag);
			save_matrix_binary(mult_noise, l_diag, l_diag, op_fp);
		} // -generate_mult_full_noise_matrix option.
		else if (t_string::compare_strings(argv[1], "-generate_constant_diagonal_matrix"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [# rows] [# cols] [Diagonal value] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			int n_row = atoi(argv[2]);
			int n_col = atoi(argv[3]);
			double diagonal_val = atof(argv[4]);
			char* op_fp = argv[5];

			// Generate the multiplicative noise that will be used to calculate inv(XtWX).
			double** mult_noise = generate_constant_diagonal_matrix(n_row, n_col, diagonal_val);
			save_matrix_binary(mult_noise, n_row, n_col, op_fp);
		} // -generate_constant_diagonal_matrix option.
		else if (t_string::compare_strings(argv[1], "-generate_constant_full_matrix"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [# rows] [# cols] [Init value] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			int n_row = atoi(argv[2]);
			int n_col = atoi(argv[3]);
			double init_val = atof(argv[4]);
			char* op_fp = argv[5];

			// Generate the multiplicative noise that will be used to calculate inv(XtWX).
			double** const_matrix = allocate_matrix(n_row, n_col, init_val);
			save_matrix_binary(const_matrix, n_row, n_col, op_fp);
		} // -generate_constant_diagonal_matrix option.

		else if (t_string::compare_strings(argv[1], "-generate_mult_full_noise_matrix"))
		{
			if (argc != 5)
			{
				fprintf(stderr, "%s %s [# rows] [# cols] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			int n_rows = atoi(argv[2]);
			int n_cols = atoi(argv[3]);
			char* op_fp = argv[4];

			// Generate the multiplicative noise that will be used to calculate inv(XtWX).
			double** mult_noise = generate_multiplicative_full_noise_matrix(n_rows, n_cols);
			save_matrix_binary(mult_noise, n_rows, n_cols, op_fp);
		} // -generate_mult_full_noise_matrix option.
		else if (t_string::compare_strings(argv[1], "-save_matrix_plain_2_bin"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [Plain file path] [Have row ids? (First column)] [Have column ids? (First row)] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* plain_matrix_fp = argv[2];
			bool have_row_ids = (argv[3][0] == '1');
			bool have_col_ids = (argv[4][0] == '1');
			char* op_fp = (argv[5]);

			int nrow = 0;
			int ncol = 0;
			double** bin_matrix = load_matrix_plain(plain_matrix_fp, nrow, ncol, have_row_ids, have_col_ids);

			fprintf(stderr, "Loaded %dx%d matrix, saving to %s\n", nrow, ncol, op_fp);
			save_matrix_binary(bin_matrix, nrow, ncol, op_fp);

			exit(0);
		} // -save_matrix_plain_2_bin option.
		else if (t_string::compare_strings(argv[1], "-assign_chisqr_pvals_per_ST_stats"))
		{
			if (argc != 4)
			{
				fprintf(stderr, "%s %s [S/T stats file path] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* chisqr_stats_fp = (argv[2]);
			char* op_fp = (argv[3]);

			vector<char*>* chisqr_stat_vs_scale_lines = buffer_file(chisqr_stats_fp);

			FILE* f_op = open_f(op_fp, "w");
			for (int i_var = 0; i_var < (int)chisqr_stat_vs_scale_lines->size(); i_var++)
			{
				double stat = 0;
				double scale = 0;
				sscanf(chisqr_stat_vs_scale_lines->at(i_var), "%lf %lf", &stat, &scale);
				double norm_chisqr_scale = stat / scale;
				double pval = get_1DOF_chisqr_pval(norm_chisqr_scale);

				fprintf(f_op, "%s\t%.4f\n", chisqr_stat_vs_scale_lines->at(i_var), pval);
			} // i_var loop.

			close_f(f_op, op_fp);
		} // -assign_chisqr_pvals_per_ST_stats option.
		else if (t_string::compare_strings(argv[1], "-diag_matrix_multiply_checks"))
		{
			double** A = allocate_matrix(5, 5);
			double** B = allocate_matrix(5, 5);
			for (int row = 0; row < 5; row++)
			{
				for (int col = 0; col < 5; col++)
				{
					A[row][col] = row + 1;
				}

				B[row][row] = row + 1;
			}

			fprintf(stderr, "Right multiplication:\n");
			save_matrix_plain(matrix_multiply(A, 5, 5, B, 5, 5, NULL), 5, 5, "stdout");
			save_matrix_plain(matrix_right_multiply_with_diag_matrix(A, 5, 5, B, 5, 5, NULL), 5, 5, "stdout");
			save_matrix_plain(matrix_right_multiply_with_diag_matrix_as_diagonal_only(A, 5, 5, get_diagonal_of_matrix(B, 5, 5, NULL), 5, 1, NULL), 5, 5, "stdout");

			fprintf(stderr, "Left multiplication:\n");
			save_matrix_plain(matrix_multiply(B, 5, 5, A, 5, 5, NULL), 5, 5, "stdout");
			save_matrix_plain(matrix_left_multiply_with_diag_matrix(B, 5, 5, A, 5, 5, NULL), 5, 5, "stdout");
			save_matrix_plain(matrix_left_multiply_with_diag_matrix_as_diagonal_only(get_diagonal_of_matrix(B, 5,5, NULL), 5, 1, A, 5, 5, NULL), 5, 5, "stdout");
		}
		else if (t_string::compare_strings(argv[1], "-dump_matrix_plain")) 
		{
			if (argc != 4)
			{
				fprintf(stderr, "USAGE: %s %s [binary matrix file] [output plain matrix file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* bin_mat_fp = argv[2];
			char* plain_mat_fp = argv[3];

			int nrow, ncol;
			double** mat = load_matrix_binary(bin_mat_fp, nrow, ncol, NULL);
			save_matrix_plain(mat, nrow, ncol, plain_mat_fp);
		} // -dump_matrix_plain
		else if (t_string::compare_strings(argv[1], "-cryptable_client_calculate_save_XtWX_XtWz"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
	[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Private working directory] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* private_working_dir = argv[10];

			char* shared_working_dir = argv[11];

			cryptable_client_calculate_save_XtWX_XtWz(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -client_calculate_save_XtWX_XtWz option.

		else if (t_string::compare_strings(argv[1], "-cryptable_client_add_mult_noise_2_XtWX"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
	[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Private working directory] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* private_working_dir = argv[10];

			char* shared_working_dir = argv[11];

			cryptable_client_add_mult_noise_2_XtWX(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_add_mult_noise_2_XtWX option.
		else if (t_string::compare_strings(argv[1], "-cryptable_client_pool_site_specific_all_site_noise_XtWX"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
	[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Private working directory] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* private_working_dir = argv[10];

			char* shared_working_dir = argv[11];

			cryptable_client_pool_site_specific_all_site_noise_XtWX(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_pool_site_specific_all_site_noise_XtWX option.
		else if (t_string::compare_strings(argv[1], "-cryptable_collaborative_decrypt_pooled_noisy_XtWX"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
		[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Private working directory] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* private_working_dir = argv[10];

			char* shared_working_dir = argv[11];

			cryptable_collaborative_decrypt_pooled_noisy_XtWX(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_collaborative_decrypt_pooled_noisy_XtWX option.
		else if (t_string::compare_strings(argv[1], "-cryptable_pool_partially_decrypted_pooled_noisy_XtWx_remove_noise"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
	[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Private working directory] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* private_working_dir = argv[10];

			char* shared_working_dir = argv[11];

			cryptable_pool_partially_decrypted_pooled_noisy_XtWx_remove_noise(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_pool_partially_decrypted_pooled_noisy_XtWx_remove_noise option.

		else if (t_string::compare_strings(argv[1], "-cryptable_client_pool_XtWX_XtWz_update_beta"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
		[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Private working directory] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* private_working_dir = argv[10];

			char* shared_working_dir = argv[11];

			cryptable_client_pool_XtWX_XtWz_update_beta(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_pool_XtWX_XtWz_update_beta option.
		else if (t_string::compare_strings(argv[1], "-cryptable_client_collaborative_decrypt_beta"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
		[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);
			char* private_working_dir = argv[10];
			char* shared_working_dir = argv[11];

			cryptable_client_collaborative_decrypt_beta(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		}
		else if (t_string::compare_strings(argv[1], "-cryptable_client_pool_partially_decrypted_beta"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
		[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);
			char* private_working_dir = argv[10];
			char* shared_working_dir = argv[11];
			cryptable_client_pool_partially_decrypted_beta(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_pool_partially_decrypted_beta option.
		else if (t_string::compare_strings(argv[1], "-cryptable_client_calculate_save_pvalue_stats"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Iter. Number] [Var block size] [GMMAT text genotype file] \
	[Features matrix (row per subject)] \
	[Observed phenotypes file] \
	[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int i_iter = atoi(argv[4]);
			int var_block_size = atoi(argv[5]);
			char* GMMAT_text_genotype_matrix_fp = argv[6];
			char* subject_per_row_feats_fp = argv[7];
			char* subject_per_row_obs_pheno_fp = argv[8];
			char* private_working_dir = argv[9];
			char* shared_working_dir = argv[10];

			cryptable_client_calculate_save_pvalue_stats(client_i, n_clients, i_iter, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_calculate_save_pvalue_stats option.
		else if (t_string::compare_strings(argv[1], "-cryptable_client_pool_pvalue_stats"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Iter. index] [Var block size] [GMMAT text genotype file] \
	[Features matrix (row per subject)] \
	[Observed phenotypes file] \
	[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int i_iter = atoi(argv[4]);
			int var_block_size = atoi(argv[5]);
			char* GMMAT_text_genotype_matrix_fp = argv[6];
			char* subject_per_row_feats_fp = argv[7];
			char* subject_per_row_obs_pheno_fp = argv[8];
			char* private_working_dir = argv[9];
			char* shared_working_dir = argv[10];

			cryptable_client_pool_pvalue_stats(client_i, n_clients, i_iter, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp,
				private_working_dir, 
				shared_working_dir);
		} // -cryptable_client_pool_pvalue_stats option.	

		else if (t_string::compare_strings(argv[1], "-cryptable_client_collaborative_decrypt_pval_stats"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Iter. index] [Var block size] [GMMAT text genotype file] \
	[Features matrix (row per subject)] \
	[Observed phenotypes file] \
	[Private working directory] \
	[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int i_iter = atoi(argv[4]);
			int var_block_size = atoi(argv[5]);
			char* GMMAT_text_genotype_matrix_fp = argv[6];
			char* subject_per_row_feats_fp = argv[7];
			char* subject_per_row_obs_pheno_fp = argv[8];
			char* private_working_dir = argv[9];
			char* shared_working_dir = argv[10];

			cryptable_client_collaborative_decrypt_pval_stats(client_i, n_clients, i_iter, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_collaborative_decrypt_pval_stats option.
		else if (t_string::compare_strings(argv[1], "-cryptable_client_pool_partially_decrypted_pval_stats"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Iter. index] [Var block size] [GMMAT text genotype file] \
	[Features matrix (row per subject)] \
	[Observed phenotypes file] \
	[Private working directory] \
	[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int i_iter = atoi(argv[4]);
			int var_block_size = atoi(argv[5]);
			char* GMMAT_text_genotype_matrix_fp = argv[6];
			char* subject_per_row_feats_fp = argv[7];
			char* subject_per_row_obs_pheno_fp = argv[8];
			char* private_working_dir = argv[9];
			char* shared_working_dir = argv[10];

			cryptable_client_pool_partially_decrypted_pval_stats(client_i, n_clients, i_iter, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp,
				private_working_dir,
				shared_working_dir);
		} // -cryptable_client_collaborative_decrypt_pval_stats option.

		// Following are the plaintext training options.
		// Following are the plaintext training options.
		// Following are the plaintext training options.
		// Following are the plaintext training options.
		else if (t_string::compare_strings(argv[1], "-client_calculate_save_XtWX_XtWz"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
	[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* shared_working_dir = argv[10];

			client_calculate_save_XtWX_XtWz(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				shared_working_dir);
		} // -client_calculate_save_XtWX_XtWz option.
		else if (t_string::compare_strings(argv[1], "-client_pool_XtWX_XtWz_update_beta"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] \
	[Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [LL epsilon for convergence] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);
			char* subject_per_row_feats_fp = argv[5];
			char* subject_per_row_pheno_fp = argv[6];
			int n_epoch = atoi(argv[7]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[8], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[8], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[8], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[9]);

			char* shared_working_dir = argv[10];

			client_pool_XtWX_XtWz_update_beta(i_iter, client_i, n_clients,
				subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				shared_working_dir);
		} // -client_pool_XtWX_XtWz_update_beta
		else if (t_string::compare_strings(argv[1], "-client_check_convergence_per_updated_beta"))
		{
			if (argc != 7)
			{
				fprintf(stderr, "%s %s [iteration index (0,1,...)] [client id (0,1...)] [# clients] [LL epsilon for convergence] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int i_iter = atoi(argv[2]);
			int client_i = atoi(argv[3]);
			int n_clients = atoi(argv[4]);

			double LL_EPSILON = atof(argv[5]);

			char* shared_working_dir = argv[6];

			// Check for convergence using all LL stats from all clients.
			client_check_convergence_per_updated_beta(i_iter, client_i, n_clients, shared_working_dir, LL_EPSILON);
		} // -client_check_convergence_per_updated_beta option.
		else if (t_string::compare_strings(argv[1], "-client_calculate_save_pvalue_stats_for_meta_analysis"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Var block size] [GMMAT text genotype file] \
	[Features matrix (row per subject)] \
	[Observed phenotypes file] [Null model phenotypes file] \
	[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int var_block_size = atoi(argv[4]);
			char* GMMAT_text_genotype_matrix_fp = argv[5];
			char* subject_per_row_feats_fp = argv[6];
			char* subject_per_row_obs_pheno_fp = argv[7];
			char* subject_per_row_null_pheno_fp = argv[8];

			char* shared_working_dir = argv[9];

			client_calculate_save_pvalue_stats_for_meta_analysis(client_i, n_clients, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp, subject_per_row_null_pheno_fp,
				shared_working_dir);
		} // -client_calculate_save_pvalue_stats option.
		else if (t_string::compare_strings(argv[1], "-client_calculate_save_pvalue_stats"))
		{
			if (argc != 10)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Var block size] [GMMAT text genotype file] \
		[Features matrix (row per subject)] \
		[Observed phenotypes file] [Null model phenotypes file] \
		[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int var_block_size = atoi(argv[4]);
			char* GMMAT_text_genotype_matrix_fp = argv[5];
			char* subject_per_row_feats_fp = argv[6];
			char* subject_per_row_obs_pheno_fp = argv[7];
			char* subject_per_row_null_pheno_fp = argv[8];

			char* shared_working_dir = argv[9];

			client_calculate_save_pvalue_stats(client_i, n_clients, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp, subject_per_row_null_pheno_fp,
				shared_working_dir);
		} // client_calculate_save_pvalue_stats option.
		else if (t_string::compare_strings(argv[1], "-client_update_pvalue_scale_stats"))
		{
			if (argc != 11)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Iter. index] [Var block size] [GMMAT text genotype file] \
			[Features matrix (row per subject)] \
			[Observed phenotypes file] [Null model phenotypes file] \
			[Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int i_iter = atoi(argv[4]);
			int var_block_size = atoi(argv[5]);
			char* GMMAT_text_genotype_matrix_fp = argv[6];
			char* subject_per_row_feats_fp = argv[7];
			char* subject_per_row_obs_pheno_fp = argv[8];
			char* subject_per_row_null_pheno_fp = argv[9];
			char* shared_working_dir = argv[10];

			client_update_pvalue_scale_stats_per_block(client_i, n_clients, i_iter, var_block_size, GMMAT_text_genotype_matrix_fp,
				subject_per_row_feats_fp,
				subject_per_row_obs_pheno_fp, subject_per_row_null_pheno_fp,
				shared_working_dir);
		} // -client_update_pvalue_scale_stats 
		else if (t_string::compare_strings(argv[1], "-client_pool_pvalue_stats"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Var block size] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int var_block_size = atoi(argv[4]);
			char* shared_working_dir = argv[5];

			client_pool_pvalue_stats(client_i, n_clients, var_block_size, shared_working_dir);
		} // -client_pool_pvalue_stats option.
		else if (t_string::compare_strings(argv[1], "-client_pool_pvalue_stats_for_meta_analysis"))
		{
			if (argc != 6)
			{
				fprintf(stderr, "%s %s [client id (0,1...)] [# clients] [Var block size] [Shared working directory]\n", argv[0], argv[1]);
				exit(0);
			}

			int client_i = atoi(argv[2]);
			int n_clients = atoi(argv[3]);
			int var_block_size = atoi(argv[4]);
			char* shared_working_dir = argv[5];

			client_pool_pvalue_stats_for_meta_analysis(client_i, n_clients, var_block_size, shared_working_dir);
		} // -client_pool_pvalue_stats option.
		else if (t_string::compare_strings(argv[1], "-test_log_sigmoid_approximations"))
		{
			if (argc != 3)
			{
				fprintf(stderr, "%s %s [max range]\n", argv[0], argv[1]);
				exit(0);
			}

			double max_range = atof(argv[2]);

			for (double val = (-1 * max_range); val <= max_range; val += 0.01)
			{
				double kim_etal_sigmoid = get_Kim_etal_poly_approx_sigmoid_per_feat_comb(val);
				double native_sigmoid = get_sigmoid_val_per_feat_comb(val);
				double tenseal_sigmoid = get_TenSeal_poly_approx_sigmoid_per_feat_comb(val);
				fprintf(stderr, "%.3f\t%.4f\t%.4f\t%.4f\n",
					val,
					kim_etal_sigmoid,
					native_sigmoid,
					tenseal_sigmoid);
			} // val loop.

			exit(0);
		} // -test_function_approximations option.
		else if (t_string::compare_strings(argv[1], "-test_LR_modified_IRLS_training_plaintext"))
		{
			if (argc != 8)
			{
				fprintf(stderr, "%s %s [Features file] [Phenotype file] [# epoch] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [Convergence epsilon] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}

			char* subject_per_row_feats_fp = argv[2];
			char* subject_per_row_pheno_fp = argv[3];
			int n_epoch = atoi(argv[4]);

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[5], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[5], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[5], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[6]);
			char* op_dir = argv[7];

			train_modified_LR_per_IRLS(subject_per_row_feats_fp, subject_per_row_pheno_fp,
				n_epoch,
				sigmoid_approx_type,
				LL_EPSILON,
				op_dir);
		} // -test_LR_modified_IRLS_training_plaintext option.
		else if (t_string::compare_strings(argv[1], "-test_LR_modified_training_plaintext"))
		{
			if (argc != 15)
			{
				fprintf(stderr, "%s %s [Features file] [Phenotype file] [Model file] [Sample size] [# covariates to generate] [Update step size] [Regularization weight] [# epoch] \
	[Phenotype noise std. dev.] [Gradient addition order option (PER_FEAT/PER_SUBJECT/PER_SAMPLE)] [Sigmoid approximation type (NATIVE/KIM_ETAL/TENSEAL)] [Convergence epsilon] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}

			char* subject_per_row_feats_fp = argv[2];
			char* subject_per_row_pheno_fp = argv[3];
			char* model_fp = argv[4];
			int sample_size = atoi(argv[5]);
			int n_covars = atoi(argv[6]);
			double update_step_size = atof(argv[7]);
			double regularization_weight = atof(argv[8]);
			int n_epoch = atoi(argv[9]);
			double residual_pheno_noise_std_dev = atof(argv[10]);
			int gradient_addition_order_option = 0;
			if (t_string::compare_strings(argv[11], "PER_FEAT"))
			{
				gradient_addition_order_option = GRADIENT_ADD_PER_FEAT;
			}
			else if (t_string::compare_strings(argv[11], "PER_SUBJECT"))
			{
				gradient_addition_order_option = GRADIENT_ADD_PER_SUBJECT;
			}
			else if (t_string::compare_strings(argv[11], "PER_SAMPLE"))
			{
				gradient_addition_order_option = GRADIENT_ADD_PER_SAMPLE;
			}
			else
			{
				fprintf(stderr, "Unknown gradient update order option; use PER_FEAT/PER_SUBJECT/PER_SAMPLE\n");
				exit(0);
			}

			int sigmoid_approx_type = 0;
			if (t_string::compare_strings(argv[12], "NATIVE"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_NATIVE;
			}
			else if (t_string::compare_strings(argv[12], "KIM_ETAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_KIM_ETAL;
			}
			else if (t_string::compare_strings(argv[12], "TENSEAL"))
			{
				sigmoid_approx_type = SIGMOID_APPROX_TENSEAL;
			}
			else
			{
				fprintf(stderr, "Unknown sigmoid approximation option; use NATIVE/KIM_ETAL/TENSEAL\n");
				exit(0);
			}

			double LL_EPSILON = atof(argv[13]);
			char* op_dir = argv[14];

			train_modified_LR(subject_per_row_feats_fp, subject_per_row_pheno_fp, model_fp,
				sample_size, n_covars,
				residual_pheno_noise_std_dev,
				update_step_size,
				regularization_weight,
				n_epoch,
				gradient_addition_order_option,
				sigmoid_approx_type,
				LL_EPSILON,
				op_dir);
		} // test_LR_modified_training_plaintext option.
		else if (t_string::compare_strings(argv[1], "-assign_ChiSqr_pvalues_per_GMMAT_text_geno"))
		{
			if (argc != 7)
			{
				fprintf(stderr, "%s %s [GMMAT text genotype file] [Features matrix (row per subject)] [Observed phenotypes file] [Null model phenotypes file] [Output file]\n", argv[0], argv[1]);
				exit(0);
			}

			char* GMMAT_text_genotype_matrix_fp = argv[2];
			char* subject_per_row_feats_fp = argv[3];
			char* obs_pheno_fp = argv[4];
			char* null_model_pheno_fp = argv[5];
			char* op_fp = argv[6];

			assign_ChiSqr_p_values_per_fit_model_pheno_per_GMMAT_text_geno(GMMAT_text_genotype_matrix_fp, subject_per_row_feats_fp, obs_pheno_fp, null_model_pheno_fp, op_fp);
		} // assign_ChiSqr_pvalues_per_GMMAT_text_geno option.
		else if (t_string::compare_strings(argv[1], "-test_LR_baseline_training_plaintext"))
		{
			if (argc != 12)
			{
				fprintf(stderr, "%s %s [Features file] [Phenotype file] [Model file] [Sample size] [# covariates to generate] [Update step size] [Regularization weight] [# epoch] [Phenotype noise std. dev.] [Output directory]\n", argv[0], argv[1]);
				exit(0);
			}

			char* subject_per_row_feats_fp = argv[2];
			char* subject_per_row_pheno_fp = argv[3];
			char* model_fp = argv[4];
			int sample_size = atoi(argv[5]);
			int n_covars = atoi(argv[6]);
			double update_step_size = atof(argv[7]);
			double regularization_weight = atof(argv[8]);
			int n_epoch = atoi(argv[9]);
			double residual_pheno_noise_std_dev = atof(argv[10]);
			char* op_dir = argv[11];

			train_LR_baseline(subject_per_row_feats_fp, subject_per_row_pheno_fp, model_fp,
				sample_size, n_covars,
				residual_pheno_noise_std_dev,
				update_step_size,
				regularization_weight,
				n_epoch, op_dir);
		} // -test_LR_training_plaintext option.
		else
		{
			fprintf(stderr, "Unknown option: %s\n", argv[1]);
			exit(0);
		}
	} // msin try block.
	catch (std::exception& e)
	{
		fprintf(stderr, "Fatal error: %s\n", e.what());
		exit(1);
	}
}
