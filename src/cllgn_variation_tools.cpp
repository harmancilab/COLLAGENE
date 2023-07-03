#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cllgn_variation_tools.h"
#include "cllgn_genomics_coords.h"
#include "cllgn_annot_region_tools.h"
#include "cllgn_genome_sequence_tools.h"
#include "cllgn_gff_utils.h"
#include "cllgn_file_utils.h"
#include "cllgn_genomics_coords.h"
#include "cllgn_ansi_string.h"
#include "cllgn_nucleotide.h"
#include <string.h>
#include <ctype.h>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cllgn_nomenclature.h"
#include "cllgn_genome_sequence_tools.h"
#include "cllgn_nucleotide.h"
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

bool __DUMP_VARIATION_TOOLS_MSGS__ = false;

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

void replace_variant_alleles_per_genotypes_signals(char* geno_regs_BED_fp, char* sample_ids_list_fp, char* new_alleles_BED_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(geno_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d genotype regions with %d samples.\n", (int)geno_regs->size(), (int)sample_ids->size());

	// Set the scores.
	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		geno_regs->at(i_reg)->score = 0;
	}

	vector<t_annot_region*>* new_allele_regs = load_BED(new_alleles_BED_fp);
	fprintf(stderr, "Loaded %d new allele regions.\n", (int)new_allele_regs->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(geno_regs, new_allele_regs, false);
	fprintf(stderr, "Processing %d intersections.\n", (int)intersects->size());

	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* geno_reg = int_info->src_reg;
		t_annot_region* new_allele_reg = int_info->dest_reg;

		vector<char*>* geno_reg_name_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(geno_reg->name, "_"));
		vector<char*>* new_allele_reg_name_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(new_allele_reg->name, "_"));

		int n_geno_toks = (int)geno_reg_name_toks->size();
		int n_new_allele_toks = (int)new_allele_reg_name_toks->size();

		if (n_geno_toks != n_new_allele_toks)
		{
			fprintf(stderr, "Sanity check failed, name formatting does not match: %s, %s\n", geno_reg->name, new_allele_reg->name);
			exit(0);
		}

		char allele_codings[10];
		allele_codings[0] = new_allele_reg_name_toks->at(n_new_allele_toks - 2)[0];
		allele_codings[1] = new_allele_reg_name_toks->at(n_new_allele_toks - 1)[0];
		char allele_values[10];
		allele_values[0] = new_allele_reg_name_toks->at(n_new_allele_toks - 2)[2];
		allele_values[1] = new_allele_reg_name_toks->at(n_new_allele_toks - 1)[2];

		if (new_allele_reg->strand == '-')
		{
			allele_values[0] = get_complementary_dna_nuc_per_dna_nuc(allele_values[0]);
			allele_values[1] = get_complementary_dna_nuc_per_dna_nuc(allele_values[1]);
		}

		char new_ref_allele = 0;
		char new_alt_allele = 0;
		if (geno_reg_name_toks->at(n_geno_toks - 2)[0] == allele_codings[0])
		{
			new_ref_allele = allele_values[0];
		}
		else if (geno_reg_name_toks->at(n_geno_toks - 2)[0] == allele_codings[1])
		{
			new_ref_allele = allele_values[1];
		}
		else
		{
			fprintf(stderr, "Could not map reference allele (%c): %s ;; %s\n", geno_reg_name_toks->at(n_geno_toks - 2)[0], geno_reg->name, new_allele_reg->name);
			exit(0);
		}

		// Map alternate allele.
		if (geno_reg_name_toks->at(n_geno_toks - 1)[0] == allele_codings[0])
		{
			new_alt_allele = allele_values[0];
		}
		else if (geno_reg_name_toks->at(n_geno_toks - 1)[0] == allele_codings[1])
		{
			new_alt_allele = allele_values[1];
		}
		else
		{
			fprintf(stderr, "Could not map alternate allele (%c): %s ;; %s\n", geno_reg_name_toks->at(n_geno_toks - 1)[0], geno_reg->name, new_allele_reg->name);
			exit(0);
		}

		if (new_ref_allele == new_alt_allele)
		{
			fprintf(stderr, "Sanity check failed, mapped ref/alt to same allele: %c, %c:: %s ;; %s\n", new_ref_allele, new_alt_allele, 
					geno_reg->name, new_allele_reg->name);

			exit(0);
		}

		// Everything looks good, replace new ref/alt alleles, then move on.
		geno_reg_name_toks->at(n_geno_toks - 2)[0] = new_ref_allele;
		geno_reg_name_toks->at(n_geno_toks - 1)[0] = new_alt_allele;

		// Concatenate to get the final identifier.
		char concat_new_geno_id[1000];
		strcpy(concat_new_geno_id, geno_reg_name_toks->at(0));

		for (int i_tok = 1; i_tok < n_geno_toks; i_tok++)
		{
			strcat(concat_new_geno_id, "_");
			strcat(concat_new_geno_id, geno_reg_name_toks->at(i_tok));
		} // i_tok loop.

		geno_reg->score = 1;
		fprintf(stderr, "New id: %s:%d::%s / %s\n", geno_reg->chrom, geno_reg->start, concat_new_geno_id, geno_reg->name);
		geno_reg->name = t_string::copy_me_str(concat_new_geno_id);
	} // i_int loop.

	vector<t_annot_region*>* reset_geno_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		if (geno_regs->at(i_reg)->score == 1)
		{
			reset_geno_regs->push_back(geno_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Saving allele recoded regions for %d/%d regions.\n", (int)reset_geno_regs->size(), (int)geno_regs->size());

	binarize_variant_genotype_signal_regions(reset_geno_regs, NULL, sample_ids, op_fp);
} // replace_variant_alleles_per_genotypes_signals 

vector<t_annot_region*>* load_variant_signal_regions_wrapper(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp)
{
	vector<t_annot_region*>* var_sig_regs = NULL;

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	if (sample_ids == NULL)
	{
		fprintf(stderr, "Could not load the sample ids from %s", sample_ids_list_fp);
		exit(0);
	}

	if (t_string::compare_strings(geno_sig_regs_BED_fp, "stdin") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed") || 
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt.gz") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed.gz"))
	{
		var_sig_regs = load_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else if (t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed.gz"))
	{
		var_sig_regs = load_binarized_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else
	{
		fprintf(stderr, "%s(%d): Could not identify variant signal regions file type: %s\n", __FILE__, __LINE__, geno_sig_regs_BED_fp);
		exit(0);
	}

	t_string::clean_string_list(sample_ids);

	return(var_sig_regs);
}

void merge_samples_per_regions_list(char* regions_list_fp,
									char* samples_list_fp,
									bool match_region_names,
									char* regions_op_fp,
									char* samples_op_fp)
{
	vector<char*>* regions_list = buffer_file(regions_list_fp);
	fprintf(stderr, "Pooling %d regions.\n", (int)regions_list->size());

	vector<char*>* samples_list = buffer_file(samples_list_fp);
	fprintf(stderr, "Pooling %d sample lists.\n", (int)samples_list->size());

	if ((int)samples_list->size() != (int)regions_list->size())
	{
		fprintf(stderr, "Samples list does not match to regions list.\n");
		exit(0);
	}

	vector<char*>* pooled_sample_ids = new vector<char*>();
	for (int i_f = 0; i_f < (int)regions_list->size(); i_f++)
	{
		vector<char*>* cur_sample_ids = buffer_file(samples_list->at(i_f));
		pooled_sample_ids->insert(pooled_sample_ids->end(), cur_sample_ids->begin(), cur_sample_ids->end());
	} // i_f loop.

	fprintf(stderr, "Loaded %d samples, pooling.\n", (int)pooled_sample_ids->size());

	vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(regions_list->at(0), samples_list->at(0));
	vector<char*>* first_samples_list = buffer_file(samples_list->at(0));
	int max_coded_geno = get_max_genotype_value(geno_regs, first_samples_list);

	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		void** cur_reg_info = (void**)(geno_regs->at(i_reg)->data);
		char* cur_geno = (char*)(cur_reg_info[0]);
		char* pooled_geno = new char[(int)pooled_sample_ids->size() + 2];
		memset(pooled_geno, 0, sizeof(char) * (int)pooled_sample_ids->size());

		delete[] cur_geno;

		// Replace the genotypes.
		cur_reg_info[0] = pooled_geno;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d regions.\n", (int)geno_regs->size());

	vector<char*>* samples_so_far = new vector<char*>();
	for (int i_f = 0; i_f < (int)regions_list->size(); i_f++)
	{
		vector<t_annot_region*>* cur_geno_regs = load_variant_signal_regions_wrapper(regions_list->at(i_f), samples_list->at(i_f));
		vector<char*>* cur_samples_list = buffer_file(samples_list->at(i_f));

		int cur_max_coded_geno = get_max_genotype_value(cur_geno_regs, cur_samples_list);

		if (cur_max_coded_geno != max_coded_geno)
		{
			fprintf(stderr, "Sanity check failed: Genocoding of variants are not the same @ %s (%d/%d).\n", regions_list->at(i_f), max_coded_geno, cur_max_coded_geno);
			exit(0);
		}

		vector<t_annot_region*>* intersects = intersect_annot_regions(cur_geno_regs, geno_regs, false);
		fprintf(stderr, "Pooling %s (%d samples) using %d/%d intersects\n", regions_list->at(i_f), (int)cur_samples_list->size(),
			(int)cur_geno_regs->size(), (int)geno_regs->size());

		for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
		{
			t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* cur_geno_reg = cur_int_info->src_reg;
			t_annot_region* cur_pooled_reg = cur_int_info->dest_reg;

			void** pooled_reg_info = (void**)(cur_pooled_reg->data);
			char* pooled_geno = (char*)(pooled_reg_info[0]);

			void** cur_reg_info = (void**)(cur_geno_reg->data);
			char* cur_geno = (char*)(cur_reg_info[0]);

			// Copy the genotypes.
			int cur_starting_i_s = (int)samples_so_far->size();
			for (int i_s = cur_starting_i_s; 
				i_s < cur_starting_i_s + (int)cur_samples_list->size();
				i_s++)
			{
				pooled_geno[i_s] = cur_geno[i_s - cur_starting_i_s];
			} // i_s loop.

			delete cur_int_info;
		} // i_int loop.
		delete_annot_regions(intersects);

		// Delete memory for the current matrix.
		for (int i_reg = 0; i_reg < (int)cur_geno_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(cur_geno_regs->at(i_reg)->data);
			char* cur_reg_geno = (char*)(cur_reg_info[0]);
			delete[] cur_reg_geno;
			delete[] cur_reg_info;
		} // i_reg loop.
		delete_annot_regions(cur_geno_regs);

		samples_so_far->insert(samples_so_far->end(), cur_samples_list->begin(), cur_samples_list->end());
	} // i_f loop.

	fprintf(stderr, "Saving %d samples (%d samples).\n", (int)samples_so_far->size(), (int)pooled_sample_ids->size());
	binarize_variant_genotype_signal_regions(geno_regs, NULL, pooled_sample_ids, regions_op_fp);

	// Save sample id's.
	FILE* f_samples_list = open_f(samples_op_fp, "w");
	for (int i_s = 0; i_s < (int)pooled_sample_ids->size(); i_s++)
	{
		fprintf(f_samples_list, "%s\n", pooled_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_samples_list);
}

void merge_samples_per_matching_genotype_signal_regions(char* sample1_regs_fp, char* sample1_list_fp, 
														char* sample2_regs_fp, char* sample2_list_fp, bool match_region_names, char* op_fp)
{
	vector<char*>* sample1_ids = buffer_file(sample1_list_fp);
	vector<t_annot_region*>* geno1_sig_regs = load_variant_signal_regions_wrapper(sample1_regs_fp, sample1_list_fp);

	int max_coded_geno1 = get_max_genotype_value(geno1_sig_regs, sample1_ids);

	vector<char*>* sample2_ids = buffer_file(sample2_list_fp);
	vector<t_annot_region*>* geno2_sig_regs = load_variant_signal_regions_wrapper(sample2_regs_fp, sample2_list_fp);

	int max_coded_geno2 = get_max_genotype_value(geno2_sig_regs, sample2_ids);

	if (max_coded_geno1 != max_coded_geno2)
	{
		fprintf(stderr, "Genocoding of variants seem to be different (%d/%d), make sure they are the same.\n", max_coded_geno1, max_coded_geno2);
		exit(0);
	}

	fprintf(stderr, "Merging %d regions of %d samples in %s and %d regions of %d samples and saving to %s; %s\n",
		(int)geno1_sig_regs->size(), (int)sample1_ids->size(), sample1_regs_fp,
		(int)geno2_sig_regs->size(), (int)sample2_ids->size(), sample2_regs_fp,
		op_fp);

	vector<t_annot_region*>* intersects = intersect_annot_regions(geno1_sig_regs, geno2_sig_regs, true);
	fprintf(stderr, "Processing %d matching regions.\n", (int)intersects->size());

	vector<t_annot_region*>* merged_sample_genotype_regions = new vector<t_annot_region*>();
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* geno1_reg = int_info->src_reg;
		t_annot_region* geno2_reg = int_info->dest_reg;

		if (geno1_reg->start == geno2_reg->start &&
			geno1_reg->end == geno2_reg->end)
		{
			if (!match_region_names ||
				t_string::compare_strings(geno1_reg->name, geno2_reg->name))
			{
				t_annot_region* new_reg = get_empty_region();
				new_reg->chrom = t_string::copy_me_str(geno1_reg->chrom);
				new_reg->start = geno1_reg->start;
				new_reg->end = geno1_reg->end;
				new_reg->strand = geno1_reg->strand;

				if (match_region_names)
				{
					new_reg->name = t_string::copy_me_str(geno1_reg->name);
				}

				void** reg1_info = (void**)(geno1_reg->data);
				char* reg1_geno_sig = (char*)(reg1_info[0]);
				void** reg2_info = (void**)(geno2_reg->data);
				char* reg2_geno_sig = (char*)(reg2_info[0]);

				char* pooled_geno_sig = new char[(int)sample1_ids->size() + (int)sample2_ids->size() + 2];
				for (int i_s = 0; i_s < (int)sample1_ids->size(); i_s++)
				{
					pooled_geno_sig[i_s] = reg1_geno_sig[i_s];
				}

				for (int i_s = 0; i_s < (int)sample2_ids->size(); i_s++)
				{
					pooled_geno_sig[i_s + (int)sample1_ids->size()] = reg2_geno_sig[i_s];
				}

				// Copy the region info.
				void** new_reg_info = new void*[2];
				new_reg_info[0] = pooled_geno_sig;
				new_reg->data = new_reg_info;

				merged_sample_genotype_regions->push_back(new_reg);
			} // region name check.
		} // coordinate check.
	} // i_int loop.

	vector<char*>* merged_sample_ids = new vector<char*>();
	merged_sample_ids->insert(merged_sample_ids->end(), sample1_ids->begin(), sample1_ids->end());
	merged_sample_ids->insert(merged_sample_ids->end(), sample2_ids->begin(), sample2_ids->end());

	fprintf(stderr, "Saving the %d regions with %d samples to %s\n", (int)merged_sample_genotype_regions->size(), (int)merged_sample_ids->size(), op_fp);

	// Write the sample ids.
	char merged_sample_ids_fp[1000];
	sprintf(merged_sample_ids_fp, "%s_samples.list", op_fp);
	FILE* f_merged_sample_ids = open_f(merged_sample_ids_fp, "w");
	for (int i_s = 0; i_s < (int)merged_sample_ids->size(); i_s++)
	{
		fprintf(f_merged_sample_ids, "%s\n", merged_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_merged_sample_ids);

	binarize_variant_genotype_signal_regions(merged_sample_genotype_regions, NULL, merged_sample_ids, op_fp);
}

void extract_genotype_signals_per_chr_ids(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* chr_ids_list_fp, char* op_fp)
{
	//vector<t_annot_region*>* roi_regs = load_BED(regions_BED_fp);
	vector<char*>* all_chr_ids_2_use = buffer_file(chr_ids_list_fp);
	vector<char*>* chr_ids_2_use = new vector<char*>();
	for (int i_chr = 0; i_chr < (int)all_chr_ids_2_use->size(); i_chr++)
	{
		char* cur_chr_id = t_string::copy_me_str(all_chr_ids_2_use->at(i_chr));
		normalize_chr_id(cur_chr_id);
		chr_ids_2_use->push_back(cur_chr_id);
		fprintf(stderr, "%s (%s)\n", cur_chr_id, all_chr_ids_2_use->at(i_chr));
	} // i_reg loop.
	fprintf(stderr, "Extracting %d chromosomes.\n", (int)chr_ids_2_use->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_sig_regs->size(), (int)sample_ids->size());

	// Count the number of ROI with signal regions.
	vector<t_annot_region*>* roi_regs_w_signals = new vector<t_annot_region*>();
	//int n_roi_w_signal = 0;
	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		if(t_string::get_i_str(chr_ids_2_use, geno_sig_regs->at(i_reg)->chrom) < (int)chr_ids_2_use->size())
		{
			roi_regs_w_signals->push_back(geno_sig_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Matched genotype signals to %d regions.\n", (int)roi_regs_w_signals->size());

	// Save the regions with signals on them.
	binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
} // extract_genotype_signals_per_chr_ids

void extract_genotype_signals_per_region_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* regions_BED_fp, char* op_fp)
{
	vector<t_annot_region*>* roi_regs = load_BED(regions_BED_fp);
	for (int i_reg = 0; i_reg < (int)roi_regs->size(); i_reg++)
	{
		roi_regs->at(i_reg)->data = NULL;
	} // i_reg loop.
	fprintf(stderr, "Extracting genotype signals for %d ROIs.\n", (int)roi_regs->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_sig_regs->size(), (int)sample_ids->size());

	fprintf(stderr, "Intersecting ROI's with genotype signal regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(roi_regs, geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", (int)intersects->size());
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_roi_reg = int_info->src_reg;
		t_annot_region* cur_geno_sig_reg = int_info->dest_reg;

		if (t_string::compare_strings(cur_roi_reg->chrom, cur_geno_sig_reg->chrom) &&
			cur_roi_reg->start == cur_geno_sig_reg->start &&
			cur_roi_reg->end == cur_geno_sig_reg->end)
		{
			if (cur_roi_reg->data != NULL)
			{
				fprintf(stderr, "***WARNING: Already assigned: %s:%d-%d ***\n", cur_geno_sig_reg->chrom, cur_geno_sig_reg->start, cur_geno_sig_reg->end);
				//exit(0);
			}

			cur_roi_reg->data = cur_geno_sig_reg->data;
		}
	} // i_int loop.

	// Count the number of ROI with signal regions.
	vector<t_annot_region*>* roi_regs_w_signals = new vector<t_annot_region*>();
	//int n_roi_w_signal = 0;
	for (int i_reg = 0; i_reg < (int)roi_regs->size(); i_reg++)
	{
		if (roi_regs->at(i_reg)->data != NULL)
		{
			roi_regs_w_signals->push_back(roi_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Matched genotype signals to %d regions.\n", (int)roi_regs_w_signals->size());

	// Save the regions with signals on them.
	binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
}

// Following function extracts the genotype signals for samples from an input VCF file. It focuses on a chromosome and a subset of regions to decrease memory usage.
void extract_genotype_signals_per_VCF_memoptimized(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	char* var_regions_BED_fp,			// Regions to focus on while extracting.
	char* chr_id_2_process,				// Chromosome to process.
	char* bin_seq_dir,					// Sequences directory.
	bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
	bool match_region_names_flag,		// This is a flag to enforce matching of region names.
	bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
	char* op_fp)						// Output file path.
{
	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(0);
	}

	vector<t_annot_region*>* all_var_regions = load_BED(var_regions_BED_fp);

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", (int)vcf_sample_ids->size());

	vector<char*>* chr_ids = get_chr_ids(all_var_regions);

	vector<t_annot_region*>* var_regions = NULL;
	int chrom_i_2_process = t_string::get_i_str(chr_ids, chr_id_2_process);
	if (chrom_i_2_process == (int)chr_ids->size())
	{
		var_regions = all_var_regions;
	}
	else
	{
		var_regions = get_regions_per_chromosome(all_var_regions, chr_ids->at(chrom_i_2_process));
	}

	chr_ids = get_chr_ids(var_regions);

	fprintf(stderr, "Extracting signals on %d regions.\n", (int)var_regions->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}
	else
	{
	}

	char** per_chr_seq = new char*[(int)chr_ids->size() + 2];

	if (match_ref_alleles_flag)
	{
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			int l_chrom = 0;
			fprintf(stderr, "Loading %s.\n", chr_ids->at(i_chr));
			char cur_chr_seq_fp[1000];
			sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));

			if (check_file(cur_chr_seq_fp))
			{
				per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
			}
			else
			{
				sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));

				if (check_file(cur_chr_seq_fp))
				{
					per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
				}
				else
				{
					fprintf(stderr, "Could not load the sequence for %s\n", chr_ids->at(i_chr));
					exit(0);
				}
			}
		} // i_chr loop.
	}
	else
	{
		fprintf(stderr, "Skipping ref genome loading.\n");
	}
	// Set the genotype signal array on all regions.
	for (int i_v = 0; i_v < (int)var_regions->size(); i_v++)
	{
		//char* pooled_var_alleles = var_regions->at(i_v)->name;

		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = NULL;
		cur_reg_info[1] = NULL;
		var_regions->at(i_v)->data = cur_reg_info;
	} // i_v loop.

	FILE* f_vcf = open_f(vcf_fp, "r");

	// Set sorting info per variant regions on each chromosome.
	t_restr_annot_region_list* restr_var_regions = restructure_annot_regions(var_regions);
	for (int i_chr = 0; i_chr < (int)restr_var_regions->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(restr_var_regions->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	char* buff = new char[100 * 1000];
	int n_processed_regions = 0;
	int n_assigned_regions = 0;
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if (n_processed_regions % 10000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region (%d).           \r", n_processed_regions, n_assigned_regions);
		}

		n_processed_regions++;

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT
		char chrom[100];
		char posn_str[100];
		char id[100];
		char ref[100];
		char alt[100];
		char qual[100];
		char filter[100];
		char info[10000];
		char format[100];

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(format, buff);

		// Check for intersect with variant regions.
		int var_i_chr = t_string::get_i_str(restr_var_regions->chr_ids, chrom);
		if (var_i_chr == (int)restr_var_regions->chr_ids->size())
		{
			delete[] cur_line;
			continue;
		}

		int cur_var_posn = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);
		int cur_var_end = cur_var_posn + t_string::string_length(ref) - 1;
		vector<t_annot_region*>* var_regs_per_cur_chr = restr_var_regions->regions_per_chrom[var_i_chr];
		int i_leftmost_reg = locate_posn_region_per_region_starts(cur_var_posn, var_regs_per_cur_chr, 0, (int)var_regs_per_cur_chr->size() - 1);

		// Go back till the cumulative end for the reg2 is to the left of reg1.
		while (i_leftmost_reg > 0 &&
				var_regs_per_cur_chr->at(i_leftmost_reg)->sort_info->cumulative_sorted_end >= cur_var_posn)
		{
			i_leftmost_reg--;
		} // i_leftmost_reg loop.

		// Check if there is any overlap. Check while we are not strictly passed over the region.
		bool all_checks_pass = false;
		int ol_start = 0;
		int ol_end = 0;
		int matching_i_reg = -1;
		while (i_leftmost_reg < (int)var_regs_per_cur_chr->size() &&
				var_regs_per_cur_chr->at(i_leftmost_reg)->start <= cur_var_end)
		{
			ol_start = MAX(var_regs_per_cur_chr->at(i_leftmost_reg)->start, cur_var_posn);
			ol_end = MIN(var_regs_per_cur_chr->at(i_leftmost_reg)->end, cur_var_posn);

			if (ol_end >= ol_start)
			{
				bool coord_check_pass = true;
				if (cur_var_posn != var_regs_per_cur_chr->at(i_leftmost_reg)->start ||
					cur_var_posn != var_regs_per_cur_chr->at(i_leftmost_reg)->end)
				{
					coord_check_pass = false;
				}

				// Do we need to have name matching?
				bool name_check_pass = true;
				if (match_region_names_flag &&
					var_regs_per_cur_chr->at(i_leftmost_reg)->name != NULL &&
					id != NULL)
				{
					int i_char = 0;

					// If there is a dot for the name, skip.
					if (!t_string::compare_strings(var_regs_per_cur_chr->at(i_leftmost_reg)->name, ".") &&
						!t_string::compare_strings(id, ".") &&
						!t_string::compare_substrings_ci(var_regs_per_cur_chr->at(i_leftmost_reg)->name, id, i_char) &&
						!t_string::compare_substrings_ci(id, var_regs_per_cur_chr->at(i_leftmost_reg)->name, i_char))
					{
						name_check_pass = false;
					}
				}

				if (name_check_pass &&
					coord_check_pass)
				{
					all_checks_pass = true;
					matching_i_reg = i_leftmost_reg;
					break;
				}
			} // overlap check.

			i_leftmost_reg++;
		} // left most position tracing loop.

		if (!all_checks_pass)
		{
			delete[] cur_line;
			continue;
		}

		if (matching_i_reg == -1)
		{
			fprintf(stderr, "Matching region index is invalid.\n");
			exit(0);
		}

		// Set the vcf region.
		t_annot_region* cur_vcf_reg = var_regs_per_cur_chr->at(matching_i_reg);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int GT_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];

		// It is important to provide the actual length of buffer not the string length.
		while (t_string::get_next_token(format, cur_tok, l_format + 2, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, "GT"))
			{
				GT_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.
		delete[] cur_tok;

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "Found GT entry @ %d\n", GT_entry_i);
		}

		if (GT_entry_i == -1)
		{
			fprintf(stderr, "Could not find GT entry in %s, skipping\n", cur_line);
			delete[] cur_line;
			continue;
		}

		if (GT_entry_i != 0)
		{
			fprintf(stderr, "Format string is not as expected: %s\n", format);
		}

		// Make sure we have the chromosome id match and we are using alleles.
		char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size() + 2];

		// Initialize everything to not set: This contains the variant signals for each individual.
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			cur_var_geno_sig[i_s] = -1;
		} // i_s loop.

		bool correctly_parsed_all_genotypes = true;
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);

			//if (buff[3] != 0)
			//{
			//	fprintf(stderr, "Failed to read the %d. genotype correctly for line:%s\n", i_s, cur_line);
			//	exit(0);
			//}

			// Check if all the variants are bi-allelic, only.
			if ((buff[0] != '.' && buff[0] != '.' && buff[0] != '0' && buff[0] != '1') ||
				(buff[1] != '|' && buff[1] != '/') ||
				(buff[2] != '0' && buff[2] != '1' && buff[2] != '.' && buff[2] != '.'))
			{
				correctly_parsed_all_genotypes = false;
				fprintf(stderr, "Failed to read the genotype correctly for entry %d in line: %s:%s: %s. sample: %s (Potentially multiallelic)\n", i_s, chrom, posn_str, id, buff);
				//exit(0);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[1] == '|')
				{
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
			}
			else
			{
				cur_var_geno_sig[i_s] = 0;

				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[0] == '1' && buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 2;
				}
				else if (buff[0] == '1' || buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 1;
				}
			}

			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(0);
		}

		if (correctly_parsed_all_genotypes)
		{
			// Update assigned regions count.
			n_assigned_regions++;

			// Intersect with the variant regions.
			int l_ref_allele = t_string::string_length(ref);

			// Don't change the name of the target region.
			//cur_vcf_reg->name = t_string::copy_me_str(id);
			
			cur_vcf_reg->strand = '+';

			// Set the signal to data.
			void** vcf_reg_info = new void*[2];

			// 1st data: signals.
			vcf_reg_info[0] = cur_var_geno_sig;

			// 2nd data: alleles.
			char** ref_alt_alleles = new char*[2];
			ref_alt_alleles[0] = t_string::copy_me_str(ref);
			ref_alt_alleles[1] = t_string::copy_me_str(alt);
			vcf_reg_info[1] = ref_alt_alleles;

			// Set the region info.
			cur_vcf_reg->data = vcf_reg_info;

			// If matching ref alleles is requested, check the alleles.
			if (l_ref_allele < 100 &&
				match_ref_alleles_flag)
			{
				int i_chr_cur_line = t_string::get_i_str(chr_ids, cur_vcf_reg->chrom);
				for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
				{
					if (toupper(per_chr_seq[i_chr_cur_line][i]) != toupper(ref[i - cur_vcf_reg->start]))
					{
						fprintf(stderr, "Could not match the reference allele to chromosome sequence @ %s:%d: %c, %c\n",
							cur_vcf_reg->chrom, cur_vcf_reg->start,
							toupper(per_chr_seq[i_chr_cur_line][i]), toupper(ref[i - cur_vcf_reg->start]));

						exit(0);
					}
				} // i loop.

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Match @ %s:%d-%d (%s): \n%s\n", cur_vcf_reg->chrom, cur_vcf_reg->start, cur_vcf_reg->end, cur_vcf_reg->name, ref);
					for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
					{
						fprintf(stderr, "%c", toupper(per_chr_seq[i_chr_cur_line][i]));
					} // i loop.
					fprintf(stderr, "\n");
				}
			} // allele check is done for variants smaller than 100 bp.
		} // genotype parsing check.

		delete[] cur_line;
	} // vcf file reading loop.

	fprintf(stderr, "Assigned to %d regions in total out of %d VCF regions.\n", n_assigned_regions, n_processed_regions);

	// Copy the current haplotype alleles: Make sure every region has some signal.
	for (int i_reg = 0; i_reg < (int)var_regions->size(); i_reg++)
	{
		void** cur_var_reg_info = (void**)(var_regions->at(i_reg)->data);

		if (cur_var_reg_info[0] == NULL)
		{
			char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size()];
			for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
			{
				cur_var_geno_sig[i_s] = -1;
			} // i_s loop.

			cur_var_reg_info[0] = cur_var_geno_sig;
		}
	} // i_reg loop.

	// Save the binarized signals.
	binarize_variant_genotype_signal_regions(var_regions, NULL, vcf_sample_ids, op_fp);
}

// Following function extracts the genotype signals for samples from an input VCF file. It focuses on a chromosome and a subset of regions to decrease memory usage.
void extract_genotype_signals_per_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
										char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
										char* var_regions_BED_fp,			// Regions to focus on while extracting.
										char* chr_id_2_process,				// Chromosome to process.
										char* bin_seq_dir,					// Sequences directory.
										bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
										bool match_region_names_flag,		// This is a flag to enforce matching of region names.
										bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
										char* op_fp)						// Output file path.
{
	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(0);
	}

	vector<t_annot_region*>* all_var_regions = load_BED(var_regions_BED_fp);

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", (int)vcf_sample_ids->size());

	vector<char*>* chr_ids = get_chr_ids(all_var_regions);

	vector<t_annot_region*>* var_regions = NULL;
	int chrom_i_2_process = t_string::get_i_str(chr_ids, chr_id_2_process);
	if (chrom_i_2_process == (int)chr_ids->size())
	{
		var_regions = all_var_regions;
	}
	else
	{
		var_regions = get_regions_per_chromosome(all_var_regions, chr_ids->at(chrom_i_2_process));
	}

	chr_ids = get_chr_ids(var_regions);

	fprintf(stderr, "Extracting signals on %d regions.\n", (int)var_regions->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}
	else
	{
	}

	char** per_chr_seq = new char*[(int)chr_ids->size() + 2];

	if (match_ref_alleles_flag)
	{
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			int l_chrom = 0;
			fprintf(stderr, "Loading %s.\n", chr_ids->at(i_chr));
			char cur_chr_seq_fp[1000];
			sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));

			if (check_file(cur_chr_seq_fp))
			{
				per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
			}
			else
			{
				sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));

				if (check_file(cur_chr_seq_fp))
				{
					per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
				}
				else
				{
					fprintf(stderr, "Could not load the sequence for %s\n", chr_ids->at(i_chr));
					exit(0);
				}
			}
		} // i_chr loop.
	}
	else
	{
		fprintf(stderr, "Skipping ref genome loading.\n");
	}
	  // Set the genotype signal array on all regions.
	for (int i_v = 0; i_v < (int)var_regions->size(); i_v++)
	{
		//char* pooled_var_alleles = var_regions->at(i_v)->name;

		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = NULL;
		cur_reg_info[1] = NULL;
		var_regions->at(i_v)->data = cur_reg_info;
	} // i_v loop.

	FILE* f_vcf = open_f(vcf_fp, "r");

	int L_BUFF = 100 * 1000;
	char* buff = new char[L_BUFF];
	char* chrom = new char[L_BUFF];
	char* posn_str = new char[L_BUFF];
	char* id = new char[L_BUFF];
	char* ref = new char[L_BUFF];
	char* alt = new char[L_BUFF];
	char* qual = new char[L_BUFF];
	char* filter = new char[L_BUFF];
	char* info = new char[L_BUFF];
	char* format = new char[L_BUFF];

	vector<t_annot_region*>* vcf_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if ((int)vcf_regs->size() % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region.           \r", (int)vcf_regs->size());
		}

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(format, buff);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int GT_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];

		// It is important to provide the actual length of buffer not the string length.
		while (t_string::get_next_token(format, cur_tok, l_format+2, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, "GT"))
			{
				GT_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.
		delete[] cur_tok;

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "Found GT entry @ %d\n", GT_entry_i);
		}

		if (GT_entry_i == -1)
		{
			fprintf(stderr, "Could not find GT entry in %s, skipping\n", cur_line);
			delete[] cur_line;
			continue;
		}

		if (GT_entry_i != 0)
		{
			fprintf(stderr, "Format string is not as expected: %s\n", format);
		}
		
		// Make sure we have the chromosome id match and we are using alleles.
		char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size() + 2];

		// Initialize everything to not set: This contains the variant signals for each individual.
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			cur_var_geno_sig[i_s] = -1;
		} // i_s loop.

		bool correctly_parsed_all_genotypes = true;
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);

			//if (buff[3] != 0)
			//{
			//	fprintf(stderr, "Failed to read the %d. genotype correctly for line:%s\n", i_s, cur_line);
			//	exit(0);
			//}

			// Check if all the variants are bi-allelic, only.
			if ((buff[0] != '.' && buff[0] != '.' && buff[0] != '0' && buff[0] != '1') ||
				(buff[1] != '|' && buff[1] != '/') ||
				(buff[2] != '0' && buff[2] != '1' && buff[2] != '.' && buff[2] != '.'))
			{
				correctly_parsed_all_genotypes = false;
				fprintf(stderr, "Failed to read the %d. genotype correctly for entry %s:%s: %s. sample: %s (Potentially multiallelic)\n", i_s, chrom, posn_str, id, buff);
				//exit(0);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[1] == '|')
				{
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
			}
			else
			{
				cur_var_geno_sig[i_s] = 0;

				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[0] == '1' && buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 2;
				}
				else if (buff[0] == '1' || buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 1;
				}
			}

			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(0);
		}

		if (correctly_parsed_all_genotypes)
		{
			// Intersect with the variant regions.
			t_annot_region* cur_vcf_reg = get_empty_region();
			cur_vcf_reg->chrom = t_string::copy_me_str(chrom);
			normalize_chr_id(cur_vcf_reg->chrom);
			cur_vcf_reg->start = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);

			int l_ref_allele = t_string::string_length(ref);
			cur_vcf_reg->end = cur_vcf_reg->start + l_ref_allele - 1;
			cur_vcf_reg->name = t_string::copy_me_str(id);
			cur_vcf_reg->strand = '+';

			// Set the signal to data.
			void** vcf_reg_info = new void*[2];

			// 1st data: signals.
			vcf_reg_info[0] = cur_var_geno_sig;

			// 2nd data: alleles.
			char** ref_alt_alleles = new char*[2];
			ref_alt_alleles[0] = t_string::copy_me_str(ref);
			ref_alt_alleles[1] = t_string::copy_me_str(alt);
			vcf_reg_info[1] = ref_alt_alleles;

			// If matching ref alleles is requested, check the alleles.
			if (l_ref_allele < 100 &&
				match_ref_alleles_flag)
			{
				int i_chr_cur_line = t_string::get_i_str(chr_ids, cur_vcf_reg->chrom);
				for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
				{
					if (toupper(per_chr_seq[i_chr_cur_line][i]) != toupper(ref[i - cur_vcf_reg->start]))
					{
						fprintf(stderr, "Could not match the reference allele to chromosome sequence @ %s:%d: %c, %c\n",
							cur_vcf_reg->chrom, cur_vcf_reg->start,
							toupper(per_chr_seq[i_chr_cur_line][i]), toupper(ref[i - cur_vcf_reg->start]));

						exit(0);
					}
				} // i loop.

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Match @ %s:%d-%d (%s): \n%s\n", cur_vcf_reg->chrom, cur_vcf_reg->start, cur_vcf_reg->end, cur_vcf_reg->name, ref);
					for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
					{
						fprintf(stderr, "%c", toupper(per_chr_seq[i_chr_cur_line][i]));
					} // i loop.
					fprintf(stderr, "\n");
				}
			} // allele check is done for variants smaller than 100 bp.

			cur_vcf_reg->data = vcf_reg_info;

			vcf_regs->push_back(cur_vcf_reg);
		} // genotype parsing check.

		delete[] cur_line;
	} // vcf file reading loop.

	vector<t_annot_region*>* intersects = intersect_annot_regions(var_regions, vcf_regs, true);
	fprintf(stderr, "\nProcessing %d intersections using coordinate matching.               \n", (int)intersects->size());

	if (match_region_names_flag)
	{
		fprintf(stderr, "Enforcing substring based name matching.\n");
	}

	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* pooled_var_reg = cur_int_info->src_reg;
		t_annot_region* vcf_reg = cur_int_info->dest_reg;

		bool coord_check_pass = true;
		if (pooled_var_reg->start != vcf_reg->start ||
			pooled_var_reg->end != vcf_reg->end)
		{
			coord_check_pass = false;
		}

		// Do we need to have name matching?
		bool name_check_pass = true;
		if (match_region_names_flag &&
			pooled_var_reg->name != NULL &&
			vcf_reg->name != NULL)
		{
			int i_char = 0;

			// If there is a dot for the name, skip.
			if (!t_string::compare_strings(pooled_var_reg->name, ".") &&
				!t_string::compare_strings(vcf_reg->name, ".") &&
				!t_string::compare_substrings_ci(pooled_var_reg->name, vcf_reg->name, i_char) &&
				!t_string::compare_substrings_ci(vcf_reg->name, pooled_var_reg->name, i_char))
			{
				name_check_pass = false;
			}
		}

		if (name_check_pass &&
			coord_check_pass)
		{
			//int i_chr = t_string::get_i_str(chr_ids, pooled_var_reg->chrom);
			//char* chr_seq = per_chr_seq[i_chr];

			// Get the genotype signal for the vcf region.
			void** vcf_reg_info = (void**)(vcf_reg->data);
			char* cur_vcf_reg_geno_sig = (char*)(vcf_reg_info[0]);

			// Set the genotype signal.
			void** cur_var_reg_info = (void**)(pooled_var_reg->data);
			cur_var_reg_info[0] = cur_vcf_reg_geno_sig;
		}
		else
		{

		}

		delete cur_int_info;
	} // i_int loop.
	delete_annot_regions(vcf_regs);

	// Copy the current haplotype alleles.
	for (int i_reg = 0; i_reg < (int)var_regions->size(); i_reg++)
	{
		void** cur_var_reg_info = (void**)(var_regions->at(i_reg)->data);

		if (cur_var_reg_info[0] == NULL)
		{
			char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size()];
			for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
			{
				cur_var_geno_sig[i_s] = -1;
			} // i_s loop.

			cur_var_reg_info[0] = cur_var_geno_sig;
		}
	} // i_reg loop.

	// Save the binarized signals.
	binarize_variant_genotype_signal_regions(var_regions, NULL, vcf_sample_ids, op_fp);

	//// Dump the plain matrix, too.
	//dump_geno_sig_regs_plain(var_regions, vcf_sample_ids, "op.bed");
	
	//// Now dump the variant signal profiles.
	//FILE* f_geno_sig = open_f("op.bed", "w");
	//for (int i_v = 0; i_v < var_regions->size(); i_v++)
	//{
	//	fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
	//		var_regions->at(i_v)->chrom, 
	//		translate_coord(var_regions->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//		translate_coord(var_regions->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
	//		var_regions->at(i_v)->name);

	//	void** cur_var_reg_info = (void**)(var_regions->at(i_v)->data);
	//	int* cur_var_signal = (int*)(cur_var_reg_info[0]);
	//	for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
	//	{
	//		fprintf(f_geno_sig, "\t%d", cur_var_signal[i_s]);
	//	} // i_s loop.

	//	fprintf(f_geno_sig, "\n");
	//} // i_v loop.
	//fclose(f_geno_sig);
}

void convert_genocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* op_VCF_fp)
{
	fprintf(stderr, "Writing the VCF from genocoded %s genotypes.\n", geno_sigs_reg_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sigs_reg_fp, sample_ids_list_fp);

	fprintf(stderr, "Loaded %d regions.\n", (int)geno_sig_regs->size());
	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(geno_sig_regs);

	FILE* f_VCF = open_f(op_VCF_fp, "w");

	// Write the header.
	fprintf(f_VCF, "##fileformat=VCFv4.2\n");
	fprintf(f_VCF, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		fprintf(f_VCF, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_VCF, "\n");

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			if (i_reg % 1000 == 0)
			{
				fprintf(stderr, "Converting %d/%d variant.             \r", i_reg, (int)cur_chr_input_geno_sig_regs->size());
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			if ((int)toks->size() < 3)
			{
				fprintf(stderr, "Could not parse alleles from %s; make sure the region names are formatted as: [Variant ID]_[Ref allele]_[Alt allele]\n", cur_chr_input_geno_sig_regs->at(i_reg)->name);
				exit(0);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_VCF, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				int geno = (int)(per_sample_haplocoded_geno[i_s]);

				if (geno == 0)
				{
					fprintf(f_VCF, "\t0/0");
				}
				else if (geno == 1)
				{
					fprintf(f_VCF, "\t0/1");
				}
				else if (geno == 2)
				{
					fprintf(f_VCF, "\t1/1");
				}
				else
				{
					fprintf(f_VCF, "\t./.");
					//fprintf(stderr, "Could not parse the genotype: %s:%d::%d\n",
					//	cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
					//	(int)(per_sample_haplocoded_geno[i_s]));

					//exit(0);
				}
			} // input genotype existence check.

			fprintf(f_VCF, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_VCF, op_VCF_fp);
}

void convert_haplocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, bool phase_op, char* op_VCF_fp)
{
	fprintf(stderr, "Writing the VCF from %s using phased=%d genotypes.\n", geno_sigs_reg_fp, phase_op);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sigs_reg_fp, sample_ids_list_fp);

	fprintf(stderr, "Loaded %d regions.\n", (int)geno_sig_regs->size());
	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(geno_sig_regs);

	FILE* f_VCF = open_f(op_VCF_fp, "w");

	// Write the header.
	fprintf(f_VCF, "##fileformat=VCFv4.2\n");
	fprintf(f_VCF, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		fprintf(f_VCF, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_VCF, "\n");

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			if (i_reg % 1000 == 0)
			{
				fprintf(stderr, "Converting %d/%d variant.             \r", i_reg, (int)cur_chr_input_geno_sig_regs->size());
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			if ((int)toks->size() < 3)
			{
				fprintf(stderr, "Could not parse alleles from %s; make sure the region names are formatted as: [Variant ID]_[Ref allele]_[Alt allele]\n", cur_chr_input_geno_sig_regs->at(i_reg)->name);
				exit(0);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_VCF, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				if (!phase_op)
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					int geno = cur_hap0 + cur_hap1;

					if (geno == 0)
					{
						fprintf(f_VCF, "\t0/0");
					}
					else if (geno == 1)
					{
						fprintf(f_VCF, "\t0/1");
					}
					else if (geno == 2)
					{
						fprintf(f_VCF, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(0);
					}
				}
				else
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					char geno_str[10];
					strcpy(geno_str, "0|0");

					if (cur_hap0 == 1)
					{
						geno_str[0] = '1';
					}

					if (cur_hap1 == 1)
					{
						geno_str[2] = '1';
					}

					fprintf(f_VCF, "\t%s", geno_str);

					int geno = cur_hap0 + cur_hap1;
					if (geno > 2)
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(0);
					}
				}
			} // input genotype existence check.

			fprintf(f_VCF, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_VCF, op_VCF_fp);
}

int get_genotype_per_haplocoded_genotype(char haplocoded_geno)
{
	int all1 = get_allele_per_haplotype(haplocoded_geno, 0);
	int all2 = get_allele_per_haplotype(haplocoded_geno, 1);
	return(all1 + all2);
}

int get_max_genotype_value(vector<t_annot_region*>* regs, vector<char*>* sample_ids)
{
	int max_geno = 0;
	for (int i_reg = 0; i_reg < (int)regs->size(); i_reg++)
	{
		void** reg_info = (void**)(regs->at(i_reg)->data);
		char* reg_geno = (char*)(reg_info[0]);
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			max_geno = (max_geno < reg_geno[i_s]) ? (reg_geno[i_s]) : (max_geno);

			if (max_geno == 3)
			{
				return(max_geno);
			}
		} // i_s loop.
	} // i_reg loop.

	return(max_geno);
}

bool is_haplocoded_allele_het(char haplocoded_geno)
{
	if (haplocoded_geno == 1 ||
		haplocoded_geno == 2)
	{
		return true;
	}

	return false;
}

int get_allele_per_haplotype(char geno, int hap_i)
{
	int allele = ((geno & (1 << hap_i)) >> hap_i);
	return(allele);
}

void dump_haplo_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, int haplo_2_extract, const char* op_fp)
{
	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < (int)geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);


		void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
		char* cur_var_signal = (char*)(cur_var_reg_info[0]);
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			char cur_sample_haplo_str[10];
			memset(cur_sample_haplo_str, 0, 10);

			if (haplo_2_extract == 2)
			{
				cur_sample_haplo_str[1] = '|';
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
				cur_sample_haplo_str[2] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}
			else if (haplo_2_extract == 3)
			{
				cur_sample_haplo_str[1] = '\t';
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
				cur_sample_haplo_str[2] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}
			else if (haplo_2_extract == 0)
			{
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
			}
			else
			{
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}

			fprintf(f_geno_sig, "\t%s", cur_sample_haplo_str);
		} // i_s loop.

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	fclose(f_geno_sig);
}

void dump_geno_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, bool dump_geno_sig_regs_plain, const char* op_fp)
{
	if (dump_geno_sig_regs_plain)
	{
		fprintf(stderr, "Saving regions only.\n");
	}

	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < (int)geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);

		if (!dump_geno_sig_regs_plain)
		{
			void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
			char* cur_var_signal = (char*)(cur_var_reg_info[0]);
			for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
			{
				fprintf(f_geno_sig, "\t%d", (int)cur_var_signal[i_s]);
			} // i_s loop.
		}

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	fclose(f_geno_sig);
}

void extract_genotype_signals_per_subsample_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* subsample_ids_list_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_sig_regs->size(), (int)sample_ids->size());

	vector<char*>* geno_subsample_ids = buffer_file(subsample_ids_list_fp);
	if (geno_subsample_ids == NULL)
	{
		fprintf(stderr, "Could not load subsample id's list from %s\n", subsample_ids_list_fp);
		exit(0);
	}
	else
	{
		fprintf(stderr, "Extracting genotype signals for %d individuals.\n", (int)geno_subsample_ids->size());
	}

	vector<int>* per_subsample_geno_sample_i = new vector<int>();
	int n_unmatched_samples = 0;
	for (int i_ss = 0; i_ss < (int)geno_subsample_ids->size(); i_ss++)
	{
		int cur_ss_sample_i = t_string::get_i_str(sample_ids, geno_subsample_ids->at(i_ss));
		per_subsample_geno_sample_i->push_back(cur_ss_sample_i);

		if (cur_ss_sample_i == (int)sample_ids->size())
		{
			n_unmatched_samples++;
		}
	} // i_sub_s loop.

	fprintf(stderr, "Found %d unmatched samples.\n", n_unmatched_samples);

	vector<t_annot_region*>* geno_sig_regs_w_subsample_signals = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Extracting %d. region            \r", i_reg);
		}

		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		char* cur_reg_geno_sigs = (char*)(cur_reg_info[0]);

		// Copy the signals.
		char* cur_copy_reg_geno_sigs = new char[(int)geno_subsample_ids->size() + 2];
		for (int i_ss = 0; i_ss < (int)geno_subsample_ids->size(); i_ss++)
		{			
			if (per_subsample_geno_sample_i->at(i_ss) < (int)sample_ids->size())
			{
				cur_copy_reg_geno_sigs[i_ss] = cur_reg_geno_sigs[per_subsample_geno_sample_i->at(i_ss)];
			}
			else
			{
				cur_copy_reg_geno_sigs[i_ss] = -1;
			}
		} // i_ss loop.

		void** cur_copy_reg_info = new void*[2];		
		cur_copy_reg_info[0] = cur_copy_reg_geno_sigs;

		t_annot_region* copy_reg = duplicate_region(geno_sig_regs->at(i_reg));
		copy_reg->data = cur_copy_reg_info;
		geno_sig_regs_w_subsample_signals->push_back(copy_reg);
	} // i_reg loop.

	// Save the regions.
	binarize_variant_genotype_signal_regions(geno_sig_regs_w_subsample_signals, NULL, geno_subsample_ids, op_fp);
}

// Following can dump the supplied regions or load them then dump.
void binarize_variant_genotype_signal_regions(vector<t_annot_region*>* genotype_signal_regions, char* variant_geno_sig_regs_BED_fp, vector<char*>* geno_sample_ids, const char* op_fp)
{
	if (genotype_signal_regions == NULL)
	{
		if (t_string::compare_strings(variant_geno_sig_regs_BED_fp, "stdin") || check_file(variant_geno_sig_regs_BED_fp))
		{
			genotype_signal_regions = load_variant_genotype_signal_regions(variant_geno_sig_regs_BED_fp, geno_sample_ids);
		}		
		else
		{
			fprintf(stderr, "Signal regions are not supplied and could not load them using %s.\n", variant_geno_sig_regs_BED_fp);
			exit(0);
		}
	}

	vector<char*>* chr_ids = get_chr_ids(genotype_signal_regions);

	// This should depend on the type of file from extension.
	FILE* f_op = open_f(op_fp, "wb");
	int n_chrs = (int)chr_ids->size();
	fwrite(&n_chrs, sizeof(int), 1, f_op);
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_chr[1000];
		strcpy(cur_chr, chr_ids->at(i_chr));
		fwrite(cur_chr, sizeof(char), 1000, f_op);
	} // i_chr loop.

	int sample_size = (int)geno_sample_ids->size();
	int n_regs = (int)genotype_signal_regions->size();

	fwrite(&sample_size, sizeof(int), 1, f_op);
	fwrite(&n_regs, sizeof(int), 1, f_op);
	for (int i_reg = 0; i_reg < (int)genotype_signal_regions->size(); i_reg++)
	{
		// Write the chromosome index.
		int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		fwrite(&i_chr, sizeof(int), 1, f_op);
		fwrite(&(reg_BED_start), sizeof(int), 1, f_op);
		fwrite(&(reg_BED_end), sizeof(int), 1, f_op);

		// Write the regions's name.
		if (genotype_signal_regions->at(i_reg)->name == NULL)
		{
			genotype_signal_regions->at(i_reg)->name = t_string::copy_me_str(".");
		}

		int l_reg_name_str = t_string::string_length(genotype_signal_regions->at(i_reg)->name);
		fwrite(&l_reg_name_str, sizeof(int), 1, f_op);
		fwrite(genotype_signal_regions->at(i_reg)->name, sizeof(char), l_reg_name_str, f_op);

		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		fwrite(cur_reg_geno_sig, sizeof(char), (int)geno_sample_ids->size(), f_op);
	} // i_reg loop.

	fprintf(stderr, "\nClosing the file.\n");
	close_f(f_op, op_fp);
	fprintf(stderr, "Finished writing. Doing a sanity check.\n");

	bool perform_check = true;

	if (!perform_check)
	{
		fprintf(stderr, "Skipping sanity check.\n");
		return;
	}

	clock_t cur_time = clock();
	vector<t_annot_region*>* loaded_sig_regs = load_binarized_variant_genotype_signal_regions(op_fp, geno_sample_ids);
	fprintf(stderr, "Loaded in %.4f seconds.\n", ((double)(clock() - cur_time)) / CLOCKS_PER_SEC);

	if (loaded_sig_regs->size() != genotype_signal_regions->size())
	{
		fprintf(stderr, "Dumping/Loading failed.\n");
		exit(0);
	}

	for (int i_reg = 0; i_reg < (int)loaded_sig_regs->size(); i_reg++)
	{
		if (loaded_sig_regs->at(i_reg)->start != genotype_signal_regions->at(i_reg)->start ||
			loaded_sig_regs->at(i_reg)->end != genotype_signal_regions->at(i_reg)->end)
		{
			fprintf(stderr, "Non-match: %s:%d-%d, %s:%d-%d\n", 
					loaded_sig_regs->at(i_reg)->chrom, loaded_sig_regs->at(i_reg)->start, loaded_sig_regs->at(i_reg)->end, 
					genotype_signal_regions->at(i_reg)->chrom, genotype_signal_regions->at(i_reg)->start, genotype_signal_regions->at(i_reg)->end);
			exit(0);
		}

		// Go over the genotype signals.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

		void** cur_loaded_reg_info = (void**)(loaded_sig_regs->at(i_reg)->data);
		char* cur_loaded_reg_geno_sig = (char*)(cur_loaded_reg_info[0]);

		// Write the signal values.
		for (int i_s = 0; i_s < (int)geno_sample_ids->size(); i_s++)
		{
			if (cur_reg_geno_sig[i_s] != cur_loaded_reg_geno_sig[i_s])
			{
				fprintf(stderr, "Non-match: %d. sample genotype: %d, %d\n", i_s,
						cur_reg_geno_sig[i_s], cur_loaded_reg_geno_sig[i_s]);
				exit(0);
			}
		} // i_s loop.
	} // i_reg loop.

	fprintf(stderr, "Check success!\n");
}

vector<t_annot_region*>* load_binarized_variant_genotype_signal_regions(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids)
{
	FILE* f_bin_geno_sig_regs = open_f(bin_geno_sig_bed_fp, "rb");

	// Load the chromosomes.
	int n_chrs = 0;
	fread(&n_chrs, sizeof(int), 1, f_bin_geno_sig_regs);
	vector<char*>* chr_ids = new vector<char*>();
	for (int i_chr = 0; i_chr < n_chrs; i_chr++)
	{
		char cur_chr[1000];
		fread(cur_chr, sizeof(char), 1000, f_bin_geno_sig_regs);
		chr_ids->push_back(t_string::copy_me_str(cur_chr));
	} // i_chr loop.

	fprintf(stderr, "Loaded %d chromosomes.\n", (int)chr_ids->size());
	int sample_size = 0;
	fread(&sample_size, sizeof(int), 1, f_bin_geno_sig_regs);
	fprintf(stderr, "Reading sample size of %d.\n", sample_size);

	if (sample_size != (int)geno_sample_ids->size())
	{
		fprintf(stderr, "Sanity check failed: Sample sizes do not match: %d, %d\n", sample_size, (int)geno_sample_ids->size());
		exit(0);
	}

	int n_regs = 0;
	fread(&n_regs, sizeof(int), 1, f_bin_geno_sig_regs);
	fprintf(stderr, "Reading %d regions.\n", n_regs);

	vector<t_annot_region*>* geno_sig_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < n_regs; i_reg++)
	{
		//int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		//int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		//int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		int i_chr = 0;
		int reg_BED_start = 0;
		int reg_BED_end = 0;
		fread(&i_chr, sizeof(int), 1, f_bin_geno_sig_regs);
		fread(&(reg_BED_start), sizeof(int), 1, f_bin_geno_sig_regs);
		fread(&(reg_BED_end), sizeof(int), 1, f_bin_geno_sig_regs);

		// Read the region's name.
		int l_reg_name_str = 0;
		fread(&l_reg_name_str, sizeof(int), 1, f_bin_geno_sig_regs);
		char* cur_reg_name = new char[l_reg_name_str + 2];
		memset(cur_reg_name, 0, sizeof(char) * (l_reg_name_str + 2));
		fread(cur_reg_name, sizeof(char), l_reg_name_str, f_bin_geno_sig_regs);

		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
		reg->start = translate_coord(reg_BED_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		reg->end = translate_coord(reg_BED_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		reg->strand = '+';
		reg->name = cur_reg_name;

		void** cur_reg_info = new void*[2];
		char* cur_reg_geno_sig = new char[sample_size + 2];
		fread(cur_reg_geno_sig, sizeof(char), sample_size, f_bin_geno_sig_regs);

		cur_reg_info[0] = cur_reg_geno_sig;
		cur_reg_info[1] = NULL;

		reg->data = cur_reg_info;

		geno_sig_regs->push_back(reg);
	} // i_reg loop.

	// Close the file.
	close_f(f_bin_geno_sig_regs, bin_geno_sig_bed_fp);

	return(geno_sig_regs);
}

vector<t_annot_region*>* load_variant_genotype_signal_regions(char* link_variant_genotype_signal_fp, vector<char*>* geno_sample_ids)
{
	vector<t_annot_region*>* link_var_regs = new vector<t_annot_region*>();
	char tok_buff[1000];
	FILE* f_link_variant_genotype_signal = open_f(link_variant_genotype_signal_fp, "r");

	if (f_link_variant_genotype_signal == NULL)
	{
		fprintf(stderr, "Could not open %s\n", link_variant_genotype_signal_fp);
		exit(0);
	}

	int i_reg = 0;
	while (1)
	{
		char* cur_reg_line = getline(f_link_variant_genotype_signal);
		if (cur_reg_line == NULL)
		{
			break;
		}
		else
		{
			i_reg++;
		}

		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Parsing genotype signal for %d. region.                  \r", i_reg);
		}

		//char* cur_reg_line = (char*)(link_var_regs->at(i_reg)->data);
		//t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
		char* cur_var_geno_signal = new char[(int)geno_sample_ids->size()];

		t_annot_region* cur_reg = get_empty_region();

		// Copy the name.
		int i_cur_char = 0;
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->chrom = t_string::copy_me_str(tok_buff);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->start = translate_coord(atoi(tok_buff), BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->end = translate_coord(atoi(tok_buff), BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_reg->strand = '+';

		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->name = t_string::copy_me_str(tok_buff);

		//fprintf(stderr, "%s\t%d\t%d: %s            \r", link_var_regs->at(i_reg)->chrom, link_var_regs->at(i_reg)->start, link_var_regs->at(i_reg)->end, link_var_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < (int)geno_sample_ids->size(); i_s++)
		{
			if (t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char) == false)
			{
				fprintf(stderr, "Could not parse the %d. sample's genotype signal entry.\n", i_s);
				exit(0);
			}

			cur_var_geno_signal[i_s] = (char)(atoi(tok_buff));
		} // i_tok loop.

		void** link_var_reg_info = new void*[3];
		link_var_reg_info[0] = cur_var_geno_signal;
		link_var_reg_info[1] = NULL;

		cur_reg->data = link_var_reg_info;

		link_var_regs->push_back(cur_reg);
	} // i_reg loop.

	close_f(f_link_variant_genotype_signal, link_variant_genotype_signal_fp);

	return(link_var_regs);
}

char* copy_nucs(char* seq)
{
	char* copy = new char[strlen(seq) + 2];
	memset(copy, 0, strlen(seq) + 2 * sizeof(char));

	for (int i = 0; i < (int)strlen(seq); i++)
	{
		copy[i] = seq[i];
	} // i loop.

	return(copy);
}

int* copy_bps(char* seq, int* bps)
{
	int* copy = new int[strlen(seq) + 2];
	memset(copy, 0, strlen(seq) + 2 * sizeof(int));

	for (int i = 0; i < (int)strlen(seq); i++)
	{
		copy[i] = bps[i];
	} // i loop.

	return(copy);
}

void parse_variants_per_info_string(t_vcf_entry* snp_vcf_entry,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	if (snp_vcf_entry->alt_alleles->size() > 1)
	{
		total_allele_count = 0;
		return;
	}

	t_string* info_string = new t_string(snp_vcf_entry->info_str);
	t_string_tokens* info_str_tokens = info_string->tokenize_by_chars("=;");

	int i_tok = 0;
	while (i_tok < (int)info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			ancestral_allele = (info_str_tokens->at(i_tok + 1)->str())[0];
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	delete info_string;
	t_string::clean_tokens(info_str_tokens);
}

void parse_variants_per_info_string(t_vcf_entry* snp_vcf_entry,
	char* ancestral_allele,
	char* derived_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	if (snp_vcf_entry->alt_alleles->size() > 1)
	{
		total_allele_count = 0;
		return;
	}

	t_string* info_string = new t_string(snp_vcf_entry->info_str);
	t_string_tokens* info_str_tokens = info_string->tokenize_by_chars("=;");

	int i_tok = 0;
	while (i_tok < (int)info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			//ancestral_allele = (info_str_tokens->at(i_tok+1)->str())[0];
			strcpy(ancestral_allele, info_str_tokens->at(i_tok + 1)->str());
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	if (t_string::compare_strings_ci(ancestral_allele, snp_vcf_entry->ref_allele->str()))
	{
		// Reference allele is the ancestral allele.
		strcpy(derived_allele, snp_vcf_entry->alt_alleles->at(0)->str());
	}
	else
	{
		// Reference allele is the derived allele.
		strcpy(derived_allele, snp_vcf_entry->ref_allele->str());
	}
}

bool verify_snp(t_vcf_entry* snp, char* ncrna_nucs, char ncrna_strand, int snp_i_seq)
{
	// If there are multiple alternate alleles, do not process this snp.
	if (snp->alt_alleles->size() > 1)
	{
		fprintf(stderr, "Invalid SNP @ %s:%d: Multiple alternate alleles.\n", snp->chrom, snp->ref_pos);
		return(false);
	}

	// Check if the reference allele matches the ncrna's nucleotide at the position.
	char ncrna_dna_nuc = get_transcribing_dna_nuc_per_rna_nuc(ncrna_nucs[snp_i_seq]);
	char pos_strand_nuc = toupper(ncrna_dna_nuc);
	if (ncrna_strand == '-')
	{
		pos_strand_nuc = get_complementary_dna_nuc_per_dna_nuc(pos_strand_nuc);
	}

	// Check whether the reference allele matches the sequence data. Note that the 
	char ref_allele_dna_nuc = (snp->ref_allele->str())[0];
	char alt_allele_dna_nuc = (snp->alt_alleles->at(0)->str())[0];

	if (!is_valid_nuc_char(ref_allele_dna_nuc) ||
		!is_valid_nuc_char(alt_allele_dna_nuc))
	{
		return(false);
	}

	if (ref_allele_dna_nuc == pos_strand_nuc)
	{
		// Check whether the ancestral allele matches one of variants or the reference allele.
		char reference_allele = (snp->ref_allele->str())[0];
		char alternate_allele = (snp->alt_alleles->at(0)->str())[0];

		char ancestral_allele = 0;
		int alt_all_cnt = 0;
		int total_all_cnt = 0;
		parse_variants_per_info_string(snp, ancestral_allele, alt_all_cnt, total_all_cnt);

		if (toupper(ancestral_allele) != toupper(reference_allele) &&
			toupper(ancestral_allele) != toupper(alternate_allele))
		{
			fprintf(stderr, "Invalid SNP @ %s:%d(%d): Ancestral allele does not match neither ref nor alternate allele. %c: %c, %c (%s).\n", snp->chrom, snp->ref_pos, snp_i_seq,
				ancestral_allele, reference_allele, alternate_allele, snp->info_str);
			return(false);
		}

		// The snp seems to be ok to process.
		return(true);
	}
	else
	{
		// The SNP's reference allele does not match the reference nucleotide in the ncrna. Something is messed up in snp entry.
		fprintf(stderr, "Invalid SNP @ %s:%d: SNP's reference allele does not match (%c, %c) the hg19 reference sequence:\n%s\n%d\n", snp->chrom, snp->ref_pos, ref_allele_dna_nuc, pos_strand_nuc, ncrna_nucs, snp_i_seq);
		return(false);
	}
}

void parse_variants_per_info_string(char* info_str,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	t_string_tokens* info_str_tokens = t_string::tokenize_by_chars(info_str, "=;");

	int i_tok = 0;
	while (i_tok < (int)info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			ancestral_allele = (info_str_tokens->at(i_tok + 1)->str())[0];
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	t_string::clean_tokens(info_str_tokens);
}

void get_mutation_matrices_per_DAF(char* snp_vcf_fp, double min_daf, double max_daf)
{
	int** mutation_cnts = new int*[4];
	for (int i = 0; i < 4; i++)
	{
		mutation_cnts[i] = new int[4];
		for (int j = 0; j < 4; j++)
		{
			mutation_cnts[i][j] = 0;
		} // j loop.
	} // i loop.

	  //t_VCF* vcf = new t_VCF(snp_vcf_fp);

	FILE* f_vcf = open_f(snp_vcf_fp, "r");

	//for(int i_snp = 0; i_snp < vcf->entries->size(); i_snp++)
	int i_line = 0;
	while (1)
	{
		if (i_line % 100000 == 0)
		{
			fprintf(stderr, "Processing %d. variant.       \r", i_line);
		}

		char* cur_line = getline(f_vcf);

		if (cur_line == NULL)
		{
			break;
		}

		//t_vcf_entry* cur_snp = vcf->entries->at(i_snp);
		// Get the reference and alternate alleles.
		// 10      68533   .       T       G       .       PASS    AA=.;AC=19;AN=120;DP=80 GT:DP:CB
		char chrom[100];
		char posn[100];
		char dot[100];
		char ref_allele[100];
		char alt_allele[100];
		char dot2[100];
		char pass_str[100];
		char info_str[100];
		if (sscanf(cur_line, "%s %s %s %s %s %s %s %s", chrom, posn, dot, ref_allele, alt_allele, dot2, pass_str, info_str) != 8)
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		char anc_all = 0;
		char der_all = 0;
		int alt_all_cnt = 0;
		int anc_all_cnt = 0;
		int total_all_cnt = 0;

		parse_variants_per_info_string(info_str,
			anc_all,
			alt_all_cnt,
			total_all_cnt);

		// Set the derived allele.
		if (anc_all == ref_allele[0])
		{
			der_all = alt_allele[0];
		}
		else
		{
			der_all = ref_allele[0];
		}

		// Valid nucleotide check on the ancestral allele.
		if (is_valid_nuc_char(anc_all))
		{
			// Make sure that the ancestral allele is either the 
			if (toupper(anc_all) == toupper(ref_allele[0]) ||
				toupper(anc_all) == toupper(alt_allele[0]))
			{
				if (toupper(anc_all) == toupper(alt_allele[0]))
				{
					der_all = toupper(ref_allele[0]);

					// Ancestral allele is the alternate allele. The derived allele is the reference allele.
					anc_all_cnt = alt_all_cnt;
				}
				else
				{
					der_all = toupper(alt_allele[0]);

					anc_all_cnt = total_all_cnt - alt_all_cnt;
				}

				if (is_valid_nuc_char(der_all))
				{
					double der_all_freq = (double)(total_all_cnt - anc_all_cnt) / total_all_cnt;

					if (max_daf > der_all_freq &&
						min_daf < der_all_freq)
					{
						//fprintf(stderr, "%s:%d, %s: %c -> %c (%lf)\n", cur_snp->chrom, cur_snp->ref_pos, cur_snp->info_str, anc_all, der_all, der_all_freq);
						mutation_cnts[nuc_2_num(anc_all)][nuc_2_num(der_all)]++;
					}
				}
			}
		} // validity check on the ancestral allele.

		delete[] cur_line;

		i_line++;
	} // i_snp loop.

	fclose(f_vcf);

	int n_total_trans = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			n_total_trans += mutation_cnts[i][j];
		} // j loop.
	} // i loop.	

	fprintf(stderr, "\n%d nucleotide transitions:\n", n_total_trans);

	fprintf(stderr, "  ");
	for (int j = 0; j < 4; j++)
	{
		fprintf(stderr, "%4c ", num_2_nuc(j));
	}
	fprintf(stderr, "\n");
	for (int i = 0; i < 4; i++)
	{
		fprintf(stderr, "%c ", num_2_nuc(i));
		for (int j = 0; j < 4; j++)
		{
			//fprintf(stderr, "%4d", mutation_cnts[i][j]);
			fprintf(stderr, "%.2lf ", (double)mutation_cnts[i][j] / n_total_trans);
		} // j loop.

		fprintf(stderr, "\n");
	} // i loop.	
}

char mutate_nuc_per_transition_transversion(char nuc)
{
	if (toupper(nuc) == 'A')
	{
		return('G');
	}
	else if (toupper(nuc) == 'C')
	{
		return('T');
	}
	else if (toupper(nuc) == 'G')
	{
		return('A');
	}
	else if (toupper(nuc) == 'T')
	{
		return('C');
	}
	else
	{
		fprintf(stderr, "WTF??\n");
		exit(0);
	}
}

char mutate_nuc(char nuc, t_rng* rng)
{
	char cur_rand_nuc = toupper(nuc);
	while (cur_rand_nuc == toupper(nuc))
	{
		int rand_i = (int)(floor(4 * rng->random_double_ran3()));
		char rand_nucs[] = "ACGT";
		cur_rand_nuc = rand_nucs[rand_i];
	}

	return(cur_rand_nuc);
}

bool mutate_GC_2_AU(char* seq)
{
	// Count GC content.
	int n_gc = 0;
	for (int i = 0; i < (int)strlen(seq); i++)
	{
		if (toupper(seq[i]) == 'G' ||
			toupper(seq[i]) == 'C')
		{
			n_gc++;
		}
	} // i loop.

	if (n_gc == 0)
	{
		return(false);
	}

	int i_rand = rand() % n_gc;

	int i_gc = 0;
	for (int i = 0; i < (int)strlen(seq); i++)
	{
		if (toupper(seq[i]) == 'G' ||
			toupper(seq[i]) == 'C')
		{
			if (i_gc == i_rand)
			{
				// Choose a nuc.
				int n_rand = rand() % 2;
				char rand_nuc = (n_rand == 0) ? ('A') : ('U');
				fprintf(stderr, "%c->%c @ %d\n", seq[i], rand_nuc, i);
				seq[i] = rand_nuc;
				break;
			}

			i_gc++;
		}
	} // i loop.

	return(true);
}

char mutate_nuc_2_GC(char nuc, t_rng* rng)
{
	if (toupper(nuc) != 'G' &&
		toupper(nuc) != 'C')
	{
		int val = (int)(floor(2 * rng->random_double_ran3()));
		if (val == 0)
		{
			return('G');
		}
		else
		{
			return('C');
		}
	}
	else
	{
		return(toupper(nuc));
	}
}

char* get_random_seq(int l)
{
	char* seq = new char[l + 2];
	memset(seq, 0, (l + 1) * sizeof(char));
	for (int i = 0; i < l; i++)
	{
		int random = rand() % 4;
		seq[i] = num_2_nuc(random);
	} // i loop

	return(seq);
}


