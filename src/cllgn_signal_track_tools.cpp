#include "cllgn_signal_track_tools.h"
#include "cllgn_annot_region_tools.h"
#include "cllgn_nomenclature.h"
#include "cllgn_genomics_coords.h"
#include "cllgn_xlog_math.h"
#include "cllgn_file_utils.h"
#include "cllgn_rng.h"
#include "cllgn_ansi_string.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

#include "cllgn_ansi_cli.h"
#include "cllgn_config.h"

#include <ctype.h>
#include <string.h>
#include "cllgn_histogram.h"
#include "cllgn_mapped_read_tools.h"
#include "cllgn_seed_manager.h"
#include <algorithm>

bool __DUMP_SIGNAL_TRACK_MSGS__ = false;

bool sort_signal_nodes_per_increasing_signal(t_signal_node* node1, t_signal_node* node2)
{
	return(node1->signal < node2->signal);
}

double* invert_1based_signal(double* signal, int l_signal)
{
	double* inverted_signal = new double[l_signal + 2];
	memset(inverted_signal, 0, sizeof(double) * (l_signal + 2));

	for (int i = 1; i <= l_signal; i++)
	{
		inverted_signal[i] = signal[l_signal - i + 1];
	} // i loop.

	return(inverted_signal);
}

void quantile_normalize_signal_matrix_4th_col_signals(char* signal_regions_BED_fp, char* op_fp)
{
	int n_loaded_samples = 0;
	fprintf(stderr, "Loading signal regions from %s starting from 4th column.\n", signal_regions_BED_fp);
	vector<t_annot_region*>* signal_node_regs = load_signal_regs_BED(signal_regions_BED_fp, n_loaded_samples, 4);
	fprintf(stderr, "Loaded %d signal regions with %d samples.\n", (int)signal_node_regs->size(), n_loaded_samples);

	// Following stores the signal level along samples.
	vector<t_signal_node*>** per_sample_signal_nodes = new vector<t_signal_node*>*[n_loaded_samples];
	for (int i_s = 0; i_s < n_loaded_samples; i_s++)
	{
		per_sample_signal_nodes[i_s] = new vector<t_signal_node*>();
	} // i_s loop.

	for (int i_reg = 0; i_reg < (int)signal_node_regs->size(); i_reg++)
	{
		double* cur_reg_per_sample_sigs = (double*)(signal_node_regs->at(i_reg)->data);
		vector<t_signal_node*>* cur_reg_per_sample_nodes = new vector<t_signal_node*>();

		for (int i_s = 0; i_s < n_loaded_samples; i_s++)
		{
			t_signal_node* cur_sample_node = new t_signal_node();
			cur_sample_node->i_reg = -1; // This is the rank that will be assigned later.
			cur_sample_node->signal = cur_reg_per_sample_sigs[i_s];

			cur_reg_per_sample_nodes->push_back(cur_sample_node);

			// Add the node to per sample signal nodes, too.
			per_sample_signal_nodes[i_s]->push_back(cur_sample_node);
		} // i_s loop.

		signal_node_regs->at(i_reg)->data = cur_reg_per_sample_nodes;

		delete[] cur_reg_per_sample_sigs;
	} // i_reg loop.

	  // Assign the ranks.
	fprintf(stderr, "Assigning ranks.\n");
	for (int i_s = 0; i_s < n_loaded_samples; i_s++)
	{
		sort(per_sample_signal_nodes[i_s]->begin(), per_sample_signal_nodes[i_s]->end(), sort_signal_nodes_per_increasing_signal);

		for (int rank_i = 0; rank_i < (int)per_sample_signal_nodes[i_s]->size(); rank_i++)
		{
			per_sample_signal_nodes[i_s]->at(rank_i)->i_reg = rank_i;
		} // rank_i loop.
	} // i_s loop.

	fprintf(stderr, "Assigning per rank signals.\n");
	double* per_rank_avg_signals = new double[(int)signal_node_regs->size() + 2];
	for (int rank_i = 0; rank_i < (int)signal_node_regs->size(); rank_i++)
	{
		double cur_rank_total_sig = 0;
		for (int i_s = 0; i_s < n_loaded_samples; i_s++)
		{
			cur_rank_total_sig += per_sample_signal_nodes[i_s]->at(rank_i)->signal;
		} // i_s loop.

		per_rank_avg_signals[rank_i] = cur_rank_total_sig / n_loaded_samples;
	} // rank_i loop.

	fprintf(stderr, "Assigning normalized signals.\n");
	for (int i_s = 0; i_s < n_loaded_samples; i_s++)
	{
		for (int rank_i = 0; rank_i < (int)per_sample_signal_nodes[i_s]->size(); rank_i++)
		{
			per_sample_signal_nodes[i_s]->at(rank_i)->signal = per_rank_avg_signals[rank_i];
		} // rank_i loop.
	} // i_s loop.

	vector<char*>* header_col_ids = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(load_header(signal_regions_BED_fp), "\t"));

	// Assign ranks to each entry.
	fprintf(stderr, "Saving normalized signals to %s.\n", op_fp);
	FILE* f_op = open_f(op_fp, "w");
	for (int col_i = 0; col_i < (int)header_col_ids->size(); col_i++)
	{
		if (col_i > 0)
		{
			fprintf(f_op, "\t%s", header_col_ids->at(col_i));
		}
		else
		{
			fprintf(f_op, "%s", header_col_ids->at(col_i));
		}
	} // col_i loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < (int)signal_node_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", signal_node_regs->at(i_reg)->chrom,
			translate_coord(signal_node_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(signal_node_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			signal_node_regs->at(i_reg)->name);

		vector<t_signal_node*>* cur_reg_per_sample_nodes = (vector<t_signal_node*>*)(signal_node_regs->at(i_reg)->data);
		for (int i_s = 0; i_s < n_loaded_samples; i_s++)
		{
			fprintf(f_op, "\t%lf", cur_reg_per_sample_nodes->at(i_s)->signal);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

void generate_per_region_sample_summarized_signal_stats(vector<char*>* signal_reg_BED_fps, char* op_fp)
{
	vector<t_annot_region*>** per_file_signal_regs = new vector<t_annot_region*>*[(int)signal_reg_BED_fps->size() + 2];
	vector<int>* per_file_n_samples = new vector<int>();
	vector<char*>* per_file_headers = new vector<char*>();
	for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
	{
		int n_loaded_samples = 0;
		per_file_signal_regs[i_f] = load_signal_regs_BED(signal_reg_BED_fps->at(i_f), n_loaded_samples);
		per_file_n_samples->push_back(n_loaded_samples);

		char* cur_file_header = load_header(signal_reg_BED_fps->at(i_f));
		per_file_headers->push_back(cur_file_header);
	} // i_f loop.

	vector<t_annot_region*>* agg_signal_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)per_file_signal_regs[0]->size(); i_reg++)
	{
		t_annot_region* cur_agg_reg = duplicate_region(per_file_signal_regs[0]->at(i_reg));
		t_annot_region** cur_agg_reg_per_file_regs = new t_annot_region*[(int)signal_reg_BED_fps->size() + 2];
		memset(cur_agg_reg_per_file_regs, 0, sizeof(t_annot_region*) * ((int)signal_reg_BED_fps->size() + 2));
		cur_agg_reg->data = cur_agg_reg_per_file_regs;
		agg_signal_regs->push_back(cur_agg_reg);
	} // i_reg loop.

	for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
	{
		fprintf(stderr, "Adding %s\n", signal_reg_BED_fps->at(i_f));
		vector<t_annot_region*>* cur_intersects = intersect_annot_regions(agg_signal_regs, per_file_signal_regs[i_f], true);

		for (int i_int = 0; i_int < (int)cur_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(cur_intersects->at(i_int)->data);
			t_annot_region* cur_agg_reg = int_info->src_reg;
			t_annot_region* cur_sig_reg = int_info->dest_reg;

			if (t_string::compare_strings(cur_agg_reg->name, cur_sig_reg->name))
			{
				t_annot_region** cur_agg_reg_per_file_regs = (t_annot_region**)(cur_agg_reg->data);
				cur_agg_reg_per_file_regs[i_f] = cur_sig_reg;
			}

			delete int_info;
		} // i_int loop.

		delete_annot_regions(cur_intersects);
	} // i_f loop

	fprintf(stderr, "Saving concatenated signal regions.\n");
	FILE* f_op = open_f(op_fp, "w");

	// Write the aggregate header.
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
	{
		fprintf(f_op, "\t%s_MEAN\t%s_VAR", signal_reg_BED_fps->at(i_f), signal_reg_BED_fps->at(i_f));
	} // i_f loop.
	fprintf(f_op, "\n");

	// Go over all the regions.
	for (int i_reg = 0; i_reg < (int)agg_signal_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", agg_signal_regs->at(i_reg)->chrom,
			translate_coord(agg_signal_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(agg_signal_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			agg_signal_regs->at(i_reg)->name);

		t_annot_region** cur_agg_reg_per_file_regs = (t_annot_region**)(agg_signal_regs->at(i_reg)->data);

		// Go over all the regions and dump.
		for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
		{
			if (cur_agg_reg_per_file_regs[i_f] == NULL)
			{
				fprintf(stderr, "Could not find sample %s for region %s:%d-%d\n", signal_reg_BED_fps->at(i_f),
					agg_signal_regs->at(i_reg)->chrom, agg_signal_regs->at(i_reg)->start, agg_signal_regs->at(i_reg)->end);
				exit(0);
			}

			double* cur_reg_signals = (double*)(cur_agg_reg_per_file_regs[i_f]->data);
			double cur_class_mean, cur_class_var;
			get_stats(cur_reg_signals, per_file_n_samples->at(i_f), cur_class_mean, cur_class_var);

			fprintf(f_op, "\t%lf\t%lf", cur_class_mean, cur_class_var);
		} // i_f loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}
void reheader_signal_regions_BED(char* signal_regions_BED_fp, char* sample_ids_fp, char* op_fp)
{
	int n_loaded_samples = 0;
	vector<t_annot_region*>* all_sig_regs = load_signal_regs_BED(signal_regions_BED_fp, n_loaded_samples);
	fprintf(stderr, "Loaded %d regions and %d samples.\n", (int)all_sig_regs->size(), n_loaded_samples);

	vector<char*>* sample_ids = buffer_file(sample_ids_fp);
	if ((int)sample_ids->size() != n_loaded_samples)
	{
		fprintf(stderr, "Could not match the number of loaded samples and the matrix columns: %d, %d\n", (int)sample_ids->size(), n_loaded_samples);
		exit(0);
	}

	fprintf(stderr, "Writing the reheadered file to %s\n", op_fp);

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < (int)all_sig_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", all_sig_regs->at(i_reg)->chrom,
			translate_coord(all_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(all_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			all_sig_regs->at(i_reg)->name);

		double* cur_reg_sigs = (double*)(all_sig_regs->at(i_reg)->data);

		for (int i_s = 0; i_s < n_loaded_samples; i_s++)
		{
			// Make sure we do not count the first 4 columns.
			fprintf(f_op, "\t%lf", cur_reg_sigs[i_s]);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

void extract_signal_regions_per_sample_list_per_4th_col_signals(char* signal_reg_BED_fp, char* selected_sample_ids_list_fp, char* op_fp)
{
	int n_loaded_samples = 0;
	vector<t_annot_region*>* all_sig_regs = load_signal_regs_BED(signal_reg_BED_fp, n_loaded_samples);
	fprintf(stderr, "Loaded %d regions and %d samples.\n", (int)all_sig_regs->size(), n_loaded_samples);

	char* header = load_header(signal_reg_BED_fp);
	vector<char*>* header_cols = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(header, "\t"));

	vector<char*>* selected_sample_ids = buffer_file(selected_sample_ids_list_fp);
	vector<char*>* config_sample_ids = new vector<char*>();
	vector<int>* per_config_sample_id_signal_reg_col_i = new vector<int>();
	for (int i_s = 0; i_s < (int)selected_sample_ids->size(); i_s++)
	{
		config_sample_ids->push_back(t_string::copy_me_str(selected_sample_ids->at(i_s)));

		int cur_sample_i_col_i = t_string::get_i_str(header_cols, selected_sample_ids->at(i_s));
		if (cur_sample_i_col_i == (int)header_cols->size())
		{
			fprintf(stderr, "Could not find %s\n", selected_sample_ids->at(i_s));
			exit(0);
		}

		per_config_sample_id_signal_reg_col_i->push_back(cur_sample_i_col_i);
	} // i_s loop.
	fprintf(stderr, "Extracting %d samples.\n", (int)per_config_sample_id_signal_reg_col_i->size());

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_s = 0; i_s < (int)per_config_sample_id_signal_reg_col_i->size(); i_s++)
	{
		fprintf(f_op, "\t%s", header_cols->at(per_config_sample_id_signal_reg_col_i->at(i_s)));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < (int)all_sig_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", all_sig_regs->at(i_reg)->chrom,
			translate_coord(all_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(all_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			all_sig_regs->at(i_reg)->name);

		double* cur_reg_sigs = (double*)(all_sig_regs->at(i_reg)->data);

		for (int i_s = 0; i_s < (int)per_config_sample_id_signal_reg_col_i->size(); i_s++)
		{
			// Make sure we do not count the first 4 columns.
			fprintf(f_op, "\t%lf", cur_reg_sigs[per_config_sample_id_signal_reg_col_i->at(i_s) - 4]);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

void concatenate_signal_regions(vector<char*>* signal_reg_BED_fps, char* op_fp)
{
	vector<t_annot_region*>** per_file_signal_regs = new vector<t_annot_region*>*[(int)signal_reg_BED_fps->size() + 2];
	vector<int>* per_file_n_samples = new vector<int>();
	vector<char*>* per_file_headers = new vector<char*>();
	for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
	{
		int n_loaded_samples = 0;
		per_file_signal_regs[i_f] = load_signal_regs_BED(signal_reg_BED_fps->at(i_f), n_loaded_samples);
		fprintf(stderr, "%s: %d regions, %d samples.\n", signal_reg_BED_fps->at(i_f), (int)per_file_signal_regs[i_f]->size(), n_loaded_samples);
		per_file_n_samples->push_back(n_loaded_samples);

		char* cur_file_header = load_header(signal_reg_BED_fps->at(i_f));
		per_file_headers->push_back(cur_file_header);
	} // i_f loop.

	vector<t_annot_region*>* agg_signal_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)per_file_signal_regs[0]->size(); i_reg++)
	{
		t_annot_region* cur_agg_reg = duplicate_region(per_file_signal_regs[0]->at(i_reg));
		t_annot_region** cur_agg_reg_per_file_regs = new t_annot_region*[(int)signal_reg_BED_fps->size() + 2];
		memset(cur_agg_reg_per_file_regs, 0, sizeof(t_annot_region*) * ((int)signal_reg_BED_fps->size() + 2));
		cur_agg_reg->data = cur_agg_reg_per_file_regs;
		agg_signal_regs->push_back(cur_agg_reg);
	} // i_reg loop.

	fprintf(stderr, "%d aggregate signal regions.\n", (int)agg_signal_regs->size());

	for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
	{
		fprintf(stderr, "Adding %s\n", signal_reg_BED_fps->at(i_f));
		vector<t_annot_region*>* cur_intersects = intersect_annot_regions(agg_signal_regs, per_file_signal_regs[i_f], true);

		for (int i_int = 0; i_int < (int)cur_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(cur_intersects->at(i_int)->data);
			t_annot_region* cur_agg_reg = int_info->src_reg;
			t_annot_region* cur_sig_reg = int_info->dest_reg;

			if (t_string::compare_strings(cur_agg_reg->name, cur_sig_reg->name))
			{
				t_annot_region** cur_agg_reg_per_file_regs = (t_annot_region**)(cur_agg_reg->data);
				cur_agg_reg_per_file_regs[i_f] = cur_sig_reg;
			}

			delete int_info;
		} // i_int loop.

		delete_annot_regions(cur_intersects);
	} // i_f loop

	fprintf(stderr, "Saving concatenated signal regions.\n");
	if (check_file(op_fp))
	{
		fprintf(stderr, "%s exists, will not overwrite.\n", op_fp);
		exit(0);
	}

	FILE* f_op = open_f(op_fp, "w");

	// Write the aggregate header.
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(per_file_headers->at(i_f), "\t");

		for (int i_t = 4; i_t < (int)toks->size(); i_t++)
		{
			fprintf(f_op, "\t%s", toks->at(i_t)->str());
		} // i_t loop.

		t_string::clean_tokens(toks);
	} // i_f loop.
	fprintf(f_op, "\n");

	// Go over all the regions.
	for (int i_reg = 0; i_reg < (int)agg_signal_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", agg_signal_regs->at(i_reg)->chrom,
			translate_coord(agg_signal_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(agg_signal_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			agg_signal_regs->at(i_reg)->name);

		t_annot_region** cur_agg_reg_per_file_regs = (t_annot_region**)(agg_signal_regs->at(i_reg)->data);

		// Go over all the regions and dump.
		for (int i_f = 0; i_f < (int)signal_reg_BED_fps->size(); i_f++)
		{
			if (cur_agg_reg_per_file_regs[i_f] == NULL)
			{
				fprintf(stderr, "Could not find sample %s for region %s:%d-%d\n", signal_reg_BED_fps->at(i_f),
					agg_signal_regs->at(i_reg)->chrom, agg_signal_regs->at(i_reg)->start, agg_signal_regs->at(i_reg)->end);
				exit(0);
			}

			double* cur_reg_signals = (double*)(cur_agg_reg_per_file_regs[i_f]->data);
			for (int i_s = 0; i_s < per_file_n_samples->at(i_f); i_s++)
			{
				fprintf(f_op, "\t%lf", cur_reg_signals[i_s]);
			} // i_s lo
		} // i_f loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

// Use this to quantify total signal summaries over a range of samples.
// "Official" RPKM/TPM/RPM computation is in expression_tools and it is much more extendable.
void extract_signal_features_per_regions(char* argv[], int argc, char* config_fp, char* op_fp)
{
	t_config* config = new t_config(config_fp, "\t");
	t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");

	char feats_bed_fp[1000];
	if (!get_cli_value(config, cli, "feature_regs_fp", "-feature_regs_fp", feats_bed_fp, true))
	{
		fprintf(stderr, "Could not read: feature_regs_fp");
		exit(0);
	}
	
	// Read the element length norm. flag.
	char length_norm_flag_str[1000];
	if (!get_cli_value(config, cli, "length_normalization", "-length_normalization", length_norm_flag_str, true))
	{
		fprintf(stderr, "Could not read: length_normalization");
		exit(0);
	}
	bool length_norm_flag = (atoi(length_norm_flag_str) == 1);

	// Read the lib. size norm. flag.
	char lib_size_norm_flag_str[1000];
	if (!get_cli_value(config, cli, "lib_size_normalization", "-lib_size_normalization", lib_size_norm_flag_str, true))
	{
		fprintf(stderr, "Could not read: lib_size_normalization");
		exit(0);
	}
	bool lib_size_norm_flag = (atoi(lib_size_norm_flag_str) == 1);

	if (lib_size_norm_flag)
	{
		fprintf(stderr, "Performing lib. size normalization.\n");
	}
	else
	{
		fprintf(stderr, "NOT performing lib. size normalization.\n");
	}

	if (length_norm_flag)
	{
		fprintf(stderr, "Performing element length normalization.\n");
	}
	else
	{
		fprintf(stderr, "NOT performing element length normalization.\n");
	}

	// Load the set the signals.
	vector<t_annot_region*>* feat_regs_w_intervals = load_Regions_as_Interval(feats_bed_fp);
	for (int i_reg = 0; i_reg < (int)feat_regs_w_intervals->size(); i_reg++)
	{
		double* cur_reg_per_sample_signals = new double[(int)sample_entries->size() + 2];
		memset(cur_reg_per_sample_signals, 0, sizeof(double) * (int)sample_entries->size());
		feat_regs_w_intervals->at(i_reg)->data = cur_reg_per_sample_signals;
	} // i_reg loop.

	vector<char*>* chr_ids = get_chr_ids(feat_regs_w_intervals);

	fprintf(stderr, "Extracting signal features on %d regions from %s.\n", (int)feat_regs_w_intervals->size(), feats_bed_fp);

	// Total signal on each sample.
	double* per_sample_total_signal = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_signal_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal_in_regs, 0, sizeof(double) * sample_entries->size());
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_feat_regs = get_regions_per_chromosome(feat_regs_w_intervals, chr_ids->at(i_chr));
		fprintf(stderr, "Processing %d regions on %s\n", (int)cur_chr_feat_regs->size(), chr_ids->at(i_chr));

		int signal_starting_col_i = 2;

		for (int i_s = 0; i_s < (int)sample_entries->size(); i_s++)
		{
			// Pool the pileups/BGRs.
			int* cur_sample_cur_chr_covg = NULL;
			int l_pileup = -1;

			// Go over all the tracks for this sample and pool them.
			for (int col_i = signal_starting_col_i; col_i < (int)sample_entries->at(i_s)->size(); col_i++)
			{
				int l_loaded_track_covg = 0;
				int* cur_sample_cur_chr_cur_track_covg = load_int_signal_covg_per_directory_chr_id(sample_entries->at(i_s)->at(col_i), chr_ids->at(i_chr), l_loaded_track_covg);

				// If the sample is not loaded, do not process.
				if (cur_sample_cur_chr_cur_track_covg == NULL)
				{
					fprintf(stderr, "Could not load signal from %s/%s; Skipping.\n", sample_entries->at(i_s)->at(col_i), chr_ids->at(i_chr));
					continue;
				}

				// Pool the current signal track.
				if (cur_sample_cur_chr_covg == NULL)
				{
					// Copy for the first track.
					l_pileup = l_loaded_track_covg + 1000;
					cur_sample_cur_chr_covg = new int[l_pileup];
					memset(cur_sample_cur_chr_covg, 0, sizeof(int) * l_pileup);
					for (int i = 1; i <= MIN(l_pileup, l_loaded_track_covg); i++)
					{
						cur_sample_cur_chr_covg[i] = cur_sample_cur_chr_cur_track_covg[i];
					} // i loop.
				}
				else
				{
					for (int i = 1; i <= MIN(l_pileup, l_loaded_track_covg); i++)
					{
						cur_sample_cur_chr_covg[i] += cur_sample_cur_chr_cur_track_covg[i];
					} // i loop.
				}

				delete[] cur_sample_cur_chr_cur_track_covg;
			} // col_i loop.

			// Update the total signal on this sample.
			for (int i = 1; i <= l_pileup; i++)
			{
				per_sample_total_signal[i_s] += cur_sample_cur_chr_covg[i];
			} // i loop.

			fprintf(stderr, "Sample %s: Total signal (chrom: %s): %.1f\n", sample_entries->at(i_s)->at(0), chr_ids->at(i_chr), per_sample_total_signal[i_s]);

			// Extract the total signal for all the regions on this chromosome.
			for (int i_reg = 0; i_reg < (int)cur_chr_feat_regs->size(); i_reg++)
			{
				double* cur_reg_per_sample_signals = (double*)(cur_chr_feat_regs->at(i_reg)->data);

				vector<t_annot_region*>* cur_reg_int_regs = (vector<t_annot_region*>*)(cur_chr_feat_regs->at(i_reg)->intervals);
				vector<t_annot_region*>* merged_cur_reg_int_regs = merge_annot_regions(cur_reg_int_regs, 0);

				// Set the coverage to the dbl score of the region.
				double cur_reg_ints_coverage = coverage(merged_cur_reg_int_regs);
				cur_chr_feat_regs->at(i_reg)->dbl_score = cur_reg_ints_coverage;

				// Update the total signal in this region's intervals for the current sample; go over all the intervals of this element and update the signal.
				cur_reg_per_sample_signals[i_s] = 0;
				for (int i_int = 0; i_int < (int)merged_cur_reg_int_regs->size(); i_int++)
				{
					for (int i = merged_cur_reg_int_regs->at(i_int)->start; i <= merged_cur_reg_int_regs->at(i_int)->end; i++)
					{
						if (i <= l_pileup)
						{
							cur_reg_per_sample_signals[i_s] += cur_sample_cur_chr_covg[i];
						}
						//else
						//{
						//	fprintf(stderr, "FATAL ERROR: Position further than chromosome length; %s: %d-%d (%d)\n",
						//		merged_cur_reg_int_regs->at(i_int)->chrom, merged_cur_reg_int_regs->at(i_int)->start, merged_cur_reg_int_regs->at(i_int)->end,
						//		l_pileup);

						//	exit(0);
						//}
					} // i loop.
				} // interval loop.

				  // Update the current total signal
				per_sample_total_signal_in_regs[i_s] += cur_reg_per_sample_signals[i_s];

				fprintf(stderr, "Sample %s: %s: %d exons (%d nucleotides in total): %.2f signal\r",
					sample_entries->at(i_s)->at(0),
					cur_chr_feat_regs->at(i_reg)->name,
					(int)merged_cur_reg_int_regs->size(),
					(int)cur_reg_ints_coverage,
					cur_reg_per_sample_signals[i_s]);
				if (i_reg % 200 == 0)
				{
					fprintf(stderr, "\n");
				}

				delete_annot_regions(merged_cur_reg_int_regs);
			} // i_reg loop.

			delete[] cur_sample_cur_chr_covg;
		} // i_s loop.	
	} // i_chr loop

	  // Dump the normalized matrix of signal levels.
	FILE* f_op = open_f(op_fp, "w");

	// Dump the header.
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME ID TYPE\tFEATURE_LENGTH");

	for (int i_s = 0; i_s < (int)sample_entries->size(); i_s++)
	{
		fprintf(f_op, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.

	fprintf(f_op, "\n");

	// Dump the signals for all the regions.
	for (int i_reg = 0; i_reg < (int)feat_regs_w_intervals->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s\t%.1f", feat_regs_w_intervals->at(i_reg)->chrom,
			translate_coord(feat_regs_w_intervals->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(feat_regs_w_intervals->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			feat_regs_w_intervals->at(i_reg)->name,
			feat_regs_w_intervals->at(i_reg)->dbl_score);

		double* cur_reg_per_sample_signal = (double*)(feat_regs_w_intervals->at(i_reg)->data);
		for (int i_s = 0; i_s < (int)sample_entries->size(); i_s++)
		{
			double cur_reg_normalized_signal = 0;

			if (lib_size_norm_flag)
			{
				cur_reg_normalized_signal = 1000 * 1000 * (cur_reg_per_sample_signal[i_s] / per_sample_total_signal[i_s]);
			}
			else
			{
				cur_reg_normalized_signal = cur_reg_per_sample_signal[i_s];
			}

			// Length normalization, if requested.
			if (length_norm_flag)
			{
				cur_reg_normalized_signal /= (feat_regs_w_intervals->at(i_reg)->dbl_score / 1000);
			}

			fprintf(f_op, "\t%lf", cur_reg_normalized_signal);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);

	// Dump the signal statistics in the regions.
	FILE* f_per_sample_sig_stats = open_f("quant_signal_stats.txt", "w");
	for (int i_s = 0; i_s < (int)sample_entries->size(); i_s++)
	{
		fprintf(f_per_sample_sig_stats, "%.2f\t%.2f\t%.3f\n",
			per_sample_total_signal_in_regs[i_s], per_sample_total_signal[i_s],
			per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
	} // i_s loop.
	fclose(f_per_sample_sig_stats);
}

vector<t_annot_region*>* load_signal_regs_BED(char* signal_regions_BED_fp, int& n_loaded_samples, int matrix_starting_col_i)
{
	vector<t_annot_region*>* regs_w_lines = load_BED_with_line_information(signal_regions_BED_fp);
	vector<t_annot_region*>* signal_regs = new vector<t_annot_region*>();

	n_loaded_samples = -1;
	for (int i_reg = 0; i_reg < (int)regs_w_lines->size(); i_reg++)
	{
		t_annot_region* cur_sig_reg = duplicate_region(regs_w_lines->at(i_reg));
		char* cur_reg_line = (char*)(regs_w_lines->at(i_reg)->data);

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");

		double* cur_sig_reg_sigs = new double[(int)toks->size() + 2];
		for (int i_tok = matrix_starting_col_i; i_tok < (int)toks->size(); i_tok++)
		{
			cur_sig_reg_sigs[i_tok - matrix_starting_col_i] = atof(toks->at(i_tok)->str());
		} // i_tok loop.
		cur_sig_reg->name = t_string::copy_me_str(toks->at(3)->str());
		cur_sig_reg->data = cur_sig_reg_sigs;

		if (n_loaded_samples == -1)
		{
			n_loaded_samples = (int)toks->size() - matrix_starting_col_i;
		}
		else if (n_loaded_samples != (int)toks->size() - matrix_starting_col_i)
		{
			fprintf(stderr, "Could not match the number of loaded samples: %d, %d\n", n_loaded_samples, (int)toks->size() - matrix_starting_col_i);
			exit(0);
		}

		signal_regs->push_back(cur_sig_reg);

		t_string::clean_tokens(toks);
		delete[] cur_reg_line;
	} // i_reg loop.

	delete_annot_regions(regs_w_lines);

	return(signal_regs);
}

vector<t_annot_region*>* load_signal_regs_BED(char* signal_regions_BED_fp, int& n_loaded_samples)
{
	vector<t_annot_region*>* regs_w_lines = load_BED_with_line_information(signal_regions_BED_fp);
	vector<t_annot_region*>* signal_regs = new vector<t_annot_region*>();

	n_loaded_samples = -1;
	for (int i_reg = 0; i_reg < (int)regs_w_lines->size(); i_reg++)
	{
		t_annot_region* cur_sig_reg = duplicate_region(regs_w_lines->at(i_reg));
		char* cur_reg_line = (char*)(regs_w_lines->at(i_reg)->data);

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");

		double* cur_sig_reg_sigs = new double[(int)toks->size() + 2];
		for (int i_tok = 4; i_tok < (int)toks->size(); i_tok++)
		{
			cur_sig_reg_sigs[i_tok - 4] = atof(toks->at(i_tok)->str());
		} // i_tok loop.
		cur_sig_reg->name = t_string::copy_me_str(toks->at(3)->str());
		cur_sig_reg->data = cur_sig_reg_sigs;

		if (n_loaded_samples == -1)
		{
			n_loaded_samples = (int)toks->size() - 4;
		}
		else if (n_loaded_samples != (int)toks->size() - 4)
		{
			fprintf(stderr, "Could not match the number of loaded samples: %d, %d\n", n_loaded_samples, (int)toks->size() - 4);
			exit(0);
		}

		signal_regs->push_back(cur_sig_reg);

		t_string::clean_tokens(toks);
		delete[] cur_reg_line;
	} // i_reg loop.

	delete_annot_regions(regs_w_lines);

	return(signal_regs);
}

void remove_duplicate_samples_from_signal_matrix(char* signal_regs_bed_fp, char* op_signal_regs_bed_fp)
{
	int n_loaded_samples = 0;
	vector<t_annot_region*>* signal_regs = load_signal_regs_BED(signal_regs_bed_fp, n_loaded_samples);
	char* header = load_header(signal_regs_bed_fp);
	t_string_tokens* header_toks = t_string::tokenize_by_chars(header, "\t");
	vector<char*>* sample_ids = t_string::copy_tokens_2_strs(header_toks, 4, -1);
	if ((int)sample_ids->size() != n_loaded_samples)
	{
		fprintf(stderr, "Loaded samples do not match the header.\n");
		exit(0);
	}

	vector<char*>* uniq_sample_ids = t_string::get_unique_entries(sample_ids);
	fprintf(stderr, "Mapping to %d unique samples out of %d samples.\n", (int)uniq_sample_ids->size(), (int)sample_ids->size());

	FILE* f_op = open_f(op_signal_regs_bed_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME ID TYPE");
	for (int i_s = 0; i_s < (int)uniq_sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", uniq_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < (int)signal_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s",
			signal_regs->at(i_reg)->chrom,
			translate_coord(signal_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(signal_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			signal_regs->at(i_reg)->name);

		double* cur_reg_signals = (double*)(signal_regs->at(i_reg)->data);

		for (int i_s = 0; i_s < (int)uniq_sample_ids->size(); i_s++)
		{
			int cur_uniq_sample_i_s = t_string::get_i_str(sample_ids, uniq_sample_ids->at(i_s));

			if (cur_uniq_sample_i_s == (int)sample_ids->size())
			{
				fprintf(stderr, "Could not find sample id %s\n", uniq_sample_ids->at(i_s));
				exit(0);
			}

			fprintf(f_op, "\t%lf", cur_reg_signals[cur_uniq_sample_i_s]);
		} // i_s loop.
		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}


void rank_normalize_signal_matrix(char* config_fp, char* signal_matrix_fp, char* op_fp)
{
	t_config* config = new t_config(config_fp, "\t");

	vector<char*>* matrix_lines = buffer_file(signal_matrix_fp);

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");

	vector<double>** per_sample_signals = new vector<double>*[(int)sample_entries->size() + 2];
	for (int i_s = 0; i_s < (int)sample_entries->size(); i_s++)
	{
		per_sample_signals[i_s] = new vector<double>();
	} // i_s loop.

	  // Parse the sample id's.
	t_string_tokens* sample_toks = t_string::tokenize_by_chars(matrix_lines->at(0), "\t");
	for (int i_t = 4; i_t < (int)sample_toks->size(); i_t++)
	{
		fprintf(stderr, "Sample %d: %s\n", i_t - 4, sample_toks->at(i_t)->str());
	}

	for (int i_l = 1; i_l < (int)matrix_lines->size(); i_l++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(matrix_lines->at(i_l), "\t");

		for (int i_t = 4; i_t < (int)toks->size(); i_t++)
		{
			int i_s = i_t - 4;

			double cur_exp = atof(toks->at(i_t)->str());
			per_sample_signals[i_s]->push_back(cur_exp);
		} // i_t loop.
	} // i_l loop.

	fprintf(stderr, "Loaded %d regions with signals.\n", (int)per_sample_signals[0]->size());

	for (int i_s = 0; i_s < (int)sample_entries->size(); i_s++)
	{
		vector<t_signal_node*>* nodes = new vector<t_signal_node*>();

		double cur_sample_total_sig = 0;
		for (int i_reg = 0; i_reg < (int)per_sample_signals[i_s]->size(); i_reg++)
		{
			t_signal_node* node = new t_signal_node();
			node->signal = per_sample_signals[i_s]->at(i_reg);
			node->i_reg = i_reg;
			nodes->push_back(node);

			cur_sample_total_sig += per_sample_signals[i_s]->at(i_reg);
		} // i_reg loop.

		  //for (int i_reg = 0; i_reg < per_sample_signals[i_s]->size(); i_reg++)
		  //{
		  //	per_sample_signals[i_s]->at(i_reg) /= cur_sample_total_sig;
		  //	per_sample_signals[i_s]->at(i_reg) *= 1000;
		  //} // i_reg loop.

		sort(nodes->begin(), nodes->end(), sort_signal_nodes_per_increasing_signal);

		for (int rank_i = 0; rank_i < (int)nodes->size(); rank_i++)
		{
			per_sample_signals[i_s]->at(nodes->at(rank_i)->i_reg) = rank_i;

			delete nodes->at(rank_i);
		} // rank_i loop.

		delete nodes;
	} // i_s loop.

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%s\n", matrix_lines->at(0));
	for (int i_l = 1; i_l < (int)matrix_lines->size(); i_l++)
	{
		int i_reg = i_l - 1;

		t_string_tokens* toks = t_string::tokenize_by_chars(matrix_lines->at(i_l), "\t");
		fprintf(f_op, "%s\t%s\t%s\t%s", toks->at(0)->str(), toks->at(1)->str(), toks->at(2)->str(), toks->at(3)->str());
		for (int i_t = 4; i_t < (int)toks->size(); i_t++)
		{
			int i_s = i_t - 4;

			fprintf(f_op, "\t%lf", per_sample_signals[i_s]->at(i_reg));
		} // i_t loop.

		fprintf(f_op, "\n");
	} // i_l loop.
	fclose(f_op);
}

double* load_signal_covg_per_signal_file(const char* cur_dat_fp,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded)
{
	double* covg_signal = NULL;

	reads_loaded = false;

	// Search for mapped reads.
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for BGR.
	if (check_file(cur_dat_fp) &&
		t_string::ends_with(cur_dat_fp, ".bgr.gz"))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	if (check_file(cur_dat_fp) &&
		t_string::ends_with(cur_dat_fp, ".bgr"))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

int* load_int_signal_covg_per_directory_chr_id(const char* dat_dir, char* chr_id, int& l_loaded_covg)
{
	int* covg_signal = NULL;

	// Search for pileup.
	char cur_dat_fp[1000];
	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for mapped reads.
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt", dat_dir, chr_id);
	if (!check_file(cur_dat_fp))
	{
		sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt.gz", dat_dir, chr_id);
	}

	int l_fragment = 0;
	if (check_file(cur_dat_fp))
	{
		//reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		double* dbl_covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			dbl_covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		int* covg_signal = new int[l_loaded_covg + 10];
		for (int i = 1; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.
		delete[] dbl_covg_signal;

		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/%s.bbgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bbgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bbgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

double* load_signal_covg_per_directory_chr_id(const char* dat_dir,
	char* chr_id,
	int l_fragment,
	char* cur_dat_fp,
	int& l_loaded_covg,
	bool& reads_loaded)
{
	if (dat_dir == NULL)
	{
		return(NULL);
	}

	double* covg_signal = NULL;

	reads_loaded = false;

	// Search for mapped reads.
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for mapped reads.
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for pileup.
	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		int* int_covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		covg_signal = new double[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (double)(int_covg_signal[i]);
		} // i loop.

		delete[] int_covg_signal;

		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		int* int_covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		covg_signal = new double[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (double)(int_covg_signal[i]);
		} // i loop.

		delete[] int_covg_signal;

		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/%s.bbgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/%s.bbgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s.bin", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_binary_profile(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_binary_profile(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_logR.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_binary_profile(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

// Following bins a profile into a 1-based array.
// It also sets the dbl_score for each quantified region such that the requested signal is the score.
double* bin_profile(double* signal_profile, int l_profile, // Signal profile that will be binned.
					vector<t_annot_region*>* target_regions, // The target regions where binning will be done, NULL means the whole profile is used as target.
					vector<t_annot_region*>* quantified_regions, // Empty list of quantified regions to be filled.
					bool flag_span_region_borders,	// This says whether we would like to go across the ends of consecutive windows while enumerating them.
					int _l_bin,						// l_bin = 0 makes per region signal computation.
					int& n_bins)					// shows the number of bins in the array. Must match the number of quantified regions.
{
	n_bins = 0;

	// Set the bin length to be used.
	int l_bin = 0;
	if(_l_bin > 0)
	{
		l_bin = _l_bin;
	}
	else
	{
		// If bin length is 0, we are asked to do per region signal, do not span borders.
		flag_span_region_borders = false;
	}
	
	// Set the target regions.
	vector<t_annot_region*>* all_target_regions = NULL;

	if(target_regions == NULL)
	{
		if (__DUMP_SIGNAL_TRACK_MSGS__)
		{
			fprintf(stderr, "Using the whole region in binning.\n");
		}

		all_target_regions = new vector<t_annot_region*>();
		t_annot_region* whole_profile_reg = get_empty_region();
		whole_profile_reg->start = 1;
		whole_profile_reg->end = l_profile;
		all_target_regions->push_back(whole_profile_reg);
	}
	else
	{
		all_target_regions = target_regions;
	}

	// Get an upper bound on the # of bins.
	unsigned int l_target_covg = (unsigned int)(coverage(all_target_regions));

	if(_l_bin == 0)
	{
		n_bins = (int)all_target_regions->size() + 10;
	}
	else
	{
		n_bins = (l_target_covg / l_bin) + 1;
	}

	if (__DUMP_SIGNAL_TRACK_MSGS__)
	{
		fprintf(stderr, "Generating at most %d bins signal.\n", n_bins);
	}

	// Initialize the binned signal to all 0's.
	double* binned_signal = new double[n_bins + 2];
	memset(binned_signal, 0, sizeof(double) * (n_bins+1));

	//// The targets are sorted, for what reason?
	//sort(all_target_regions->begin(), all_target_regions->end(), sort_regions);

	// Following are not dependent on which region we are at.
	int bin_i = 0; // This is the index which points to the bin that we are filling.
	double cur_bin_signal = 0; // The signal in the current bin.
	int cur_bin_posn = 0; // The index that points to the position that we are at in the bin.

	// The idea is to go over all the regions and trace each bin.
	char* cur_bin_chrom = all_target_regions->at(0)->chrom;
	int cur_bin_start = all_target_regions->at(0)->start;
	for(int i_reg = 0; i_reg < (int)all_target_regions->size(); i_reg++)
	{
		if (__DUMP_SIGNAL_TRACK_MSGS__)
		{
			fprintf(stderr, "%d. (%d) region                 \r", i_reg, (int)all_target_regions->size());
		}

		// 0 bin length corresponds to per target total signal computes.
		if(_l_bin == 0)
		{
			flag_span_region_borders = false;
			l_bin = all_target_regions->at(i_reg)->end - all_target_regions->at(i_reg)->start + 1;
		}

		// If we do not want to span the region borders, reset the bin position and signal.
		if(!flag_span_region_borders)
		{
			/*
			if(cur_bin_signal > 0)
			{
				fprintf(stderr, "The region border signal reset!\n");
			}
			*/

			cur_bin_signal = 0;
			cur_bin_posn = 0;
			cur_bin_chrom = all_target_regions->at(i_reg)->chrom;
			cur_bin_start = all_target_regions->at(i_reg)->start;
		}

if(__DUMP_SIGNAL_TRACK_MSGS__)
{
		fprintf(stderr, "Processing bins in %s:%d-%d (%d)\n", all_target_regions->at(i_reg)->chrom, 
				all_target_regions->at(i_reg)->start, all_target_regions->at(i_reg)->end, cur_bin_posn);
}

		// Go over all the bins over this region.
		for(int posn = all_target_regions->at(i_reg)->start; posn <= all_target_regions->at(i_reg)->end; posn++)
		{
			if(posn <= l_profile)
			{
				cur_bin_signal += signal_profile[posn];
				cur_bin_posn++;
			}

			// If we finished the current bin, set the binned signal and reset the signal and bin position.
			// Otherwise we hit the end of the region, which broke the above loop. Do not assign, yet.
			if(cur_bin_posn == l_bin)
			{
				// This is the region whose signal is just added; make sure this is done before resetting cur_bin_signal and cur_bin_posn
				if(quantified_regions != NULL)
				{
					t_annot_region* new_quant_reg = get_empty_region();
					new_quant_reg->chrom = t_string::copy_me_str(cur_bin_chrom);
					new_quant_reg->start = cur_bin_start;
					new_quant_reg->end = cur_bin_start + l_bin - 1;
					new_quant_reg->strand = '+';
					new_quant_reg->dbl_score = cur_bin_signal;

					quantified_regions->push_back(new_quant_reg);

if(__DUMP_SIGNAL_TRACK_MSGS__)
{
					fprintf(stderr, "Added %s:%d-%d\n", cur_bin_chrom, cur_bin_start, cur_bin_start + l_bin - 1);getc(stdin);
}
				}

				// Set the current value.
				bin_i++; // Update the bin index; also points to the number of bins.
				binned_signal[bin_i] = cur_bin_signal;

				// Note that this reset must be done after adding the region; do not reset before adding the region.
				cur_bin_signal = 0;
				cur_bin_posn = 0;

				if(posn != cur_bin_start + l_bin - 1)
				{
					fprintf(stderr, "Sanity check failed; %d, %d\n", 
						posn, cur_bin_start + l_bin - 1);
					exit(0);
				}

				// This is the new bin's start.
				cur_bin_chrom = all_target_regions->at(i_reg)->chrom;
				cur_bin_start = posn + 1; // Note that the new start is right after this position. posn is updated in the for loop.

				//if(binned_signal[bin_i] < 0)
				//{
				//	fprintf(stderr, "Binned signal smaller than 0 signal.\n");
				//	exit(0);
				//}
			} // bin ending check.

			if(bin_i > n_bins)
			{
				fprintf(stderr, "bin_i (%d) is larger than the max number of bins expected (%d).\n", bin_i, n_bins);
				exit(0);
			}
		} // posn loop.
	} // i_reg loop.

	// Add the final final value to the array, if generated; in case we are spanning regions.
	// This is only added if we are spanning region borders; this is by choice, we do not have to treat it like this.
	if(flag_span_region_borders &&
		cur_bin_posn > 0)
	{
		// Add the final region.
		if(quantified_regions != NULL)
		{
			t_annot_region* new_quant_reg = get_empty_region();
			new_quant_reg->chrom = t_string::copy_me_str(cur_bin_chrom);
			new_quant_reg->start = cur_bin_start;
			new_quant_reg->end = cur_bin_start + l_bin - 1;
			new_quant_reg->strand = '+';
			new_quant_reg->dbl_score = cur_bin_signal;
			quantified_regions->push_back(new_quant_reg);
		}

		bin_i++;
		binned_signal[bin_i] = cur_bin_signal;
	}

	n_bins = bin_i;

	if(quantified_regions != NULL &&
		(int)quantified_regions->size() != n_bins)
	{
		fprintf(stderr, "The sanity check failed: Quantified region size does not match counted bin number: %d, %d\n", 
			(int)quantified_regions->size(), n_bins);
		exit(0);
	}

	return(binned_signal);
}

double* concat_profile2_2_profile1(double* prof1, int l_prof1, 
									double* prof2, int l_prof2, 
									int& l_new_track)
{
	double* conc_prof = NULL;
	if(prof2 == NULL)
	{
		l_new_track = l_prof1;

		conc_prof = new double[l_new_track + 1];

		for(int i = 1; i <= l_prof1; i++)
		{
			conc_prof[i] = prof1[i];
		} // i loop.
	}
	else if(prof1 == NULL)
	{
		l_new_track = l_prof2;

		conc_prof = new double[l_new_track + 1];

		for(int i = 1; i <= l_prof2; i++)
		{
			conc_prof[i] = prof2[i];
		} // i loop.

	}
	else
	{
		l_new_track = l_prof1 + l_prof2;

		conc_prof = new double[l_new_track + 1];
		for(int i = 1; i <= l_prof1; i++)
		{
			conc_prof[i] = prof1[i];
		} // i loop.

		for(int i = 1; i <= l_prof2; i++)
		{
			conc_prof[i + l_prof1] = prof2[i];
		} // i loop.
	}

	return(conc_prof);
}

// This function is a wrapper around bin_profile to normalize and log the signal.
double* bin_and_normalize_target_profile(double* cur_chr_signal_buff, 
										vector<t_annot_region*>* cur_chr_target_regs, 
										vector<t_annot_region*>* cur_chr_quantified_regions,
										bool flag_span_target_ends,
										double log_base,
										int l_track, 
										int l_bin, 
										int& l_binned_track)
{
	// Bin the profile, delete the buffer.
	l_binned_track = 0;

	double* cur_binned_profile = bin_profile(cur_chr_signal_buff, l_track,
											cur_chr_target_regs,
											cur_chr_quantified_regions,
											flag_span_target_ends,
											l_bin, 
											l_binned_track);

	if(log_base > 0)
	{
		fprintf(stderr, "Converting binned signal to logarithm.\n");
		for(int bin_i = 1; bin_i <= l_binned_track; bin_i++)
		{
			// We can have negative scores now; for example, 
			//if(cur_binned_profile[bin_i] < 0)
			//{
			//	fprintf(stderr, "%d. bin value is negative.\n", bin_i);
			//	exit(0);
			//}

			// Must be careful with the data transformations.
			cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i]+1) / xlog(log_base);
		} // bin_i loop.

		if(cur_chr_quantified_regions != NULL)
		{
			for(int i_reg = 0; i_reg < (int)cur_chr_quantified_regions->size(); i_reg++)
			{
				cur_chr_quantified_regions->at(i_reg)->dbl_score = xlog(cur_chr_quantified_regions->at(i_reg)->dbl_score+1) / xlog(log_base);
			} // i_reg loop.
		} // target regions existence check.
	}
	else
	{
		fprintf(stderr, "Using linear signals.\n");
	}

	return(cur_binned_profile);
}

/****************************************************************************************************************************************************
This function generates a multitrack profile and conveniently returns it in a matrix. It handles regions, signal tracks, reads. 
For base-resolution signal, use l_bin=1; target_regions=NULL.
For gw binned signal, use l_bin>0; target_regions=NULL.
For targeted binned signal, use l_bin>0; target_regions~=NULL.
****************************************************************************************************************************************************/
double** generate_multitrack_signal_profiles(vector<char*>* preprocessed_reads_dirs,
											vector<char*>* bin_track_dirs,
											vector<char*>* region_fps,
											vector<char*>* chr_ids_2_use, 											
											vector<t_annot_region*>* _target_regions, // null means genomewide signal; i.e., uses the length of the first preprocessed reads track.
											vector<t_annot_region*>* quantified_regions, // This is the list of regions for which signals are generated. The size must be the same as l_track.
											bool flag_span_target_region_ends, 
											int l_bin, // <=0 for no binning and per target signal.
											double log_base, // 0 means linear signal return.
											bool use_random_track, // Adds a random track to the end if this is true.
											int& l_multitrack)	// Length of the track.
{
	vector<t_annot_region*>* target_regions = NULL;
	if(_target_regions != NULL)
	{
		target_regions = merge_annot_regions(_target_regions, 1);
		fprintf(stderr, "Merged targets into %d regions.\n", (int)target_regions->size());
	}

	double target_covg = coverage(target_regions);
	fprintf(stderr, "Generating multitrack signal profiles using %d read dirs, %d binary tracks, %d regions over %d target regions (%d bps) with bin length of %d and log base of %.4f.\n", 
		(int)preprocessed_reads_dirs->size(), (int)bin_track_dirs->size(), (int)region_fps->size(), (int)target_regions->size(), (int)target_covg, l_bin, log_base);

	int n_tracks = (int)(preprocessed_reads_dirs->size() + (int)bin_track_dirs->size() + (int)region_fps->size());

	// Load the chromosome id's.
	vector<char*>* chr_ids = NULL;
	if(chr_ids_2_use == NULL)
	{
		if(target_regions != NULL)
		{
			chr_ids = get_chr_ids(target_regions);
		}
		else if(preprocessed_reads_dirs != NULL &&
			(int)preprocessed_reads_dirs->size() > 0)
		{
			// Load chromosome id's from first preprocessed reads directory.
			char chr_ids_fp[1000];
			sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dirs->at(0));
			chr_ids = buffer_file(chr_ids_fp);
		}
		else if(region_fps != NULL)
		{
			vector<t_annot_region*>* regs = load_BED(region_fps->at(0));
			chr_ids = get_chr_ids(regs);
			delete_annot_regions(regs);
		}
		else
		{
			fprintf(stderr, "Could not load chromosome id's from nowhere, check the inputs.\n");
			return(NULL);
		}
	}
	else
	{
		chr_ids = chr_ids_2_use;
	}

	if(use_random_track)
	{
		fprintf(stderr, "Using a randomized control track.\n");
		n_tracks++;
	}

	fprintf(stderr, "Loading %d preprocessed dirs, %d binary track dirs, %d region files.\n", (int)preprocessed_reads_dirs->size(),
		(int)bin_track_dirs->size(), (int)region_fps->size());

	if(l_bin == 0)
	{
		fprintf(stderr, "Generating whole genome signals for %d chromosomes.\n", (int)chr_ids_2_use->size());
	}

	// Do not do tag extension.
	int l_ext_tag = 0;

	double** gw_tracks = new double*[n_tracks+2];
	memset(gw_tracks, 0, sizeof(double*) * n_tracks);
	int* l_tracks = new int[n_tracks];

	// Process all the chromosomes.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int min_n_bins = -1;

		fprintf(stderr, "Processing %s\n", chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_target_regs = NULL;
		if(target_regions != NULL)
		{
			cur_chr_target_regs = get_regions_per_chromosome(target_regions, chr_ids->at(i_chr));
			fprintf(stderr, "%d target regions.\n", (int)cur_chr_target_regs->size());
		}

		double** cur_chr_tracks = new double*[n_tracks+2];
		int track_i = 0;

		// Allocate and reset per track regions.
		vector<t_annot_region*>** per_track_quantified_regions = new vector<t_annot_region*>*[n_tracks];
		for(int i_t = 0; i_t < n_tracks; i_t++)
		{
			per_track_quantified_regions[i_t] = NULL;
		} // i_t loop

		// First, process the list of processed reads.
		enum{DOM_LINEAR_AVG, DOM_LINEAR_TOTAL, DOM_LOG};
		for(int i_dir = 0; i_dir < (int)preprocessed_reads_dirs->size(); i_dir++)
		{
			//int domain_type = -1;
			enum{DOM_LINEAR_AVG, DOM_LINEAR_TOTAL, DOM_LOG};

			fprintf(stderr, "Loading signal profile from %s\n", preprocessed_reads_dirs->at(i_dir));
			int l_buff = 300*1000*1000;
			double* cur_chr_signal_buff = new double[l_buff + 2];
			int l_track = 0;

			//bool cur_track_linear = false;
			//bool cur_track_log = false;

			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dirs->at(i_dir), chr_ids->at(i_chr));
			if(check_file(cur_chr_reads_fp))
			{
				buffer_per_nucleotide_profile_no_buffer(cur_chr_reads_fp, l_ext_tag, cur_chr_signal_buff, NULL, NULL, l_buff, l_track);
			}
			else
			{
				fprintf(stderr, "Could not find %s\n", cur_chr_reads_fp);
				break;
			}

			int l_binned_track = 0;
			double* cur_binned_profile = NULL;
			vector<t_annot_region*>* cur_track_quantified_regions = new vector<t_annot_region*>();
			cur_binned_profile = bin_and_normalize_target_profile(cur_chr_signal_buff, 
																	cur_chr_target_regs, 
																	cur_track_quantified_regions,
																	flag_span_target_region_ends,
																	log_base,
																	l_track, 
																	l_bin, 
																	l_binned_track);

			fprintf(stderr, "%d quantified regions.\n", (int)cur_track_quantified_regions->size());

			// Set the quantified regions.
			sort(cur_track_quantified_regions->begin(), cur_track_quantified_regions->end(), sort_regions);
			per_track_quantified_regions[track_i] = cur_track_quantified_regions;

			//char cur_quant_fp[1000];
			//sprintf(cur_quant_fp, "%d_quant.bed", track_i);
			//dump_BED(cur_quant_fp, cur_track_quantified_regions);

#undef __MANUAL_CHECK__
#ifdef __MANUAL_CHECK__
			// Do a sanity check on the computed signals.
			// We are guaranteed to have the quantified regions.
			for(int i_r = 0; i_r < (int)per_track_quantified_regions[track_i]->size(); i_r++)
			{
				fprintf(stderr, "Checking %s:%d-%d: %lf                      \r", per_track_quantified_regions[track_i]->at(i_r)->chrom, 
						per_track_quantified_regions[track_i]->at(i_r)->start, per_track_quantified_regions[track_i]->at(i_r)->end,
						cur_binned_profile[i_r+1]);

				if(per_track_quantified_regions[track_i]->at(i_r)->start > per_track_quantified_regions[track_i]->at(i_r)->end)
				{
					fprintf(stderr, "\nThere is a problem with this region, skipping.\n");
					continue;
				}

				int i = per_track_quantified_regions[track_i]->at(i_r)->start;
				int l_cur_bin = l_bin;
				if(l_bin == 0)
				{
					l_cur_bin = (per_track_quantified_regions[track_i]->at(i_r)->end - per_track_quantified_regions[track_i]->at(i_r)->start);
				}

				while(1)
				{
					if(i + l_cur_bin > per_track_quantified_regions[track_i]->at(i_r)->end)
					{
						break;
					}

					double cur_total_sig = 0;
					int cur_bin_end = i+l_cur_bin;
					for(int pos = i; pos <= cur_bin_end; pos++)
					{
						cur_total_sig += cur_chr_signal_buff[pos];
					} // loop.
				
					cur_total_sig = (log(cur_total_sig+1) / log(log_base));

					if(MAX(cur_binned_profile[i_r+1], cur_total_sig) > 0.01 &&
						fabs(cur_binned_profile[i_r+1] - cur_total_sig) / MAX(cur_binned_profile[i_r+1], cur_total_sig) > 0.01)
					{
						fprintf(stderr, "%s:%d-%d does not match: %lf, %lf\n", 
								per_track_quantified_regions[track_i]->at(i_r)->chrom, per_track_quantified_regions[track_i]->at(i_r)->start, per_track_quantified_regions[track_i]->at(i_r)->end,
								cur_binned_profile[i_r+1], cur_total_sig);

						exit(0);
					}

					// Update the current index.
					i += l_cur_bin;
				} // i loop.
			} // i_r loop.
#endif // __MANUAL_CHECK__

			// Set the quantified tracks.
			cur_chr_tracks[track_i] = cur_binned_profile;
			track_i++;

			// Set the consistent length of tracks.
			if(min_n_bins < 0 ||
				min_n_bins > l_binned_track)
			{
				min_n_bins = l_binned_track;
			}

			delete [] cur_chr_signal_buff;
		} // i_dir loop.

		// Process the list of profiles in the config next.
		for(int i_p_dir = 0; i_p_dir < (int)bin_track_dirs->size(); i_p_dir++)
		{
			fprintf(stderr, "Loading signal profile from %s\n", bin_track_dirs->at(i_p_dir));

			// Load the profile.
			char cur_chr_profile_fp[1000];
			sprintf(cur_chr_profile_fp, "%s/%s.bin", bin_track_dirs->at(i_p_dir), chr_ids->at(i_chr));
			if(check_file(cur_chr_profile_fp))
			{

			}
			else
			{
				break;
			}

			int l_track = 0;
			double* cur_chr_profile = load_per_nucleotide_binary_profile(cur_chr_profile_fp, l_track);

			int l_binned_track = 0;
			double* cur_binned_profile = NULL;
			vector<t_annot_region*>* cur_track_quantified_regions = new vector<t_annot_region*>();
			cur_binned_profile = bin_and_normalize_target_profile(cur_chr_profile, 
																	cur_chr_target_regs, 
																	cur_track_quantified_regions,
																	flag_span_target_region_ends,
																	log_base,
																	l_track, 
																	l_bin, 
																	l_binned_track);

			// Set the quantified regions.
			per_track_quantified_regions[track_i] = cur_track_quantified_regions;
			sort(cur_track_quantified_regions->begin(), cur_track_quantified_regions->end(), sort_regions);
			fprintf(stderr, "%d quantified regions.\n", (int)cur_track_quantified_regions->size());

			// Set the quantified tracks.
			cur_chr_tracks[track_i] = cur_binned_profile;
			track_i++;

			// Set the consistent length of tracks.
			if(min_n_bins < 0 ||
				min_n_bins > l_binned_track)
			{
				min_n_bins = l_binned_track;
			}

			delete [] cur_chr_profile;
		} // i_p_dir

		// Go over all the region files and add them as tracks.
		for(int i_reg_fp = 0; i_reg_fp < (int)region_fps->size(); i_reg_fp++)
		{
			// Load the regions on this chromosome, and generate the binary profile.
			fprintf(stderr, "Loading signal profile from regions file %s\n", region_fps->at(i_reg_fp));

			// Load the profile.
			vector<t_annot_region*>* all_regions = load_BED(region_fps->at(i_reg_fp));
			vector<t_annot_region*>* cur_chr_regions = get_regions_per_chromosome(all_regions, chr_ids->at(i_chr));

			vector<t_annot_region*>* cur_chr_target_regions = get_regions_per_chromosome(target_regions, chr_ids->at(i_chr));

			// If there are no regions, set a 
			int l_track = 0;
			double* cur_chr_profile = NULL;
			if(cur_chr_regions->size() == 0)
			{
				sort(cur_chr_target_regions->begin(), cur_chr_target_regions->end(), sort_regions_per_ends);
				l_track = cur_chr_target_regions->back()->end + 1000;

				cur_chr_profile = new double[l_track];
				memset(cur_chr_profile, 0, sizeof(double) * l_track);
			}
			else
			{
				sort(cur_chr_regions->begin(), cur_chr_regions->end(), sort_regions_per_ends);
				sort(cur_chr_target_regions->begin(), cur_chr_target_regions->end(), sort_regions_per_ends);

				l_track = MAX(cur_chr_regions->back()->end, cur_chr_target_regions->back()->end) + 1000;

				cur_chr_profile = new double[l_track];
				memset(cur_chr_profile, 0, sizeof(double) * l_track);
				for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
				{
					for(int i = cur_chr_regions->at(i_reg)->start; i <= cur_chr_regions->at(i_reg)->end; i++)
					{
						cur_chr_profile[i] = 1;
					} // i loop.
				} // i_reg loop.
			}

			int l_binned_track = 0;
			vector<t_annot_region*>* cur_track_quantified_regions = new vector<t_annot_region*>();
			double* cur_binned_profile = bin_and_normalize_target_profile(cur_chr_profile, 
																			cur_chr_target_regs, 
																			cur_track_quantified_regions,
																			flag_span_target_region_ends,
																			log_base,
																			l_track, 
																			l_bin, 
																			l_binned_track);

			// Set the quantified regions.
			per_track_quantified_regions[track_i] = cur_track_quantified_regions;
			fprintf(stderr, "%d quantified regions.\n", (int)cur_track_quantified_regions->size());

			// Set the quantified tracks.
			cur_chr_tracks[track_i] = cur_binned_profile;
			track_i++;

			// Set the consistent length of tracks.
			if(min_n_bins < 0 ||
				min_n_bins > l_binned_track)
			{
				min_n_bins = l_binned_track;
			}

			delete [] cur_chr_profile;
		} // i_reg_fp loop.

		// Add the random signal track: This is random selection of signals from any one of the tracks. It should not have any predictability.
		if(use_random_track)
		{
			fprintf(stderr, "Adding a randomized control track.\n");

			// Generate the random track.
			double* cur_chr_random_track = new double[min_n_bins+2];
				
			// Copy the randomized track's values by selecting randomly from other tracks.
			for(int bin_i = 1; bin_i <= min_n_bins; bin_i++)
			{
				// Track based randomization.
				//int rand_track_i = ((int)(rng->random_double_ran3() * n_tracks)) % n_tracks;
				////fprintf(stderr, "Selected %d/%d track\n", rand_track_i, n_tracks);
				//cur_chr_random_track[bin_i] = cur_chr_tracks[rand_track_i][bin_i];

				// Total uniform binary randomization.
				cur_chr_random_track[bin_i] = floor(rng->random_double_ran3() * 2);
			} // bin_i loop.

			// Update the number of tracks and track ids, too. Note that the random track does not obey log values.
			cur_chr_tracks[track_i] = cur_chr_random_track;
			track_i++;
		} // random track check.

		// Make sure that all the tracks are set.
		bool all_tracks_set = true;
		for(int i_t = 0; i_t < n_tracks; i_t++)
		{
			if(per_track_quantified_regions[i_t] == NULL)
			{
				all_tracks_set = false;
				break;
			}
		} // i_t loop.

		if(all_tracks_set == false)
		{
			fprintf(stderr, "Skipping processing %s\n", chr_ids->at(i_chr));
			continue;
		}

		// Reprocess the quantified regions and merge their signals.
		vector<t_annot_region*>* cur_chr_quantified_regions = new vector<t_annot_region*>();
		
		// Following loop copies the per region signal levels to a single vector of regions. Can delete them after copy.
		for(int i_reg = 0; i_reg < min_n_bins; i_reg++)
		{
			t_annot_region* cur_quant_reg = duplicate_region(per_track_quantified_regions[0]->at(i_reg));
			cur_chr_quantified_regions->push_back(cur_quant_reg);

			double* cur_reg_signals = new double[n_tracks];
			for(int i_t = 0; i_t < n_tracks; i_t++)
			{
				if(per_track_quantified_regions[i_t]->at(i_reg)->start != per_track_quantified_regions[0]->at(i_reg)->start ||
					per_track_quantified_regions[i_t]->at(i_reg)->end != per_track_quantified_regions[0]->at(i_reg)->end)
				{
					fprintf(stderr, "The regions are not consistent: @ %d; %s:%d-%d, @ %d; %s:%d-%d\n", 
						0,
						per_track_quantified_regions[0]->at(i_reg)->chrom, per_track_quantified_regions[0]->at(i_reg)->start, per_track_quantified_regions[0]->at(i_reg)->end,
						i_t,
						per_track_quantified_regions[i_t]->at(i_reg)->chrom, per_track_quantified_regions[i_t]->at(i_reg)->start, per_track_quantified_regions[i_t]->at(i_reg)->end);
					exit(0);
				}

				cur_reg_signals[i_t] = per_track_quantified_regions[i_t]->at(i_reg)->dbl_score;
			} // i_t loop.

			// Set the data of the current region with joint signals.
			cur_quant_reg->data = cur_reg_signals;
		} // i_reg loop.

		// Free memory.
		for(int i_t = 0; i_t < n_tracks; i_t++)
		{
			delete_annot_regions(per_track_quantified_regions[i_t]);
		} // i_t loop.

		// Insert the quantified regions.
		quantified_regions->insert(quantified_regions->end(), cur_chr_quantified_regions->begin(), cur_chr_quantified_regions->end());

		delete cur_chr_quantified_regions;

		// Concatenate the currently generated tracks to the genomewide tracks.
		for(track_i = 0; track_i < n_tracks; track_i++)
		{
			int l_new_track = 0;
			double* new_gw_track = concat_profile2_2_profile1(gw_tracks[track_i], l_tracks[track_i], cur_chr_tracks[track_i], min_n_bins, l_new_track);

			// If there is already a track at this gw track, free its memory.
			if(gw_tracks[track_i] != NULL)
			{
				delete [] gw_tracks[track_i];
			}

			// Free the current chromosome track, which is not required any more.
			delete [] cur_chr_tracks[track_i];
			
			// Update the buffer and the length for the tracks.
			gw_tracks[track_i] = new_gw_track;
			l_tracks[track_i] = l_new_track;
			fprintf(stderr, "Concatted genomewide track length is %d (%d quant regions) @ %s (%d. track).\n", l_new_track, (int)quantified_regions->size(),
					chr_ids->at(i_chr), track_i);
		} // track_i loop.
	} // i_chr loop.

	// Assign the array index of the quantified regions.
	fprintf(stderr, "Setting region array index and comparing the per region signals with array signals.\n");
	for(int i_reg = 0; i_reg < (int)quantified_regions->size(); i_reg++)
	{
		quantified_regions->at(i_reg)->score = i_reg;

		double* cur_reg_signals = (double*)(quantified_regions->at(i_reg)->data);

		// Following is a sanity check for the computed signals so that we make sure the array really corresponds to the region signals.
		for(int i_t = 0; i_t < n_tracks; i_t++)
		{
			if(MAX(gw_tracks[i_t][i_reg+1], cur_reg_signals[i_t]) > 0.01 &&
				fabs(gw_tracks[i_t][i_reg+1] - cur_reg_signals[i_t]) / MAX(gw_tracks[i_t][i_reg+1], cur_reg_signals[i_t]) > 0.01)
			{
				fprintf(stderr, "Failed; %lf, %lf\n", gw_tracks[i_t][i_reg], cur_reg_signals[i_t]);
				exit(0);
			}
			else
			{
				//fprintf(stderr, "@ %d, check.               \r", i_reg);
			}
		} // i_t loop.
	} // i_reg loop.

	// Set the number and length of tracks.
	l_multitrack = l_tracks[0];

	return(gw_tracks);
}

// This is a function that generates multitrack signal matrix using a config file.
double** generate_multitrack_signal_profiles_per_config(char* config_fp, 
														bool flag_span_target_region_ends,
														vector<char*>* chr_ids_2_use, 
														int& n_gw_tracks, 
														int& l_gw_track, 
														vector<char*>* track_ids)
{
	// Get the joint entropy of the multiple tracks in the config file.
	t_config* config = new t_config(config_fp, "\t ");

	int use_per_target_signal = 0;
	if(config->get_int_val("use_per_target_signal", use_per_target_signal))
	{
		if(use_per_target_signal == 1)
		{
			use_per_target_signal = 0;
			fprintf(stderr, "Using per target signal.\n");
		}
	}

	// Set the bin size.
	int bin_size;
	if(!config->get_int_val("bin_size", bin_size))
	{
		if(use_per_target_signal == 0)
		{
			fprintf(stderr, "Need a bin_size or use_per_target_signal entry in the config file.\n");
			exit(0);
		}
		else
		{
			// No bin size but per target signal is requested.
		}
	}
	else
	{
		// There is a bin size, check if there is an option to use per target signal.
		if(use_per_target_signal == 1)
		{
			fprintf(stderr, "Using per target signal, although there is a binning option.\n");
		}
	}

	vector<char*>* chr_ids = NULL;
	if(chr_ids_2_use == NULL)
	{
		char chr_ids_fp[1000];
		if(!config->get_str_val("chr_ids_fp", chr_ids_fp))
		{
			fprintf(stderr, "Need the chr_ids_fp entry in the config file.\n");
			exit(0);
		}

		chr_ids = buffer_file(chr_ids_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load %s\n", chr_ids_fp);
			exit(0);
		}
	}
	else
	{
		chr_ids = chr_ids_2_use;
	}

	// This tunes the amount of quantization level in the signal.
	double log_base = exp(1.0);
	config->get_double_val("log_base", log_base);

	fprintf(stderr, "Set the log base to %lf for quantization.\n", log_base);

	char target_regions_fp[1000];
	if(!config->get_str_val("target_regions", target_regions_fp))
	{
		fprintf(stderr, "Need the target_regions entry in the config file.\n");
		exit(0);
	}

	// Load the target regions.
	vector<t_annot_region*>* target_regs = NULL;	
	if(check_file(target_regions_fp))
	{
		target_regs = load_BED(target_regions_fp);
		fprintf(stderr, "Loaded %d target regions.\n", (int)target_regs->size());
	}	
	else
	{
		fprintf(stderr, "No target region.\n");
	}

	// Track profiles from reads.
	vector<vector<char*>*>* reads_dirs = config->get_all_entries_per_id("reads_dir");

	// Track profiles from binarized profiles.
	vector<vector<char*>*>* track_bin_dirs = config->get_all_entries_per_id("profile_dir");

	// Track profiles from point events; i.e., TF peaks, SNPs, etc.
	vector<vector<char*>*>* region_fps = config->get_all_entries_per_id("region_fp");

	// Set the number of tracks: We need this to set the genomewide track buffer.
	int n_tracks = (int)reads_dirs->size() + (int)track_bin_dirs->size() + (int)region_fps->size();
	int use_random_track_flag = 0;
	if(config->get_int_val("use_random_track", use_random_track_flag))
	{
		if(use_random_track_flag > 0)
		{
			fprintf(stderr, "Adding a randomized control track.\n");
			n_tracks++;
		}
	} // random track usage flag.

	// Process all the tracks.
	fprintf(stderr, "Generating all the %d tracks.\n", n_tracks);

	// Do not do tag extension.
	int l_ext_tag = 0;

	double** gw_tracks = new double*[n_tracks+2];
	memset(gw_tracks, 0, sizeof(double*) * n_tracks);
	int* l_tracks = new int[n_tracks];

	// Process all the chromosomes.
	//FILE* f_per_chrom_stats = open_f("per_chrom_stats.txt", "w");
	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int min_n_bins = -1;

		fprintf(stderr, "Processing %s\n", chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_target_regs = NULL;
		if(target_regs != NULL)
		{
			cur_chr_target_regs = get_regions_per_chromosome(target_regs, chr_ids->at(i_chr));
			fprintf(stderr, "%d target regions.\n", (int)cur_chr_target_regs->size());
		}		

		double** cur_chr_tracks = new double*[n_tracks+2];
		//track_ids = new vector<char*>();
		track_ids->clear();
		int track_i = 0;

		// First, process the list of processed reads.
		enum{DOM_LINEAR_AVG, DOM_LINEAR_TOTAL, DOM_LOG};
		for(int i_dir = 0; i_dir < (int)reads_dirs->size(); i_dir++)
		{
			if(reads_dirs->at(i_dir)->size() != 3)
			{
				fprintf(stderr, "Could not read all the column for one of the entries.\n");
				exit(0);
			}

			int domain_type = -1;
			enum{DOM_LINEAR_AVG, DOM_LINEAR_TOTAL, DOM_LOG};
			if(t_string::compare_strings(reads_dirs->at(i_dir)->at(2), "linear_average"))
			{
				domain_type = DOM_LINEAR_AVG;
			}
			else if(t_string::compare_strings(reads_dirs->at(i_dir)->at(2), "linear_total"))
			{
				domain_type = DOM_LINEAR_TOTAL;
			}
			else if(t_string::compare_strings(reads_dirs->at(i_dir)->at(2), "log"))
			{
				domain_type = DOM_LOG;
			}
			else
			{
				fprintf(stderr, "Could not understand domain.\n");
				exit(0);
			}

			fprintf(stderr, "Loading signal profile from %s\n", reads_dirs->at(i_dir)->at(1));
			int l_buff = 300*1000*1000;
			double* cur_chr_signal_buff = new double[l_buff];
			int l_track = 0;

			track_ids->push_back(t_string::copy_me_str(reads_dirs->at(i_dir)->at(0)));

			//bool cur_track_linear = false;
			//bool cur_track_log = false;

			//if(t_string::compare_strings(reads_dirs->at(i_dir)->at(2), "linear"))
			//{
			//	cur_track_linear = true;
			//}			
			//else if(t_string::compare_strings(reads_dirs->at(i_dir)->at(2), "log"))
			//{
			//	cur_track_log = true;
			//}
			//else
			//{
			//	fprintf(stderr, "Could not determine the domain of signal for %s\n", reads_dirs->at(i_dir)->at(2));
			//	exit(0);
			//}

			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", reads_dirs->at(i_dir)->at(1), chr_ids->at(i_chr));
			buffer_per_nucleotide_profile_no_buffer(cur_chr_reads_fp, l_ext_tag, cur_chr_signal_buff, NULL, NULL, l_buff, l_track);

			// Bin the profile, delete the buffer.
			int n_bins = 0;
			double* cur_binned_profile = bin_profile(cur_chr_signal_buff, l_track,
													cur_chr_target_regs,
													NULL,
													flag_span_target_region_ends,
													bin_size,
													n_bins);
			
			if(domain_type == DOM_LINEAR_TOTAL)
			{
				// No need to do anything.
				fprintf(stderr, "Using linear total signal.\n");
			}
			else if(domain_type == DOM_LOG)
			{
				fprintf(stderr, "Converting binned signal to logarithm.\n");
				for(int bin_i = 1; bin_i <= n_bins; bin_i++)
				{
					if(cur_binned_profile[bin_i] < 0)
					{
						fprintf(stderr, "%d. bin value is negative.\n", bin_i);
						exit(0);
					}

					// Must be careful with the data transformations.
					//cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i] + 1) / xlog(log_base);
					cur_binned_profile[bin_i] = floor(xlog(cur_binned_profile[bin_i]+1) / xlog(log_base));
				} // bin_i loop.
			}
			else if(domain_type == DOM_LINEAR_AVG)
			{
				fprintf(stderr, "Converting binned signal to average.\n");
				for(int bin_i = 1; bin_i <= n_bins; bin_i++)
				{
					if(cur_binned_profile[bin_i] < 0)
					{
						fprintf(stderr, "%d. bin value is negative.\n", bin_i);
						exit(0);
					}

					// Must be careful with the data transformations.
					//cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i] + 1) / xlog(log_base);
					cur_binned_profile[bin_i] = floor(cur_binned_profile[bin_i] / bin_size);
				} // bin_i loop.
			}

			// Add to the list of tracks.
			cur_chr_tracks[track_i] = cur_binned_profile;
			track_i++;

			// Set the consistent length of tracks.
			if(min_n_bins < 0 ||
				min_n_bins > n_bins)
			{
				min_n_bins = n_bins;
			}

			delete [] cur_chr_signal_buff;
		} // i_dir loop.

		// Process the list of profiles in the config next.
		for(int i_p_dir = 0; i_p_dir < (int)track_bin_dirs->size(); i_p_dir++)
		{
			if(track_bin_dirs->at(i_p_dir)->size() != 3)
			{
				fprintf(stderr, "Could not read all the column for one of the entries.\n");
				exit(0);
			}

			fprintf(stderr, "Loading signal profile from %s\n", track_bin_dirs->at(i_p_dir)->at(1));

			int domain_type = -1;
			enum{DOM_LINEAR_AVG, DOM_LINEAR_TOTAL, DOM_LOG};
			if(t_string::compare_strings(track_bin_dirs->at(i_p_dir)->at(2), "linear_average"))
			{
				domain_type = DOM_LINEAR_AVG;
			}
			else if(t_string::compare_strings(track_bin_dirs->at(i_p_dir)->at(2), "linear_total"))
			{
				domain_type = DOM_LINEAR_TOTAL;
			}
			else if(t_string::compare_strings(track_bin_dirs->at(i_p_dir)->at(2), "log"))
			{
				domain_type = DOM_LOG;
			}
			else
			{
				fprintf(stderr, "Could not understand domain.\n");
				exit(0);
			}

			// Add the current profile's id to the list of ids.
			track_ids->push_back(t_string::copy_me_str(track_bin_dirs->at(i_p_dir)->at(0)));

			//bool cur_track_linear = false;
			//bool cur_track_log = false;
			//if(t_string::compare_strings(track_bin_dirs->at(i_p_dir)->at(2), "linear"))
			//{
			//	cur_track_linear = true;
			//}			
			//else if(t_string::compare_strings(track_bin_dirs->at(i_p_dir)->at(2), "log"))
			//{
			//	cur_track_log = true;
			//}
			//else
			//{
			//	fprintf(stderr, "Could not determine the domain of signal for %s\n", track_bin_dirs->at(i_p_dir)->at(0));
			//	exit(0);
			//}

			// Load the profile.
			char cur_chr_profile_fp[1000];
			sprintf(cur_chr_profile_fp, "%s/%s.bin", track_bin_dirs->at(i_p_dir)->at(1), chr_ids->at(i_chr));
			int l_profile = 0;
			double* cur_chr_profile = load_per_nucleotide_binary_profile(cur_chr_profile_fp, l_profile);

			int n_bins = 0;
			double* cur_binned_profile = bin_profile(cur_chr_profile, l_profile,
													cur_chr_target_regs,
													NULL,
													flag_span_target_region_ends,
													bin_size,
													n_bins);

			if(domain_type == DOM_LINEAR_TOTAL)
			{
				// No need to do anything.
				fprintf(stderr, "Using linear total signal.\n");
			}
			else if(domain_type == DOM_LOG)
			{
				fprintf(stderr, "Converting binned signal to logarithm.\n");
				for(int bin_i = 1; bin_i <= n_bins; bin_i++)
				{
					if(cur_binned_profile[bin_i] < 0)
					{
						fprintf(stderr, "%d. bin value is negative.\n", bin_i);
						exit(0);
					}

					// Must be careful with the data transformations.
					//cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i] + 1) / xlog(log_base);
					cur_binned_profile[bin_i] = floor(xlog(cur_binned_profile[bin_i]+1) / xlog(log_base));
				} // bin_i loop.
			}
			else if(domain_type == DOM_LINEAR_AVG)
			{
				fprintf(stderr, "Converting binned signal to average.\n");
				for(int bin_i = 1; bin_i <= n_bins; bin_i++)
				{
					if(cur_binned_profile[bin_i] < 0)
					{
						fprintf(stderr, "%d. bin value is negative.\n", bin_i);
						exit(0);
					}

					// Must be careful with the data transformations.
					//cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i] + 1) / xlog(log_base);
					cur_binned_profile[bin_i] = floor(cur_binned_profile[bin_i] / bin_size);
				} // bin_i loop.
			}

			// Add to the list of tracks.
			cur_chr_tracks[track_i] = cur_binned_profile;
			track_i++;

			if(min_n_bins < 0 ||
				min_n_bins > n_bins)
			{
				min_n_bins = n_bins;
			}

			delete [] cur_chr_profile;
		} // i_p_dir

		// Go over all the region files and add them as tracks.
		for(int i_reg_fp = 0; i_reg_fp < (int)region_fps->size(); i_reg_fp++)
		{
			if(region_fps->at(i_reg_fp)->size() != 3)
			{
				fprintf(stderr, "Could not read all the column for one of the entries.\n");
				exit(0);
			}

			// Load the regions on this chromosome, and generate the binary profile.
			fprintf(stderr, "Loading signal profile from regions file %s\n", region_fps->at(i_reg_fp)->at(1));

			int domain_type = -1;
			if(t_string::compare_strings(region_fps->at(i_reg_fp)->at(2), "linear_average"))
			{
				domain_type = DOM_LINEAR_AVG;
			}
			else if(t_string::compare_strings(region_fps->at(i_reg_fp)->at(2), "linear_total"))
			{
				domain_type = DOM_LINEAR_TOTAL;
			}
			else if(t_string::compare_strings(region_fps->at(i_reg_fp)->at(2), "log"))
			{
				domain_type = DOM_LOG;
			}
			else
			{
				fprintf(stderr, "Could not understand domain.\n");
				exit(0);
			}

			// Add the current profile's id to the list of ids.
			track_ids->push_back(t_string::copy_me_str(region_fps->at(i_reg_fp)->at(0)));

			//bool cur_track_linear = false;
			//bool cur_track_log = false;
			//if(t_string::compare_strings(region_fps->at(i_reg_fp)->at(2), "linear"))
			//{
			//	cur_track_linear = true;
			//}			
			//else if(t_string::compare_strings(region_fps->at(i_reg_fp)->at(2), "log"))
			//{
			//	cur_track_log = true;
			//}
			//else
			//{
			//	fprintf(stderr, "Could not determine the domain of signal for %s\n", region_fps->at(i_reg_fp)->at(0));
			//	exit(0);
			//}

			// Load the profile.
			vector<t_annot_region*>* all_regions = load_BED(region_fps->at(i_reg_fp)->at(1));
			vector<t_annot_region*>* cur_chr_regions = get_regions_per_chromosome(all_regions, chr_ids->at(i_chr));
			sort(cur_chr_regions->begin(), cur_chr_regions->end(), sort_regions_per_ends);
			int l_profile = cur_chr_regions->back()->end + 1000;
			double* cur_chr_profile = new double[l_profile];
			memset(cur_chr_profile, 0, sizeof(double) * l_profile);
			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
			{
				for(int i = cur_chr_regions->at(i_reg)->start; i <= cur_chr_regions->at(i_reg)->end; i++)
				{
					cur_chr_profile[i] = 1;
				} // i loop.
			} // i_reg loop.

			int n_bins = 0;
			double* cur_binned_profile = bin_profile(cur_chr_profile, l_profile,
													cur_chr_target_regs,
													NULL,
													flag_span_target_region_ends,
													(int)bin_size,
													n_bins);

			if(domain_type == DOM_LINEAR_TOTAL)
			{
				// No need to do anything.
				fprintf(stderr, "Using linear total signal.\n");
			}
			else if(domain_type == DOM_LOG)
			{
				fprintf(stderr, "Converting binned signal to logarithm.\n");
				for(int bin_i = 1; bin_i <= n_bins; bin_i++)
				{
					if(cur_binned_profile[bin_i] < 0)
					{
						fprintf(stderr, "%d. bin value is negative.\n", bin_i);
						exit(0);
					}

					// Must be careful with the data transformations.
					//cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i] + 1) / xlog(log_base);
					cur_binned_profile[bin_i] = floor(xlog(cur_binned_profile[bin_i]+1) / xlog(log_base));
				} // bin_i loop.
			}
			else if(domain_type == DOM_LINEAR_AVG)
			{
				fprintf(stderr, "Converting binned signal to average.\n");
				for(int bin_i = 1; bin_i <= n_bins; bin_i++)
				{
					if(cur_binned_profile[bin_i] < 0)
					{
						fprintf(stderr, "%d. bin value is negative.\n", bin_i);
						exit(0);
					}

					// Must be careful with the data transformations.
					//cur_binned_profile[bin_i] = xlog(cur_binned_profile[bin_i] + 1) / xlog(log_base);
					cur_binned_profile[bin_i] = floor(cur_binned_profile[bin_i] / bin_size);
				} // bin_i loop.
			}

			// Add to the list of tracks.
			cur_chr_tracks[track_i] = cur_binned_profile;
			track_i++;

			if(min_n_bins < 0 ||
				min_n_bins > n_bins)
			{
				min_n_bins = n_bins;
			}

			delete [] cur_chr_profile;
		} // i_reg_fp loop.

		// Add the random signal track: This is random selection of signals from any one of the tracks. It should not have any predictability.
		if(use_random_track_flag == 1)
		{
			fprintf(stderr, "Adding a randomized control track.\n");

			// Generate the random track.
			double* cur_chr_random_track = new double[min_n_bins+2];
				
			// Copy the randomized track's values by selecting randomly from other tracks.
			for(int bin_i = 1; bin_i <= min_n_bins; bin_i++)
			{
				// Track based randomization.
				//int rand_track_i = ((int)(rng->random_double_ran3() * n_tracks)) % n_tracks;
				////fprintf(stderr, "Selected %d/%d track\n", rand_track_i, n_tracks);
				//cur_chr_random_track[bin_i] = cur_chr_tracks[rand_track_i][bin_i];

				// Total uniform binary randomization.
				cur_chr_random_track[bin_i] = floor(rng->random_double_ran3() * 2);
			} // bin_i loop.

			// Update the number of tracks and track ids, too.
			track_ids->push_back(t_string::copy_me_str("Random_Control"));
			cur_chr_tracks[track_i] = cur_chr_random_track;
			track_i++;
		}

		// All the tracks must be processed at this point.
		// Sanity check on the number of tracks read and enumerated.
		if(track_i != n_tracks)
		{
			fprintf(stderr, "Could not match the track number count from enumeration: %d, %d\n", track_i, n_tracks);
			exit(0);
		}

		// Concatenate the currently generated tracks to the genomewide tracks.
		for(track_i = 0; track_i < n_tracks; track_i++)
		{
			int l_new_track = 0;
			double* new_gw_track = concat_profile2_2_profile1(gw_tracks[track_i], l_tracks[track_i], cur_chr_tracks[track_i], min_n_bins, l_new_track);

			// If there is already a track at this gw track, free its memory.
			if(gw_tracks[track_i] != NULL)
			{
				delete [] gw_tracks[track_i];
			}
			
			// Update the buffer and the length for the tracks.
			gw_tracks[track_i] = new_gw_track;
			l_tracks[track_i] = l_new_track;
			fprintf(stderr, "Concatted genomewide track length is %d @ %s (%d. track).\n", l_new_track, track_ids->at(track_i), track_i);
		} // track_i loop.
	} // i_chr loop.

	if((int)track_ids->size() != n_tracks)
	{
		fprintf(stderr, "The number of tracks is not consistent: %d ids, Read %d tracks.\n", (int)track_ids->size(), n_tracks);
		exit(0);
	}

	// Set the number and length of tracks.
	n_gw_tracks = n_tracks;
	l_gw_track = l_tracks[0];

	return(gw_tracks);
} // generate_multitrack_signal_profiles_per_config option

//#define MAX(x,y) ((x)>(y)?(x):(y))
//#define MIN(x,y) ((x)<(y)?(x):(y))
double* get_random_circular_shifted_profile(double* profile_data, t_rng* rng, int l_profile)
{
	double* shifted_track_data = new double[l_profile+2];

	int l_shift = (int)(l_profile * rng->random_double_ran3());
	while(l_shift < 0)
	{
		l_shift += l_profile;
	}

	for(int i = 1; i <= l_profile; i++)
	{
		int i_shifted = (i + l_shift);
		while(i_shifted > l_profile)
		{
			i_shifted = (i_shifted - l_profile);
		}

		shifted_track_data[i_shifted] = profile_data[i];
	} // i loop.

	return(shifted_track_data);
}

void dump_per_nucleotide_uchar_binary_profile(unsigned char* signal_profile, int l_profile, char* op_fp)
{
	// Dump the per nucleotide signal profile.
	FILE* f_op = open_f(op_fp, "wb");

	// Write data length to first couple bytes.
	fwrite(&(l_profile), sizeof(int), 1, f_op);

	// Dump the data: Dump 0 based data.
	fwrite(&(signal_profile[1]), sizeof(unsigned char), l_profile+1, f_op);

	fclose(f_op);
}

unsigned char* load_per_nucleotide_binary_uchar_profile(char* binary_per_nucleotide_profile_fp, int& l_profile)
{
	FILE* f_prof = open_f(binary_per_nucleotide_profile_fp, "rb");

	// Read the profile length.
	int l_data = 0;
	fread(&l_data, sizeof(int), 1, f_prof);
	l_profile = l_data;

	// Read the data.
	unsigned char* signal_profile_buffer = new unsigned char[l_profile+2];

if(__DUMP_SIGNAL_TRACK_MSGS__)
{
	fprintf(stderr, "Loading %d data values.\n", l_profile);
}

	// Following is to use the codebase indexing: 1 based indexing.
	fread(&(signal_profile_buffer[1]), sizeof(char), l_profile+1, f_prof);
	
	fclose(f_prof);

	return(signal_profile_buffer);
}

double* get_block_permute_profile(double* profile_data, t_rng* rng, int l_profile, int& l_permuted, int l_block)
{
	//int l_cur_profile;
	//double* cur_profile_data = load_per_nucleotide_binary_profile(profile_fp, l_cur_profile);

	//  Count the blocks.
	int n_blocks = (int)(l_profile / l_block);
	fprintf(stderr, "There are %d blocks in total.\n", n_blocks);

	// Permute the block ordering.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	//vector<int>* permuted_block_i = rng->permute_indices(n_blocks, n_blocks);
	vector<int>* permuted_block_i = rng->fast_permute_indices(0, n_blocks);

	double* permuted_profile_data = new double[l_profile + 2];
	l_permuted = 1;
	for(int block_i = 0; block_i < (int)permuted_block_i->size(); block_i++)
	{
		int cur_block_start = (permuted_block_i->at(block_i) * l_block)+1;
		int cur_block_end = (l_profile > cur_block_start+l_block)?(cur_block_start+l_block):(l_profile);

		for(int i = cur_block_start; i < cur_block_end; i++)
		{
			permuted_profile_data[l_permuted] = profile_data[i];
			l_permuted++;
		} // i loop.
	} // block_i loop.

	// Dump the profile.
	return(permuted_profile_data);
}

void get_profile_extrema(double* profile, int l_profile, double& prof_min, double& prof_max)
{
	prof_min = 1000*1000;
	prof_max = -1000*1000;
	for(int i = 1; i <= l_profile; i++)
	{
		if(profile[i] > prof_max)
		{
			prof_max = profile[i];
		}

		if(profile[i] < prof_min)
		{
			prof_min = profile[i];
		}
	} // i loop.
}

double* get_zero_indexed_per_one_indexed_data(double* one_indexed_data, int l_profile)
{
	double* zero_indexed_data = new double[l_profile];
	for(int i = 1; i <= l_profile; i++)
	{
		zero_indexed_data[i-1] = one_indexed_data[i];
	} // i loop.

	return(zero_indexed_data);
}

double* get_one_indexed_per_zero_indexed_data(double* zero_indexed_data, int l_profile)
{
	double* one_indexed_data = new double[l_profile+2];
	for(int i = 0; i < l_profile; i++)
	{
		one_indexed_data[i+1] = zero_indexed_data[i];
	} // i loop.

	return(one_indexed_data);
}

double* quantize_per_nucleotide_profiles(double* profile, int l_profile, vector<double>* thresholds, vector<double>* quantized_vals)
{
	//int l_profile = 0;
	//double* profile = load_per_nucleotide_binary_profile(prof_fp, l_profile);

	//vector<char*>* thresh_val_lines = buffer_file(thresh_val_fp);

	//vector<double>* quantized_vals = new vector<double>();
	//vector<double>* thresholds = new vector<double>();
	//for(int i_l = 0; i_l < thresh_val_lines->size(); i_l++)
	//{
	//	double cur_quantized_val;
	//	double cur_thresh;
	//	if(sscanf(thresh_val_lines->at(i_l), "%lf %lf", &cur_thresh, &cur_quantized_val) != 2)
	//	{
	//		fprintf(stderr, "Could not parse the line: %s\n", thresh_val_lines->at(i_l));
	//		exit(0);
	//	}

	//	quantized_vals->push_back(cur_quantized_val);
	//	thresholds->push_back(cur_thresh);
	//	fprintf(stderr, "%lf -> %lf\n", cur_thresh, cur_quantized_val);
	//} // i_l loop.

	if(thresholds->size() != quantized_vals->size())
	{
		fprintf(stderr, "The size of thresholds is not the same as size of quantized values.\n");
		exit(0);
	}

	double* quantized_profile = new double[l_profile+2];
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_prof_val = profile[i];

		// Find the largest threshold that this value is greater than or equal to.
		bool quantized_current_val = false;
		for(int i_th = (int)thresholds->size()-1;
			!quantized_current_val && i_th >= 0; 
			i_th--)
		{
			if(cur_prof_val >= thresholds->at(i_th))
			{
				quantized_current_val = true;
				quantized_profile[i] = quantized_vals->at(i_th);
			}
		} // i_th loop.

		if(!quantized_current_val)
		{
			fprintf(stderr, "Could not quantize the current value: %lf\n", profile[i]);
			exit(0);
		}
	} // i loop.

	return(quantized_profile);
}

double** joint_prune_multiple_profiles_per_thresholds(double** profiles, int n_profiles, int l_profile, double* threshold_per_profile, int& l_pruned_profile)
{
	// Go over all the profiles, find the regions.	
	double** pruned_profiles = new double*[n_profiles];
	for(int i_p = 0; i_p < n_profiles; i_p++)
	{
		pruned_profiles[i_p] = new double[l_profile+2];
	} // i_p loop.

	// Prune the columns.
	l_pruned_profile = 1;
	for(int i_sig = 1; i_sig <= l_profile; i_sig++)
	{
		bool include_cur_val = false;
		for(int i_p = 0; i_p < n_profiles; i_p++)
		{
			if(threshold_per_profile[i_p] <= -1000)
			{
				//fprintf(stderr, "Not using %d. profile for joint pruning.\n");
			}
			else if(profiles[i_p][i_sig] > threshold_per_profile[i_p])
			{
				include_cur_val = true;
			}
		} // i_p loop.

		if(include_cur_val)
		{
			// Add this columns to all the pruned profiles.
			for(int i_p = 0; i_p < n_profiles; i_p++)
			{
				// Copy the value.
				pruned_profiles[i_p][l_pruned_profile] = profiles[i_p][i_sig];
			} // i_p loop.	

			l_pruned_profile++;
		}
	} // i_sig loop.

	// Pruned profiles are indexed 1 based.
	l_pruned_profile--;

	return(pruned_profiles);
}

double* copy_profile(double* signal_profile, int l_profile)
{
	double* copy_prof = new double[l_profile+2];

	for(int i_sig = 1; i_sig <= l_profile; i_sig++)
	{
		copy_prof[i_sig] = signal_profile[i_sig];
	} // i_sig loop.

	return(copy_prof);
}

double* extract_one_indexed_profile_per_profile_stranded(double* signal_profile_buffer, int l_profile, int start, int end, char strand, int& l_extracted_profile)
{
	double* extracted_prof = new double[end - start + 3];
	for (int i_sig = start; i_sig <= end; i_sig++)
	{
		extracted_prof[i_sig - start] = 0;
	} // i_sig loop.

	// Count the extracted profile.
	l_extracted_profile = 1;

	if (strand == '+' ||
		strand == 'F')
	{
		for (int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
		{
			extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
			l_extracted_profile++;
		} // i_sig loop.
	}
	else
	{
		//for (int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
		for (int i_sig = MIN(l_profile, end); i_sig >= start; i_sig--)
		{
			extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
			l_extracted_profile++;
		} // i_sig loop.
	}
	// Profiles are one based.
	l_extracted_profile--;

	return(extracted_prof);
}


double* extract_one_indexed_profile_per_profile(double* signal_profile_buffer, int l_profile, int start, int end, int& l_extracted_profile)
{
	double* extracted_prof = new double[end - start + 3];	
	for(int i_sig = start; i_sig <= end; i_sig++)
	{
		extracted_prof[i_sig - start] = 0;
	} // i_sig loop.
	
	// Count the extracted profile.
	l_extracted_profile = 1;
	for(int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
	{
		extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
		l_extracted_profile++;
	} // i_sig loop.

	// Profiles are one based.
	l_extracted_profile--;

	return(extracted_prof);
}

// The signal profile is 1 based, consistent with the codebase indexing.
void dump_bedGraph_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* op_fp)
{
	FILE* f_op = NULL;
	//bool concatting = false;
	if(check_file(op_fp))
	{
		fprintf(stderr, "%s exists, concatting.\n", op_fp);
		f_op = open_f(op_fp, "a");
		//concatting = true;
	}
	else
	{
		f_op = open_f(op_fp, "w");
	}

	// Get the bedgraph for the current profile.
	int i_nuc = 1; 
	double cur_height = signal_profile_buffer[i_nuc];
	int cur_block_start_i = i_nuc;
	i_nuc++;
	while(1)
	{
		// Find the point where the height changes: The end goes till it's equal to the profile length since the profile is 1-based.
		while(i_nuc <= l_profile)
		{
			// Wait till there is a change in the height, which marks the start of a new block.
			if(cur_height != signal_profile_buffer[i_nuc])
			{
				break;
			}

			i_nuc++;
		} // i_nuc loop.

		// At this point, either this is the end of the profile, or there was a change in the height, either way, this was the end of the current block. Definitely dump it.
		if(cur_height != signal_profile_buffer[i_nuc])
		{
			// Dump the current block.
			fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom, 
				translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(i_nuc-1, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
				cur_height);

			// Update the new height and new start.
			cur_height = signal_profile_buffer[i_nuc];
			
			// Current position starts the next block.
			cur_block_start_i = i_nuc; 
		}

		// If the above block end was the end of the whole profile, we are done, otherwise continue to the next block.
		if(i_nuc > l_profile)
		{
			break;
		}

		//i_nuc++;
	} // i_nuc loop.

	//fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom, 
	//			translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//			translate_coord(i_nuc-1, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
	//			cur_height);

	close_f(f_op, op_fp);
}

// The signal profile is 1 based, consistent with the codebase indexing: This is the latest version of bBGR dumping.
void dump_bBGR_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* bbgr_op_fp)
{
	FILE* f_bbgr_op = open_f(bbgr_op_fp, "wb");

	// Get the bedgraph for the current profile.
	int i_nuc = 1;
	double cur_height = signal_profile_buffer[i_nuc];
	int cur_block_start_i = i_nuc;
	i_nuc++;
	while (1)
	{
		// Find the point where the height changes: The end goes till it's equal to the profile length since the profile is 1-based.
		while (i_nuc <= l_profile)
		{
			// Wait till there is a change in the height, which marks the start of a new block.
			if (cur_height != signal_profile_buffer[i_nuc])
			{
				break;
			}

			i_nuc++;
		} // i_nuc loop.

		  // At this point, either this is the end of the profile, or there was a change in the height, either way, this was the end of the current block. Definitely dump it.
		if (cur_height != signal_profile_buffer[i_nuc])
		{
			int start = translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			int end = translate_coord(i_nuc - 1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			//// Dump the current block.
			//fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom,
			//	translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			//	translate_coord(i_nuc - 1, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			//	cur_height);

			// Dump the current block.
			fwrite(&start, sizeof(int), 1, f_bbgr_op);
			fwrite(&end, sizeof(int), 1, f_bbgr_op);
			fwrite(&cur_height, sizeof(double), 1, f_bbgr_op);

			// Update the new height and new start.
			cur_height = signal_profile_buffer[i_nuc];

			// Current position starts the next block.
			cur_block_start_i = i_nuc;
		}

		// If the above block end was the end of the whole profile, we are done, otherwise continue to the next block.
		if (i_nuc > l_profile)
		{
			break;
		}

		//i_nuc++;
	} // i_nuc loop.

	fclose(f_bbgr_op);

	char comp_bbgr_op_fp[1000];
	sprintf(comp_bbgr_op_fp, "%s.gz", bbgr_op_fp);
	compressFile(bbgr_op_fp, comp_bbgr_op_fp);
	delete_file(bbgr_op_fp);
}

void dump_bedGraph_per_bBGR_v1(char* bbgr_fp, char* chrom, char* op_fp)
{
	int l_profile = 0;
	double* profile = load_per_nucleotide_bBGR_v1_track(bbgr_fp, l_profile);

	dump_bedGraph_per_per_nucleotide_binary_profile(profile, l_profile, chrom, op_fp);
}

void dump_bedGraph_per_bBGR(char* bbgr_fp, char* chrom, char* op_fp)
{
	int l_profile = 0;
	double* profile = load_per_nucleotide_bBGR_track(bbgr_fp, l_profile);

	dump_bedGraph_per_per_nucleotide_binary_profile(profile, l_profile, chrom, op_fp);
}

void dump_bBGR_per_per_bedGraph(char* bgr_fp, char* bbgr_op_fp)
{
	FILE* f_bgr = open_f(bgr_fp, "r");

	FILE* f_bbgr_op = open_f(bbgr_op_fp, "wb");
	while (1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if (cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if (sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		//fwrite(cur_chrom, 1, 20, f_bbgr_op);
		fwrite(&start, sizeof(int), 1, f_bbgr_op);
		fwrite(&end, sizeof(int), 1, f_bbgr_op);
		fwrite(&cur_sig, sizeof(double), 1, f_bbgr_op);
	} // bgr reading loop.

	close_f(f_bbgr_op, bbgr_op_fp);
	close_f(f_bgr, bgr_fp);
}

double* load_per_nucleotide_bBGR_track(char* bgr_fp, int& l_profile)
{
	FILE* f_bbgr = open_f(bgr_fp, "rb");

	// Initialize the signal profile.
	int l_max = 300 * 1000 * 1000;
	double* signal_profile = new double[l_max + 1];
	memset(signal_profile, 0, sizeof(double) * l_max);

	// Go over all the input file and process all the lines.
	while (1)
	{
		//int i_chr = 0;
		int start = 0;
		int end = 0;
		double cur_sig = 0;
		if (fread(&start, sizeof(int), 1, f_bbgr) == 0)
		{
			break;
		}

		fread(&end, sizeof(int), 1, f_bbgr);
		fread(&cur_sig, sizeof(double), 1, f_bbgr);

		l_profile = end + 10;

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for (int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if (i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d\n", i_nuc);
				exit(0);
			}
			signal_profile[i_nuc] += cur_sig;
		} // i_nuc loop.
	} // file reading loop.	

	// Close the bedgraph file.
	close_f(f_bbgr, bgr_fp);

	return(signal_profile);
}

double* load_per_nucleotide_bBGR_v1_track(char* bgr_fp, int& l_profile)
{
	FILE* f_bbgr = open_f(bgr_fp, "rb");

	// Initialize the signal profile.
	int l_max = 300 * 1000 * 1000;
	double* signal_profile = new double[l_max + 1];
	memset(signal_profile, 0, sizeof(double) * l_max);

	// Go over all the input file and process all the lines.
	while (1)
	{
		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if (fread(cur_chrom, 1, 20, f_bbgr) == 0)
		{
			break;
		}

		fread(&start, sizeof(int), 1, f_bbgr);
		fread(&end, sizeof(int), 1, f_bbgr);
		fread(&cur_sig, sizeof(double), 1, f_bbgr);

		l_profile = end + 10;

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for (int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if (i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d\n", i_nuc);
				exit(0);
			}
			signal_profile[i_nuc] = cur_sig;
		} // i_nuc loop.
	} // file reading loop.	

	  // Close the bedgraph file.
	close_f(f_bbgr, bgr_fp);

	return(signal_profile);
}

double* load_per_nucleotide_BGR_track(const char* bgr_fp, int& l_profile)
{
	FILE* f_bgr = open_f(bgr_fp, "r");

	// Initialize the signal profile.
	int l_max = 300 * 1000 * 1000;
	double* signal_profile = new double[l_max + 1];
	memset(signal_profile, 0, sizeof(double) * l_max);

	// Go over all the input file and process all the lines.
	while (1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if (cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if (sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		l_profile = end + 10;

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for (int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if (i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d (%s)\n", i_nuc, cur_bgr_line);
				exit(0);
			}

			// Note that this pools if there are overlapping positions in the bedgraph file.
			signal_profile[i_nuc] += cur_sig;
		} // i_nuc loop.

		delete[] cur_bgr_line;
	} // file reading loop.	

	// Close the bedgraph file.
	close_f(f_bgr, bgr_fp);

	return(signal_profile);
}

void dump_per_nucleotide_binary_profile_per_bedgraph(const char* bgr_fp, bool dump_binary, const char* op_fp)
{
	FILE* f_bgr = open_f(bgr_fp, "r");

	// Initialize the signal profile.
	int l_max = 300*1000*1000;
	double* signal_profile = new double[l_max+1];
	for(int i_nuc = 0; i_nuc <= l_max; i_nuc++)
	{
		signal_profile[i_nuc] = 0.0;
	} // i_nuc loop.

	// Go over all the input file and process all the lines.
	fprintf(stderr, "Dumping the profile to %s.\n", op_fp);
	while(1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if(cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if(sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for(int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if(i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d (%s)\n", i_nuc, cur_bgr_line);
				exit(0);
			}
			signal_profile[i_nuc] = cur_sig;
		} // i_nuc loop.

		delete [] cur_bgr_line;
	} // file reading loop.	

	// Close the bedgraph file.
	close_f(f_bgr, bgr_fp);

	// Get the end of the signal.
	int l_data = l_max;
	while(signal_profile[l_data] == 0.0)
	{
		l_data--;
	} // i_nuc loop.

	fprintf(stderr, "Signal length is %d, dumping the per nucleotide profile.\n", l_data);

	if(dump_binary)
	{
		dump_per_nucleotide_binary_profile(signal_profile, l_data, op_fp);
		return;
	}

	// Dump the per nucleotide signal profile.
	FILE* f_op = NULL;
	
	fprintf(stderr, "Saving Plain + Compressing.\n");
	f_op = open_f(op_fp, "w");
	
	fprintf(f_op, "%d\n",  l_data);
	
	for(int i_nuc = 1; i_nuc <= l_data; i_nuc++)
	{
		if(dump_binary)
		{
			fwrite(&(signal_profile[i_nuc]), sizeof(double), 1, f_op);
		}
		else
		{
			fprintf(f_op, "%lf ",  signal_profile[i_nuc]);
		}
	} // i_nuc loop.

	close_f(f_op, op_fp);
}

void dump_per_nucleotide_binary_profile(double* signal_profile, int l_profile, const char* op_fp)
{
	// Dump the per nucleotide signal profile.
	FILE* f_op = open_f(op_fp, "wb");

	// Write data length to first couple bytes.
	fwrite(&(l_profile), sizeof(int), 1, f_op);

	// Dump the data: Dump 0 based data.
	fwrite(&(signal_profile[1]), sizeof(double), l_profile, f_op);

	close_f(f_op, op_fp);
}

double* load_per_nucleotide_binary_profile(char* binary_per_nucleotide_profile_fp, int& l_profile)
{
	FILE* f_prof = open_f(binary_per_nucleotide_profile_fp, "rb");

	if (f_prof == NULL)
	{
		fprintf(stderr, "Could not open %s, unexpected extension other than bin/bin.gz?\n", binary_per_nucleotide_profile_fp);
		exit(0);
	}

	// Read the profile length.
	int l_data = 0;
	fread(&l_data, sizeof(int), 1, f_prof);
	l_profile = l_data;

	// Read the data.
	double* signal_profile_buffer = new double[l_profile+2];
	fprintf(stderr, "Loading %d data values.\n", l_profile);

	// Following is to use the codebase indexing: 1 based indexing.
	fread(&(signal_profile_buffer[1]), sizeof(double), l_profile, f_prof);
	
	close_f(f_prof, binary_per_nucleotide_profile_fp);

	return(signal_profile_buffer);
}

void exclude_regions_from_signal_profiles(double* signal_profile, int l_profile, vector<t_annot_region*>* regions_2_exclude, double* pruned_signal_profile, int& l_pruned_profile)
{
	sort(regions_2_exclude->begin(), regions_2_exclude->end(), sort_regions);

	bool* exclusion_profile = new bool[l_profile+2];
	for(int i = 1; i <= l_profile; i++)
	{
		exclusion_profile[i] = false;
	} // i loop.

	for(int i_reg = 0; i_reg < (int)regions_2_exclude->size(); i_reg++)
	{
		for(int i = regions_2_exclude->at(i_reg)->start; i <= regions_2_exclude->at(i_reg)->end; i++)
		{
			if(i <= l_profile)
			{
				exclusion_profile[i] = true;
			}
		} // i loop.
	} // i_reg loop.

	l_pruned_profile = 1;

	for(int i = 1; i <= l_profile; i++)
	{
		if(exclusion_profile[i])
		{
		}
		else
		{
			pruned_signal_profile[l_pruned_profile] = signal_profile[i];
			l_pruned_profile++;
		}
	} // i loop.

	// Must move the pruned profile index back by one.
	l_pruned_profile--;

	fprintf(stderr, "Pruned the signal of %d values to %d values.\n", l_profile, l_pruned_profile);
	delete [] exclusion_profile;
}

void get_log_plus_one_profile(double* signal_profile, double base, int l_profile)
{
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_log_value = xlog(signal_profile[i]+1.0)/xlog(base);
		signal_profile[i] = cur_log_value;
	} // i loop.
}

void floorize_profile(double* signal_profile, int l_profile)
{
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_floor_val = floor(signal_profile[i]);
		signal_profile[i] = cur_floor_val;
	} // i loop.
}

double* get_zero_profile(int l_profile)
{
	double* cur_profile = new double[l_profile+2];
	for(int i = 0; i <= l_profile; i++)
	{
		cur_profile[i] = 0.0;
	} // i loop.

	return(cur_profile);
}

vector<t_annot_region*>** get_peaks_per_per_nucleotide_signal_profile(double* signal_profile, char* chrom, int l_data, double min_thresh, double max_thresh)
{
	// Allocate the states of peaks at each position.
	vector<t_annot_region*>** peaks_per_threshes = new vector<t_annot_region*>*[(int)max_thresh-(int)min_thresh+1];
	int* peak_starts = new int[(int)max_thresh-(int)min_thresh+1];
	for(int thresh = (int)min_thresh; thresh <= (int)max_thresh; thresh++)
	{
		peaks_per_threshes[thresh-(int)min_thresh] = new vector<t_annot_region*>();
		peak_starts[thresh-(int)min_thresh] = 0;
	} // thresh loop.

	for(int i = 1; i <= l_data; i++)
	{
		for(int thresh = (int)min_thresh; thresh <= (int)max_thresh; thresh++)
		{
			if(signal_profile[i] >= thresh)
			{
				// The signal is below the threshold.
				if(peak_starts[thresh-(int)min_thresh] == 0)
				{
					// Set a new peak start.
					peak_starts[thresh-(int)min_thresh] = i;
				}
				else
				{
					// This was a peak already.
				}
			}
			else
			{
				// The signal is below the threshold.
				if(peak_starts[thresh-(int)min_thresh] != 0)
				{
					// Add this as a peak.
					t_annot_region* new_peak = get_empty_region();
					new_peak->chrom = t_string::copy_me_str(chrom);
					new_peak->start = peak_starts[thresh-(int)min_thresh];
					new_peak->end = i-1;
					new_peak->strand = '+';

					// Set the current start at this threshold to 0; no peak.
					peak_starts[thresh-(int)min_thresh] = 0;

					// Add the new peak.
					peaks_per_threshes[thresh-(int)min_thresh]->push_back(new_peak);
				}
				else
				{
					// This was not a peak already.
				}
			}
		} // thresh loop.
	} // i loop.

	delete [] peak_starts;

	return(peaks_per_threshes);
}

