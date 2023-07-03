#include <stdio.h>
#include <stdlib.h>
#include "cllgn_ansi_string.h"
#include "cllgn_file_utils.h"
#include "cllgn_multicolumn_processing.h"
#include <vector>
#include <algorithm>
using namespace std;

// Define the static arrays.
vector<char*>* t_g_Column_IDs::column_ids[N_COL_IDS];

bool sort_multi_col_line_info(char** line1_info, char** line2_info)
{
	return(t_string::sort_strings_per_prefix(line1_info[0], line2_info[0]));
}

void extract_rows_per_query_column_preserve_query_order(char* queries_col_fp, int col_i, char* multicol_fp, bool write_all_matches, char* op_fp)
{
	vector<char*>* query_entries = buffer_file(queries_col_fp);
	fprintf(stderr, "Searching for entries of %d query values at %d. column.\n", (int)query_entries->size(), col_i);

	vector<char*>* multicol_lines = buffer_file(multicol_fp);
	vector<char**>* sorted_multicol_line_info = new vector<char**>();
	for (int i_l = 0; i_l < (int)multicol_lines->size(); i_l++)
	{
		vector<char*>* toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(multicol_lines->at(i_l), "\t", true));
		char** cur_line_info = new char*[2];
		cur_line_info[0] = toks->at(col_i);
		cur_line_info[1] = multicol_lines->at(i_l);

		sorted_multicol_line_info->push_back(cur_line_info);
	} // i_l loop.

	// Sort the multicolumn lines per selected column.
	vector<char*>* sorted_multicol_values = new vector<char*>();
	sort(sorted_multicol_line_info->begin(), sorted_multicol_line_info->end(), sort_multi_col_line_info);
	for (int i_l = 0; i_l < (int)sorted_multicol_line_info->size(); i_l++)
	{
		sorted_multicol_values->push_back(sorted_multicol_line_info->at(i_l)[0]);
	} // i_l loop.

	int n_matches = 0;
	FILE* f_op = open_f(op_fp, "w");
	for (int i_q = 0; i_q < (int)query_entries->size(); i_q++)
	{
		int row_i = t_string::fast_search_string_per_prefix(query_entries->at(i_q), sorted_multicol_values, 0, (int)sorted_multicol_values->size());

		while (row_i < (int)sorted_multicol_values->size() &&
			row_i > 0 &&
			(t_string::sort_strings_per_prefix(query_entries->at(i_q), sorted_multicol_values->at(row_i)) ||
				t_string::compare_strings(query_entries->at(i_q), sorted_multicol_values->at(row_i))))
		{
			row_i--;
		}

		bool found_match = false;
		while (row_i < (int)sorted_multicol_values->size() &&
			(t_string::sort_strings_per_prefix(sorted_multicol_values->at(row_i), query_entries->at(i_q)) ||
				t_string::compare_strings(sorted_multicol_values->at(row_i), query_entries->at(i_q))))
		{
			if (t_string::compare_strings(query_entries->at(i_q), sorted_multicol_values->at(row_i)))
			{
				fprintf(f_op, "%s\n", sorted_multicol_line_info->at(row_i)[1]);

				if (!found_match)
				{
					found_match = true;
					n_matches++;
				}

				// If not writing all matches, break out of the search.
				if (!write_all_matches)
				{
					break;
				}
			}

			row_i++;
		}
	} // i_l loop.
	close_f(f_op, op_fp);

	fprintf(stderr, "Matched %d queries to %d rows.\n", n_matches, (int)query_entries->size());
}

void extract_rows_per_query_column_preserve_multicol_order(char* queries_col_fp, int col_i, char* multicol_fp, bool write_all_matches, char* op_fp)
{
	vector<char*>* query_entries = buffer_file(queries_col_fp);
	fprintf(stderr, "Searching for entries of %d query values at %d. column.\n", (int)query_entries->size(), col_i);

	sort(query_entries->begin(), query_entries->end(), t_string::sort_strings_per_prefix);

	vector<char*>* multicol_lines = buffer_file(multicol_fp);
	int n_matches = 0;
	FILE* f_op = open_f(op_fp, "w");
	for (int i_l = 0; i_l < (int)multicol_lines->size(); i_l++)
	{
		vector<char*>* toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(multicol_lines->at(i_l), "\t", true));
		int query_i = t_string::fast_search_string_per_prefix(toks->at(col_i), query_entries, 0, (int)query_entries->size());

		while (query_i < (int)query_entries->size() &&
			query_i > 0 &&
			(t_string::sort_strings_per_prefix(toks->at(col_i), query_entries->at(query_i)) ||
				t_string::compare_strings(toks->at(col_i), query_entries->at(query_i))))
		{
			query_i--;
		}

		bool found_match = false;
		while (query_i < (int)query_entries->size() &&
			(t_string::sort_strings_per_prefix(query_entries->at(query_i), toks->at(col_i)) ||
				t_string::compare_strings(toks->at(col_i), query_entries->at(query_i))))
		{
			if (t_string::compare_strings(toks->at(col_i), query_entries->at(query_i)))
			{
				fprintf(f_op, "%s\n", multicol_lines->at(i_l));

				if (!found_match)
				{
					found_match = true;
					n_matches++;
				}

				// If not writing all matches, break out of the search.
				if (!write_all_matches)
				{
					break;
				}
			}

			query_i++;
		}
	} // i_l loop.
	close_f(f_op, op_fp);

	fprintf(stderr, "Matched %d queries to %d rows.\n", n_matches, (int)query_entries->size());
}

void sort_columns_per_external_col_ids(char* multicolumn_fp, char* sorting_col_ids_fp, char* op_fp)
{
	vector<char*>* sorting_col_ids = buffer_file(sorting_col_ids_fp);
	fprintf(stderr, "Loaded %d sorting columns.\n", (int)sorting_col_ids->size());

	vector<char*>* multicol_lines = buffer_file(multicolumn_fp);
	fprintf(stderr, "Loaded %d lines\n", (int)multicol_lines->size());

	vector<char*>* multicol_col_ids = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(multicol_lines->at(0), "\t", true));

	FILE* f_op = open_f(op_fp, "w");
	for (int i_l = 0; i_l < (int)multicol_lines->size(); i_l++)
	{
		vector<char*>* cur_line_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(multicol_lines->at(i_l), "\t", true));
		if (cur_line_toks->size() != multicol_col_ids->size())
		{
			fprintf(stderr, "Could not find %d columns in: %s\n", (int)multicol_col_ids->size(), multicol_lines->at(i_l));
			exit(0);
		}

		for (int i_sorting_col = 0; i_sorting_col < (int)sorting_col_ids->size(); i_sorting_col++)
		{
			if (i_sorting_col > 0)
			{
				fprintf(f_op, "\t");
			}

			int file_col_i = t_string::get_i_str(multicol_col_ids, sorting_col_ids->at(i_sorting_col));
			if (file_col_i == (int)multicol_col_ids->size())
			{
				fprintf(f_op, ".");
			}
			else
			{
				fprintf(f_op, "%s", cur_line_toks->at(file_col_i));
			}
		} // i_sorting_col loop.

		// This line is processed.
		fprintf(f_op, "\n");
	} // i_l loop.
	close_f(f_op, op_fp);
}

int link_multi_column_files(char* file1_path, int file1_col_i, char* file2_path, int file2_col_i, char* op_fp)
{
	fprintf(stderr, "Linking %s and %s @ %d. and %d. columns.\n", file1_path, file2_path, file1_col_i, file2_col_i);

	vector<char*>* file1_lines = buffer_file(file1_path);
	vector<char*>* file2_lines = buffer_file(file2_path);
	fprintf(stderr, "Loaded %d and %d lines.\n", (int)file1_lines->size(), (int)file2_lines->size());

	vector<t_str_node*>* file1_cols = new vector<t_str_node*>();
	for (int il1 = 0; il1 < (int)file1_lines->size(); il1++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(file1_lines->at(il1), "\t", true);

		if (file1_col_i < (int)toks->size())
		{
			t_str_node* node = new t_str_node();
			node->str = t_string::copy_me_str(toks->at(file1_col_i)->str());
			node->dat = file1_lines->at(il1);
			node->i = il1;
			file1_cols->push_back(node);

			if ((int)file1_cols->size() < 5)
			{
				fprintf(stderr, "%s: %s (%d)\n", file1_path, toks->at(file1_col_i)->str(), (int)toks->size());
			}
		}

		t_string::clean_tokens(toks);
	} // il1 loop.

	vector<t_str_node*>* file2_cols = new vector<t_str_node*>();
	for (int il2 = 0; il2 < (int)file2_lines->size(); il2++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(file2_lines->at(il2), "\t", true);

		if (file2_col_i < (int)toks->size())
		{
			t_str_node* node = new t_str_node();
			node->str = t_string::copy_me_str(toks->at(file2_col_i)->str());
			node->dat = file2_lines->at(il2);
			node->i = il2;
			file2_cols->push_back(node);

			if ((int)file2_cols->size() < 5)
			{
				fprintf(stderr, "%s: %s (%d)\n", file2_path, toks->at(file2_col_i)->str(), (int)toks->size());
			}
		}

		t_string::clean_tokens(toks);
	} // il1 loop.

	int n_matched_cols = 0;
	FILE* f_op = open_f(op_fp, "w");
	for (int il1 = 0; il1 < (int)file1_cols->size(); il1++)
	{
		for (int il2 = 0; il2 < (int)file2_cols->size(); il2++)
		{
			if (t_string::compare_strings(file1_cols->at(il1)->str, file2_cols->at(il2)->str))
			{
				char* cur_line1 = (char*)(file1_cols->at(il1)->dat);
				char* cur_line2 = (char*)(file2_cols->at(il2)->dat);
				fprintf(f_op, "%s\t%s\n", cur_line1, cur_line2);
				n_matched_cols++;

				//if (!link_all)
				//{
				break;
				//}				
			}
		} // il2 loop.
	} // il1 loop.
	fclose(f_op);

	fprintf(stderr, "%d columns matched.\n", n_matched_cols);
	return(0);
}

int link_multi_column_files_per_header_col_name(char* file1_path, char* file1_col_id, char* file2_path, char* file2_col_id, char* op_fp)
{
	char* header1 = load_header(file1_path);
	char* header2 = load_header(file2_path);

	vector<char*>* header1_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(header1, "\t", true));
	vector<char*>* header2_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(header2, "\t", true));

	fprintf(stderr, "%d columns in %s and %d columns in %s\n", (int)header1_toks->size(), file1_path, (int)header2_toks->size(), file2_path);

	int file1_col_i = t_string::get_i_str(header1_toks, file1_col_id);
	int file2_col_i = t_string::get_i_str(header2_toks, file2_col_id);

	if (file1_col_i == (int)header1_toks->size() ||
		file2_col_i == (int)header2_toks->size())
	{
		fprintf(stderr, "Could not find column id's in the header.\n");
		exit(0);
	}

	fprintf(stderr, "Found columns @ %s::%d (%s, %s)\n%s::%d (%s, %s)\n",
		file1_path, file1_col_i, header1_toks->at(file1_col_i), file1_col_id,
		file2_path, file2_col_i, header2_toks->at(file2_col_i), file2_col_id);

	link_multi_column_files(file1_path, file1_col_i, file2_path, file2_col_i, op_fp);

	// Add the header.
	char* new_header = new char[t_string::string_length(header1) + t_string::string_length(header2) + 2];
	sprintf(new_header, "%s\t%s", header1, header2);

	vector<char*>* pasted_file_lines = buffer_file(op_fp);

	fprintf(stderr, "Saving %d pasted lines to %s.\n", (int)pasted_file_lines->size(), op_fp);
	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%s\n", new_header);
	for (int i_l = 0; i_l < (int)pasted_file_lines->size(); i_l++)
	{
		fprintf(f_op, "%s\n", pasted_file_lines->at(i_l));
	} // i_l loop.
	close_f(f_op, op_fp);

	return(0);
}

vector<vector<char*>*>* load_multicolumn_file(char* multicolumn_fp, const char* delims, int& n_loaded_columns)
{
	vector<vector<char*>*>* per_line_per_col_entries = new vector<vector<char*>*>();

	n_loaded_columns = 0;
	FILE* f_dat = open_f(multicolumn_fp, "r");
;	while (1)
	{
		char* cur_line = getline(f_dat);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, delims, true);
		if (n_loaded_columns == 0)
		{
			n_loaded_columns = (int)toks->size();
		}
		else if(n_loaded_columns != (int)toks->size())
		{
			fprintf(stderr, "The loaded cols are not consistent!\n");
			return(NULL);
		}
		
		per_line_per_col_entries->push_back(t_string::copy_tokens_2_strs(toks));
		t_string::clean_tokens(toks);
	} // file reading loop.

	close_f(f_dat, multicolumn_fp);

	return(per_line_per_col_entries);
}

void extract_subcolumn_list_from_multicolumn_file(char* multicol_fp, char* col_ids_list_fp, char* op_fp)
{
	// The multicolumn file.
	vector<char*>* multi_col_lines = buffer_file(multicol_fp);
	fprintf(stderr, "Loaded %d lines from %s.\n", (int)multi_col_lines->size(), multicol_fp);

	char* header = multi_col_lines->at(0);
	t_string_tokens* header_toks = t_string::tokenize_by_chars(header, "\t", true);
	vector<char*>* header_col_strs = t_string::copy_tokens_2_strs(header_toks);
	fprintf(stderr, "Found %d columns in the header.\n", (int)header_col_strs->size());

	// Match the column id's.
	vector<char*>* col_ids_2_extract = buffer_file(col_ids_list_fp);
	fprintf(stderr, "Looking for %d column id's from %s.\n", (int)col_ids_2_extract->size(), col_ids_list_fp);

	vector<int>* matched_col_ids = new vector<int>();
	for (int i_col = 0; i_col < (int)col_ids_2_extract->size(); i_col++)
	{
		int matched_col_i = t_string::get_i_str(header_col_strs, col_ids_2_extract->at(i_col));
		if (matched_col_i < (int)header_col_strs->size())
		{
			fprintf(stderr, "Matched %s (%d)\n", col_ids_2_extract->at(i_col), matched_col_i);
			matched_col_ids->push_back(matched_col_i);
		}
		else
		{
			fprintf(stderr, "Could not match %s\n", col_ids_2_extract->at(i_col));
		}
	} // i_col loop.

	if ((int)matched_col_ids->size() == 0)
	{
		fprintf(stderr, "Could not match any column id's.\n");
		exit(0);
	}

	// Sort the matched column id's.
	sort(matched_col_ids->begin(), matched_col_ids->end());

	fprintf(stderr, "Extracting %d columns.\n", (int)matched_col_ids->size());

	// We do not have to explicitly write the header.
	FILE* f_op = open_f(op_fp, "w");
	char cur_tok_buffer[1000];
	for (int i_l = 0; i_l < (int)multi_col_lines->size(); i_l++)
	{
		char* cur_multicol_line = multi_col_lines->at(i_l);
		int str_char_i = 0;
		
		int cur_i_tok = -1;
		for (int i_col = 0; i_col < (int)matched_col_ids->size(); i_col++)
		{
			// Find and read the next column id.
			for (int i_tok = cur_i_tok; i_tok < matched_col_ids->at(i_col); i_tok++)
			{
				t_string::get_next_token(cur_multicol_line, cur_tok_buffer, 1000, "\t", str_char_i);
			} // Read the next token.

			cur_i_tok = matched_col_ids->at(i_col);

			if (i_col > 0)
			{
				fprintf(f_op, "\t");
			}

			fprintf(f_op, "%s", cur_tok_buffer);
		} // i_col loop.

		fprintf(f_op, "\n");
	} // i_l loop.
	fclose(f_op);
}

// Copy all the column id's.
void t_g_Column_IDs::init_global_IDs()
{
	t_g_Column_IDs::column_ids[COL_ID_CHR] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_CHR]->push_back(t_string::copy_me_str("chrom"));
	t_g_Column_IDs::column_ids[COL_ID_CHR]->push_back(t_string::copy_me_str("chr"));

	t_g_Column_IDs::column_ids[COL_ID_CONVERTED_VCF_N_ALT_SUPP_READS] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_CONVERTED_VCF_N_ALT_SUPP_READS]->push_back(t_string::copy_me_str("n_alt_supp_reads"));

	t_g_Column_IDs::column_ids[COL_ID_CONVERTED_VCF_N_REF_SUPP_READS] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_CONVERTED_VCF_N_REF_SUPP_READS]->push_back(t_string::copy_me_str("n_ref_supp_reads"));

	t_g_Column_IDs::column_ids[COL_ID_CONVERTED_VCF_N_TOTAL_READS] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_CONVERTED_VCF_N_TOTAL_READS]->push_back(t_string::copy_me_str("n_total_supp_reads"));

	// BED column ids.
	t_g_Column_IDs::column_ids[COL_ID_BED_START] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_BED_START]->push_back(t_string::copy_me_str("start"));

	t_g_Column_IDs::column_ids[COL_ID_BED_END] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_BED_END]->push_back(t_string::copy_me_str("end"));

	t_g_Column_IDs::column_ids[COL_ID_BED_NAME] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_BED_NAME]->push_back(t_string::copy_me_str("name"));

	t_g_Column_IDs::column_ids[COL_ID_BED_SCORE] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_BED_SCORE]->push_back(t_string::copy_me_str("score"));

	t_g_Column_IDs::column_ids[COL_ID_BED_STRAND] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_BED_STRAND]->push_back(t_string::copy_me_str("strand"));

	// VCF column id's
	t_g_Column_IDs::column_ids[COL_ID_VCF_CHR] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_CHR]->push_back(t_string::copy_me_str("chr"));
	t_g_Column_IDs::column_ids[COL_ID_VCF_CHR]->push_back(t_string::copy_me_str("chrom"));
	t_g_Column_IDs::column_ids[COL_ID_VCF_CHR]->push_back(t_string::copy_me_str("chromosome"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_POSN] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_POSN]->push_back(t_string::copy_me_str("posn"));
	t_g_Column_IDs::column_ids[COL_ID_VCF_POSN]->push_back(t_string::copy_me_str("pos"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_REF_ALL] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_REF_ALL]->push_back(t_string::copy_me_str("ref_allele"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_ALT_ALL] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_ALT_ALL]->push_back(t_string::copy_me_str("alt_allele"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_QUAL] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_QUAL]->push_back(t_string::copy_me_str("qual"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_FILTER] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_FILTER]->push_back(t_string::copy_me_str("filter"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_VARIANT_ID] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_VARIANT_ID]->push_back(t_string::copy_me_str("variant_id"));
	t_g_Column_IDs::column_ids[COL_ID_VCF_VARIANT_ID]->push_back(t_string::copy_me_str("id"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_FORMAT] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_FORMAT]->push_back(t_string::copy_me_str("format"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_GT] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_GT]->push_back(t_string::copy_me_str("genotype"));
	t_g_Column_IDs::column_ids[COL_ID_VCF_GT]->push_back(t_string::copy_me_str("gt"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_FILTER] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_FILTER]->push_back(t_string::copy_me_str("filter"));

	t_g_Column_IDs::column_ids[COL_ID_VCF_INFO] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_VCF_INFO]->push_back(t_string::copy_me_str("info"));

	// Annotated variant column id's.
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_TYPE] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_TYPE]->push_back(t_string::copy_me_str("type"));

	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_EFFECT] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_EFFECT]->push_back(t_string::copy_me_str("effect"));
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_EFFECT]->push_back(t_string::copy_me_str("variant_effect"));

	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_GENE_SYMBOL] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_GENE_SYMBOL]->push_back(t_string::copy_me_str("gene_symbol"));
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_GENE_SYMBOL]->push_back(t_string::copy_me_str("variant_gene_symbol"));

	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_ENSEMBL_TRANSCRIPT] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_ENSEMBL_TRANSCRIPT]->push_back(t_string::copy_me_str("ensembl_transcript"));
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_ENSEMBL_TRANSCRIPT]->push_back(t_string::copy_me_str("variant_ensembl_transcript"));

	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_FEAT_LENGTH] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_ANNO_VAR_ADDED_FEAT_LENGTH]->push_back(t_string::copy_me_str("feat_length"));

	// KB column id's.
	t_g_Column_IDs::column_ids[COL_ID_KB_dbSNP_COL_ID] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_dbSNP_COL_ID]->push_back(t_string::copy_me_str("dbSNP"));

	// COSMIC column id's; two of them currently.
	t_g_Column_IDs::column_ids[COL_ID_KB_COSMIC_MATCH_VARIANT_COL_ID] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_COSMIC_MATCH_VARIANT_COL_ID]->push_back(t_string::copy_me_str("COSMIC_variant_match"));

	t_g_Column_IDs::column_ids[COL_ID_KB_COSMIC_MATCH_TRANSCRIPT_COL_ID] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_COSMIC_MATCH_TRANSCRIPT_COL_ID]->push_back(t_string::copy_me_str("COSMIC_transcript_match"));
	
	t_g_Column_IDs::column_ids[COL_ID_KB_GWAS] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_GWAS]->push_back(t_string::copy_me_str("GWAS"));

	t_g_Column_IDs::column_ids[COL_ID_KB_REACTOME] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_REACTOME]->push_back(t_string::copy_me_str("REACTOME"));

	t_g_Column_IDs::column_ids[COL_ID_KB_BIOGRID] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_BIOGRID]->push_back(t_string::copy_me_str("BIOGRID"));

	t_g_Column_IDs::column_ids[COL_ID_KB_EXPRESSION_ATLAS] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_EXPRESSION_ATLAS]->push_back(t_string::copy_me_str("Expression_Atlas"));

	t_g_Column_IDs::column_ids[COL_ID_KB_CLINVAR] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_CLINVAR]->push_back(t_string::copy_me_str("ClinVar"));

	t_g_Column_IDs::column_ids[COL_ID_KB_CLINICAL_TRIALS] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_CLINICAL_TRIALS]->push_back(t_string::copy_me_str("Clinical_Trials"));

	t_g_Column_IDs::column_ids[COL_ID_KB_DRUGBANK] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_DRUGBANK]->push_back(t_string::copy_me_str("Drugbank"));

	t_g_Column_IDs::column_ids[COL_ID_KB_OMIM] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_OMIM]->push_back(t_string::copy_me_str("OMIM"));

	t_g_Column_IDs::column_ids[COL_ID_KB_1kG] = new vector<char*>();
	t_g_Column_IDs::column_ids[COL_ID_KB_1kG]->push_back(t_string::copy_me_str("1kG"));
}

int t_g_Column_IDs::get_col_i_per_identifier(vector<char*>* col_ids, int col_identifier)
{
	vector<char*>* col_aliases = t_g_Column_IDs::column_ids[col_identifier];
	for(int i_alias = 0; i_alias < (int)col_aliases->size(); i_alias++)
	{
		int cur_alias_col_i = t_string::get_i_str(col_ids, col_aliases->at(i_alias));
		if(cur_alias_col_i < (int)col_ids->size())
		{
			return(cur_alias_col_i);
		}
	} // i_alias loop.

	return((int)col_ids->size());
}

// Add the column id's as a hdeader.
void add_col_ids_list_2_header(char* col_ids_list_fp, char* multicol_fp, char* op_fp)
{
	vector<char*>* col_ids = buffer_file(col_ids_list_fp);
	if(col_ids == NULL)
	{
		fprintf(stderr, "Could not load column id's from %s\n", col_ids_list_fp);
		exit(0);
	}

	FILE* f_op = open_f(op_fp, "w");
	FILE* f_multicol = open_f(multicol_fp, "r");

	fprintf(f_op, "#%s", col_ids->at(0));
	for(int col_id_i = 1; col_id_i < (int)col_ids->size(); col_id_i++)
	{
		fprintf(f_op, "\t%s", col_ids->at(col_id_i));
	} // col_id_i loop.
	fprintf(f_op, "\n");

	while(1)
	{
		char* cur_line = getline(f_multicol);
		if(cur_line == NULL)
		{
			break;
		}
		else
		{
			fprintf(f_op, "%s\n", cur_line);
			delete [] cur_line;
		}
	} // file reading loop.

	fclose(f_multicol);
	fclose(f_op);
}

t_g_Column_IDs::t_g_Column_IDs()
{
	t_g_Column_IDs::init_global_IDs();
}


t_g_Column_IDs::~t_g_Column_IDs()
{
}