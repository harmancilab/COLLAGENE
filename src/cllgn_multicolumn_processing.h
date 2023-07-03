#ifndef _MULTICOLUMN_PROCESSING_
#define _MULTICOLUMN_PROCESSING_

#include <vector>
using namespace std;

// This is an attempt to standardize the column id's used by the codebase. 
// When the multicolumn files are dumped by the codebase tools where the column id's are not standardized, the column id processing enables handling the 
// variable columns without hardcoding the column id's multiple times. They are hardcoded once in this class.

enum
{
	COL_ID_CHR,
	COL_ID_CONVERTED_VCF_N_ALT_SUPP_READS,
	COL_ID_CONVERTED_VCF_N_REF_SUPP_READS,
	COL_ID_CONVERTED_VCF_N_TOTAL_READS,

	// BED column ids.
	COL_ID_BED_START, 
	COL_ID_BED_END, 
	COL_ID_BED_NAME, 
	COL_ID_BED_SCORE, 
	COL_ID_BED_STRAND, 

	// VCF column id's
	COL_ID_VCF_CHR, 
	COL_ID_VCF_POSN, 
	COL_ID_VCF_REF_ALL,
	COL_ID_VCF_ALT_ALL,
	COL_ID_VCF_QUAL,
	COL_ID_VCF_FILTER,
	COL_ID_VCF_INFO,
	COL_ID_VCF_VARIANT_ID,
	COL_ID_VCF_FORMAT,
	COL_ID_VCF_GT,

	// Annotated variant column id's.
	COL_ID_ANNO_VAR_ADDED_TYPE,
	COL_ID_ANNO_VAR_ADDED_EFFECT,
	COL_ID_ANNO_VAR_ADDED_GENE_SYMBOL,
	COL_ID_ANNO_VAR_ADDED_ENSEMBL_TRANSCRIPT,
	COL_ID_ANNO_VAR_ADDED_FEAT_LENGTH,

	// KB column id's.
	COL_ID_KB_dbSNP_COL_ID,

	// COSMIC column id's; two of them currently.
	COL_ID_KB_COSMIC_MATCH_VARIANT_COL_ID,
	COL_ID_KB_COSMIC_MATCH_TRANSCRIPT_COL_ID,

	COL_ID_KB_GWAS,
	COL_ID_KB_REACTOME,
	COL_ID_KB_BIOGRID,
	COL_ID_KB_EXPRESSION_ATLAS,
	COL_ID_KB_CLINVAR,
	COL_ID_KB_CLINICAL_TRIALS,
	COL_ID_KB_DRUGBANK,
	COL_ID_KB_OMIM,
	COL_ID_KB_1kG,

	N_COL_IDS
};

void extract_subcolumn_list_from_multicolumn_file(char* multicol_fp, char* col_ids_list_fp, char* op_fp);

vector<vector<char*>*>* load_multicolumn_file(char* multicolumn_fp, const char* delims, int& n_loaded_columns);

int link_multi_column_files(char* file1_path, int file1_col_i, char* file2_path, int file2_col_i, char* op_fp);

bool sort_multi_col_line_info(char** line1_info, char** line2_info);

void extract_rows_per_query_column_preserve_query_order(char* queries_col_fp, int col_i, char* multicol_fp, bool write_all_matches, char* op_fp);

void extract_rows_per_query_column_preserve_multicol_order(char* queries_col_fp, int col_i, char* multicol_fp, bool write_all_matches, char* op_fp);

void sort_columns_per_external_col_ids(char* multicolumn_fp, char* sorting_col_ids_fp, char* op_fp);

int link_multi_column_files_per_header_col_name(char* file1_path, char* file1_col_id, char* file2_path, char* file2_col_id, char* op_fp);

class t_g_Column_IDs
{
public:
	static vector<char*>* column_ids[N_COL_IDS];

	// Initialize and set the column ids.
	t_g_Column_IDs();
	~t_g_Column_IDs();

	static void init_global_IDs();

	static void validate_columns_per_data_type();

	static int get_col_i_per_identifier(vector<char*>* col_ids, int col_i);
};

#endif // _MULTICOLUMN_PROCESSING_