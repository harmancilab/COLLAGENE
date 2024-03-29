#ifndef __MAPPED_READ_FILE_INTERFACE__
#define __MAPPED_READ_FILE_INTERFACE__

#include <vector>
using namespace std;

//struct t_profile_site;
struct t_mapped_fragment;
struct t_frag_end;

const int MEG_BASE = 1000 * 1000;
const int K_BASE = 1000;

class t_string;
struct t_annot_region;

// MApping information for a mapped read.
struct t_mapping_info
{
	char* chrom;
	int base_posn;
	int sequenced_length;
	char strand;
	char* mapping_quality_str; // CIGAR string.
};

struct t_variant_per_cell_stats
{
	char* ref_allele;
	char* alt_allele;

	int* ref_allele_cnt_per_cell;
	int* alt_allele_cnt_per_cell;

	// Comes from variant_tools.h
	int var_type;
};

// This is a sequenced tag w/o any mapping information: Does not have any information about mapping, initially. In case there is
// mapping information, it is stored under the member named mapping_info.
struct t_sequenced_read
{
	char* id;
	t_mapping_info* mapping_info;
	char* quality_str;
	char* nucs;
};

// This is a mapped fragment: It can be a part of a read.
struct t_mapped_fragment
{
	int base_index;
	int sequenced_fragment_length; // Necessary when building enrichment profile. This is the fragment length coming from ChIP-Seq data.
	char strand_char;
};

struct t_mapped_read
{
	char* mapping_str;

	// The start of the read: 5' position.
	int base_index;

	// Left side posn.
	//int left_base_posn;

	// This is the total span of the read.
	int span; 
	char strand;
};

struct t_per_bin_10x_discordancy_stat
{
	int cur_bin_start;
	int n_discordant_reads;
	double avg_multimapp;
	vector<double>* cis_discordant_linked_read_barcodes;
	int n_cis_discordant_reads;
	vector<double>* trans_discordant_linked_read_barcodes;
	int n_trans_discordant_reads;
};

struct t_10x_linked_read_info
{
	int posn;
	int flag;
};

// Following packs some values to store space.
struct t_10X_BX_barcode_group_entry
{
	double BX_barcode;
	vector<int>* self_chr_w_mate_chr_w_flag;
	vector<int>* posn;
	vector<int>* mate_posn;
};

void compress_nucleotide_seq(char* fragment, int l_frag, char* compressed_seq, int& l_comp_seq);
void decompress_nucleotide_seq(char* compressed_seq, int l_comp_seq, char* fragment_buffer, int l_frag);
void convert_bSAM_2_SAM_tagged(char* bSAM_fp, char* SAM_op_fp);
void convert_SAM_2_bSAM_tagged(char* sam_fp, char* op_fp);
void convert_SAM_2_bSAM_untagged(char* sam_fp, char* op_fp);
void convert_bSAM_2_SAM_untagged(char* bSAM_fp, char* SAM_op_fp);
void get_empty_qual(int l_read, char* qual_buffer, char base_qual);

void get_unique_10X_CBZ_tags(char* tenX_sam_fp, int min_mapQ, char* op_fp);
void get_unique_10X_BX_tags(char* tenX_sam_fp, int min_mapQ, char* op_fp);
void parse_10x_linked_reads_per_BX_tag(char* tenX_sam_fp, vector<char*>* chr_ids, char* GEM_barcodes_list_fp, int min_mapQ, char* op_fp);

void compute_discordancy_stats_per_10x_linked_reads(char* parsed_linked_reads_fp, 
													char* chr_ids_lengths_list_fp,
													int l_bin,
													char* multimapp_dir,
													char* op_fp);

void get_differential_discordancy_stats_per_10x_linked_reads(vector<char*>* chr_ids,																
																char* tumor_discordancy_fp, 
																char* normal_discordancy_fp,
																int l_bin,
																char* op_fp);

void compute_single_cell_allelic_stats_per_10X_SAM(char* per_cell_barcodes_fp, 
													char* SAM_fp, 
													char* chr_info_list_fp, 
													char* genome_seq_dir, 
													char* candidate_snvs_bed_fp, char* candidate_indels_bed_fp, char* op_fp);

double* load_normalized_multimappability_profile(char* mapability_signal_profile_fp, int& l_mapability_profile);

#define MAX_N_PAIRS (10)
struct t_mapped_PE_read
{
	char* mapping_str;
	int base_index;

	int span;
	char strand;

	// The list of pairs that are matching in id with this pair.
	t_mapped_PE_read** matching_id_pairs;
	int n_matching_id_pairs;
};

// 4 column preprocessed reads w id processing.
struct t_read_line_w_id
{
	char* id;
	char pair_position; // 0=first, 1=last.
	char* read_line;
};

bool get_SAM_read_tag_entry(char* sam_read_line, char* tag_entry_buffer, const char* tag_prefix);

bool set_preprocessed_read_file_path_per_dir_chr(char* preprocessed_reads_dir, char* chrom, char* preprocessed_reads_fp);

void dump_pileup_SNV_candidates(char* chr_info_fp, char* pileup_dir, char* bin_seq_dir, int min_covg, int min_alternate_covg, double min_alternate_freq, char* op_fp);

void dump_pileup_SNV_candidates_per_stranded_pileup(char* chr_info_fp,
	char* pileup_strand_0_dir, char* pileup_strand_1_dir,
	char* bin_seq_dir, double min_covg, double min_alternate_covg, double min_alternate_freq, double max_strand_imbalance,
	char* op_fp);

vector<int>* count_preprocessed_reads(char* preprocessed_reads_dir, vector<char*>* chr_ids);

void get_signal_statistics_in_regions_per_pileups_dir(char* chrom_info_fp, char* pileups_dir_path, char* regs_bed_fp);

bool sort_read_line_entries_per_id(t_read_line_w_id* read1, t_read_line_w_id* read2);
vector<char*>* sort_read_lines_per_id(vector<char*>* cur_chr_read_lines);
void sort_read_lines_per_id_in_place(vector<char*>* cur_chr_read_lines);

// Following represents one connection.
struct t_connection_info_entry
{
	t_annot_region* connecting_region;
	vector<char*>* read_ids;

	// Strands for the read.
	vector<char>* first_read_strands;
	vector<char>* last_read_strands;

	// Starts and ends of the read starts.
	vector<int>* first_read_starts;
	vector<int>* first_read_i_chrs;
	vector<int>* last_read_starts;
	vector<int>* last_read_i_chrs;
};

struct t_per_bin_connection_info
{
	double avg_mappability;
	double avg_read_mappability;
	vector<t_connection_info_entry*>* connecting_bins;
};

void get_insert_stats_per_concordant_PE_reads(char* sorted_pe_reads_fp, 
	int l_read,
	int max_concordant_mapping_separation, 
	char valid_first_mapper_strand,
	char valid_last_mapper_strand);

void label_merge_external_sort_PE_reads_per_id(char* first_mapped_reads_dir,
	char* last_mapped_reads_dir,
	char* temp_sorted_op_dir,
	char* sorted_op_fp);

void label_PE_reads_file(char* reads_fp, char side_char, char* chr_id, char* op_fp);

void sort_PE_reads_file_per_id_in_memory(char* mapped_reads_fp,
	char* sorted_op_fp);

void external_sort_PE_reads_per_file_list(char* read_fps, char* sorted_op_fp);

void label_merge_external_sort_SE_reads_per_id(char* mapped_reads_dir,
	char* temp_sorted_op_dir,
	char* sorted_op_fp);

//void generate_discordant_PE_read_connections_per_sorted_PE_reads_file(char* chr_ids_fp, char* sorted_pe_reads_fp, int l_bin, int min_l_same_chr_separation, int max_concordant_mapping_separation,
//	char valid_first_mapper_strand,
//	char valid_last_mapper_strand);

vector<int>* get_chromosome_lengths_per_mapped_reads(char* mapped_reads_dir);

void generate_discordant_PE_read_connections_per_sorted_PE_reads_file(char* chr_ids_lengths_fp, char* sorted_pe_reads_fp, 
	char* mapability_dir,
	int l_bin, 
	int l_step,
	int min_l_same_chr_separation, 
	int max_concordant_mapping_separation,
	char valid_first_mapper_strand,
	char valid_last_mapper_strand);

void generate_discordant_read_connections_per_id_sorted_SE_mapped_reads_file(char* chr_ids_fp, char* sorted_se_mapped_reads_fp, 
	int l_bin, int min_l_same_chr_separation, int max_concordant_mapping_separation,
	char valid_first_mapper_strand,
	char valid_last_mapper_strand);

void get_differential_discordant_PE_pair_connections(char* sample_fp, char* control_fp, char* op_fp);

bool check_genome_index_update_per_CIGAR_entry(char entry_char);
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char);

// Load sequenced reads directly from the fastq file.
// Loads and pools the reads.
void load_sequenced_reads_per_fastq(char* fastq_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void dump_phred_quality_distribution(vector<t_sequenced_read*>* sequenced_reads, char* op_fp);

// Subsample the reads to a desired size.
void subsample_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads, 
	vector<t_sequenced_read*>* subsampled_sequenced_reads, 
	int n_subsample_size);

void subsample_mapped_reads(vector<t_mapped_read*>* mapped_reads, 
	vector<t_mapped_read*>* subsampled_mapped_reads, 
	int n_subsample_size);

void dump_fastq(vector<t_sequenced_read*>* sequenced_reads, char* op_fastq_fp);
void delete_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads);

// Load the mapped read files and preprocess them, from which the sequenced fragments can be loaded.
void load_mapped_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_ELAND(char* eland_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_tagAlign(char* tagalign_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_bowtie(char* bowtie_fp, vector<t_sequenced_read*>* sequenced_reads);

void count_mapped_reads_per_file(char* mrf_fp, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	double n_total_reads);

// Following are the preprocessing functions for preprocessing the mapped reads file 
void preprocess_mapped_PE_SAM_file(char* pe_sam_fp,
	int min_mapping_qual, 
	char* first_reads_dir, char* last_reads_dir);

void preprocess_PE_fragments_per_pos_sorted_SAM_file(char* mrf_fp, char* parsed_reads_op_dir, 
	void (preprocess_mapped_read_line)(char* cur_line,
		char* read_id,
		char* chrom,
		int& chr_index, int& sequenced_length,
		char& strand_char,
		char* mapping_quality_str),
	int max_l_fragment,
	bool dump_read_id);

void preprocess_mapped_reads_file_single_file_buffering(char* mrf_fp, char* parsed_reads_op_dir,
	void (preprocess_mapped_read_line)(char* cur_line,
		char* read_id,
		char* chrom,
		int& chr_index, int& sequenced_length,
		char& strand_char,
		char* mapping_quality_str),
	bool dump_read_id);

void preprocess_mapped_reads_file(char* mrf_fp,
	char* parsed_reads_op_dir,
	vector<char*>* preset_chr_ids,
	void (preprocess_mapped_read_line)(char* cur_line,
		char* read_id,
		char* chrom,
		int& chr_index, int& sequenced_length,
		char& strand_char,
		char* mapping_quality_str),
	bool dump_read_id);

void preprocess_tagAlign_read_line(char* cur_line,
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);

void preprocess_SAM_read_line_positional_info(char* cur_line,
	char* chrom,
	int& chr_index,
	int& mapQ,
	char& strand_char,
	char* cigar_str);

void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_preprocessed_LH_GFF3_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED5_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED6_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str);
void preprocess_BED456_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str);

void preprocess_PE_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	bool& first_segment_in_template,
	bool& last_segment_in_template,
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	int& mapping_quality,
	char* cigar_str);

void count_10X_reads_per_cell_per_barcodes_list(char* TenX_SAM_file, char* per_cell_barcode_list_fp, char* op_fp);
bool get_10X_cell_barcode_per_SAM_read(char* TenX_SAM_read_line, char* barcode_buffer);

bool sort_regions_coords_first_names_second(t_annot_region* reg1, t_annot_region* reg2);
void extract_summarize_indel_containing_read_blocks_per_SAM(char* SAM_fp, char* chrom_info_fp, char* genome_seq_dir, char* op_fp);

void scan_indels_per_summarized_indel_blocks(char* chr_ids_lengths_fp, char* indel_blocks_dir, char* genome_seq_dir,
	char* postv_pileup_dir, char* negtv_pileup_dir,
	int min_covg_per_indel, int min_alternate_covg_per_indel, double min_alternate_freq_per_indel,
	int l_indel_scanning_window,
	char* op_fp);

//void dump_nucleotide_pileup_per_SAM_file(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int& n_processed_reads);
void dump_nucleotide_pileup_per_SAM_file(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, unsigned long long& n_processed_reads);
void dump_nucleotide_pileup_per_SAM_file_phred_partitioning(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, vector<int>* phred_qual_partitions, unsigned long long& n_processed_reads);
void compress_nucleotide_pileup_track(unsigned short** pileup, int l_sig, char* op_fp);
unsigned short** allocate_pileup(int l_sig);
void delete_pileup(unsigned short** loaded_pileup);
int* load_coverage_per_compressed_pileup_file(char* comp_pileup_fp, int& l_sig);
int* load_coverage_per_pileups(unsigned short** pileups, int l_sig);
unsigned short** load_compressed_pileups(char* cur_comp_allele_fp, int& l_pileup);

// Partitioned pileup computes.
void compress_partitioned_nucleotide_pileup_track(unsigned short*** pileup, int n_partitions, int l_sig, char* op_fp);
int* load_coverage_per_compressed_partitioned_pileup_file(char* comp_pileup_fp, int& l_sig);
unsigned short*** load_partitioned_compressed_pileups(char* cur_comp_allele_fp, int& n_partitions, int& l_pileup);

int get_l_signal_per_reads(char* reads_fp, int l_ext_tag);

//void buffer_per_nucleotide_profile_no_buffer(char* sorted_read_fp, int l_extended_tag, double* signal_profile_buffer, int l_buffer, int& l_data);
void buffer_per_nucleotide_profile_no_buffer(const char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data);

void buffer_per_nucleotide_preprocessed_read_profile_no_buffer(char* sorted_read_fp,
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int max_l_read, // This is the length of the longest spliced read to be processed.
	int l_buffer, int& l_data);

double get_n_mapped_nucs(vector<t_mapped_fragment*>* fragments);

struct t_read_line_sorting_info
{
	int start;
	char* read_line;
};

bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2);

vector<char*>* sort_bucket_read_lines(char* bucket_fp);


void get_per_strand_read_stats(vector<t_annot_region*>* regions,
	char* preprocessed_reads_dir);




// BINARY INTERFACE IS OBSOLETE
void load_mapped_reads_file(char* read_file_path, vector<char*>* chr_ids, 
														vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
														void (preprocess_mapped_read_line)(char* cur_line, 
														char* chrom, 
														int& chr_index, int& sequenced_length, 
														char& strand_char, 
														char* mapping_quality_str));

// Load reads reads/fragments.
void load_fragments_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
	vector<vector<t_mapped_fragment*>*>* fore_strand_fragments_per_chr, vector<vector<t_mapped_fragment*>*>* rev_strand_fragments_per_chr, 
	int max_n_pcr_amplified_reads);

void load_fragments(char* mapped_reads_fp, 
	vector<t_mapped_fragment*>* fore_strand_fragments, vector<t_mapped_fragment*>* rev_strand_frag, 
	int max_n_pcr_amplified_reads = 10000);

void load_reads_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
	vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
	int max_n_pcr_amplified_reads);

void load_reads(char* mapped_reads_fp, 
	vector<t_mapped_read*>* fore_strand_reads, vector<t_mapped_read*>* rev_strand_reads, 
	int max_n_pcr_amplified_reads);

void get_mapped_fragments_per_mapped_reads(vector<t_mapped_read*>* mapped_reads, vector<t_mapped_fragment*>* mapped_fragments);
void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments);

// Following are for doing binary searches over the fragments.
bool sort_mapped_reads_per_3p(t_mapped_read* read1, t_mapped_read* read2);
bool sort_mapped_reads_per_5p(t_mapped_read* read1, t_mapped_read* read2);
int read_5p_accessor(void* obj_ptr);
//int read_3p_accessor(void* obj_ptr);

void delete_mapped_reads(vector<t_mapped_read*>* mapped_read);
void delete_mapped_read(t_mapped_read* mapped_read);

int fragment_5p_accessor(void* obj_ptr);
int fragment_3p_accessor(void* obj_ptr);

void get_read_statistics_per_region(vector<t_mapped_fragment*>* fragments, 
									vector<int>* all_fragment_3p_posns, 
									t_annot_region* region, 
									double& n_mapped_reads, 
									double& n_mapped_nucs);

bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced);
void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										int& l_cur_entry,
										char& entry_type_char);

void delete_fragments(vector<t_mapped_fragment*>* fragment_list);
void delete_fragments(t_mapped_fragment** fragment_list);

bool check_fragment_quality(char* fragment, char* chr_subseq, char* quality_str);
bool sort_mapped_fragments(t_mapped_fragment* frag1, t_mapped_fragment* frag2);
bool sort_mapped_fragments_per_3p(t_mapped_fragment* frag1, t_mapped_fragment* frag2);

void preprocessed_read_file_iterator(char* mapped_reads_fp, 
	void (per_read_callback)(char*, char, int, void*), 
	void (per_fragment_callback)(char*, char, int, void*),
	void* per_read_callback_param,
	void* per_fragment_callback_param);

// Prune fragments/reads
void prune_reads(vector<t_mapped_read*>* mapped_reads, int n_max_reps_per_posn, 
				vector<t_mapped_read*>* pruned_forward_reads, 
				vector<t_mapped_read*>* pruned_reverse_reads);
//vector<t_mapped_fragment*>* prune_fragment_list_by_strand(vector<t_mapped_fragment*>* chr_fragments, char strand_2_prune, int n_max_reps);

/*
Following function is the only function that does tag extension with enrichment_mapped_fragment_length parameter. This should not be handled anywhere else.
*/
vector<t_mapped_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_mapped_fragment*>* fore_frag_list, vector<t_mapped_fragment*>* rev_frag_list, int enrichment_mapped_fragment_length);

int* get_n_reads_per_window(int n_wins, vector<t_mapped_fragment*>* frag_list);

void exclude_reads_per_regions(vector<char*>* read_chr_ids, 
	vector<vector<t_mapped_read*>*>* reads_per_chrs, 
	vector<vector<t_mapped_read*>*>* no_overlap_reads_per_chrs,
	vector<t_annot_region*>* regions_2_exclude);

void exclude_reads_per_regions_per_chr(char* read_chr_id, 
	vector<t_mapped_read*>* cur_chr_reads,
	vector<t_mapped_read*>* no_overlap_reads,
	vector<t_annot_region*>* regions_2_exclude);

#endif // __MAPPED_READ_FILE_INTERFACE__

