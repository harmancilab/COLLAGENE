#ifndef __SIGNAL_TRACK_FILE_INTERFACE__
#define __SIGNAL_TRACK_FILE_INTERFACE__

#include <vector>

using namespace std;

struct t_annot_region;
class t_rng;

double* load_signal_covg_per_directory_chr_id(const char* dat_dir,
	char* chr_id,
	int l_fragment,
	char* cur_dat_fp,
	int& l_loaded_covg,
	bool& reads_loaded);

double* load_signal_covg_per_signal_file(const char* cur_dat_fp,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded);

double* bin_profile(double* signal_profile, int l_profile,
					vector<t_annot_region*>* target_regions,
					vector<t_annot_region*>* quantified_regions,
					bool flag_span_region_borders,	// This says whether we would like to go across the ends of consecutive windows while enumerating them.
					int l_bin,
					int& n_bins);

double* invert_1based_signal(double* signal, int l_signal);

double* concat_profile2_2_profile1(double* prof1, int l_prof1, 
									double* prof2, int l_prof2, 
									int& l_new_track);

double** generate_multitrack_signal_profiles_per_config(char* config_fp, bool flag_span_target_region_ends, vector<char*>* chr_ids_2_use, int& n_gw_tracks, int& l_gw_track, vector<char*>* track_ids);

double** generate_multitrack_signal_profiles(vector<char*>* preprocessed_reads_dirs,
											vector<char*>* bin_track_dirs,
											vector<char*>* region_fps,
											vector<char*>* chr_ids_2_use, 											
											vector<t_annot_region*>* target_regions, // null means genomewide signal; i.e., uses the length of the first preprocessed reads track.
											vector<t_annot_region*>* quantified_regions, // This is the list of regions for which signals are generated. The size must be the same as l_track.
											bool flag_span_target_region_ends,
											int l_bin, // <=0 for no binning and per target signal.
											double log_base, // 0 means linear signal return.
											bool use_random_track, // Adds a random track to the end if this is true.
											int& l_track);	// Length of the track.

double* load_per_nucleotide_binary_profile(char* binary_per_nucleotide_profile_fp, int& l_profile);
void dump_per_nucleotide_binary_profile_per_bedgraph(const char* bgr_fp, bool dump_binary, const char* op_dir);
void dump_per_nucleotide_binary_profile(double* profile, int l_profile, const char* op_fp);
void dump_bedGraph_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* op_fp);
void dump_bBGR_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* bbgr_op_fp);
//void convert_BGR_2_bBGR(char* bgr_fp, char* bbgr_op_fp);
void dump_bBGR_per_per_bedGraph(char* bgr_fp, char* bbgr_op_fp);
void dump_bedGraph_per_bBGR(char* bbgr_fp, char* chrom, char* op_fp);
void dump_bedGraph_per_bBGR_v1(char* bbgr_fp, char* chrom_id_2_use, char* bgr_op_fp);

struct t_signal_node
{
	double signal;
	int i_reg;
};

bool sort_signal_nodes_per_increasing_signal(t_signal_node* node1, t_signal_node* node2);


int* load_int_signal_covg_per_directory_chr_id(const char* dat_dir, char* chr_id, int& l_loaded_covg);
void remove_duplicate_samples_from_signal_matrix(char* signal_regs_bed_fp, char* op_signal_regs_bed_fp);
void extract_signal_features_per_regions(char* argv[], int argc, char* config_fp, char* op_fp);
void quantile_normalize_signal_matrix_4th_col_signals(char* signal_regions_BED_fp, char* op_fp);
void rank_normalize_signal_matrix(char* config_fp, char* signal_matrix_fp, char* op_fp);
void extract_signal_regions_per_sample_list_per_4th_col_signals(char* signal_reg_BED_fp, char* selected_sample_ids_list_fp, char* op_fp);
void reheader_signal_regions_BED(char* signal_regions_BED_fp, char* sample_ids_fp, char* op_fp);
void concatenate_signal_regions(vector<char*>* signal_reg_BED_fps, char* op_fp);
void generate_per_region_sample_summarized_signal_stats(vector<char*>* signal_reg_BED_fps, char* op_fp);
vector<t_annot_region*>* load_signal_regs_BED(char* signal_regions_BED_fp, int& n_loaded_samples);
vector<t_annot_region*>* load_signal_regs_BED(char* signal_regions_BED_fp, int& n_loaded_samples, int matrix_starting_0_based_col_i);

double* load_per_nucleotide_bBGR_track(char* bgr_fp, int& l_profile);
double* load_per_nucleotide_bBGR_v1_track(char* bgr_fp, int& l_profile);

double* load_per_nucleotide_BGR_track(const char* bgr_fp, int& l_profile);

double* get_zero_indexed_per_one_indexed_data(double* one_indexed_data, int l_profile);
double* get_one_indexed_per_zero_indexed_data(double* zero_indexed_data, int l_profile);

double* quantize_per_nucleotide_profiles(double* profile, int l_profile, vector<double>* thresholds, vector<double>* quantized_vals);

double** joint_prune_multiple_profiles_per_thresholds(double** profiles, int n_profiles, int l_profile, double* threshold_per_profile, int& l_pruned_profile);

unsigned char* load_per_nucleotide_binary_uchar_profile(char* binary_per_nucleotide_profile_fp, int& l_profile);
void dump_per_nucleotide_uchar_binary_profile(unsigned char* signal_profile, int l_profile, char* op_fp);

double* get_block_permute_profile(double* profile_data, t_rng* rng, int l_profile, int& l_permuted, int l_block);
double* get_random_circular_shifted_profile(double* profile_data, t_rng* rng, int l_profile);

void get_profile_extrema(double* profile, int l_profile, double& prof_min, double& prof_max);

double* extract_one_indexed_profile_per_profile_stranded(double* signal_profile_buffer, int l_profile, int start, int end, char strand, int& l_extracted_profile);
double* extract_one_indexed_profile_per_profile(double* signal_profile_buffer, int l_profile, int start, int end, int& l_extracted_profile);
double* copy_profile(double* signal_profile, int l_profile);

vector<t_annot_region*>** get_peaks_per_per_nucleotide_signal_profile(double* signal_profile, char* chrom, int l_data, double min_thresh, double max_thresh);
vector<t_annot_region*>** get_valleys_per_per_nucleotide_signal_profile(double* signal_profile, char* chrom, int l_data, double min_thresh, double max_thresh);

void exclude_regions_from_signal_profiles(double* signal_profile, int l_profile, vector<t_annot_region*>* regions_2_exclude, double* pruned_signal_profile, int& l_pruned_profile);

void floorize_profile(double* signal_profile, int l_profile);
void get_log_plus_one_profile(double* signal_profile, double base, int l_profile);
double* get_zero_profile(int l_profile);

#endif // __SIGNAL_TRACK_FILE_INTERFACE__

