#ifndef __HUMAN_DATA_PROCESSING__
#define __HUMAN_DATA_PROCESSING__

#include <vector>
using namespace std;

struct t_transcript_info;
class t_config;
struct t_annot_region;

void parse_GENCODE_gff_grp_strs(vector<t_annot_region*>* gff_entry);

// Load the human trancripts, then restructure them to get the transcript information and the exons.
vector<t_transcript_info*>* get_human_transcript_info_per_GFF(t_config* config);

void extract_human_transcript_sequences(char* human_genome_dir,
	vector<t_transcript_info*>* transcript_info_list);

bool compare_human_GENCODE_gene_ids(char* gene_id1, char* gene_id2);

#endif // __HUMAN_DATA_PROCESSING__

