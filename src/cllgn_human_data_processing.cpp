#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cllgn_human_data_processing.h"
#include "cllgn_ansi_string.h"
#include "cllgn_annot_region_tools.h"
#include "cllgn_gff_utils.h"

//char* get_human_gff_transcript_id(t_annot_region* gff_entry);
//vector<char*>* get_human_gff_parent_transcript_ids(t_annot_region* gff_entry);
//
//// Load the human trancripts, then restructure them to get the transcript information and the exons.
//vector<t_transcript_info*>* get_human_transcript_info_per_GFF(t_config* config);
//
//void extract_human_transcript_sequences(char* human_genome_dir,
//	vector<t_transcript_info*>* transcript_info_list);

bool compare_human_GENCODE_gene_ids(char* gene_id1, char* gene_id2)
{
	if(strcmp(gene_id1, gene_id2) == 0)
	{
		return(true);
	}

	int next_char1_i;
	int next_char2_i;
	vector<t_string*>* id1_f_token = t_string::get_first_n_tokens(gene_id1, 1, ".", next_char1_i);
	vector<t_string*>* id2_f_token = t_string::get_first_n_tokens(gene_id2, 1, ".", next_char2_i);

	if(t_string::compare_strings(id1_f_token->at(0), id2_f_token->at(0)))
	{
		return(true);
	}

	t_string::clean_tokens(id1_f_token);
	t_string::clean_tokens(id2_f_token);

	return(false);
}

void parse_GENCODE_gff_grp_strs(vector<t_annot_region*>* all_regions)
{
	for(int i_r = 0; i_r < (int)all_regions->size(); i_r++)
	{
		if(i_r % 1000 == 0)
		{
			fprintf(stderr, "Parsing grp string: %d (%d)      \r", i_r, (int)all_regions->size());
		}

		t_gff_info* cur_gff_info = (t_gff_info*)(all_regions->at(i_r)->data);
		char* grp_str = cur_gff_info->group_str;

		// Allocate the property id's and values for this entry.
		cur_gff_info->prop_ids = new vector<char*>();
		cur_gff_info->prop_vals = new vector<char*>();

		// Parse the parent group string.
		t_string* grp_str_str = new t_string(grp_str);
		
		// The group string has semi-colon separated entries, each entry is separated by '=' or ' '.
		t_string_tokens* grp_str_tokens = grp_str_str->tokenize_by_chars(";");

		for(int i_tok = 0; i_tok < (int)grp_str_tokens->size(); i_tok++)
		{
			t_string_tokens* cur_grp_str_entry_tokens = grp_str_tokens->at(i_tok)->tokenize_by_chars("= \"");

			// chr1    HAVANA  CDS     874420  874509  .       +       2       ID=CDS:ENST00000420190.1:6;Parent=ENST00000420190.1;gene_id=ENSG00000187634.6;transcript_id=ENST00000420190.1;gene_type=protein_coding;gene_status=KNOWN;gene_name=SAMD11;transcript_type=protein_coding;transcript_status=KNOWN;transcript_name=SAMD11-011;exon_number=6;exon_id=ENSE00002686739.1;level=2;protein_id=ENSP00000411579.1;havana_gene=OTTHUMG00000040719.8;havana_transcript=OTTHUMT00000316521.1;tag=alternative_5_UTR,mRNA_end_NF,cds_end_NF
			if(cur_grp_str_entry_tokens->size() == 2)
			{
				cur_gff_info->prop_ids->push_back(t_string::copy_me_str(cur_grp_str_entry_tokens->at(0)->str()));
				cur_gff_info->prop_vals->push_back(t_string::copy_me_str(cur_grp_str_entry_tokens->at(1)->str()));
			}
			else
			{
			}

			// Free group string token memory.
			t_string::clean_tokens(cur_grp_str_entry_tokens);
		} // i_tok loop.

		//t_string_tokens* grp_str_tokens = grp_str_str->tokenize_by_chars(";\" ");

		//int i = 0; 
		//while(i < grp_str_tokens->size())
		//{
		//	if(i+1 >= grp_str_tokens->size())
		//	{
		//		fprintf(stderr, "Could not parse the GENCODE gff group string: %s:%d-%d\n%s\n", 
		//			all_regions->at(i_r)->chrom, all_regions->at(i_r)->start, all_regions->at(i_r)->end, grp_str);
		//		exit(0);
		//	}

		//	//fprintf(stderr, "Prop %s: %s\n", grp_str_tokens->at(i)->str(), grp_str_tokens->at(i+1)->str());

		//	cur_gff_info->prop_ids->push_back(t_string::copy_me_str(grp_str_tokens->at(i)->str()));
		//	i++;
		//	cur_gff_info->prop_vals->push_back(t_string::copy_me_str(grp_str_tokens->at(i)->str()));
		//	i++;
		//} // i loop.

		t_string::clean_tokens(grp_str_tokens);
		delete grp_str_str;
	} // i_r loop.
}

