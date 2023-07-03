#include <stdio.h>
#include <stdlib.h>
#include "cllgn_mapped_read_tools.h"
#include <vector>
#include <ctype.h>
#include <math.h>
#include "cllgn_signal_track_tools.h"
#include "cllgn_annot_region_tools.h"
#include "cllgn_genome_sequence_tools.h"
#include "cllgn_variation_tools.h"
#include "cllgn_genomics_coords.h"
#include "cllgn_file_utils.h"
#include "cllgn_nomenclature.h"
#include "cllgn_rng.h"
#include "cllgn_seed_manager.h"
#include "cllgn_nucleotide.h"
#include <string.h>
#include <algorithm>
#include "cllgn_ansi_string.h"
#include "cllgn_annot_region_tools.h"

using namespace std; 

bool __DUMP_MAPPED_READ_TOOLS_MSGS__ = false;

FILE* get_processed_reads_ptr_wrapper(char* cur_chr_reads_fp)
{
	FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
//	if (t_string::compare_strings(cur_chr_reads_fp, "stdin"))
//	{
//		f_cur_chr_reads = stdin;
//	}
//	else if (t_string::ends_with(cur_chr_reads_fp, ".gz"))
//	{
//		char ungzip_cmd[1000];
//		sprintf(ungzip_cmd, "gzip -cd %s", cur_chr_reads_fp);
//#ifdef _WIN32
//		f_cur_chr_reads = _popen(ungzip_cmd, "r");
//#else 
//		f_cur_chr_reads = popen(ungzip_cmd, "r");
//#endif
//	}
//	else
//	{
//		f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
//	}

	return f_cur_chr_reads;
}

void close_processed_reads_ptr_wrapper(FILE* f_cur_chr_reads, char* cur_chr_reads_fp)
{
	if (t_string::compare_strings(cur_chr_reads_fp, "stdin"))
	{
	}
	else if (t_string::ends_with(cur_chr_reads_fp, ".gz"))
	{
#ifdef _WIN32
		_pclose(f_cur_chr_reads);
#else 
		pclose(f_cur_chr_reads);
#endif
	}
	else
	{
		fclose(f_cur_chr_reads);
	}
}

#define __UCHAR_MAPPABILITY__

bool sort_read_lines(char* read1, char* read2)
{
	return(t_string::sort_strings(read1, read2));
}

#define MAX(x, y) (((x)>(y))?(x):(y))
#define MIN(x, y) (((x)<(y))?(x):(y))

#define L_CHROM (250*1000*1000)

vector<int>* get_chromosome_lengths_per_mapped_reads(char* mapped_reads_dir)
{
	vector<int>* chr_lengths = new vector<int>();

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_cur_chr = 0;
		char cur_line[1000];
		char cur_mapped_reads_fp[1000];
		sprintf(cur_mapped_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, chr_ids->at(i_chr));
		FILE* f_mapped_reads = get_processed_reads_ptr_wrapper(cur_mapped_reads_fp);

		while(1)
		{
			if(fgets(cur_line, 1000, f_mapped_reads) == NULL)
			{
				break;
			}

			int cur_pos = 0;
			sscanf(cur_line, "%*s %*s %d", &cur_pos);

			if(cur_pos > l_cur_chr)
			{
				l_cur_chr = cur_pos + 1000;
			}
		} // file reading loop.
		//fclose(f_mapped_reads);
		close_processed_reads_ptr_wrapper(f_mapped_reads, cur_mapped_reads_fp);

		chr_lengths->push_back(l_cur_chr);
	} // i_chr loop.

	return(chr_lengths);
}

vector<char*>* sort_bucket_read_lines(char* bucket_fp)
{
	// Load the reads.
	vector<char*>* bucket_read_lines = buffer_file(bucket_fp);	
	vector<int>* read_starts = new vector<int>();
	vector<t_read_line_sorting_info*>* sorting_info_list = new vector<t_read_line_sorting_info*>();
	for(int i_read = 0; i_read < (int)bucket_read_lines->size(); i_read++)
	{
		int cur_read_start = 0;
		sscanf(bucket_read_lines->at(i_read), "%*s %*s %d", &cur_read_start);

		t_read_line_sorting_info* cur_line_info = new t_read_line_sorting_info();
		cur_line_info->start = cur_read_start;
		cur_line_info->read_line = bucket_read_lines->at(i_read);

		sorting_info_list->push_back(cur_line_info);
	} // i_read loop.

	sort(sorting_info_list->begin(), sorting_info_list->end(), sort_read_line_info);
	vector<char*>* sorted_bucket_read_lines = new vector<char*>();

	for(int i_read = 0; i_read < (int)sorting_info_list->size(); i_read++)
	{
		sorted_bucket_read_lines->push_back(sorting_info_list->at(i_read)->read_line);

		delete sorting_info_list->at(i_read);
	} // i_read loop.

	delete(sorting_info_list);
	delete(read_starts);

	delete(bucket_read_lines);

	return(sorted_bucket_read_lines);
}

// Following is for sorting the mapped reads offline.
bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2)
{
	return(info1->start < info2->start);
}

// This function open the last file; optimal for sorted data.
void preprocess_mapped_reads_file_single_file_buffering(char* mrf_fp, char* parsed_reads_op_dir, 
	void (preprocess_mapped_read_line)(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str),
	bool dump_read_id)
{
	FILE* f_mrf = open_f(mrf_fp, "r");
	if (f_mrf == NULL)
	{
		fprintf(stderr, "Could not open %s\n", mrf_fp);
		return;
	}

	//char cur_line[100000];
	int n_frags = 0;
	//int n_total_frags = 0;

	FILE* cur_preprocessed_read_f_ptr = NULL;
	char* cur_preprocessed_read_fp = NULL;
	char* cur_preprocessed_chr = NULL;
	vector<char*>* chr_ids = new vector<char*>();

	while (1)
	{
		char* cur_line = getline(f_mrf);
		if (cur_line == NULL)
		{
			break;
		}

		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char;
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);

		// Make sure that the line is valid.
		if (chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if (i_chr == (int)chr_ids->size())
			{
				// Add the chromosome id.
				chr_ids->push_back(t_string::copy_me_str(chrom));

				fprintf(stderr, "Added %s\n", chrom);
			}

			// Select the file pointer.
			if (cur_preprocessed_chr == NULL ||
				(cur_preprocessed_chr != NULL && !t_string::compare_strings(cur_preprocessed_chr, chrom)))
			{
				// Copy the current chromosome.
				if (cur_preprocessed_chr != NULL)
				{
					delete[] cur_preprocessed_chr;
				}

				cur_preprocessed_chr = t_string::copy_me_str(chrom);

				// If the last file is open, close it.
				if (cur_preprocessed_read_f_ptr != NULL)
				{
					// Close the file.
					fprintf(stderr, "Closing %s\n", cur_preprocessed_read_fp);
					close_f(cur_preprocessed_read_f_ptr, cur_preprocessed_read_fp);
				}

				// Re-open the file.
				if (cur_preprocessed_read_fp == NULL)
				{
					cur_preprocessed_read_fp = new char[1000];
				}

				sprintf(cur_preprocessed_read_fp, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chrom);
				cur_preprocessed_read_f_ptr = open_f(cur_preprocessed_read_fp, "a");
			}

			FILE* cur_frag_file = cur_preprocessed_read_f_ptr;

			if (cur_frag_file == NULL)
			{
				//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
			}
			else
			{
				if (dump_read_id)
				{
					fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
				}
				else
				{
					fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
				}
				n_frags++;
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete[] cur_line;
	} // file reading loop.

	// Unload/close the mapped read file.
	close_f(f_mrf, mrf_fp);

	if (cur_preprocessed_read_f_ptr != NULL)
	{
		close_f(cur_preprocessed_read_f_ptr, cur_preprocessed_read_fp);
	}

	// (Re)Dump the chromosome id list.
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for (int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);
}

// Generic preprocessing function for mapped read files.
void preprocess_mapped_reads_file(char* mrf_fp, 
	char* parsed_reads_op_dir, 
	vector<char*>* preset_chr_ids,
	void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	bool dump_read_id)
{
    FILE* f_mrf = open_f(mrf_fp, "r");
	if(f_mrf == NULL)
	{
		fprintf(stderr, "Could not open %s\n", mrf_fp);
		return;
	}

    //char cur_line[100000];
    int n_frags = 0;
    //int n_total_frags = 0;

    vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
	vector<char*>* frag_fps = new vector<char*>();

	// Check chromosome id's list file.
	char chr_ids_fp[100000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);

	vector<char*>* chr_ids = NULL;
	if (preset_chr_ids == NULL)
	{
		if (check_file(chr_ids_fp))
		{
			chr_ids = buffer_file(chr_ids_fp);

			fprintf(stderr, "Found chromosome id's @ %s, pooling.\n", chr_ids_fp);

			// Open the files for appending.
			for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
			{
				char new_fn[1000];
				sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chr_ids->at(i_chr));
				frag_f_ptrs->push_back(open_f(new_fn, "a"));
				frag_fps->push_back(t_string::copy_me_str(new_fn));
			} // i_chr loop.
		}
		else
		{
			// The chromosomes will be added now.
			chr_ids = new vector<char*>();
		}
	} // preset_chr_ids check.
	else
	{
		fprintf(stderr, "Using preset %d chromosomes.\n", (int)preset_chr_ids->size());

		// Preset chromosomes exist.
		chr_ids = new vector<char*>();
		for (int i_chr = 0; i_chr < (int)preset_chr_ids->size(); i_chr++)
		{
			char* cur_chr_id = t_string::copy_me_str(preset_chr_ids->at(i_chr));
			normalize_chr_id(cur_chr_id);

			chr_ids->push_back(cur_chr_id);

			char new_fn[1000];
			sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, cur_chr_id);

			fprintf(stderr, "Opening %s for writing.\n", new_fn);
			frag_f_ptrs->push_back(open_f(new_fn, "a"));
			frag_fps->push_back(t_string::copy_me_str(new_fn));
		} // i_chr loop.
	} // preset chromosome check.

	if (frag_fps->size() != frag_f_ptrs->size())
	{
		fprintf(stderr, "Sanity check failed: Number of files do not match to pointers.\n");
		exit(0);
	}

	while(1)
	{
		char* cur_line = getline(f_mrf);
		if(cur_line == NULL)
		{
			break;
		}
		
		n_frags++;

		if (n_frags % (1000 * 1000) == 0)
		{
			fprintf(stderr, "@ %d reads.                  \r", n_frags);
		}

		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char; 
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
									read_id,
									chrom, 
									chr_index, sequenced_length, 
									strand_char, 
									mapping_quality_str);

		// Make sure that the line is valid.
		if(chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			bool chromosome_preset_confirm = true;
			if(i_chr == (int)chr_ids->size())
			{
				// If there were no preset chromosome id's, update the list.
				if (preset_chr_ids == NULL)
				{
					// Add the chromosome id.
					chr_ids->push_back(t_string::copy_me_str(chrom));
					i_chr = t_string::get_i_str(chr_ids, chrom);

					char new_fn[10000];
					sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chrom);

					// Re-open this file from scratch; even if it was there, it will be overwritten.
					frag_f_ptrs->push_back(open_f(new_fn, "w"));
					frag_fps->push_back(t_string::copy_me_str(new_fn));

					fprintf(stderr, "Added %s\n", chrom);
				}
				else
				{
					// Skip adding this chromosome; it is not in the preset list.
					chromosome_preset_confirm = false;

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Skipping processing %s               \n", cur_line);
					}
				}
			}

			if (chromosome_preset_confirm)
			{
				FILE* cur_frag_file = frag_f_ptrs->at(i_chr);

				if (cur_frag_file == NULL)
				{
					//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
				}
				else
				{
					if (dump_read_id)
					{
						fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
					}
					else
					{
						fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
					}
				}
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete [] cur_line;
	} // file reading loop.

	// (Re)Dump the chromosome id list.
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for(int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);

	// Close fragment file pointers.
	for(int i_f = 0; i_f < (int)frag_f_ptrs->size(); i_f++)
	{
		close_f(frag_f_ptrs->at(i_f), frag_fps->at(i_f));

		// Compress, if it is already not compressed.
		if (!t_string::ends_with(frag_fps->at(i_f), "gz"))
		{
			char comp_frag_fp[1000];
			sprintf(comp_frag_fp, "%s.gz", frag_fps->at(i_f));
			fprintf(stderr, "Compressing to %s\n", comp_frag_fp);
			compressFile(frag_fps->at(i_f), comp_frag_fp);
			delete_file(frag_fps->at(i_f));
		}

		delete[] frag_fps->at(i_f);
	} // i_f loop.

	delete frag_fps;
	delete frag_f_ptrs;

	// Unload/close the mapped read file.	
	close_f(f_mrf, mrf_fp);
}

void count_mapped_reads_per_file(char* mrf_fp, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	double n_total_reads)
{
    // Divide SAM output with respect to chromosomes.
	FILE* f_mrf = open_f(mrf_fp, "r");

	if(f_mrf == NULL)
	{
		fprintf(stderr, "mapped read file pointer and file buffer are both NULL for %s\n", mrf_fp);
		return;
	}

    //char cur_line[2000];
    //int n_frags = 0;
    //int n_total_frags = 0;
	n_total_reads = 0;

	fprintf(stderr, "Counting the mapped reads from %s\n", mrf_fp);
	while(1)
	{
		//char* cur_line = getline(f_mrf);
		char* cur_line = getline(f_mrf);

		if(cur_line == NULL)
		{
			break;
		}
		
		// Load the mapping info based on the file type.
		//char chrom[1000];
		//int chr_index;
		//int sequenced_length;
		//char strand_char; 
		//char mapping_quality_str[1000];
		//preprocess_mapped_read_line(cur_line, 
		//							chrom, 
		//							chr_index, sequenced_length, 
		//							strand_char, 
		//							mapping_quality_str);
		char phred_quality_str[1000];
		int flag;
		if(cur_line[0] == '@')
		{
		}
		else
		{
			flag = 0;
			//if(sscanf(cur_line, "%*s %*d %*s %*d %*d %*s %*s %*d %*d %*s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 1)
			if(sscanf(cur_line, "%*s %d %*s %*d %*d %*s %*s %*d %*d %*s %s", &flag, phred_quality_str) == 2)
			{
				// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
				//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
				//chr_index = translate_coord(chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

				// Check the flag and determine the strand.
				/*strand_char = 'F';
				if(flag & 0x10)
				{
					strand_char = 'R';
				}*/

				// Sanity check. Is this fragment mapped?
				if(flag & 0x04)
				{
					// The read is not mapping.
					//chrom[0] = 0;
				}
				else
				{
					n_total_reads++;
				}
			}
			else
			{
				// Could not parse the line.
			}
		}

		if( (((int)n_total_reads) % 1000000) == 0)
		{
			fprintf(stderr, "At %lf. read.\n", n_total_reads);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "%lf reads.\n", n_total_reads);

	// Unload/close the mapped read file.
	close_f(f_mrf, mrf_fp);
}

void preprocess_tagAlign_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	int chr_start_index;
	int chr_end_index;
	char strand_sign;

	if(sscanf(cur_line, "%s %d %d %*s %*d %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
	{
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		//chr_start_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		//chr_end_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		chr_start_index = translate_coord(chr_start_index, TAGALIGN_COORDS::start_base, CODEBASE_COORDS::start_base);
		chr_end_index = translate_coord(chr_end_index, TAGALIGN_COORDS::end_base, CODEBASE_COORDS::end_base);

		// Set quality to all matches.
		sprintf(cigar_str, "%dM", chr_end_index-chr_start_index+1);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = chr_end_index-chr_start_index+1;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_preprocessed_LH_GFF3_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int _chr_index;
	char strand;

	//if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 6)
	// X       14705460        35M     +
	if(sscanf(cur_line, "%s %d %s %c", chrom, &_chr_index, mapping_quality_str, &strand) == 4)
	{
		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		chr_index = _chr_index;

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand == '-')
		{
			strand_char = 'R';
		}

			chr_index = _chr_index;
			sequenced_length = 0;
	}
	else
	{
		chrom[0] = 0;
	}
}

// TODO: Make an option to load only a certain phed partition from the file, so that we conserve memory.
unsigned short*** load_partitioned_compressed_pileups(char* cur_comp_allele_fp, int& n_partitions, int& l_pileup)
{
	FILE* f_comp = open_f(cur_comp_allele_fp, "rb");

	// Load the partitions.
	int n_read_partitions = 0;
	fread(&n_read_partitions, sizeof(int), 1, f_comp);
	n_partitions = n_read_partitions;

	// Load the partitions.
	int l_read_sig = 0;
	fread(&l_read_sig, sizeof(int), 1, f_comp);
	int l_sig = l_read_sig;

	l_pileup = l_sig;

	fprintf(stderr, "Loading pileup of length %d over %d phred partitions\n", l_sig, n_partitions);

	unsigned short*** loaded_pileup = new unsigned short**[n_partitions];

	for(int part_i = 0; part_i < n_partitions; part_i++)
	{
		loaded_pileup[part_i] = allocate_pileup(l_sig);
	} // part_i loop.

	// Go over all the positions.
	int i = 1;
	while(i <= l_sig)
	{
		if(i % 1000000 == 0)
		{
			fprintf(stderr, "Loading %d            \r", i);
		}

		// Read the existence flag.
		unsigned char existence_flag = 0;
		fread(&existence_flag, sizeof(unsigned char), 1, f_comp);

		if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
		{
			fprintf(stderr, "Loading %d. position, existence flag: 0x%x\n", i, existence_flag);
		}

		// Check an RLE case.
		if(existence_flag == 0xFF)
		{
			unsigned int l_RLE = 0;
			fread(&l_RLE, sizeof(unsigned int), 1, f_comp);

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Loading RLE of length %d @ %d\n", l_RLE, i);
			}

			// When we add this, we move upto the location where a 0-run ends.
			i += l_RLE;
		} // RLE check.
		else
		{
			// Read the first byte, parse the data range and which alleles exist.
			unsigned char drange_flag = ((existence_flag & (1 << 5)) >> 5);

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Loading val @ %d: drange_flag: %d\n", i, drange_flag);
			}

			for(int part_i = 0; part_i < n_partitions; part_i++)
			{
				for(int allele_i = 0; allele_i < 5; allele_i++)
				{
					// Do we have an entry here?
					if((existence_flag & (1 << allele_i)) > 0)
					{
						if(drange_flag == 1)
						{	
							unsigned short current_count_short = 0;
							fread(&(current_count_short), sizeof(unsigned short), 1, f_comp);
							loaded_pileup[part_i][allele_i][i] = current_count_short;
						}
						else
						{
							unsigned char current_count_char = 0;
							fread(&(current_count_char), sizeof(unsigned char), 1, f_comp);
							loaded_pileup[part_i][allele_i][i] = (unsigned short)current_count_char;
						}

						//fprintf(stderr, "%d, ", loaded_pileup[part_i][allele_i][i]);
					} // count check for current position.
				} // allele_i loop.
			} // part_i loop.

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "\n");
			}

			// Update the position.
			i++;
		} // non-RLE check.
	} // i loop.
	close_f(f_comp, cur_comp_allele_fp);

	// Free memory.
	return(loaded_pileup);
}

unsigned short** allocate_pileup(int l_sig)
{
	unsigned short** loaded_pileup = new unsigned short*[5];
	for (int i_allele = 0; i_allele < 5; i_allele++)
	{
		loaded_pileup[i_allele] = new unsigned short[l_sig + 2];
		memset(loaded_pileup[i_allele], 0, sizeof(unsigned short) * (l_sig + 2));
	} // i_allele loop.

	return loaded_pileup;
}

void delete_pileup(unsigned short** loaded_pileup)
{
	//unsigned short** loaded_pileup = new unsigned short*[5];
	for (int i_allele = 0; i_allele < 5; i_allele++)
	{
		delete [] loaded_pileup[i_allele];
	} // i_allele loop.

	delete[] loaded_pileup;
}

// This code should be eventually combined with the pileup generation code.
void extract_summarize_indel_containing_read_blocks_per_SAM(char* SAM_fp, char* chrom_info_fp, char* genome_seq_dir, char* op_dir)
{
	fprintf(stderr, "Summarizing and extracting indel supporting reads from %s and saving to %s.\n", SAM_fp, op_dir);

	// Load the chromosome information.
	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chrom_info_fp, chr_ids, chr_lengths);

	vector<t_annot_region*>** per_chr_indel_containing_read_blocks = new vector<t_annot_region*>*[(int)chr_ids->size() + 2];
	char** per_chrom_seq = new char*[(int)chr_ids->size() + 2];
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		per_chr_indel_containing_read_blocks[i_chr] = new vector<t_annot_region*>();

		int l_loaded_seq = 0;
		char cur_chrom_seq_fp[1000];
		sprintf(cur_chrom_seq_fp, "%s/%s.bin", genome_seq_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chrom_seq_fp))
		{
			sprintf(cur_chrom_seq_fp, "%s/%s.bin.gz", genome_seq_dir, chr_ids->at(i_chr));

			if (!check_file(cur_chrom_seq_fp))
			{
				fprintf(stderr, "Could not find the sequence file %s\n", cur_chrom_seq_fp);
				exit(0);
			}
		}

		per_chrom_seq[i_chr] = load_binary_sequence_file(cur_chrom_seq_fp, l_loaded_seq);
		fprintf(stderr, "Loaded %s (%d, %d)                   \r", chr_ids->at(i_chr), chr_lengths->at(i_chr), l_loaded_seq);
	} // i_chr loop.

	// Enter file reading loop.
	char cur_fragment[1000];
	char phred_quality_str[100000];

	unsigned long long n_total_processed_reads = 0;
	unsigned long long n_indel_containing_reads = 0;
	unsigned long long n_unmapped_reads = 0;

	double n_matching_matching_nucs = 0;
	double n_total_matching_nucs = 0;

	FILE* f_sam = open_f(SAM_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_sam);
		if (cur_line == NULL)
		{
			break;
		}

		// This is the global parameter, update this before checking for thread outof/which check.
		n_total_processed_reads++;

		// Report only with the lead thread.
		if (n_total_processed_reads % (1000 * 1000) == 0)
		{
			fprintf(stderr, "@ %llu. read: %llu indel containing; %llu unmapped; Ref match/unmatch million nucs: %.1f/%.1f.               \r",
					n_total_processed_reads,
					n_indel_containing_reads,
					n_unmapped_reads,
					n_matching_matching_nucs / (1000 * 1000),
					n_total_matching_nucs / (1000 * 1000));
		}

		// If this is a comment line, skip it.
		if (cur_line[0] == '@')
		{
			delete[] cur_line;
			continue;
		}

		char read_id[1000];
		char chrom[100];
		int chr_index = 0;
		int sequenced_length;
		char strand_char = 0;
		char cigar_str[100];

		preprocess_SAM_read_line(cur_line,
								read_id,
								chrom,
								chr_index, sequenced_length,
								strand_char,
								cigar_str);

		// Copy the cigar string.
		char* mapping_map_str = cigar_str;

		// Make sure that the line is valid.
		if (chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int chr_i = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if (chr_i == (int)chr_ids->size())
			{
				// This read is on a chromosome we do not care about. Count it then throw away?
			}
			else
			{
				// Load the sequence and quality data.
				if (sscanf(cur_line, "%*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", 
							cur_fragment, phred_quality_str) == 2)
				{

				}
				else
				{
					fprintf(stderr, "Could not parse %s\n", cur_line);
					exit(0);
				}

				char strand_plus_min = '+';
				if (strand_char == 'R')
				{
					strand_plus_min = '-';
				}

				// Parse the cigar string.
				int i_mapp_map = 0;
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				while (mapping_map_str_valid &&
						mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map,
														is_matching,
														l_cur_entry,
														entry_type_char);

					// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
					if (is_matching)
					{
						int cur_read_i = read_nuc_index;
						for (int cur_genome_i = chr_index;
							cur_genome_i <= chr_index + l_cur_entry - 1;
							cur_genome_i++)
						{
							if (cur_genome_i <= chr_lengths->at(chr_i))
							{
								if (nuc_2_num(cur_fragment[cur_read_i]) == nuc_2_num(per_chrom_seq[chr_i][cur_genome_i]))
								{
									n_matching_matching_nucs++;
								}
								else
								{
									n_total_matching_nucs++;
								}
							}
							else
							{
								fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
									chr_ids->at(chr_i), cur_genome_i, chr_lengths->at(chr_i));

								exit(0);
							}

							cur_read_i++;
						} // cur_genome_i loop.
					} // check if this block is matching.
					else if (entry_type_char == 'D')
					{
						n_indel_containing_reads++;

						t_annot_region* new_del_reg = get_empty_region();
						new_del_reg->chrom = t_string::copy_me_str(chrom);
						new_del_reg->start = translate_coord(chr_index, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_del_reg->end = translate_coord(chr_index + l_cur_entry - 1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_del_reg->strand = strand_plus_min;
						new_del_reg->data = NULL;
						new_del_reg->score = VAR_TYPE_DELETION;

						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						int l_ref_allele = (chr_index + l_cur_entry - 1) - chr_index + 1;
						//int l_alt_allele = 1;
						char* cur_del_ref_allele = new char[l_ref_allele + 5];
						memset(cur_del_ref_allele, 0, sizeof(char) * (l_ref_allele + 2));
						for (int cur_genome_i = chr_index; cur_genome_i <= chr_index + l_cur_entry - 1; cur_genome_i++)
						{
							if (cur_genome_i <= chr_lengths->at(chr_i))
							{
								cur_del_ref_allele[cur_genome_i - chr_index] = per_chrom_seq[chr_i][cur_genome_i];
							}
							else
							{
								fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
										chr_ids->at(chr_i), cur_genome_i, chr_lengths->at(chr_i));

								exit(0);
							}
						} // cur_genome_i loop.

						// Set the ref allele to the data for deletions.
						new_del_reg->data = cur_del_ref_allele;

						if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
						{
							fprintf(stderr, "DELETION @ %d (Ref. %s): %s\n", chr_index, cur_del_ref_allele, cur_line);
							getc(stdin);
						}

						// Add the del block region.
						per_chr_indel_containing_read_blocks[chr_i]->push_back(new_del_reg);
					} // deletion check.
					else if (entry_type_char == 'I')
					{
						n_indel_containing_reads++;

						// This points to the position right before the insertion.
						t_annot_region* new_ins_reg = get_empty_region();
						new_ins_reg->chrom = t_string::copy_me_str(chrom);
						new_ins_reg->start = translate_coord(chr_index-1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_ins_reg->end = translate_coord(chr_index, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_ins_reg->strand = strand_plus_min;
						new_ins_reg->data = NULL;
						new_ins_reg->score = VAR_TYPE_INSERTION;

						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						//int l_ref_allele = 1;
						int l_alt_allele = (chr_index + l_cur_entry - 1) - chr_index + 1;
						char* cur_ins_alt_allele = new char[l_alt_allele + 5];
						memset(cur_ins_alt_allele, 0, sizeof(char) * (l_alt_allele + 2));

						// Insertion to the reference: This is included to one position.
						for (int cur_read_i = read_nuc_index;
							cur_read_i <= read_nuc_index + l_cur_entry - 1;
							cur_read_i++)
						{
							cur_ins_alt_allele[cur_read_i - read_nuc_index] = cur_fragment[cur_read_i];
						} // cur_genome_i loop.
						// Set the ref allele to the data for deletions.
						new_ins_reg->data = cur_ins_alt_allele;

						if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
						{
							fprintf(stderr, "INSERTION @ %d (Alt %s): %s\n", chr_index, cur_ins_alt_allele, cur_line);
							getc(stdin);
						}

						// Add the current insertion region to the list of regions.
						per_chr_indel_containing_read_blocks[chr_i]->push_back(new_ins_reg);
					} // insert block check.

					// Update the base for the current entry.
					if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if (check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.
			} // chromosome check.
		} // chr_index check.
		else
		{
			n_unmapped_reads++;
		}

		delete[] cur_line;
	} // sam file reading loop.

	close_f(f_sam, SAM_fp);

	// Sort each chromosome, then save.
	fprintf(stderr, "Finished processing reads, saving to %s\n", op_dir);
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "%s: %d indel blocks.\n", chr_ids->at(i_chr), (int)per_chr_indel_containing_read_blocks[i_chr]->size());
		char cur_indel_blocks_fp[1000];
		sprintf(cur_indel_blocks_fp, "%s/%s_indel_blocks.bed.gz", op_dir, chr_ids->at(i_chr));

		FILE* f_indel_blocks = NULL;
		if (check_file(cur_indel_blocks_fp))
		{
			fprintf(stderr, "Opening %s for pooling.\n", cur_indel_blocks_fp);
			f_indel_blocks = open_f(cur_indel_blocks_fp, "a");
		}
		else
		{
			fprintf(stderr, "Opening %s for saving.\n", cur_indel_blocks_fp);
			f_indel_blocks = open_f(cur_indel_blocks_fp, "w");
		}

		for (int i_block = 0; i_block < (int)per_chr_indel_containing_read_blocks[i_chr]->size(); i_block++)
		{
			char empty_str[] = ".";
			char* ref_allele = empty_str;
			char* alt_allele = empty_str;

			// Check for insert/delete.
			if (per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->score == VAR_TYPE_INSERTION)
			{
				alt_allele = (char*)(per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->data);
			} // insertion check.
			else if(per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->score == VAR_TYPE_DELETION)
			{
				ref_allele = (char*)(per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->data);
			} // deletion check.
			else
			{
				fprintf(stderr, "Sanity check failed, the block size does not make sense: %s:%d-%d\n", per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->chrom,
						per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->start,
						per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->end);
			}

			fprintf(f_indel_blocks, "%s\t%d\t%d\t%s %s\t%d\t%c\n", 
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->chrom, 
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->start,
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->end,
				ref_allele,
				alt_allele,
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->score,
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->strand);
		} // i_block loop.
		close_f(f_indel_blocks, cur_indel_blocks_fp);
	} // i_chr loop.
}

bool sort_regions_coords_first_names_second(t_annot_region* reg1, t_annot_region* reg2)
{
	if (reg1->start != reg2->start)
	{
		return(reg1->start < reg2->start);
	}
	else if (reg1->end != reg2->end)
	{
		return(reg1->end < reg2->end);
	}
	else
	{
		// Both starts and ends are the same; check the alleles.
		return t_string::sort_strings(reg1->name, reg2->name);
	}
}


unsigned short** load_compressed_pileups(char* cur_comp_allele_fp, int& l_pileup)
{
	FILE* f_comp = open_f(cur_comp_allele_fp, "rb");

	if(f_comp == NULL)
	{
		fprintf(stderr, "Could not open %s\n", cur_comp_allele_fp);
		exit(0);
	}

	int l_read_sig = 0;
	fread(&l_read_sig, sizeof(int), 1, f_comp);
	int l_sig = l_read_sig;

	l_pileup = l_sig;

	fprintf(stderr, "Loading pileup of length %d\n", l_sig);

	unsigned short** loaded_pileup = allocate_pileup(l_sig);

	// Go over all the positions.
	int i = 1;
	while(i <= l_sig)
	{
		// Read the existence flag.
		unsigned char existence_flag = 0;
		int n_read = fread(&existence_flag, sizeof(unsigned char), 1, f_comp);

		// We have reached the EOF before reading whole file.
		if (n_read == 0)
		{
			break;
		}

		// Check an RLE case.
		if(existence_flag == 0xFF)
		{
			unsigned int l_RLE = 0;
			fread(&l_RLE, sizeof(unsigned int), 1, f_comp);

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Loading RLE of length %d @ %d\n", l_RLE, i);
			}

			// When we add this, we move upto the location where a 0-run ends.
			i += l_RLE;
		} // RLE check.
		else
		{
			// Read the first byte, parse the data range and which alleles exist.
			unsigned char drange_flag = ((existence_flag & (1 << 5)) >> 5);

			for(int allele_i = 0; allele_i < 5; allele_i++)
			{
				// Do we have an entry here?
				if((existence_flag & (1 << allele_i)) > 0)
				{
					if(drange_flag == 1)
					{	
						unsigned short current_count_short = 0;
						fread(&(current_count_short), sizeof(unsigned short), 1, f_comp);
						loaded_pileup[allele_i][i] = current_count_short;
					}
					else
					{
						unsigned char current_count_char = 0;
						fread(&(current_count_char), sizeof(unsigned char), 1, f_comp);
						loaded_pileup[allele_i][i] = (unsigned short)current_count_char;
					}
				} // count check for current position.
			} // allele_i loop.

			i++;
		}
	} // i loop.

	close_f(f_comp, cur_comp_allele_fp);

	// Free memory.
	return(loaded_pileup);
}

int* load_coverage_per_pileups(unsigned short** pileups, int l_sig)
{
	int* covg_signal = new int[l_sig+2];
	memset(covg_signal, 0, sizeof(int)*(l_sig+2));

	for(int i = 1; i <= l_sig; i++)
	{
		int cur_posn_covg = 0;
		for(int allele_i = 0; allele_i < 5; allele_i++)
		{
			cur_posn_covg += pileups[allele_i][i];
		} // allele_i loop.

		covg_signal[i] = cur_posn_covg;
	} // i loop.

	return(covg_signal);
}

int* load_coverage_per_compressed_partitioned_pileup_file(char* comp_pileup_fp, int& l_sig)
{
	fprintf(stderr, "Loading coverage from compressed pileup file %s.\n", comp_pileup_fp);

	// Load the compressed pileup.
	int n_parts = 0;
	unsigned short*** pileups = load_partitioned_compressed_pileups(comp_pileup_fp, n_parts, l_sig);

	int* covg_profile = new int[l_sig + 5];
	memset(covg_profile, 0, sizeof(int) * (l_sig+2));
	for(int i = 1; i <= l_sig; i++)
	{
		for(int part_i = 0; part_i < n_parts; part_i++)
		{
			for(int all_i = 0; all_i < 5; all_i++)
			{
				covg_profile[i] += pileups[part_i][all_i][i];
			} // all_i loop.
		} // part_i loop.
	} // i loop.

	// Free the pileups, since we do not need them anymore.
	for(int part_i = 0; part_i < n_parts; part_i++)
	{
		for(int all_i = 0; all_i < 5; all_i++)
		{
			delete [] pileups[part_i][all_i];
		} // all_i loop.

		delete [] pileups[part_i];
	} // part_i loop.

	delete [] pileups;

	return(covg_profile);
}

int* load_coverage_per_compressed_pileup_file(char* comp_pileup_fp, int& l_sig)
{
	fprintf(stderr, "Loading coverage from compressed pileup file %s.\n", comp_pileup_fp);

	// Load the compressed pileup.
	unsigned short** pileups = load_compressed_pileups(comp_pileup_fp, l_sig);

	// Get the coverage profile
	int* covg_profile = load_coverage_per_pileups(pileups, l_sig);

	// Free the pileups, since we do not need them anymore.
	for(int all_i = 0; all_i < 5; all_i++)
	{
		delete [] pileups[all_i];
	} // all_i loop.

	delete [] pileups;

	return(covg_profile);
}

void compress_partitioned_nucleotide_pileup_track(unsigned short*** pileup, int n_partitions, int l_sig, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	
	// Write the # of phred partitions
	fwrite(&n_partitions, sizeof(int), 1, f_op);

	// Write the length of pileup
	fwrite(&l_sig, sizeof(int), 1, f_op);	
	int n_4bit_posns = 0;
	int n_8bit_posns = 0;
	int n_12bit_posns = 0;
	int n_14bit_posns = 0;
	int n_16bit_posns = 0;
	int n_0_signal = 0;

	int cur_0_run_start = -1;
	for(int i = 1; i <= l_sig; i++)
	{
		// Dump the current 5 levels: 		
		// Generate the existence flag.
		//if(i % 1000 == 0)
		//{
		//	fprintf(stderr, "Position %d:\n", i);
		//}

		unsigned char existence_flag = 0;
		unsigned char drange_flag = 0;
		unsigned short max_val = 0;
		unsigned int total_val = 0;
		for(int part_i = 0; part_i < n_partitions; part_i++)
		{
			for(int allele_i = 0; allele_i < 5; allele_i++)
			{
				if(pileup[part_i][allele_i][i] > 0)
				{	
					total_val += (unsigned int)(pileup[part_i][allele_i][i]);

					if(max_val < pileup[part_i][allele_i][i])
					{
						max_val = pileup[part_i][allele_i][i];
					}

					existence_flag = existence_flag | (1 << allele_i);
					//fprintf(stderr, "Nuc %d: %d\n", allele_i, pileup[allele_i][i]);
					//getc(stdin);
				} // signal check.
			} // allele_i loop.
		} // part_i loop.

		// Following is an RLE check.
		if(total_val == 0)
		{
			if(cur_0_run_start == -1)
			{
				// Initiate the 0 run.
				cur_0_run_start = i;

				if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
				{
					fprintf(stderr, "Initiated a 0-run @ %d\n", i);
				}
				
			}
			else
			{
				// We are already in a zero run.
			}
		} // total_val 0 check.
		else
		{
			// Check if we will do a RLE dump.
			// Make sure there was a run of zeros before this position.
			if(cur_0_run_start != -1)
			{
				bool do_RLE_dump = ((i - cur_0_run_start) > 5);

				if(do_RLE_dump)
				{
					//fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

					// Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
					// Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
					unsigned char RLE_indicator = 0xFF;
					fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Writing: %x\n", RLE_indicator);
					}

					// Write the length of 0 length region: Note that i does not have a 0, it must not be included.
					unsigned int l_RLE = (i - cur_0_run_start);
					fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
					}
				} // RLE dump check.
				else
				{
					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);
					}

					// Save the position with normal dump where we dump 0's.
					for(int j = cur_0_run_start; j < i; j++)
					{
						// Write 0 to each position.
						unsigned char RLE_indicator = 0;
						fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
					} // j loop.
				} // non-RLE check.

				// Reset the RLE start position.
				cur_0_run_start = -1;
			} // 0-run check.	

			if(max_val == 0)
			{
				if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
				{
					fprintf(stderr, "Sanity check failed: Max value is 0 on a non-zero dumping path.\n");
				}
				exit(0);
			}

			// Continue dumping the value for the current position.
			// Set the data range flag.
			if(max_val >= (1 << 8)-1)
			{
				drange_flag = 1;
			}

			// If val is not char, flag it to make sure we dump the right number of bytes.
			existence_flag = existence_flag | (drange_flag << 5);

			if(existence_flag == 0xFF)
			{
				if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
				{
					fprintf(stderr, "Sanity check failed: Existence flag is 0xFF.\n");
				}
				exit(0);
			}

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Normal dumping value %d @ %d, Existence flag: 0x%x, drange_flag: %d\n", total_val, i, existence_flag, drange_flag);
			}

			// Do a normal dump.
			// Write existence flag.
			fwrite(&existence_flag, sizeof(unsigned char), 1, f_op);

			// Dump the values.
			for(int part_i = 0; part_i < n_partitions; part_i++)
			{
				for(int allele_i = 0; allele_i < 5; allele_i++)
				{
					// Following ensures that we write only the positions that exist.
					//if(pileup[part_i][allele_i][i] > 0)
					if((existence_flag & (1 << allele_i)) > 0)
					{
						if(drange_flag == 1)
						{	
							unsigned short current_count_short = (unsigned short)(pileup[part_i][allele_i][i]);
							fwrite(&(current_count_short), sizeof(unsigned short), 1, f_op);
						}
						else
						{
							unsigned char current_count_char = (unsigned char)(pileup[part_i][allele_i][i]);
							fwrite(&(current_count_char), sizeof(unsigned char), 1, f_op);
						}
						//fprintf(stderr, "%d, ", pileup[part_i][allele_i][i]);
					} // count check for current position.
				} // allele_i loop.
			} // part_i loop.

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "\n");
			}
		} // total value check.		

		// Get some simple stats.
		if(max_val == 0)
		{
			n_0_signal++;
		}
		else if(max_val <= ((2<<4)-1))
		{
			n_4bit_posns++;
		}
		else if(max_val <= ((2<<8)-1))
		{
			n_8bit_posns++;
		}
		else if(max_val <= ((2<<12)-1))
		{
			n_12bit_posns++;
		}
		else if(max_val <= ((2<<14)-1))
		{
			n_14bit_posns++;
		}
		else if(max_val <= ((2<<16)-1))
		{
			n_16bit_posns++;
		}
	} // i loop.

	// If there was a 0-run at the end of the file, dump the final RLE.
	if(cur_0_run_start != -1)
	{
		int i = l_sig+1;
		bool do_RLE_dump = ((i - cur_0_run_start) > 5);

		if(do_RLE_dump)
		{
			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);
			}

			// Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
			// Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
			unsigned char RLE_indicator = 0xFF;
			fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Writing: %x\n", RLE_indicator);
			}

			// Write the length of 0 length region: Note that i does not have a 0, it must not be included.
			unsigned int l_RLE = (i - cur_0_run_start);
			fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
			}
		} // RLE dump check.
		else
		{
			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);
			}

			// Save the position with normal dump where we dump 0's.
			for(int j = cur_0_run_start; j < i; j++)
			{
				// Write 0 to each position.
				unsigned char RLE_indicator = 0;
				fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
			} // j loop.
		} // non-RLE check.
	} // RLE check for the end of file.

	close_f(f_op, op_fp);
	fprintf(stderr, "%d 0 signal, %d 4-bit, %d 8-bit, %d 12-bit, %d 14-bit, %d 16-bit positions.\n", n_0_signal, n_4bit_posns, n_8bit_posns, n_12bit_posns, n_14bit_posns, n_16bit_posns);

	// Now, load and test if we received the correct values.
	
	bool loading_check = false;
	if(loading_check)
	{
		fprintf(stderr, "Loading and comparing.\n");

		// Load the compressed pileup; then do sanity check.
		int l_loaded_pileup = 0;
		int n_loaded_partitions = 0;
		unsigned short*** loaded_phred_partitioned_pileup = load_partitioned_compressed_pileups(op_fp, n_loaded_partitions, l_loaded_pileup);

		// Compare the loaded and the actual pileups.
		for(int part_i = 0; part_i < n_loaded_partitions; part_i++)
		{
			for(int i_allele = 0; i_allele < 5; i_allele++)
			{
				for(int i = 1; i <= l_sig; i++)
				{
					if(loaded_phred_partitioned_pileup[part_i][i_allele][i] != pileup[part_i][i_allele][i])
					{
						fprintf(stderr, "Sanity check failed: partition: %d; posn: %d; allele: %d: Loaded: %d, Original: %d\n", 
							part_i, i, i_allele, loaded_phred_partitioned_pileup[part_i][i_allele][i], pileup[part_i][i_allele][i]);
						exit(0);
					}
				} // i loop.				
			} // i_allele loop.
		} // part_i loop.

		// Free memory for both loaded and dumped pileups.
		for(int part_i = 0; part_i < n_partitions; part_i++)
		{
			for(int i_allele = 0; i_allele < 5; i_allele++)
			{
				delete [] loaded_phred_partitioned_pileup[part_i][i_allele];
			} // i_allele loop.
			delete [] loaded_phred_partitioned_pileup[part_i];
		} // part_i loop.
	}
	else
	{
		fprintf(stderr, "Skipping loading check.\n");
	}
}

void compress_nucleotide_pileup_track(unsigned short** pileup, int l_sig, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	
	// Write the length of pileup
	fwrite(&l_sig, sizeof(int), 1, f_op);
	int n_4bit_posns = 0;
	int n_8bit_posns = 0;
	int n_12bit_posns = 0;
	int n_14bit_posns = 0;
	int n_16bit_posns = 0;
	int n_0_signal = 0;

	// This is the 0 RLE start position.
	int cur_0_run_start = -1;
	for(int i = 1; i <= l_sig; i++)
	{
		// Dump the current 5 levels: 		
		// Generate the existence flag.
		//fprintf(stderr, "Position %d:\n", i);
		unsigned char existence_flag = 0;
		unsigned char drange_flag = 0;
		unsigned short max_val = 0;
		unsigned int total_val = 0;
		for(int allele_i = 0; allele_i < 5; allele_i++)
		{
			if(pileup[allele_i][i] > 0)
			{
				total_val += (unsigned int)(pileup[allele_i][i]);

				if(max_val < pileup[allele_i][i])
				{
					max_val = pileup[allele_i][i];
				}

				existence_flag = existence_flag | (1 << allele_i);
				//fprintf(stderr, "Nuc %d: %d\n", allele_i, pileup[allele_i][i]);
				//getc(stdin);
			} // signal check.
		} // allele_i loop.

		// Write the 1 byte phred-based partitioning flag.
		// Following is an RLE check.
		if(total_val == 0)
		{
			if(cur_0_run_start == -1)
			{
				// Initiate the 0 run.
				cur_0_run_start = i;

				if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
				{
					fprintf(stderr, "Initiated a 0-run @ %d\n", i);
				}
			}
			else
			{
				// We are already in a zero run.
			}
		} // total_val 0 check.
		else
		{
			// Check if we will do a RLE dump.
			// Make sure there was a run of zeros before this position.
			if(cur_0_run_start != -1)
			{
				// Check the RLE length and dump if requested.
				bool do_RLE_dump = ((i - cur_0_run_start) > 5);

				if(do_RLE_dump)
				{
					//fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

					// Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
					// Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
					unsigned char RLE_indicator = 0xFF;
					fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Writing: %x\n", RLE_indicator);
					}

					// Write the length of 0 length region: Note that i does not have a 0, it must not be included.
					unsigned int l_RLE = (i - cur_0_run_start);
					fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
					}
				} // RLE dump check.
				else
				{
					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);
					}

					// Save the position with normal dump where we dump 0's.
					for(int j = cur_0_run_start; j < i; j++)
					{
						// Write 0 to each position.
						unsigned char RLE_indicator = 0;
						fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
					} // j loop.
				} // non-RLE check.

				// Reset the RLE start position.
				cur_0_run_start = -1;
			} // 0-run check.	

			if(max_val == 0)
			{
				if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
				{
					fprintf(stderr, "Sanity check failed: Max value is 0 on a non-zero dumping path.\n");
				}
				exit(0);
			}
			
			// Set the data range flag: 1 for short, 0 for char.
			if(max_val >= (1 << 8)-1)
			{
				drange_flag = 1;
			}

			// If val is not char, flag it to make sure we dump the right number of bytes.
			existence_flag = existence_flag | (drange_flag << 5);

			// Write existence flag.
			fwrite(&existence_flag, sizeof(unsigned char), 1, f_op);

			// Dump the numbers.
			for(int allele_i = 0; allele_i < 5; allele_i++)
			{
				if(pileup[allele_i][i] > 0)
				{
					if(drange_flag == 1)
					{	
						unsigned short current_count_short = (unsigned short)(pileup[allele_i][i]);
						fwrite(&(current_count_short), sizeof(unsigned short), 1, f_op);
					}
					else
					{
						unsigned char current_count_char = (unsigned char)(pileup[allele_i][i]);
						fwrite(&(current_count_char), sizeof(unsigned char), 1, f_op);
					}
				} // count check for current position.
			} // allele_i loop.

		} // total value check.

		// Get some simple stats.
		if(max_val == 0)
		{
			n_0_signal++;
		}
		else if(max_val <= ((2<<4)-1))
		{
			n_4bit_posns++;
		}
		else if(max_val <= ((2<<8)-1))
		{
			n_8bit_posns++;
		}
		else if(max_val <= ((2<<12)-1))
		{
			n_12bit_posns++;
		}
		else if(max_val <= ((2<<14)-1))
		{
			n_14bit_posns++;
		}
		else if(max_val <= ((2<<16)-1))
		{
			n_16bit_posns++;
		}
	} // i loop.
	close_f(f_op, op_fp);
	fprintf(stderr, "%d 0 signal, %d 4-bit, %d 8-bit, %d 12-bit, %d 14-bit, %d 16-bit positions.\n", n_0_signal, n_4bit_posns, n_8bit_posns, n_12bit_posns, n_14bit_posns, n_16bit_posns);

	// Now, load and test if we received the correct values.
	bool loading_check = false;
	if(loading_check)
	{
		fprintf(stderr, "Loading and comparing.\n");

		// Load the compressed pileup; then do sanity check.
		int l_loaded_pileup = 0;
		unsigned short** loaded_pileup = load_compressed_pileups(op_fp, l_loaded_pileup);

		// Compare the loaded and the actual pileups.
		for(int i_allele = 0; i_allele < 5; i_allele++)
		{
			for(int i = 1; i <= l_sig; i++)
			{
				if(loaded_pileup[i_allele][i] != pileup[i_allele][i])
				{
					fprintf(stderr, "Sanity check failed: posn: %d; allele: %d: %d,%d\n", 
						i, i_allele, loaded_pileup[i_allele][i], pileup[i_allele][i]);
					exit(0);
				}
			} // i loop.				
		} // i_allele loop.

		// Free memory for both loaded and dumped pileups.
		for(int i_allele = 0; i_allele < 5; i_allele++)
		{
			delete [] loaded_pileup[i_allele];
		} // i_allele loop.
	}
	else
	{
		fprintf(stderr, "Skipping loading check.\n");
	}
}

void dump_nucleotide_pileup_per_SAM_file_orig(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, int& n_processed_reads)
{
	//// Init the file.
	//char summary_fp[1000];
	//sprintf(summary_fp, "%s/%s", op_dir, t_config_params::OP_filenames[OP_PILEUP_SUMMARY_FN]);
	//FILE* f_summary = open_f(summary_fp, "w");
	//fclose(f_summary);

	// Initialize the number of processed readss.
	n_processed_reads = 0;

	// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
	bool allele_counts_are_there = true;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
		if(!check_file(cur_allele_fp))
		{
			// This file does not exist, not all the allele counts are there.
			allele_counts_are_there = false;
			break;
		}
	} // i_chr loop.

	if(allele_counts_are_there)
	{
		fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
		return;
	}

	fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

	fprintf(stderr, "Allocating pileup memory.\n");
	unsigned short*** per_chrom_nuc_count_per_allele = new unsigned short**[(int)chr_ids->size()];
	//int** per_chrom_coverage = new int*[(int)chr_ids->size()];
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_sig = chr_lengths->at(i_chr);

		//per_chrom_coverage[i_chr] = new int[l_sig + 2];
		//memset(per_chrom_coverage[i_chr], 0, sizeof(int) * l_sig);

		per_chrom_nuc_count_per_allele[i_chr] = allocate_pileup(l_sig);
	} // i_chr loop.
	fprintf(stderr, "Done.\n");

	int max_n_alleles_per_posn = 1<<16;
	fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

	// Start reading the file.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = fopen(sam_fp, "r");
	}

	int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
	//for(int i = 0; i < 256; i++)
	//{
	//	per_nuc_val[i] = 4;
	//} // i loop.
	//per_nuc_val[(int)'A'] = 0;
	//per_nuc_val[(int)'C'] = 1;
	//per_nuc_val[(int)'G'] = 2;
	//per_nuc_val[(int)'T'] = 3;
	//per_nuc_val[(int)'U'] = 3;
	//per_nuc_val[(int)'a'] = 0;
	//per_nuc_val[(int)'c'] = 1;
	//per_nuc_val[(int)'g'] = 2;
	//per_nuc_val[(int)'t'] = 3;
	//per_nuc_val[(int)'u'] = 3;

	// Enter file reading loop.
	char read_id[1000];
	//int flag;
	char chrom[100];
	int chr_index;
	char mapping_map_str[10000];
	char cur_fragment[1000];
	char flag_str[100];
	char _chr_index_str[100];
	char phred_quality_str[100000];
	char mapp_quality_str[100];

	int n_unmapped_reads = 0;
	int n_low_quality_reads = 0;

	//int phred_score_base = -123;

	while(1)
	{
		char* cur_line = getline(f_sam);

		if(cur_line == NULL)
		{
			break;
		}

		// If this is a comment line, skip it.
		if(cur_line[0] == '@')
		{
			delete [] cur_line;
			continue;
		}

		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
		{
			//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
			// If the quality is not adequate, do not use this read.
			if(atoi(mapp_quality_str) < min_mapp_qual)
			{
				n_low_quality_reads++;
				delete [] cur_line;
				continue;
			}

			// Make sure the normalized chromosome ids match.
			normalize_chr_id(chrom);

			int flag = atoi(flag_str);

			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If we do not have the chromosome, do not process.
			if(i_chr == (int)chr_ids->size())
			{
				n_unmapped_reads++;
				delete [] cur_line;
				continue;
			}

			int _chr_index = atoi(_chr_index_str);

			// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
			chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

			// Sanity check. Is this fragment mapped?
			if(flag & 0x04)
			{
				// The read is not mapping.
			}
			else
			{
				n_processed_reads++;
					
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %d. read             \r", n_processed_reads);
				}

				int i_mapp_map = 0;
				//t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
				while(mapping_map_str_valid && 
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map, 
														is_matching,
														l_cur_entry,
														entry_type_char);

					// If this block is matching, update the pileup.
					if(is_matching)
					{
						// Update the counts for the nucleotides.
						int cur_read_nuc_i = read_nuc_index;
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							//int cur_nuc_num = nuc_2_num(cur_fragment[cur_read_nuc_i]);
							int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

							//if(numerized_sequence_signal[cur_genome_i] > 0 &&
							//	nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
							//{
							//	fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
							//}

							if(cur_nuc_num < 4)
							{
								// Update the count: The allelic counts must be checked for bounds.
								if(per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
								{
									per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i]++;
								}

								//per_chrom_coverage[i_chr][cur_genome_i]++;
								cur_read_nuc_i++;

								//fprintf(stderr, "%d: %c, %c\n", i_cur, num_2_nuc(cur_nuc_num), num_2_nuc(numerized_sequence_signal[i_cur]));
							}
						} // i_cur loop.
					}
					else if(entry_type_char == 'D')
					{
						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							per_chrom_nuc_count_per_allele[i_chr][4][cur_genome_i]++;
						} // cur_genome_i loop.
					}
					else if(entry_type_char == 'I')
					{
						// Insertion to the reference: This is included to one position.
						per_chrom_nuc_count_per_allele[i_chr][4][chr_index]++;
					}
					// Update the base for the current entry.
					if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				//delete(cur_entry_length_str);
			} // mapping check for the current read.
		} // SAM read line parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "Finished reading.\n");

	fclose(f_sam);

	// Dump the counts per position.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));

		// compress and dump the current file.
		compress_nucleotide_pileup_track(per_chrom_nuc_count_per_allele[i_chr], chr_lengths->at(i_chr), cur_allele_fp);
	} // i_chr loop.

	//char analysis_summary_str[1000];
	//sprintf(analysis_summary_str, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality, n_unmapped_reads);
	//t_config_params::copy_analysis_summary_string(analysis_summary_str);

	//// Write the summary file.
	//f_summary = open_f(summary_fp, "w");
	//fprintf(f_summary, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality_reads, n_unmapped_reads);
	//fclose(f_summary);
}

void dump_nucleotide_pileup_per_SAM_file(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, unsigned long long& n_processed_reads)
{
	fprintf(stderr, "Generating pileup from SAM file with min mapp qual %d, min base qual %d\n", min_mapp_qual, min_phred_qual);

	//// Init the file.
	//char summary_fp[1000];
	//sprintf(summary_fp, "%s/%s", op_dir, t_config_params::OP_filenames[OP_PILEUP_SUMMARY_FN]);
	//FILE* f_summary = open_f(summary_fp, "w");
	//fclose(f_summary);

	// Initialize the number of processed readss.
	n_processed_reads = 0;

	//// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
	//bool allele_counts_are_there = true;
	//for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	//{
	//	char cur_allele_fp[1000];
	//	sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
	//	if(!check_file(cur_allele_fp))
	//	{
	//		// This file does not exist, not all the allele counts are there.
	//		allele_counts_are_there = false;
	//		break;
	//	}
	//} // i_chr loop.

	//if(allele_counts_are_there)
	//{
	//	fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
	//	return;
	//}

	fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

	fprintf(stderr, "Allocating pileup memory.\n");
	unsigned short*** per_chrom_nuc_count_per_allele = new unsigned short**[(int)chr_ids->size()];
	//int** per_chrom_coverage = new int*[(int)chr_ids->size()];
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_sig = chr_lengths->at(i_chr);

		// Check the existing pileup file, if it is there, add to it.
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
		if(check_file(cur_allele_fp))
		{
			int l_loaded_pileup = 0;
			per_chrom_nuc_count_per_allele[i_chr] = load_compressed_pileups(cur_allele_fp, l_loaded_pileup);
			if(l_loaded_pileup != chr_lengths->at(i_chr))
			{
				fprintf(stderr, "Loaded pileup is not the same length as the chromosome info: %s: %d, %d\n", chr_ids->at(i_chr), chr_lengths->at(i_chr), l_loaded_pileup);
				exit(0);
			}
			else
			{
				fprintf(stderr, "Loaded existing pileup from %s.\n", cur_allele_fp);
			}
		}
		else
		{
			per_chrom_nuc_count_per_allele[i_chr] = allocate_pileup(l_sig);
		}
	} // i_chr loop.
	fprintf(stderr, "Done.\n");

	int max_n_alleles_per_posn = 1<<16;
	fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

	// Start reading the file.
	FILE* f_sam = open_f(sam_fp, "r");

	int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
	//for(int i = 0; i < 256; i++)
	//{
	//	per_nuc_val[i] = 4;
	//} // i loop.
	//per_nuc_val[(int)'A'] = 0;
	//per_nuc_val[(int)'C'] = 1;
	//per_nuc_val[(int)'G'] = 2;
	//per_nuc_val[(int)'T'] = 3;
	//per_nuc_val[(int)'U'] = 3;
	//per_nuc_val[(int)'a'] = 0;
	//per_nuc_val[(int)'c'] = 1;
	//per_nuc_val[(int)'g'] = 2;
	//per_nuc_val[(int)'t'] = 3;
	//per_nuc_val[(int)'u'] = 3;

	// Enter file reading loop.
	char read_id[1000];
	//int flag;
	char chrom[100];
	int chr_index;
	char mapping_map_str[10000];
	char cur_fragment[1000];
	char flag_str[100];
	char _chr_index_str[100];
	char phred_quality_str[100000];
	char mapp_quality_str[100];

	int n_unmapped_reads = 0;
	int n_low_quality_reads = 0;

	int phred_score_base = -123;

	while(1)
	{
		char* cur_line = getline(f_sam);

		if(cur_line == NULL)
		{
			break;
		}

		// If this is a comment line, skip it.
		if(cur_line[0] == '@')
		{
			delete [] cur_line;
			continue;
		}

		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
		{
			bool phred_exists = true;
			if (t_string::compare_strings(phred_quality_str, "*") ||
				t_string::string_length(phred_quality_str) != t_string::string_length(cur_fragment))
			{
				phred_exists = false;
			}

			// Note that this can be updated for any read that contains a valid quality score.
			if(phred_exists &&
				phred_score_base == -123)
			{
				int l_read = strlen(phred_quality_str);
				for(int i = 0; i < l_read; i++)
				{
					if(phred_quality_str[i] > 'J')
					{
						fprintf(stderr, "Phred+64 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = ';';
						break;
					}
					else if(phred_quality_str[i] < ';')
					{
						fprintf(stderr, "Phred+33 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = '!';
						break;
					}
				} // i loop.
			}
			
			//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
			// If the quality is not adequate, do not use this read.
			if(atoi(mapp_quality_str) < min_mapp_qual)
			{
				n_low_quality_reads++;
				delete [] cur_line;
				continue;
			}

			// Make sure the normalized chromosome ids match.
			normalize_chr_id(chrom);

			int flag = atoi(flag_str);

			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If we do not have the chromosome, do not process.
			if(i_chr == (int)chr_ids->size())
			{
				n_unmapped_reads++;
				delete [] cur_line;
				continue;
			}

			int _chr_index = atoi(_chr_index_str);

			// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
			chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

			// Sanity check. Is this fragment mapped?
			if(flag & 0x04)
			{
				// The read is not mapping.
			}
			else
			{
				n_processed_reads++;
					
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %llu. read             \r", n_processed_reads);
				}

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				fprintf(stderr, "Processing:\n%s\n%s\n", cur_fragment, phred_quality_str);
}

				int i_mapp_map = 0;
				//t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
				while(mapping_map_str_valid && 
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map, 
														is_matching,
														l_cur_entry,
														entry_type_char);

					// If this block is matching, update the pileup.
					if(is_matching)
					{
						// Update the counts for the nucleotides.
						int cur_read_nuc_i = read_nuc_index;
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

							//if(numerized_sequence_signal[cur_genome_i] > 0 &&
							//	nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
							//{
							//	fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
							//}

							if(cur_nuc_num < 4)
							{
								// Update the count: The allelic counts must be checked for bounds.
								if(per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
								{
									// If phred is missing, process all of them, if not, check the quality.
									bool phred_check = true;
									if (phred_exists)
									{
										phred_check = ((phred_quality_str[cur_read_nuc_i] - phred_score_base) > min_phred_qual);
									}

									// Does phred check hold?
									if(phred_check)
									{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
										fprintf(stderr, "Adding %d: %c, %c (%d)\n", 
												cur_read_nuc_i, 
												cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
												(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}

										// Phred quality holds, update the pileup position.
										if (cur_genome_i <= chr_lengths->at(i_chr))
										{
											per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i]++;
										}
										else
										{
											fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
												chr_ids->at(i_chr), cur_genome_i, chr_lengths->at(i_chr));
										}
									} // phred qual pass check.
									else
									{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
										fprintf(stderr, "Skipping %d: %c, %c (%d)\n", 
											cur_read_nuc_i, 
											cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
											(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}
									} // phred qual nopass check.
								} // max_n_alleles_per_posn check.
							} // char < 4 check.

							// Update the read's nucleotide index.
							cur_read_nuc_i++;
						} // i_cur loop.
					} // genomic match check.
					else if(entry_type_char == 'D')
					{
						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							if (cur_genome_i <= chr_lengths->at(i_chr))
							{
								per_chrom_nuc_count_per_allele[i_chr][4][cur_genome_i]++;
							}
							else
							{
								fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
									chr_ids->at(i_chr), cur_genome_i, chr_lengths->at(i_chr));
							}
						} // cur_genome_i loop.
					}
					else if(entry_type_char == 'I')
					{
						// Insertion to the reference: This is included to one position.
						if (chr_index < chr_lengths->at(i_chr))
						{
							per_chrom_nuc_count_per_allele[i_chr][4][chr_index]++;
						}
					}
					// Update the base for the current entry.
					if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				//delete(cur_entry_length_str);
			} // mapping check for the current read.
		} // SAM read line parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "Finished reading.\n");

	close_f(f_sam, sam_fp);

	// Dump the counts per position.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin.gz", op_dir, chr_ids->at(i_chr));

		// compress and dump the current file.
		compress_nucleotide_pileup_track(per_chrom_nuc_count_per_allele[i_chr], chr_lengths->at(i_chr), cur_allele_fp);
	} // i_chr loop.

	//char analysis_summary_str[1000];
	//sprintf(analysis_summary_str, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality, n_unmapped_reads);
	//t_config_params::copy_analysis_summary_string(analysis_summary_str);

	//// Write the summary file.
	//f_summary = open_f(summary_fp, "w");
	//fprintf(f_summary, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality_reads, n_unmapped_reads);
	//fclose(f_summary);
}

void dump_nucleotide_pileup_per_SAM_file_phred_partitioning(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, vector<int>* phred_qual_partitions, unsigned long long& n_processed_reads)
{
	// Sort the quality partition boundary values.
	int max_qual = 10000;
	phred_qual_partitions->push_back(0);
	phred_qual_partitions->push_back(max_qual);
	sort(phred_qual_partitions->begin(), phred_qual_partitions->end());
	fprintf(stderr, "Generating pileup from SAM file with min mapp qual %d, with %d phred partitions @:\n", min_mapp_qual, (int)phred_qual_partitions->size()-1);

	// Check the quality partitions and merge the ones that overlap.
	vector<int>* uniq_phred_qual_partitions = new vector<int>();
	for(int i_p = 0; i_p < (int)phred_qual_partitions->size(); i_p++)
	{
		if(i_p+1 < (int)phred_qual_partitions->size() &&
			phred_qual_partitions->at(i_p) == phred_qual_partitions->at(i_p+1))
		{
			fprintf(stderr, "Merging partitions that match at qualities of %d\n", phred_qual_partitions->at(i_p));
		}
		else
		{
			uniq_phred_qual_partitions->push_back(phred_qual_partitions->at(i_p));
		}
	} // i_p loop.

	// Replace the phred quality partitions.
	phred_qual_partitions->clear();
	phred_qual_partitions->insert(phred_qual_partitions->end(), uniq_phred_qual_partitions->begin(), uniq_phred_qual_partitions->end());

	int* part_i_per_phred_qual = new int[max_qual+1];
	for(int i_p = 1; i_p < (int)phred_qual_partitions->size(); i_p++)
	{
		fprintf(stderr, "Partition %d: [%d-%d]\n", i_p-1, phred_qual_partitions->at(i_p-1)+1, phred_qual_partitions->at(i_p));
		for(int i = phred_qual_partitions->at(i_p-1)+1; i <= phred_qual_partitions->at(i_p); i++)
		{
			part_i_per_phred_qual[i] = i_p - 1;
		} // i loop.
	} // i_p loop.
	part_i_per_phred_qual[0] = 0;

	int n_partitions = (int)phred_qual_partitions->size()-1;

	// Initialize the number of processed readss.
	n_processed_reads = 0;

	// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
	bool allele_counts_are_there = true;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
		if(!check_file(cur_allele_fp))
		{
			// This file does not exist, not all the allele counts are there.
			allele_counts_are_there = false;
			break;
		}
	} // i_chr loop.

	if(allele_counts_are_there)
	{
		fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
		return;
	}

	fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

	fprintf(stderr, "Allocating pileup memory.\n");
	unsigned short**** per_chrom_per_part_nuc_count_per_allele = new unsigned short***[(int)chr_ids->size()];
	//double total_mem = 0;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_sig = chr_lengths->at(i_chr);

		per_chrom_per_part_nuc_count_per_allele[i_chr] = new unsigned short**[n_partitions];

		for(int part_i = 0; part_i < n_partitions; part_i++)
		{
			per_chrom_per_part_nuc_count_per_allele[i_chr][part_i] = allocate_pileup(l_sig);
		} // part_i loop.
	} // i_chr loop.
	fprintf(stderr, "\nDone.\n");

	int max_n_alleles_per_posn = 1<<16;
	fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

	// Start reading the file.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = fopen(sam_fp, "r");
	}

	int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
	//for(int i = 0; i < 256; i++)
	//{
	//	per_nuc_val[i] = 4;
	//} // i loop.
	//per_nuc_val[(int)'A'] = 0;
	//per_nuc_val[(int)'C'] = 1;
	//per_nuc_val[(int)'G'] = 2;
	//per_nuc_val[(int)'T'] = 3;
	//per_nuc_val[(int)'U'] = 3;
	//per_nuc_val[(int)'a'] = 0;
	//per_nuc_val[(int)'c'] = 1;
	//per_nuc_val[(int)'g'] = 2;
	//per_nuc_val[(int)'t'] = 3;
	//per_nuc_val[(int)'u'] = 3;

	// Enter file reading loop.
	char read_id[1000];
	//int flag;
	char chrom[100];
	int chr_index;
	char mapping_map_str[10000];
	char cur_fragment[1000];
	char flag_str[100];
	char _chr_index_str[100];
	char phred_quality_str[100000];
	char mapp_quality_str[100];

	int n_unmapped_reads = 0;
	int n_low_quality_reads = 0;

	int phred_score_base = -123;

	while(1)
	{
		char* cur_line = getline(f_sam);

		if(cur_line == NULL)
		{
			break;
		}

		// If this is a comment line, skip it.
		if(cur_line[0] == '@')
		{
			delete [] cur_line;
			continue;
		}

		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
		{
			if(phred_score_base == -123)
			{
				int l_read = strlen(phred_quality_str);
				for(int i = 0; i < l_read; i++)
				{
					if(phred_quality_str[i] > 'J')
					{
						fprintf(stderr, "Phred+64 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = ';';
						break;
					}
					else if(phred_quality_str[i] < ';')
					{
						fprintf(stderr, "Phred+33 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = '!';
						break;
					}
				} // i loop.
			}

			//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
			// If the quality is not adequate, do not use this read.
			if(atoi(mapp_quality_str) < min_mapp_qual)
			{
				n_low_quality_reads++;
				delete [] cur_line;
				continue;
			}

			// Make sure the normalized chromosome ids match.
			normalize_chr_id(chrom);

			int flag = atoi(flag_str);

			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If we do not have the chromosome, do not process.
			if(i_chr == (int)chr_ids->size())
			{
				n_unmapped_reads++;
				delete [] cur_line;
				continue;
			}

			int _chr_index = atoi(_chr_index_str);

			// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
			chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

			// Sanity check. Is this fragment mapped?
			if(flag & 0x04)
			{
				// The read is not mapping.
			}
			else
			{
				n_processed_reads++;
					
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %llu. read             \r", n_processed_reads);
				}

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				fprintf(stderr, "Processing:\n%s\n%s\n", cur_fragment, phred_quality_str);
}

				int i_mapp_map = 0;
				//t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
				while(mapping_map_str_valid && 
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map, 
														is_matching,
														l_cur_entry,
														entry_type_char);

					// If this block is matching, update the pileup.
					if(is_matching)
					{
						// Update the counts for the nucleotides.
						int cur_read_nuc_i = read_nuc_index;
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

							//if(numerized_sequence_signal[cur_genome_i] > 0 &&
							//	nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
							//{
							//	fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
							//}

							if(cur_nuc_num < 4)
							{
								int cur_nuc_phred_qual = phred_quality_str[cur_read_nuc_i] - phred_score_base;

								// Get the partition index.
								int cur_read_phred_part_i = part_i_per_phred_qual[cur_nuc_phred_qual];

								// Update the count: The allelic counts must be checked for bounds.
								if(per_chrom_per_part_nuc_count_per_allele[i_chr][cur_read_phred_part_i][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
								{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
									fprintf(stderr, "Adding %d: %c, %c (%d)\n", 
											cur_read_nuc_i, 
											cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
											(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}

										// Phred quality holds, update the pileup position.
										per_chrom_per_part_nuc_count_per_allele[i_chr][cur_read_phred_part_i][cur_nuc_num][cur_genome_i]++;
//									} // phred qual pass check.
//									else
//									{
//if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
//{
//										fprintf(stderr, "Skipping %d: %c, %c (%d)\n", 
//											cur_read_nuc_i, 
//											cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
//											(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
//}
//									} // phred qual nopass check.
								} // max_n_alleles_per_posn check.
							} // char < 4 check.

							// Update the read's nucleotide index.
							cur_read_nuc_i++;
						} // i_cur loop.
					} // genomic match check.
					else if(entry_type_char == 'D')
					{
						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							per_chrom_per_part_nuc_count_per_allele[i_chr][n_partitions-1][4][cur_genome_i]++;
						} // cur_genome_i loop.
					}
					else if(entry_type_char == 'I')
					{
						// Insertion to the reference: This is included to one position.
						per_chrom_per_part_nuc_count_per_allele[i_chr][n_partitions-1][4][chr_index]++;
					}
					// Update the base for the current entry.
					if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				//delete(cur_entry_length_str);
			} // mapping check for the current read.
		} // SAM read line parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "Finished reading.\n");

	fclose(f_sam);

	// Dump the counts per position.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin.gz", op_dir, chr_ids->at(i_chr));

		// compress and dump the current file.
		compress_partitioned_nucleotide_pileup_track(per_chrom_per_part_nuc_count_per_allele[i_chr], n_partitions, chr_lengths->at(i_chr), cur_allele_fp);
	} // i_chr loop.
}

void preprocess_SAM_read_line_positional_info(char* cur_line,
	char* chrom,
	int& chr_index, 
	int& mapQ,
	char& strand_char,
	char* cigar_str)
{
	// Skip the comment and headers.
	if (cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag = 0;
	char flag_str[100];
	int _chr_index = 0;
	char _chr_index_str[100];
	char mapQ_str[1000];

	int char_i = 0;
	if (!t_string::get_next_token(cur_line, NULL, 0, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	if (!t_string::get_next_token(cur_line, flag_str, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	if (!t_string::get_next_token(cur_line, chrom, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	if (!t_string::get_next_token(cur_line, _chr_index_str, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}
	
	if (!t_string::get_next_token(cur_line, mapQ_str, 1000, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}	

	if(!t_string::get_next_token(cur_line, cigar_str, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	// Is this fragment mapped?
	if (flag & 0x04)
	{
		// The read is not mapping.
		chrom[0] = 0;
		mapQ = 0;
		chr_index = 0;
		flag = 0;
		strand_char = 0;
	}
	else
	{
		// Tramslate the coordinate and strand.
		_chr_index = atoi(_chr_index_str);
		mapQ = atoi(mapQ_str);
		chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);
		flag = atoi(flag_str);

		strand_char = 'F';
		if (flag & 0x10)
		{
			strand_char = 'R';
		}
	}

	//t_string::clean_tokens(cur_tokens);
}

void get_empty_qual(int l_read, char* qual_buffer, char base_qual)
{
	for (int i = 0; i < l_read; i++)
	{
		qual_buffer[i] = base_qual;
	} // i loop.

	qual_buffer[l_read] = 0;
}

// 3 bits per nuc.
void compress_nucleotide_seq(char* fragment, int l_frag, char* compressed_seq, int& l_comp_seq)
{
	int nuc_i = 0;
	int comp_sym_i = 0;
	while (nuc_i < l_frag)
	{
		char cur_nuc = fragment[nuc_i];

		char next_nuc = 'N';
		if ((nuc_i + 1) < l_frag)
		{
			next_nuc = fragment[nuc_i + 1];
		}
			
		char next_next_nuc = 'N';
		if ((nuc_i + 2) < l_frag)
		{
			next_next_nuc = fragment[nuc_i + 2];
		}

		// Pack each 3 nucleotides into one byte.
		char cur_sym = nuc_2_num(cur_nuc) + 
						nuc_2_num(next_nuc) * 5 + 
						nuc_2_num(next_next_nuc) * 25;
		compressed_seq[comp_sym_i] = cur_sym;
		comp_sym_i++;
		nuc_i += 3;
	} // nuc_i loop.

	l_comp_seq = comp_sym_i;
}

void decompress_nucleotide_seq(char* compressed_seq, int l_comp_seq, char* fragment_buffer, int l_frag)
{
	int nuc_i = 0;
	int comp_sym_i = 0;
	while (comp_sym_i < l_comp_seq)
	{
		char cur_sym = compressed_seq[comp_sym_i];

		// Save current nucleotide.
		char cur_nuc = num_2_nuc(cur_sym % 5);
		cur_sym /= 5;
		if (nuc_i < l_frag)
		{
			fragment_buffer[nuc_i] = cur_nuc;
		}

		// Save next nucleotide.
		nuc_i++;
		char next_nuc = num_2_nuc(cur_sym % 5);
		cur_sym /= 5;
		if (nuc_i < l_frag)
		{
			fragment_buffer[nuc_i] = next_nuc;
		}

		// Save next next nucleotide.
		nuc_i++;
		char next_next_nuc = num_2_nuc(cur_sym % 5);
		cur_sym /= 5;
		if (nuc_i < l_frag)
		{
			fragment_buffer[nuc_i] = next_next_nuc;
		}

		if (cur_sym % 5 != 0)
		{
			fprintf(stderr, "Decompression error.\n");
			exit(0);
		}

		// Move to next nucleotide and next symbol.
		nuc_i++;
		comp_sym_i++;
	} // nuc_i loop.

	fragment_buffer[l_frag] = 0;

	//if (nuc_i != l_frag)
	//{
	//	fprintf(stderr, "Sanity check failed while decompressing: %d/%d\n", nuc_i, l_frag);
	//	exit(0);
	//}
}

void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	char flag_str[100];
	int _chr_index;
	char _chr_index_str[100];
	char fragment[100000];
	char phred_quality_str[100000];

	//t_string_tokens* cur_tokens = t_string::tokenize_by_chars(cur_line, "\t");
	//if(sscanf(cur_line, "%s %d %s %d %*s %s %*s %*s %*s %s %s", read_id, &flag, chrom, &_chr_index, cigar_str, fragment, phred_quality_str) == 7)
	//if(cur_tokens->size() >= 11)
	if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, cigar_str, fragment, phred_quality_str) == 7)
	{
		//t_string::copy(read_id, cur_tokens->at(0)->str());
		//flag = atoi(cur_tokens->at(1)->str());
		//t_string::copy(chrom, cur_tokens->at(2)->str());
		//_chr_index = atoi(cur_tokens->at(3)->str());
		//t_string::copy(cigar_str, cur_tokens->at(5)->str());
		//t_string::copy(fragment, cur_tokens->at(9)->str());
		//t_string::copy(phred_quality_str, cur_tokens->at(10)->str());

		_chr_index = atoi(_chr_index_str);
		flag = atoi(flag_str);

		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(flag & 0x10)
		{
			strand_char = 'R';
		}

		// Sanity check. Is this fragment mapped?
		if(flag & 0x04)
		{
			// The read is not mapping.
			chrom[0] = 0;
		}
		else
		{
			chr_index = _chr_index;
			sequenced_length = strlen(fragment);
		}
	}
	else
	{
		chrom[0] = 0;
	}

	//t_string::clean_tokens(cur_tokens);
}

void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char cur_fragment[100];
	char quality_str[100];
	int _chr_index;
	char _strand_char;

	if(sscanf(cur_line, "%s %s %s %*d %*d %*d %s %d %c", read_id, cur_fragment, quality_str, chrom, &_chr_index, &_strand_char) == 6)
	{
		chr_index = _chr_index;
		sequenced_length = strlen(cur_fragment);
		sprintf(mapping_quality_str, "%dM", sequenced_length);

		strand_char = 'F';
		if(_strand_char == '-')
		{
			strand_char = 'R';
		}
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	char nucs[1000];
	if(sscanf(cur_line, "%s %c %s %d %s", read_id, &strand_sign, chrom, &chr_start_index, nucs) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		chr_start_index = translate_coord(chr_start_index, BOWTIE_COORDS::start_base, CODEBASE_COORDS::start_base);

		sprintf(mapping_quality_str, "%dM", (int)strlen(nucs));

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = strlen(nucs);
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;

	const int l_buff = 1000;
	char per_entry_buff[4][l_buff];
	int i_cur_char = 0;

	if (!t_string::get_next_token(cur_line, per_entry_buff[0], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[1], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[2], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[3], 1000, "\t", i_cur_char))
	{
		chrom[0] = 0;
		return;
	}
	else
	{
		strcpy(chrom, per_entry_buff[0]);
		chr_start_index = atoi(per_entry_buff[1]);
		chr_end_index = atoi(per_entry_buff[2]);
		strand_sign = per_entry_buff[3][0];
	}

	// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
	sprintf(mapping_quality_str, "%dM", chr_end_index - chr_start_index);
	sequenced_length = chr_end_index - chr_start_index;

	chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

	// Check the flag and determine the strand.
	strand_char = 'F';
	if (strand_sign == '-')
	{
		strand_char = 'R';
	}

	chr_index = chr_start_index;
}

void preprocess_BED5_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;

	const int l_buff = 1000;
	char per_entry_buff[5][l_buff];
	int i_cur_char = 0;

	if (!t_string::get_next_token(cur_line, per_entry_buff[0], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[1], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[2], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[3], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[4], 1000, "\t", i_cur_char))
	{
		chrom[0] = 0;
		return;
	}
	else
	{
		strcpy(chrom, per_entry_buff[0]);
		chr_start_index = atoi(per_entry_buff[1]);
		chr_end_index = atoi(per_entry_buff[2]);
		strand_sign = per_entry_buff[4][0];
	}

	// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
	sprintf(mapping_quality_str, "%dM", chr_end_index - chr_start_index);
	sequenced_length = chr_end_index - chr_start_index;

	chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

	// Check the flag and determine the strand.
	strand_char = 'F';
	if (strand_sign == '-')
	{
		strand_char = 'R';
	}

	chr_index = chr_start_index;
}

void preprocess_BED6_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;

	const int l_buff = 1000;
	char per_entry_buff[6][l_buff];
	int i_cur_char = 0;

	if (!t_string::get_next_token(cur_line, per_entry_buff[0], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[1], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[2], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[3], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[4], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[5], 1000, "\t", i_cur_char))
	{
		chrom[0] = 0;
		return;
	}
	else
	{
		strcpy(chrom, per_entry_buff[0]);
		chr_start_index = atoi(per_entry_buff[1]);
		chr_end_index = atoi(per_entry_buff[2]);
		strand_sign = per_entry_buff[5][0];
	}

	// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
	sprintf(mapping_quality_str, "%dM", chr_end_index - chr_start_index);
	sequenced_length = chr_end_index - chr_start_index;

	chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

	// Check the flag and determine the strand.
	strand_char = 'F';
	if (strand_sign == '-')
	{
		strand_char = 'R';
	}

	chr_index = chr_start_index;
}

void preprocess_BED456_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str)
{
	int n_cols = t_string::fast_count_tokens(cur_line, true, "\t");
	
	/*fprintf(stderr, "%d columns.\r", n_cols);*/

	if (n_cols == 6)
	{
		preprocess_BED6_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);
	}
	else if(n_cols == 5)
	{
		preprocess_BED5_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);
	}
	else if (n_cols == 4)
	{
		preprocess_BED4_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);
	}
}

double get_n_mapped_nucs(vector<t_mapped_fragment*>* fragments)
{
	double n_mapped_nucs = 0.0;

	for(int i_frag = 0; i_frag < (int)fragments->size(); i_frag++)
	{
		n_mapped_nucs += fragments->at(i_frag)->sequenced_fragment_length;
	} // i_frag loop.

	return(n_mapped_nucs);
}

void buffer_per_nucleotide_preprocessed_read_profile_no_buffer(char* sorted_read_fp,
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int max_l_read,
	int l_buffer, int& l_data)
{
	for(int i = 0; i <= l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initialize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	char strand_char = 0;
	//char cur_fragment[10000];
	char mapping_map_str[10000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = open_f(sorted_read_fp, "r");
	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
		}	

		n_processed_reads++;

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		int read_start = 0;
		int read_end = 0;

		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			if(is_matching)
			{
				if(read_start == 0)
				{
					read_start = chr_index;
				}
			} // CIGAR string matching check

			read_end = chr_index + l_cur_entry - 1;

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Increase the height: The indexing is already 1 based, there is no conversion needed.
		if(signal_profile_buffer != NULL &&
			max_l_read >= (read_end-read_start+1))
		{
			for(int i_cur = read_start; i_cur <= read_end; i_cur++)
			{
				signal_profile_buffer[i_cur]++;
			} // i_cur loop.
		}

		// Update the strand signals, if requested.
		if(forward_strand_signal != NULL &&
			max_l_read >= (read_end-read_start+1))
		{
			if(strand_char == 'F')
			{
				// Update the forward strand signal.
				for(int i_cur = read_start; i_cur <= read_end; i_cur++)
				{
					forward_strand_signal[i_cur]++;
				} // i_cur loop.
			}
			else
			{
				// Update the reverse strand signal.
				for(int i_cur = read_start; i_cur <= read_end; i_cur++)
				{
					reverse_strand_signal[i_cur]++;
				} // i_cur loop.
			}
		} // strand signal check.

		delete(cur_entry_length_str);
		delete [] cur_read;
	} // file reading loop.
	fclose(f_sorted_reads);

	// Get the length of the data.
	l_data = l_buffer;
	while(l_data > 0)
	{
		if(signal_profile_buffer != NULL && signal_profile_buffer[l_data] > 0.0)
		{
			break;
		}
		else if(forward_strand_signal != NULL && forward_strand_signal[l_data] > 0.0)
		{
			break;
		}
		else
		{
			l_data--;
		}
	} // l_data setting loop.
}

int get_l_signal_per_reads(char* reads_fp, int l_ext_tag)
{
	char strand_char = 0;
	//char cur_fragment[10000];
	char mapping_map_str[10000];
	int chr_index;
	int n_processed_reads = 0;
	int l_profile = 0;
	FILE* f_sorted_reads = open_f(reads_fp, "r");
	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
		}	

		n_processed_reads++;

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int cur_l_profile = chr_index + 1000 + l_ext_tag;
		if(cur_l_profile > l_profile)
		{
			l_profile = cur_l_profile;
		}
	} // file loading loop.

	fclose(f_sorted_reads);

	return(l_profile);
} // get_l_signal_per_reads option.

// In order to ignore extended tag length, set it to non-positive value.
void buffer_per_nucleotide_profile_no_buffer(const char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data)
{
	for(int i = 0; i <= l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initilize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	// Non-positive tag extension lengths are ignored and falls back to using the length in the CIGAR string entry.
	if(l_extended_tag <= 0)
	{
		fprintf(stderr, "Ignoring tag extension.\n");
	}

	char strand_char = 0;
	//char cur_fragment[100000];
	char mapping_map_str[100000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = open_f(sorted_read_fp, "r");

	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
}
		}	

		n_processed_reads++;

		// We need a check on the current read line to make sure what we have is a valid: index must be strictly numbers; mapping map string is validated
		// below.
		char chr_index_str[1000];

		if(sscanf(cur_read, "%s %c %s", mapping_map_str, &strand_char, chr_index_str) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		// Validate the string read as chromosome index.
		int char_i = 0;
		while(chr_index_str[char_i] != 0)
		{
			bool char_is_a_num = (chr_index_str[char_i] >= '0' && chr_index_str[char_i] <= '9');
			if(!char_is_a_num)
			{
				fprintf(stderr, "Chromosome index must be a number: %s\n", cur_read);
				exit(0);
			}

			char_i++;
		} // char_i loop.

		// This may be pointing to an error in read preprocessing.
		chr_index = atoi(chr_index_str);
		if(chr_index > l_buffer)
		{
			fprintf(stderr, "%s: The read mapped out of buffer.\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		int right_most_match_pos = 0;
		int left_most_match_pos = 1000*1000*1000;

		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			if(is_matching)
			{
				// Should we ignore the tag extension?
				if(l_extended_tag <= 0)
				{
					// Increase the height: The indexing is already 1 based, there is no conversion needed.
					if(signal_profile_buffer != NULL)
					{
						for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
						{
							signal_profile_buffer[i_cur]++;
						} // i_cur loop.
					}

					// Update the strand signals, if requested.
					if(forward_strand_signal != NULL)
					{
						if(strand_char == 'F')
						{
							// Update the forward strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								forward_strand_signal[i_cur]++;
							} // i_cur loop.
						}
						else
						{
							// Update the reverse strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								reverse_strand_signal[i_cur]++;
							} // i_cur loop.
						}
					} // strand signal check.
				}
				else // tag extension validation check.
				{
					left_most_match_pos = MIN(left_most_match_pos, chr_index);
					right_most_match_pos = (chr_index + l_cur_entry - 1);
				} // extension length check.
			} // match block check.

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// If tag extension is requested, utilize the leftmost and rightmost matching position for the mapped read.
		if(l_extended_tag > 0)
		{
			int ext_start = 0;
			int ext_end = 0;
			if(strand_char == 'F')
			{
				ext_start = left_most_match_pos;
			}
			else
			{
				//ext_start = (chr_index + l_cur_entry - 1) - (l_extended_tag - 1);
				ext_start = right_most_match_pos - (l_extended_tag - 1);
			}

			// Check for negative starts.
			if(ext_start < 0)
			{
				ext_start = 1;
			}

			ext_end = ext_start + l_extended_tag - 1;

			// Update profiles for the strands.
			if(forward_strand_signal != NULL)
			{
				if(strand_char == 'F')
				{
					// Update the forward strand signal.
					for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
					{
						forward_strand_signal[i_cur]++;
					} // i_cur loop.
				}
				else
				{
					// Update the reverse strand signal.
					for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
					{
						reverse_strand_signal[i_cur]++;
					} // i_cur loop.
				}
			} // strand signal check.

			// Update the profile.
			if(signal_profile_buffer != NULL)
			{
				for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
				{
					signal_profile_buffer[i_cur]++;
				} // i_cur loop.
			} // profile existence check.
		} // signal update check.

		delete [] cur_read;
	} // file reading loop.

	close_f(f_sorted_reads, sorted_read_fp);

	// Get the length of the data.
	l_data = l_buffer;
	while(l_data > 0)
	{
		if(signal_profile_buffer != NULL && signal_profile_buffer[l_data] > 0.0)
		{
			break;
		}
		else if(forward_strand_signal != NULL && forward_strand_signal[l_data] > 0.0)
		{
			break;
		}
		else
		{
			l_data--;
		}
	} // l_data setting loop.
}

enum{VAL, TYPE};
bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced)
{
	int i = 0;

	is_read_spliced = false;

	int state = VAL;
	while(mapping_map_str[i] != 0)
	{		
		if(state == VAL)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 0);
			// MIDNSHPX=
			if(mapping_map_str[i] == 'M' ||
				mapping_map_str[i] == 'I' ||
				mapping_map_str[i] == 'D' ||
				mapping_map_str[i] == 'N' ||
				mapping_map_str[i] == 'S' ||		
				mapping_map_str[i] == 'H' ||
				mapping_map_str[i] == 'P' ||
				mapping_map_str[i] == 'X' ||
				mapping_map_str[i] == '=')
			{
				state = TYPE;

				//if(mapping_map_str[i] != 'M')
				// If the state is N, we assume that this is a spliced read: Page 7, CIGAR definition @ https://samtools.github.io/hts-specs/SAMv1.pdf 
				if (mapping_map_str[i] == 'N')
				{
					is_read_spliced = true;
				}
			}
			else if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				// State is still VAL.
			}
			else
			{
				return(false);
			}
		}
		else if(state == TYPE)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 1);
			// A number is expected.
			if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				state = VAL;
			}
			else
			{
				return(false);
			}
		}

		// Move to next character.
		i++;
	}

	return(true);
}

void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										//t_string* cur_entry_length_str,
										int& l_cur_entry,
										char& entry_type_char)
{	
	// Clean the length string.
	//cur_entry_length_str->empty();
	l_cur_entry = 0;

	// Get the next entry in the cigar string.
	while(mapping_map_str[i_mapp_map] != 0)
	{
		if(mapping_map_str[i_mapp_map] < '0' || mapping_map_str[i_mapp_map] > '9')
		{
			break;
		}
		//cur_entry_length_str->concat_char(mapping_map_str[i_mapp_map]);
		l_cur_entry = l_cur_entry*10 + (int)(mapping_map_str[i_mapp_map]-'0');
		i_mapp_map++;
	}

	is_matching = false;
	if(mapping_map_str[i_mapp_map] == 'M')
	{
		//fprintf(stderr, "Adding matching length of %d\n", l_cur_entry);
		is_matching = true;
	}
	else
	{
		//fprintf(stderr, "Adding some other length of %d\n", l_cur_entry);
	}	

	entry_type_char = mapping_map_str[i_mapp_map];

	// Move over the current entry identifier.
	i_mapp_map++;
}

// This is a generic iterator over the processed reads file, useful for doing online processing of the reads files, for example, when they are too large to load into memory.
void preprocessed_read_file_iterator(char* mapped_reads_fp,
	void (per_read_callback)(char*, char, int, void*), 
	void (per_fragment_callback)(char*, char, int, void*),
	void* per_read_callback_param,
	void* per_fragment_callback_param)
{
	printf("Loading mapped-reads from %s.\n", mapped_reads_fp);

	FILE* f_mrf = open_f(mapped_reads_fp, "r");

	//char cur_fragment[10000];
	char mapping_map_str[10000];
	char strand_char;
	int chr_index;

	// Read and validate the mapped reads in the file.
	//while(fscanf(f_mapped_reads, "%s %s %c %d", cur_fragment, quality_str, &strand_char, &chr_index) == 4)
	while(1)
	{
		char* cur_line = getline(f_mrf);

		if(cur_line == NULL)
		{
			break;
		}

		if(sscanf(cur_line, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse fragment line: %s\n", cur_line);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		//fprintf(stderr, "Processing cigar string: %s\n", quality_str);
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		//int left_posn = chr_index;

		// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
			if(is_matching && per_fragment_callback != NULL)
			{
				// Call the fragment callback.
				per_fragment_callback(mapping_map_str, strand_char, chr_index, per_fragment_callback_param);
			}

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Call the read callback.
		per_read_callback(mapping_map_str, strand_char, chr_index, per_read_callback_param);

		//fprintf(stderr, "%s: %s, %d (%d), %c\n", cur_line, new_mapped_read->mapping_str, new_mapped_read->base_index, new_mapped_read->span, new_mapped_read->strand);
		//getc(stdin);

		delete(cur_entry_length_str);
		delete [] cur_line;
	} // curent fragment data reading loop.

	close_f(f_mrf, mapped_reads_fp);
}

void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments)
{
	int i_mapp_map = 0;
	t_string* cur_entry_length_str = new t_string();
	bool is_matching = false;
	char entry_type_char;
	char strand_char = mapped_read->strand;
	//int chr_index = (mapped_read->strand=='F')?(mapped_read->base_index):(mapped_read->base_index-mapped_read->span+1);
	int chr_index = (mapped_read->base_index);
	char* mapping_map_str = mapped_read->mapping_str;

	//fprintf(stderr, "Processing cigar string: %s (%d, %c)\n", mapping_map_str, chr_index, strand_char);
	bool is_read_spliced = false;
	bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

	// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
	while(mapping_map_str_valid && 
		mapping_map_str[i_mapp_map] != 0)
	{
		int l_cur_entry = 0;
		get_next_entry_per_mapp_map_string(mapping_map_str,
											i_mapp_map, 
											is_matching,
											l_cur_entry,
											entry_type_char);

		// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.		

		//int l_fragment = strlen(cur_fragment);
		if(is_matching)
		{
			//int l_fragment = get_l_fragment_per_cigar(quality_str);
			if(strand_char == 'F')
			{
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;
		
				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
			}
			else if(strand_char == 'R')
			{
				// Allocate and initialize a fragment and add it to the reverse strand fragment list.			
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				//new_fragment->base_index = chr_index + l_cur_entry - 1;
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;

				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
				//rev_strand_frags->push_back(new_fragment);
			} // reverse strand check.
		} // maching check.

		// Update the base for the current entry.
		// Must check whether to update the mapping posn: Update only for D and M entries.
		//if(entry_type_char == 'D' || 
		//	entry_type_char == 'M' ||
		//	entry_type_char == 'N' ||
		//	entry_type_char == 'H')
		if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
		{
			chr_index += l_cur_entry;
		}
	} // mapping map string processing loop.

	delete cur_entry_length_str;
}

// Following should be current with SAM specs.
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char)
{
	if(entry_char == '=' ||
		entry_char == 'X' ||
		entry_char == 'S' ||
		entry_char == 'M' ||
		entry_char == 'I')
	{
		return(true);
	}

	return(false);
}

// Following should be current with SAM specs.
bool check_genome_index_update_per_CIGAR_entry(char entry_char)
{
	if(entry_char == 'D' || 
		entry_char == 'M' ||
		entry_char == 'N' ||
		entry_char == '=' ||
		entry_char == 'X')
	{
		return(true);
	}

	return(false);
}

int read_5p_accessor(void* obj_ptr)
{
	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;

	return(frag_obj_ptr->base_index);
}

