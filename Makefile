all: COLLAGENE_Release

CC = g++
exec_name = bin/COLLAGENE_Release
lib_flags = -lz -lpthread -lgsl -lgslcblas
LIB_DIR = src

# Use local github library.
SEAL_INCLUDE_DIR = src/SEAL-main/install/include/SEAL-4.0
SEAL_LIB_DIR = src/SEAL-main/install/lib64/libseal-4.0.a

comp_flags = -c -O3 -Wall -std=c++17 -I${SEAL_INCLUDE_DIR}

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/cllgn_main.o \
${LIB_DIR}/cllgn_seal_genomics_key_sharing_utils.o \
${LIB_DIR}/cllgn_seal_genomics_matrix_utils.o \
${LIB_DIR}/cllgn_LR_model_stats.o \
${LIB_DIR}/cllgn_LR_utils.o \
${LIB_DIR}/cllgn_LR_utils_IRLS.o \
${LIB_DIR}/cllgn_fedLR_utils.o \
${LIB_DIR}/cllgn_fedLR_meta_type_utils.o \
${LIB_DIR}/cllgn_fedLR_secure_convertible_protocol_utils.o \
${LIB_DIR}/cllgn_features_weight_utils.o \
${LIB_DIR}/cllgn_x_chisqr.o \
${LIB_DIR}/cllgn_x_sigmoid.o \
${LIB_DIR}/cllgn_ansi_cli.o \
${LIB_DIR}/cllgn_file_utils.o \
${LIB_DIR}/cllgn_ansi_thread.o \
${LIB_DIR}/cllgn_config.o \
${LIB_DIR}/cllgn_xlog_math.o \
${LIB_DIR}/cllgn_matrix_linalg_utils.o \
${LIB_DIR}/cllgn_ansi_string.o \
${LIB_DIR}/cllgn_exception_obj.o \
${LIB_DIR}/cllgn_annot_region_tools.o \
${LIB_DIR}/cllgn_multicolumn_processing.o \
${LIB_DIR}/cllgn_variation_tools.o \
${LIB_DIR}/cllgn_gff_utils.o \
${LIB_DIR}/cllgn_human_data_processing.o \
${LIB_DIR}/cllgn_signal_track_tools.o \
${LIB_DIR}/cllgn_genome_sequence_tools.o \
${LIB_DIR}/cllgn_mapped_read_tools.o \
${LIB_DIR}/cllgn_nucleotide.o \
${LIB_DIR}/cllgn_histogram.o \
${LIB_DIR}/cllgn_nomenclature.o \
${LIB_DIR}/cllgn_genomics_coords.o \
${LIB_DIR}/cllgn_rng.o \
${LIB_DIR}/cllgn_seed_manager.o 

COLLAGENE_Release: ${objs}
	@echo Linking executable $@
	@${CC} -O3 ${lib_flags} -I ${SEAL_INCLUDE_DIR} -o ${exec_name} ${objs} ${SEAL_LIB_DIR}

clean:
	@echo Cleaning..
	@rm -f ${objs} ${exec_name} 
