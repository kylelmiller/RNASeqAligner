all: FilterFusionAlignmentsByFilteredFusions filterremappedfusion FilterFusionAlignmentsByFilteredFusions MPSSam2fq collectstats SetUnmappedBitFlag sam2fq recover_fusion_alignments_order SepSamUnmapped swap_dRanger_and_MPS_matched comp_fusiondb_offset sepdRangerfusion sepMPSfusion generate_fusiongene_convert_coordinate_trim generate_fusiongene_convert_coordinate_trim_dRanger junc2bed SeparateNormalFromFusionJunc search_fusion_gene AddFusionStrandConsistent gtf2genetab FilterFusionByNormalPaired load_fusion_chrom_seq_std matchfusion2normal filteroriginalfusion Convert2FusionAlignment generate_combined_sequence DNA2StdRegion search_unmapped_reads reads2unmappedsam junc_db_fusion bsb4 mapsplice_multi_thread bowtie check_reads_format alignmenthandler_multi parseCluster cluster junc_db filterjuncbyROCarguNonCanonical filter_1hits newsam2junc fusionsam2junc_filteranchor_newfmt read_chromo_size find_mate_sam_fq alignmenthandler SepSam RemovePairNo check_index_consistency samtools

OPTFLAGS = -O3

# LDFLAGS = -static -static-libgcc 

CFLAGS += $(OPTFLAGS)

# CFLAGS += $(LDFLAGS)

ERRLOG = 2>log

FilterFusionAlignmentsByFilteredFusions:
	g++ $(CFLAGS) -o bin/FilterFusionAlignmentsByFilteredFusions src/FilterFusionAlignmentsByFilteredFusions/FilterFusionAlignmentsByFilteredFusions.cpp $(ERRLOG)

filterremappedfusion:
	g++ $(CFLAGS) -o bin/filterremappedfusion src/filterremappedfusion/filterremappedfusion.cpp $(ERRLOG)

MPSSam2fq:
	g++ $(CFLAGS) -o bin/MPSSam2fq src/MPSSam2fq/MPSSam2fq.cpp src/MPSSam2fq/JunctionSeed.cpp src/MPSSam2fq/SamRec.cpp src/MPSSam2fq/sharedlib.cpp src/MPSSam2fq/SpliceWay.cpp $(ERRLOG)

collectstats:
	g++ $(CFLAGS) -o bin/collectstats src/collectstats/collectstats.cpp $(ERRLOG)

check_index_consistency:
	g++ $(CFLAGS) -o bin/check_index_consistency src/check_index_consistency/check_index_consistency.cpp $(ERRLOG)

SetUnmappedBitFlag:
	g++ $(CFLAGS) -o bin/SetUnmappedBitFlag src/SetUnmappedBitFlag/SetUnmappedBitFlag.cpp $(ERRLOG)

sam2fq:
	g++ $(CFLAGS) -o bin/sam2fq src/sam2fq/sam2fq.cpp $(ERRLOG)

recover_fusion_alignments_order:
	g++ $(CFLAGS) -o bin/recover_fusion_alignments_order src/recover_fusion_alignments_order/recover_fusion_alignments_order.cpp src/recover_fusion_alignments_order/JunctionSeed.cpp src/recover_fusion_alignments_order/SamRec.cpp src/recover_fusion_alignments_order/sharedlib.cpp src/recover_fusion_alignments_order/SpliceWay.cpp $(ERRLOG)

SepSamUnmapped:
	g++ $(CFLAGS) -o bin/SepSamUnmapped src/SepSamUnmapped/SepSamUnmapped.cpp $(ERRLOG)

swap_dRanger_and_MPS_matched:
	g++ $(CFLAGS) -o bin/swap_dRanger_and_MPS_matched src/swap_dRanger_and_MPS_matched/swap_dRanger_and_MPS_matched.cpp $(ERRLOG)

comp_fusiondb_offset:
	g++ $(CFLAGS) -o bin/comp_fusiondb_offset src/comp_fusiondb_offset/comp_fusiondb_offset.cpp $(ERRLOG)

sepdRangerfusion:
	g++ $(CFLAGS) -o bin/sepdRangerfusion src/sepdRangerfusion/sepdRangerfusion.cpp $(ERRLOG)

sepMPSfusion:
	g++ $(CFLAGS) -o bin/sepMPSfusion src/sepMPSfusion/sepMPSfusion.cpp $(ERRLOG)

generate_fusiongene_convert_coordinate_trim:

	g++ -O3 -o bin/generate_fusiongene_convert_coordinate_trim src/generate_fusiongene_convert_coordinate/generate_fusiongene_convert_coordinate_trim.cpp src/generate_fusiongene_convert_coordinate/JunctionSeed.cpp src/generate_fusiongene_convert_coordinate/SamRec.cpp src/generate_fusiongene_convert_coordinate/sharedlib.cpp src/generate_fusiongene_convert_coordinate/SpliceWay.cpp $(ERRLOG)

generate_fusiongene_convert_coordinate_trim_dRanger:

	g++ -O3 -o bin/generate_fusiongene_convert_coordinate_trim_dRanger src/generate_fusiongene_convert_coordinate/generate_fusiongene_convert_coordinate_trim_dRanger.cpp src/generate_fusiongene_convert_coordinate/JunctionSeed.cpp src/generate_fusiongene_convert_coordinate/SamRec.cpp src/generate_fusiongene_convert_coordinate/sharedlib.cpp src/generate_fusiongene_convert_coordinate/SpliceWay.cpp $(ERRLOG)

junc2bed:
	g++ $(CFLAGS) -o bin/junc2bed src/junc2bed/junc2bed.cpp $(ERRLOG)

SeparateNormalFromFusionJunc:
	g++ $(CFLAGS) -o bin/SeparateNormalFromFusionJunc src/SeparateNormalFromFusionJunc/SeparateNormalFromFusionJunc.cpp $(ERRLOG)

search_fusion_gene:
	g++ $(CFLAGS) -o bin/search_fusion_gene src/search_fusion_gene/search_fusion_gene.cpp $(ERRLOG)

AddFusionStrandConsistent:
	g++ $(CFLAGS) -o bin/AddFusionStrandConsistent src/AddFusionStrandConsistent/AddFusionStrandConsistent.cpp $(ERRLOG)

gtf2genetab:
	g++ $(CFLAGS) -o bin/gtf2genetab src/gtf2genetab/gtf2genetab.cpp $(ERRLOG)

bsb4:
	g++ $(CFLAGS) -o bin/bsb4 src/bsb4/bsb4.cpp $(ERRLOG)

samtools:
	cd ./samtools-0.1.9;make

	cp ./samtools-0.1.9/samtools ./bin/

mapsplice_multi_thread:
	cd ./src/MapSplice;make

	cp ./src/MapSplice/bowtie ./bin/mapsplice_multi_thread

bowtie:
	cd ./src/bowtie; make

	cp ./src/bowtie/bowtie ./bin/bowtie

	cp ./src/bowtie/bowtie-build ./bin/bowtie-build

	cp ./src/bowtie/bowtie-inspect ./bin/bowtie-inspect

reads2unmappedsam:
	g++ $(CFLAGS) -o bin/reads2unmappedsam src/reads2unmappedsam/reads2unmappedsam.cpp $(ERRLOG)

DNA2StdRegion:
	g++ $(CFLAGS) -o bin/DNA2StdRegion src/DNA2StdRegion/DNA2StdRegion.cpp $(ERRLOG)

generate_combined_sequence:
	g++ $(CFLAGS) -o bin/generate_combined_sequence src/generate_combined_sequence/generate_combined_sequence.cpp $(ERRLOG)

Convert2FusionAlignment:
	g++ $(CFLAGS) -o bin/Convert2FusionAlignment src/Convert2FusionAlignment/Convert2FusionAlignment.cpp src/Convert2FusionAlignment/sharedlib.cpp src/Convert2FusionAlignment/SamRec.cpp src/Convert2FusionAlignment/JunctionSeed.cpp $(ERRLOG)

check_reads_format:
	g++ $(CFLAGS) -o bin/check_reads_format src/check_reads_format/check_reads_format.cpp $(ERRLOG)

search_unmapped_reads:
	g++ $(CFLAGS) -o bin/search_unmapped_reads src/search_unmapped_reads/search_unmapped_reads.cpp $(ERRLOG)

parseCluster:
	g++ $(CFLAGS) -o bin/parseCluster src/cluster/preparseFusion.cpp $(ERRLOG)

cluster:
	g++ $(CFLAGS) -o bin/cluster src/cluster/PE_match_junction.cpp $(ERRLOG)

junc_db:
	g++ $(CFLAGS) -o bin/junc_db src/alignmenthandler/junc_db.cpp src/alignmenthandler/JunctionHandler.cpp src/alignmenthandler/SamRec.cpp src/alignmenthandler/sharedlib.cpp src/alignmenthandler/JunctionSeed.cpp $(ERRLOG)

junc_db_fusion:
	g++ $(CFLAGS) -o bin/junc_db_fusion src/junction_database_fusion/junction_seq_construction_fusion.cpp $(ERRLOG)

filterjuncbyROCarguNonCanonical:
	g++ $(CFLAGS) -o bin/filterjuncbyROCarguNonCanonical src/filterjuncbyROCarguNonCanonical/filterjuncbyROCarguNonCanonical.cpp $(ERRLOG)

filter_1hits:
	g++ $(CFLAGS) -o bin/filter_1hits src/filter_1hits/filter_1hits.cpp $(ERRLOG)

newsam2junc:
	g++ $(CFLAGS) -o bin/newsam2junc src/newsam2junc/newsam2junc.cpp $(ERRLOG)

fusionsam2junc_filteranchor_newfmt:
	g++ $(CFLAGS) -o bin/fusionsam2junc_filteranchor_newfmt src/fusionsam2junc_filteranchor_newfmt/fusionsam2junc_filteranchor_newfmt.cpp $(ERRLOG)

read_chromo_size:
	g++ $(CFLAGS) -o bin/read_chromo_size src/read_chromo_size/read_chromo_size.cpp $(ERRLOG)

find_mate_sam_fq:
	g++ $(CFLAGS) -o bin/find_mate_sam_fq src/find_mate_sam_fq/find_mate_sam_fq.cpp $(ERRLOG)

alignmenthandler:
	g++ $(CFLAGS) -o bin/alignment_handler src/alignmenthandler/test_handler.cpp src/alignmenthandler/JunctionHandler.cpp src/alignmenthandler/SamRec.cpp src/alignmenthandler/sharedlib.cpp src/alignmenthandler/AlignmentHandler.cpp src/alignmenthandler/JunctionSeed.cpp src/alignmenthandler/FusionSamRec.cpp src/alignmenthandler/UnionExpressedRegions.cpp src/alignmenthandler/disjointset.cpp $(ERRLOG)

alignmenthandler_multi:
	g++ $(CFLAGS) -o bin/alignment_handler_multi src/alignmenthandler_multi/test_handler.cpp src/alignmenthandler_multi/JunctionHandler.cpp src/alignmenthandler_multi/SamRec.cpp src/alignmenthandler_multi/sharedlib.cpp src/alignmenthandler_multi/AlignmentHandler.cpp src/alignmenthandler_multi/JunctionSeed.cpp src/alignmenthandler_multi/SpliceWay.cpp src/alignmenthandler_multi/UnionExpressedRegions.cpp src/alignmenthandler_multi/disjointset.cpp -lpthread $(ERRLOG)

SepSam:
	g++ $(CFLAGS) -o bin/SepSam src/SepSam/SepSam.cpp $(ERRLOG)

filteroriginalfusion:
	g++ $(CFLAGS) -o bin/filteroriginalfusion src/filteroriginalfusion/filteroriginalfusion.cpp $(ERRLOG)

matchfusion2normal:
	g++ $(CFLAGS) -o bin/matchfusion2normal src/matchfusion2normal/matchfusion2normal.cpp $(ERRLOG)

load_fusion_chrom_seq_std:
	g++ $(CFLAGS) -o bin/load_fusion_chrom_seq_std src/load_fusion_chrom_seq_std/load_fusion_chrom_seq_std.cpp $(ERRLOG)

FilterFusionByNormalPaired:
	g++ $(CFLAGS) -o bin/FilterFusionByNormalPaired src/FilterFusionByNormalPaired/FilterFusionByNormalPaired.cpp $(ERRLOG)

RemovePairNo:
	g++ $(CFLAGS) -o bin/RemovePairNo src/RemovePairNo/remove_pair_no.cpp $(ERRLOG)

clean:
	cd ./samtools-0.1.9 ; make clean-recur
	cd ./src/MapSplice ; make clean
	cd ./src/bowtie ; make clean
	rm -f bin/FilterFusionAlignmentsByFilteredFusions bin/filterremappedfusion bin/MPSSam2fq bin/comp_fusiondb_offset bin/SeparateNormalFromFusionJunc bin/collectstats bin/SetUnmappedBitFlag bin/sam2fq bin/recover_fusion_alignments_order bin/SepSamUnmapped bin/swap_dRanger_and_MPS_matched bin/sepdRangerfusion bin/sepMPSfusion bin/generate_fusiongene_convert_coordinate_trim_dRanger bin/generate_fusiongene_convert_coordinate_trim bin/junc2bed bin/search_fusion_gene bin/AddFusionStrandConsistent bin/gtf2genetab bin/FilterFusionByNormalPaired bin/load_fusion_chrom_seq_std bin/matchfusion2normal bin/filteroriginalfusion bin/DNA2StdRegion bin/generate_combined_sequence bin/Convert2FusionAlignment bin/search_unmapped_reads bin/reads2unmappedsam bin/junc_db_fusion bin/bsb4 bin/fusionsam2junc_filteranchor_newfmt bin/mapsplice_multi_thread bin/bowtie-build bin/bowtie bin/bowtie-inspect bin/check_reads_format bin/samtools bin/alignment_handler_multi bin/parseCluster bin/cluster bin/junc_db bin/filterjuncbyROCarguNonCanonical bin/filter_1hits bin/newsam2junc bin/read_chromo_size bin/find_mate_sam_fq bin/alignment_handler bin/SepSam bin/RemovePairNo bin/check_index_consistency log