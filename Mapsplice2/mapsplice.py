 #    MapSplice
 #
 #    Copyright (C) 2012-2013 University of Kentucky 
 #
 #    Developers: Zheng Zeng and Kai Wang 
 #    PIs : Jinze Liu and Jan F. Prins
 #    This work was supported by 
 #      The US National Science Foundation EF-0850237 to J.L.  and J.F.P. 
 #      The US National Science Foundation Career award IIS-1054631 to J.L.
 #      The US National Institutes of Health R01-HG006272 to J.F.P and J.L.
import sys
import getopt
import subprocess
import errno
import os, glob
import tempfile
import warnings
import shutil
import math  
import locale

from datetime import datetime, date, time

use_message = '''

Usage:
    python mapsplice.py [options] -c <Reference_Sequence>
           -x <Bowtie_Index> -1 <Read_List1> -2 <Read_List2>

Required Arguments:
    -c/--chromosome-dir            <string>     reference sequence directory                       
    -x                             <string>     path and prefix of bowtie index                     
    -1/                            <string>     end 1 reads / single end reads                      
    -2/                            <string>     end 2 reads                                         

Optional Arguments:
  Input/Output and Performance options:
    -p/--threads                   <int>        number of threads, default: 1                       
    -o/--output                    <string>     output directory, default: ./mapsplice_out
    --qual-scale                   <string>     quality scale, phred64(default) or phred33 or solexa64      
    --bam                                       output alignment in BAM format, default: SAM format.           
    --keep-tmp                                  keep intermediate files, default: off     

  Alignment options:               
    -s/--seglen                    <int>        segment length, default: 25                          
    --min-map-len                  <int>        minimum alignment length, default: 50                
    -k/max-hits                    <int>        maximum alignments per read, default: 4              
    -i/--min-intron                <int>        minimum intron length, default: 50         
    -I/--max-intron                <int>        maximum intron length, default: 300,000    
    --non-canonical                             also search for non-canonical junction, default: off      
    --filtering                    <int>        filtering stringency level, default: 2     
    -m/--splice-mis                <int>        Maximum number of mismatches that are allowed in the first/last segment
                                                crossing a splice junction in the range of [0, 2]. Default is 1.
                                                (Maximum number of mismatches that are allowed in the middle segment 
                                                crossing a splice junction is always fixed at 2.)                                 
    --max-append-mis               <int>        Maximum number of mismatches allowed to append a high error exonic segment 
                                                next to an adjacent low error segment. Default is 3.
    --ins                          <int>        maximum insertion length, default: 6, must in range [0, 10]                                                
    --del                          <int>        maximum deletion length, default: 6                                 
    --fusion | --fusion-non-canonical           also search for fusion junction, default: off
    --filtering                    <int>        The stringency level of filtering splice junctions in the range of [1, 2]. 
                                                Default is 2.
                                                1: Less stringent filtering, with higher sensitivity of splice junction detection.
                                                2: Standard filtering.
        
Other Arguments:    
    -h/--help                                   print the usage message                    
    -v/--version                                print the version of MapSplice             

For more detailed manual, please visit MapSplice 2.0 website:
http://www.netlab.uky.edu/p/bioinfo/MapSplice2UserGuide

'''

ver_message = "MapSplice v2.1.6"

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class Ver(Exception):
    def __init__(self, msg):
        self.msg = msg

build_bowtie_index_dir = ""
output_dir = "mapsplice_out/"
logging_dir = output_dir + "logs/"
temp_dir = output_dir + "tmp/"
original_dir = temp_dir + "original/"
filteredbest_dir = temp_dir + "best/"
remap_dir = temp_dir + "remap/"
fusion_dir = temp_dir + "fusion/"
cluster_dir = temp_dir + "cluster/"
cluster_result_dir = cluster_dir + "result/"
cluster_data_dir = cluster_dir + "data/"
cluster_data_parsedPER_dir = cluster_data_dir + "parsedPER/"
formated_chrom_dir = temp_dir + "formated_chrom/"
formated_reads_dir = temp_dir + "formated_reads/"
filter_repeats_dir = temp_dir + "filter_repeats/"
bin_dir = sys.path[0] + "/bin/"
rerun_all = 1
DEBUG = 0
max_chromosome = 1000
fail_str = "[MapSplice Running Failed]\n"
FASTA_file_extension = "fa"
input_reads_1 = ""
input_reads_2 = ""
segment_mismatches = 1
splice_mismatches_double = 2
splice_mismatches_single = 1
min_intron_length = 50
max_intron_length = 300000
read_width = 0
min_read_len = 25
flank_case = 5
fusion_flank_case = 5
chromosome_files_dir = ""
qual_scale = ""
bwt_idx_prefix = ""
bowtie_threads = 1
max_insert = 6
max_delete = 6
is_paired = 0
output_bam = 0
filtering_level = 2
do_filter_fusion_by_repeat = 1
max_hits = 4
format_flag = ""
is_fastq = 0
do_fusion = 0
rm_temp = 1
run_mapper = 0
collect_stat = 1
seg_len = 25
seg_len_overided = False
min_map_len = 50
gene_gtf_file = ""
junction_tobe_syn = ""
append_mismatch = 3
remap_mismatch = 2
boundary = 36
unmapped_reads_sams = []
min_entropy_repeats = 1
min_entropy = -0.0001
low_support_thresh = 1


def get_version():
    return ver_message


def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def check_executables():
    if os.listdir(bin_dir)== []:
        print >> sys.stderr, fail_str, "Error: Could not find executable files in bin directory "
        print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"               
        exit(1)

def prepare_output_dir():
    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)

    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)

    if os.path.exists(temp_dir):
        pass
    else:        
        os.mkdir(temp_dir)

    if os.path.exists(original_dir):
        pass
    else:        
        os.mkdir(original_dir)

    if os.path.exists(filteredbest_dir):
        pass
    else:        
        os.mkdir(filteredbest_dir)

    if os.path.exists(remap_dir):
        pass
    else:        
        os.mkdir(remap_dir)

    if os.path.exists(fusion_dir):
        pass
    else:        
        os.mkdir(fusion_dir)
    
    if os.path.exists(cluster_dir):
        pass
    else:        
        os.mkdir(cluster_dir)

    if os.path.exists(cluster_result_dir):
        pass
    else:        
        os.mkdir(cluster_result_dir)

    if os.path.exists(cluster_data_dir):
        pass
    else:        
        os.mkdir(cluster_data_dir)

    if os.path.exists(cluster_data_parsedPER_dir):
        pass
    else:        
        os.mkdir(cluster_data_parsedPER_dir)

    if os.path.exists(filter_repeats_dir):
        pass
    else:
        os.mkdir(filter_repeats_dir)
        
    if build_bowtie_index_dir == "" or os.path.exists(build_bowtie_index_dir):
        pass
    else:
        os.mkdir(build_bowtie_index_dir)

def read_dir_by_suffix(path_name, suffix):
    suffix = '*.' + suffix
    sequences = []
    flag = 0
    if os.path.isdir(path_name) == False:
        sequences = sequences + [path_name]
        return sequences
    for infile in glob.glob( os.path.join(path_name, suffix) ):
        if flag == 1:
            sequences = sequences + [infile]
        else:
            sequences = sequences + [infile]
            flag = 1
    return sequences


def read_sequence_by_suffix(path_name, suffix):
    suffix = '*.' + suffix
    sequences = ''
    flag = 0
    chromosome_count = 0
    if os.path.isdir(path_name) == False:
        return path_name
    for infile in glob.glob( os.path.join(path_name, suffix) ):
        if flag == 1:
            sequences += ',' + infile
            chromosome_count = chromosome_count + 1
        else:
            sequences = infile
            flag = 1
    """
    if chromosome_count > max_chromosome:
        all_chromos_path = read_dir_by_suffix(chromosome_files_dir, suffix)
        cat_files(all_chromos_path, formated_chrom_dir + "combined_chromosomes.fa")
        return formated_chrom_dir + "combined_chromosomes.fa"
    """
    return sequences


def check_file_existence(reads_files_or_dir):
    print >> sys.stderr, "[%s] Checking files or directory: %s" % (right_now(), reads_files_or_dir)
    if os.path.exists(reads_files_or_dir):
        return reads_files_or_dir
    else:
        print >> sys.stderr, "Error: Could not find files or directory " + reads_files_or_dir
        exit(1)


def call_mapsplice_multithreads(segment_mis,
                                splice_mis_double,
                                splice_mis_single,
                                read_format,
                                maximum_hits,
                                threads,
                                segment_len,
                                chromosome_size,
                                ref_seq_path,
                                flankcase,
                                max_insertion,
                                max_deletion,
                                minimum_map_len,
                                mapsplice_out,
                                bowtie_index,
                                input_reads_file_1,
                                input_reads_file_2,
                                unmapped_reads,
                                bowtie_output,
                                juncdb_index,
                                no_original_junc,
                                do_fusion_and_is_paired,
                                max_intron,
                                append_mis,
                                minimum_read_len,
                                min_intron,
                                quality_scale,
                                log_file):

    start_time = datetime.now()
    print >> sys.stderr, "[%s] Running MapSplice multi-thread" % start_time.strftime("%c")
    if os.path.exists(mapsplice_out + ".1") and rerun_all == 0:
        return mapsplice_out
    splice_cmd = [bin_dir + "mapsplice_multi_thread",
                            read_format,
                            "--min_len", str(minimum_read_len),
                            "--seg_len", str(segment_len),
                            "--min_intron", str(min_intron),
                            "--max_intron_single", str(max_intron),
                            "--max_intron_double", str(max_intron),                            
                            "-v", str(segment_mis),
                            "--max_double_splice_mis", str(splice_mis_double),
                            "--max_single_splice_mis", str(splice_mis_single),                            
                            "--max_append_mis", str(append_mis),
                            "--max_ins", str(max_insertion),
                            "--max_del", str(max_deletion),
                            "-k", str(maximum_hits), 
                            "-m", str(maximum_hits), 
                            "-p", str(threads), 
                            "--chrom_tab", chromosome_size,
                            "--ref_seq_path", ref_seq_path,
                            "--mapsplice_out", mapsplice_out,
                            "--check_read", logging_dir + "check_reads_format.log"]
    if read_format == "-q":
        splice_cmd.append("--qual-scale")
        splice_cmd.append(quality_scale)
    if juncdb_index == "":
        splice_cmd.append("--splice_only")
    else:
        splice_cmd.append("--optimize_repeats")
        if no_original_junc == False:
            splice_cmd.append("--juncdb_index")
            splice_cmd.append(juncdb_index)
    if unmapped_reads != "":
        if do_fusion_and_is_paired:
            splice_cmd.append("--output_unmapped_pe")
        else:
            splice_cmd.append("--output_unmapped")
        splice_cmd.append(unmapped_reads)     
    if flankcase < 5:
        splice_cmd.append("--double_anchor_noncanon")
    #if minimum_map_len >= 0:
    splice_cmd.append("--min_map_len")
    splice_cmd.append(str(minimum_map_len))
    splice_cmd.append(bowtie_index)
    if input_reads_file_2 != "":
        splice_cmd.append("-1")
        splice_cmd.append(input_reads_file_1)
        splice_cmd.append("-2")
        splice_cmd.append(input_reads_file_2)
    else:
        splice_cmd.append(input_reads_file_1)
    splice_cmd.append(bowtie_output)

    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % splice_cmd
    mapsplice_log = open(log_file, "w")
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Running MapSplice multi-thread failed: ", retcode
            exit(1)    
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_multi_thread not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"            
        exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    return mapsplice_out


def call_mapsplice_multithreads_fusion(segment_mis,
                                       splice_mis_double,
                                       splice_mis_single,
                                       read_format,
                                       maximum_hits,
                                       threads,
                                       segment_len,
                                       chromosome_size,
                                       ref_seq_path,
                                       flankcase,
                                       fusion_flankcase,
                                       max_insertion,
                                       max_deletion,
                                       minimum_map_len,
                                       mapsplice_out,
                                       bowtie_index,
                                       input_reads_file_1,
                                       input_reads_file_2,
                                       unmapped_reads,
                                       bowtie_output,
                                       juncdb_index,
                                       fusion_juncdb_index,
                                       cluster,
                                       fusion,
                                       no_original_junc,
                                       no_original_fusion_junc,
                                       max_intron,
                                       append_mis,
                                       minimum_read_len,
                                       min_intron,
                                       quality_scale,
                                       log_file):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Running MapSplice multi-thread fusion" % start_time.strftime("%c")
    if os.path.exists(fusion) and rerun_all == 0:
        return fusion
    splice_cmd = [bin_dir + "mapsplice_multi_thread",
                            read_format,
                            "--min_len", str(minimum_read_len),
                            "--seg_len", str(segment_len), 
                            "--min_intron", str(min_intron),
                            "--max_intron_single", str(max_intron),
                            "--max_intron_double", str(max_intron),
                            "-v", str(segment_mis),
                            "--max_double_splice_mis", str(splice_mis_double),
                            "--max_single_splice_mis", str(splice_mis_single),
                            "--max_append_mis", str(append_mis),
                            "--max_ins", str(max_insertion),
                            "--max_del", str(max_deletion),
                            "-k", str(maximum_hits),
                            "-m", str(maximum_hits), 
                            "-p", str(threads),
                            "--chrom_tab", chromosome_size, 
                            "--ref_seq_path", ref_seq_path,
                            "--mapsplice_out", mapsplice_out,
                            "--check_read", logging_dir + "check_reads_format.log",
                            "--optimize_repeats",
                            "--fusion", fusion]
    if read_format == "-q":
        splice_cmd.append("--qual-scale")
        splice_cmd.append(quality_scale)
    if no_original_junc == False:
        splice_cmd.append("--juncdb_index")    
        splice_cmd.append(juncdb_index)
    if flankcase < 5:
        splice_cmd.append("--double_anchor_noncanon")
    if fusion_flankcase < 5:
        splice_cmd.append("--fusion_double_anchor_noncanon")
        splice_cmd.append("--fusion_single_anchor_noncanon")
    #if minimum_map_len >= 0:
    splice_cmd.append("--min_map_len")
    splice_cmd.append(str(minimum_map_len))
    if fusion_juncdb_index != "":    
        splice_cmd.append("--fusiondb_index")
        splice_cmd.append(fusion_juncdb_index)
    if unmapped_reads != "":    
        splice_cmd.append("--output_unmapped")
        splice_cmd.append(unmapped_reads)
    if input_reads_file_2 != "":
        splice_cmd.append("--cluster")
        splice_cmd.append(cluster)
    splice_cmd.append(bowtie_index)
    if input_reads_file_2 != "":
        splice_cmd.append("-1")
        splice_cmd.append(input_reads_file_1)
        splice_cmd.append("-2")
        splice_cmd.append(input_reads_file_2)
    else:
        splice_cmd.append(input_reads_file_1)
    splice_cmd.append(bowtie_output)
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % splice_cmd
    mapsplice_log = open(log_file, "w")
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Running MapSplice multi-thread fusion failed: ", retcode
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_multi_thread not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1) 
    finish_time = datetime.now()
    duration = finish_time - start_time
    return fusion


def run_alignment_handler_multi(junction_file, 
                   pair_end,
                   max_mate_dist,
                   maximum_hits,
                   filtered_alignment_base,
                   filtered_alignment_append,
                   max_read_length,
                   chromosome_dir,
                   filter_flag,
                   max_insertion,
                   max_deletion,
                   min_anchor,
                   min_junction_anchor,
                   min_mismatch,
                   add_soft_clip,
                   mate_dist_sd,
                   intron_dist_sd,
                   max_anchor_diff,
                   chromosome_size_file,
                   encompassing_fusion_region_extension,
                   number_of_threads,
                   min_coverage,
                   fragment_length,
                   fragment_length_sd,
                   avearge_fragment_length,
                   bound,
                   min_isoform_length,
                   min_encompass_count,
                   minimum_entropy,
                   low_support_threshold,
                   input_sam_file,                   
                   log_file,
                   err_file):
    
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Running alignment handler" % start_time.strftime("%c")
    splice_cmd = ""
    filtered_alignment_file = input_sam_file + filtered_alignment_append
    if os.path.exists(filtered_alignment_base + ".fil.junc") and rerun_all == 0:
        return filtered_alignment_file
    mapsplice_log = open(log_file, "w")
    err_log = open(err_file, "w")
    splice_cmd = [bin_dir + "alignment_handler_multi",
                  junction_file,
                  str(pair_end),
                  str(max_mate_dist),
                  str(maximum_hits),
                  filtered_alignment_base,
                  str(max_read_length),
                  chromosome_dir,
                  str(filter_flag),
                  str(max_insertion),
                  str(max_deletion),
                  str(min_anchor),
                  str(min_junction_anchor),
                  str(min_mismatch),
                  str(add_soft_clip),
                  str(mate_dist_sd),
                  str(max_anchor_diff),
                  chromosome_size_file,
                  str(encompassing_fusion_region_extension), 
                  str(number_of_threads),
                  str(intron_dist_sd),
                  str(min_coverage),
                  str(fragment_length),
                  str(fragment_length_sd),
                  str(avearge_fragment_length),
                  str(bound),
                  str(min_isoform_length),
                  str(min_encompass_count),
                  str(minimum_entropy),
                  str(low_support_threshold),
                  input_sam_file #tmp dir
                  ]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % splice_cmd
    try:
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log, stderr=err_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter alignment failed"
            exit(1)    
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: alignment_handler not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    finish_time = datetime.now()
    duration = finish_time - start_time
    return filtered_alignment_file


def filterbyrepeats(in_fusion,
                    in_normal,
                    out_fusion,
                    chromosomes_dir,
                    segment_len,
                    max_repeat_length,
                    max_anchor_diff,
                    min_seg_len,
                    minimum_entropy,
                    threads,
                    bowidx,
                    chromosome_size_file,
                    logging_dir):
    
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Filtering fusion by repeats" % start_time.strftime("%c")
    if segment_len < 25:
        segment_len = 25
    if os.path.exists(out_fusion) and rerun_all == 0:
        return out_fusion

    ##############match fusion to normal junction#################
    
    match_junc_log = open(logging_dir + "matchfusion2normal.log", "w")
    match_junc_cmd = [bin_dir + "matchfusion2normal",
                      in_normal,
                      in_fusion,
                      filter_repeats_dir + "fusions_matched_normal",
                      str(10)]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % match_junc_cmd
    try:
        retcode = subprocess.call(match_junc_cmd, stdout=match_junc_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: match fusion to normal failed"
            exit(1)    
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: matchfusion2normal not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    ##############load fusion chromosome seq##################
    
    load_fusion_chrom_seq_std_log = open(logging_dir + "load_fusion_chrom_seq_std.log", "w")
    load_fusion_chrom_seq_std_cmd = [bin_dir + "load_fusion_chrom_seq_std",
                      chromosomes_dir,
                      str(segment_len),
                      str(max_anchor_diff),
                      str(min_seg_len),
                      filter_repeats_dir + "fusions_matched_normal.matched_junc",
                      filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq",
                      str(1),
                      str(minimum_entropy)]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % load_fusion_chrom_seq_std_cmd
    try:
        retcode = subprocess.call(load_fusion_chrom_seq_std_cmd, stdout=load_fusion_chrom_seq_std_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: add chromosome sequence failed"
            exit(1)  
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: load_fusion_chrom_seq_std not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    ##############align anchor reads by bowtie##################
    
    anchor_bowtie_log = open(logging_dir + "anchor_bowtie.log", "w")
    anchor_bowtie_err = open(logging_dir + "anchor_bowtie.err", "w")
    readslist = filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq.1.fa" + "," + filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq.2.fa"
    anchor_bowtie_cmd = [bin_dir + "bowtie",
                         "-S",
                         "--sam-nohead",
                         "-f",
                         "--threads", str(threads),
                         "-X", str(50000),
                         "--un", filter_repeats_dir + "unmapped_segments",
                         '-a', #"-k", str(40000),
                         "-m", str(40001),
                         "-v", str(3),
                         "--max", filter_repeats_dir + "repeat_segments",
                         bowidx,
                         readslist,
                         filter_repeats_dir + "repeats.sam"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % anchor_bowtie_cmd
    try:
        retcode = subprocess.call(anchor_bowtie_cmd, stdout=anchor_bowtie_log, stderr=anchor_bowtie_err)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: map repeat reads failed"
            exit(1) 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    ##############sort repeat alignments##################
    
    sort_repeat_log = open(logging_dir + "sort_repeat.log", "w")
    sort_repeat_cmd = ["sort",
                       "-k1,1n",
                       "-t:",
                       "-o", filter_repeats_dir + "repeats_sorted.sam",
                       filter_repeats_dir + "repeats.sam"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % sort_repeat_cmd
    try:
        retcode = subprocess.call(sort_repeat_cmd, stdout=sort_repeat_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort repeat reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    
    ##############load fusion chromosome long seq##################
    
    load_fusion_chrom_seq_std_long_log = open(logging_dir + "load_fusion_chrom_seq_std_long_seq.log", "w")
    load_fusion_chrom_seq_std_long_cmd = [bin_dir + "load_fusion_chrom_seq_std",
                      chromosomes_dir,
                      str(max_repeat_length),
                      str(max_anchor_diff),
                      str(min_seg_len),
                      filter_repeats_dir + "fusions_matched_normal.matched_junc",
                      filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq_long",
                      str(1),
                      str(minimum_entropy),
                      str(1)]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % load_fusion_chrom_seq_std_long_cmd
    try:
        retcode = subprocess.call(load_fusion_chrom_seq_std_long_cmd, stdout=load_fusion_chrom_seq_std_long_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: add chromosome long sequence failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: load_fusion_chrom_seq_std not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)


    ##############align long fusion junction sequence by bowtie##################
    
    long_anchor_bowtie_log = open(logging_dir + "long_fusion_seq_bowtie.log", "w")
    long_anchor_bowtie_err = open(logging_dir + "long_fusion_seq_bowtie.err", "w")
    readslonglist = filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq_long.1.fa" + "," + filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq_long.2.fa"
    long_anchor_bowtie_cmd = [bin_dir + "bowtie",
                         "-f",
                         "--threads", str(threads),
                         "--un", filter_repeats_dir + "unmapped_segments_long",
                         "-k", str(1),
                         "-m", str(1),
                         "-v", str(2),
                         "--max", filter_repeats_dir + "repeat_segments_long",
                         bowidx,
                         readslonglist,
                         filter_repeats_dir + "repeats_long.sam"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % long_anchor_bowtie_cmd
    try:
        retcode = subprocess.call(long_anchor_bowtie_cmd, stdout=long_anchor_bowtie_log, stderr=long_anchor_bowtie_err)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: map long fusion sequence failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    ##############pair repeat reads##################
    
    mapsplice_log = open(logging_dir + "pair_repeat_reads.log", "w")
    err_log = open(logging_dir + "pair_repeat_reads.err", "w")
    splice_cmd = [bin_dir + "alignment_handler_multi",
                  "",
                  str(1), # Anchor length
                  str(50000), # Mismatches allowed in extension
                  str(4000), # Maxmimum intron length
                  filter_repeats_dir + "_filtered_normal_alignments_fix_pair", # Minimum intron length
                  str(segment_len), # Seed size for reads
                  "", # read width for reads
                  str(2348), # islands extension
                  str(6), # block size for reading chromosome
                  str(10), # nothing important
                  str(0), #if is 1, only output flank string that is not case 0
                  str(10), #number of anchors
                  str(5), #number of segments
                  str(1), #segment length
                  str(100), #segment length
                  str(50), #segment length
                  chromosome_size_file, #segment length
                  str(50000), 
                  str(threads),
                  str(500),
                  str(2),
                  str(400),
                  str(100),
                  str(225),
                  str(36),
                  str(100),
                  str(1),
                  str(0),
                  filter_repeats_dir + "repeats_sorted.sam" #tmp dir
                  ]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % splice_cmd
    try:
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log, stderr=err_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pair repeat reads failed"
            exit(1)    
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: alignment_handler not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
       
    ##############sort repeat alignments##################
    
    FilterFusionByNormalPaired_log = open(logging_dir + "FilterFusionByNormalPaired.log", "w")
    FilterFusionByNormalPaired_cmd = [bin_dir + "FilterFusionByNormalPaired",
                                      filter_repeats_dir + "repeats_sorted.sam_filtered_normal_alignments_fix_pair.bothunspliced.paired",
                                      filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq",
                                      out_fusion,
                                      str(1),
                                      filter_repeats_dir + "repeat_segments",
                                      filter_repeats_dir + "repeat_segments_long",
                                      filter_repeats_dir + "fusion_repeats_filtered.txt"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % FilterFusionByNormalPaired_cmd
    try:
        retcode = subprocess.call(FilterFusionByNormalPaired_cmd, stdout=FilterFusionByNormalPaired_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filtering fusion by normalPaired failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterFusionByNormalPaired not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    finish_time = datetime.now()
    duration = finish_time - start_time
    return filter_repeats_dir + "repeats_sorted.sam_filtered_normal_alignments_fix_pair"


def filterbyannotation(in_fusion,
                       out_fusion,
                       out_fusion_not_well_annotated,
                       circular_RNAs,                       
                       gtf_file,
                       logging_dir):
    
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Filtering fusion by annotation" % start_time.strftime("%c")
    if os.path.exists(out_fusion) and rerun_all == 0:
        return out_fusion

    ##############convert gene gtf file to annotation format#################
    
    gtf2genetab_log = open(logging_dir + "gtf2genetab.log", "w")
    gtf2genetab_cmd = [bin_dir + "gtf2genetab",
                       gtf_file,
                       temp_dir + "gene_exons.txt"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % gtf2genetab_cmd
    try:
        retcode = subprocess.call(gtf2genetab_cmd, stdout=gtf2genetab_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert gene gtf file to annotation failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: gtf2genetab not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    ##############annotate fusion junctions##################
    
    annotate_fusion_log = open(logging_dir + "annotate_fusion.log", "w")
    annotate_fusion_cmd = [bin_dir + "search_fusion_gene",
                           "-g", temp_dir + "gene_exons.txt",
                           "-f", in_fusion,
                           "-o", filter_repeats_dir + "filter_repeats.txt.annot"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % annotate_fusion_cmd
    try:
        retcode = subprocess.call(annotate_fusion_cmd, stdout=annotate_fusion_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: annotate fusion junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: search_fusion_gene not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    ##############add fusion strand consistent##################
    
    add_fusion_strand_consistent_log = open(logging_dir + "add_fusion_strand_consistent.log", "w")
    add_fusion_strand_consistent_err = open(logging_dir + "add_fusion_strand_consistent.err", "w")
    add_fusion_strand_consistent_cmd = [bin_dir + "AddFusionStrandConsistent",
                                        filter_repeats_dir + "filter_repeats.txt.annot",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % add_fusion_strand_consistent_cmd
    try:
        retcode = subprocess.call(add_fusion_strand_consistent_cmd, stdout=add_fusion_strand_consistent_log, stderr=add_fusion_strand_consistent_err)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: add fusion strand consistent failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: AddFusionStrandConsistent not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)


    ##############separate matched strand fusion from fusion junction##################
    
    SeparateNormalFromFusionJunc_log = open(logging_dir + "SeparateNormalFromFusionJunc.log", "w")
    SeparateNormalFromFusionJunc_cmd = [bin_dir + "SeparateNormalFromFusionJunc",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.match",
                                        out_fusion,#output_dir + "fusion.from.fusion",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.match.normal",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.match.circularRNAs"]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % SeparateNormalFromFusionJunc_cmd
    try:
        retcode = subprocess.call(SeparateNormalFromFusionJunc_cmd, stdout=SeparateNormalFromFusionJunc_log)     
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate matched strand fusion from fusion junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SeparateNormalFromFusionJunc not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)


    ##############separate normal junction and circular RNAs from fusion junction##################
    
    SeparateNormalFromFusionJunc_log = open(logging_dir + "SeparateNormalFromFusionJunc.log", "w")
    SeparateNormalFromFusionJunc_cmd = [bin_dir + "SeparateNormalFromFusionJunc",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.notmatch",
                                        out_fusion_not_well_annotated,#output_dir + "fusion.from.fusion",
                                        filter_repeats_dir + "normal.from.fusion",
                                        circular_RNAs]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % SeparateNormalFromFusionJunc_cmd
    try:
        retcode = subprocess.call(SeparateNormalFromFusionJunc_cmd, stdout=SeparateNormalFromFusionJunc_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate normal junction and circular RNAs from fusion junction failed"
            exit(1)    
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SeparateNormalFromFusionJunc not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
        
    finish_time = datetime.now()
    duration = finish_time - start_time
    return output_dir + "fusion.from.fusion"


def extract_maxlen(log_file):
    fh = open(log_file,"r")
    igot = fh.readlines()
    maxlen = 0
    num_read = 0
    low_support_junction_threshold = 1
    read_format = ""
    quality_scale = ""
    for line in igot:
        if line[0] == '\n' or line[0] == '#':
            comments = line
        else:
            line = line.rstrip()
            command = line.split(' ')
            if command[0] == "maxlen":
                maxlen = int(command[1])
            if command[0] == "read_format":
                read_format = command[1]
            if command[0] == "qual_scale":
                quality_scale = command[1]
            if command[0] == "num_read":
                num_read = int(command[1])
                low_support_junction_threshold = int(math.floor(math.log(max(2, (num_read / 10000000)), 2)))
    if maxlen != 0 and read_format != "":
        return (maxlen, read_format, quality_scale, low_support_junction_threshold)
    else:
        print >> sys.stderr, fail_str, "Error: No max read length or read_format found"
    exit(1)


def parseCluster(merge_pair_sam, output_dir, log_file):
    print >> sys.stderr, "[%s] Parsing cluster regions" % (right_now())
    resortmapreads_cmd = ["sed", 
                         "-e",
                         "s/^/1~/",
                         merge_pair_sam]
    sed_log = open(fusion_dir + "sed.fusion.sam", "w")
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=sed_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sed failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sed not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    
    mapsplice_log = open(log_file, "w")
    resortmapreads_cmd = [bin_dir + "parseCluster", 
                         fusion_dir + "sed.fusion.sam",
                         output_dir]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Parseing cluster regions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: parseCluster not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)


def cluster(cluster_dir, log_file):
    print >> sys.stderr, "[%s] Clustering regions" % (right_now())
    mapsplice_log = open(log_file, "w")
    if os.path.exists(cluster_dir + "result/cluster.txt") and rerun_all == 0:
        return cluster_dir + "result/cluster.txt"
    resortmapreads_cmd = [bin_dir + "cluster", 
                         cluster_dir]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Clustering regions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cluster not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)


def sam2juncarray(sam_file_array,
                  converted_junc,
                  chromosomes_file_dir,
                  max_read_width,
                  min_intron,
                  max_intron,
                  cur_output_dir,
                  min_anchor_length,
                  log_file):
    print >> sys.stderr, "[%s] Generating junctions from sam file" % (right_now())
    output_dir_converted_junc = cur_output_dir + converted_junc
    sam2junc_cmd = [bin_dir + "newsam2junc", 
                    output_dir_converted_junc,
                    str(max_read_width),
                    chromosomes_file_dir,
                    str(min_intron),str(max_intron),
                    str(min_anchor_length)]
    for sam_file in sam_file_array:
        sam2junc_cmd = sam2junc_cmd + [sam_file]
    mapsplice_log = open(log_file, "w")
    if os.path.exists(output_dir_converted_junc) and rerun_all == 0:
        if os.path.exists(output_dir_converted_junc + ".sepdir"):
            shutil.rmtree(output_dir_converted_junc + ".sepdir")
        return output_dir_converted_junc
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % sam2junc_cmd
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2junc not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)

    if os.path.exists(output_dir_converted_junc + ".sepdir"):
        shutil.rmtree(output_dir_converted_junc + ".sepdir")
    return output_dir_converted_junc


def filter_junc_by_min_mis_lpq(junction_file,
                               remained_junc,
                               filtered_out_junc,
                               min_mismatch,
                               min_lpq,
                               log_file):
    print >> sys.stderr, "[%s] Filtering junction by min mis and min lpq" % (right_now())
    if os.path.exists(filtered_out_junc) and os.path.exists(remained_junc) and rerun_all == 0:
        return (remained_junc, filtered_out_junc)
    mapsplice_log = open(log_file, "w")
    if min_lpq > 0.5:
        min_lpq = 0.5
    sam2junc_cmd = [bin_dir + "filter_1hits",
                    junction_file,
                    str(min_mismatch),
                    str(min_lpq),
                    remained_junc,
                    filtered_out_junc]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % sam2junc_cmd
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junction by min mis and min lpq failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filter_1hits not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return (remained_junc, filtered_out_junc)


def filteroriginalfusion(junction_file, remained_junc, filtered_out_junc, min_mismatch, min_lpq, log_file):
    print >> sys.stderr, "[%s] Filtering original fusion junction" % (right_now())
    if os.path.exists(filtered_out_junc) and os.path.exists(remained_junc) and rerun_all == 0:
        return (remained_junc, filtered_out_junc)
    mapsplice_log = open(log_file, "w")
    if min_lpq > 0.5:
        min_lpq = 0.5
    sam2junc_cmd = [bin_dir + "filteroriginalfusion",
                    junction_file,
                    str(min_mismatch),
                    str(min_lpq),
                    remained_junc,
                    filtered_out_junc]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % sam2junc_cmd
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
        if retcode == 100:
            print >> sys.stderr, "Warning: No original fusions found, skip build index step"
            return ("", "")
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter original fusion junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filteroriginalfusion not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return (remained_junc, filtered_out_junc)


def filterremappedfusion(junction_file,
                         remained_junc,
                         filtered_out_junc,
                         min_mismatch, min_hits,
                         min_lpq, minimum_entropy,
                         log_file):
    print >> sys.stderr, "[%s] Filtering remapped fusion junction" % (right_now())
    if os.path.exists(filtered_out_junc) and os.path.exists(remained_junc) and rerun_all == 0:
        return (remained_junc, filtered_out_junc)
    mapsplice_log = open(log_file, "w")
    if min_lpq > 0.5:
        min_lpq = 0.5
    sam2junc_cmd = [bin_dir + "filterremappedfusion",
                    junction_file,
                    str(min_mismatch),
                    str(min_hits),
                    str(min_lpq),
                    str(minimum_entropy),
                    remained_junc,
                    filtered_out_junc]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % sam2junc_cmd
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
        if retcode == 100:
            print >> sys.stderr, "Warning: No remapped fusions found"
            return ("", "")
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter remapped fusion junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterremappedfusion not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return (remained_junc, filtered_out_junc)


def filterjuncbyROCarguNoncanon(junction_file, remained_junc, filtered_out_junc, entropy_weight, lpq_weight, ave_mis_weight, min_score, log_file):
    print >> sys.stderr, "[%s] Filtering junction by ROC argu noncanonical" % (right_now())
    if os.path.exists(filtered_out_junc) and os.path.exists(remained_junc) and rerun_all == 0:
        return (remained_junc, filtered_out_junc)
    mapsplice_log = open(log_file, "w")
    sam2junc_cmd = [bin_dir + "filterjuncbyROCarguNonCanonical",
                    junction_file,
                    str(entropy_weight),
                    str(lpq_weight),
                    str(ave_mis_weight),
                    '0',
                    '0',
                    str(min_score),
                    '5',
                    remained_junc,
                    filtered_out_junc]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
        if retcode == 100:
            print >> sys.stderr, "Waring: No original junctions found, skip build index step"
            return ("", "")
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junc by ROC argu failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterjuncbyROCargu not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return (remained_junc, filtered_out_junc)


def fusionsam2junc_filteranchor_newfmt(mps_unique_mapped,
                                       converted_junc, 
                                       max_read_width,
                                       cur_output_dir,
                                       min_anchor, 
                                       chromosomes_file_dir,
                                       log_file):
    print >> sys.stderr, "[%s] Generating fusion junctions from sam file and filter by anchor" % (right_now())
    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped
    output_dir_converted_junc = cur_output_dir + converted_junc
    if os.path.exists(output_dir_converted_junc) and rerun_all == 0:
        return output_dir_converted_junc
    mapsplice_log = open(log_file, "w")
    sam2junc_cmd = [bin_dir + "fusionsam2junc_filteranchor_newfmt", 
                    output_dir_converted_junc,
                    str(max_read_width),
                    str(min_anchor),
                    chromosomes_file_dir,
                    output_dir_mps_unique_mapped]    
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % sam2junc_cmd
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert new sam file to junctions and filter by anchor failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fusionsam2junc_filteranchor_newfmt not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return output_dir_converted_junc


def build_chromosome_index(user_splice_fasta, user_splices_out_prefix, build_target):
    print >> sys.stderr, "[%s] Building Bowtie index for %s" % (right_now(), build_target)
    bowtie_build_log = open(logging_dir + "bowtie_build.log", "w")
    bowtie_build_cmd = [bin_dir + "bowtie-build", 
                        user_splice_fasta,
                        user_splices_out_prefix]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] bowtie-build commandline" % (bowtie_build_cmd)
    try:    
        retcode = subprocess.call(bowtie_build_cmd, stdout=bowtie_build_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Splice sequence indexing failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-build not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return user_splices_out_prefix


def check_bowtie_index(idx_prefix, chromo_dir, rerun, file_extension, build_target):
    print >> sys.stderr, "[%s] Checking Bowtie index files" % right_now()
    idx_fwd_1 = idx_prefix + ".1.ebwt"
    idx_fwd_2 = idx_prefix + ".2.ebwt"
    idx_fwd_3 = idx_prefix + ".3.ebwt"
    idx_fwd_4 = idx_prefix + ".4.ebwt"
    idx_rev_1 = idx_prefix + ".rev.1.ebwt"
    idx_rev_2 = idx_prefix + ".rev.2.ebwt"
    if os.path.exists(idx_fwd_1) and \
       os.path.exists(idx_fwd_2) and \
       os.path.exists(idx_fwd_3) and \
       os.path.exists(idx_fwd_4) and \
       os.path.exists(idx_rev_1) and \
       os.path.exists(idx_rev_2) and \
       (rerun == 0 or build_target == "reference sequence"):
        return 
    else:
        #print >> sys.stderr, "Warning: Could not find Bowtie index specified by -x: " + idx_prefix + ".*"
        chromo_sequences = read_sequence_by_suffix(chromo_dir, file_extension)
        idx_prefix = build_chromosome_index(chromo_sequences, idx_prefix, build_target)
        idx_fwd_1 = idx_prefix + ".1.ebwt"
        idx_fwd_2 = idx_prefix + ".2.ebwt"
        idx_fwd_3 = idx_prefix + ".3.ebwt"
        idx_fwd_4 = idx_prefix + ".4.ebwt"
        idx_rev_1 = idx_prefix + ".rev.1.ebwt"
        idx_rev_2 = idx_prefix + ".rev.2.ebwt"
        if os.path.exists(idx_fwd_1) and \
           os.path.exists(idx_fwd_2) and \
           os.path.exists(idx_fwd_3) and \
           os.path.exists(idx_fwd_4) and \
           os.path.exists(idx_rev_1) and \
           os.path.exists(idx_rev_2):
            return 
        else:
            print >> sys.stderr, "Error: Building Bowtie index files failed " + idx_prefix + ".*"
            exit(1)

def inspect_bowtie_index(idx_prefix, bowtie_refnames):
    print >> sys.stderr, "[%s] Inspecting Bowtie index files" % right_now()
    bowtie_inspect_cmd = [bin_dir + "bowtie-inspect", 
                          "-s",
                          idx_prefix]   
    bowtie_refnames_fs = open(bowtie_refnames, "w")
    if DEBUG == 1:
        print >> sys.stderr, "[%s] bowtie-inspect command" % (bowtie_inspect_cmd)
    try:    
        retcode = subprocess.call(bowtie_inspect_cmd, stdout=bowtie_refnames_fs)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Inspect bowtie index failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-inspect not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return

def check_index_consistency(bowtie_refnames, refseq_refnames):
    print >> sys.stderr, "[%s] Checking consistency of Bowtie index and reference sequence" % right_now()
    check_index_cmd = [bin_dir + "check_index_consistency", 
                       bowtie_refnames,
                       refseq_refnames]   
    if DEBUG == 1:
        print >> sys.stderr, "[%s] check_index_consistency command" % (check_index_cmd)
    try:    
        retcode = subprocess.call(check_index_cmd)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Checking consistency of Bowtie index and reference sequence failed"
            print >> sys.stderr, "Please check if Bowtie Index and Reference Sequence parameters are set correctly and they comply with MapSplice requirements"
            print >> sys.stderr, "Visit MapSplice 2.0 online manual for details:"
            print >> sys.stderr, "http://www.netlab.uky.edu/p/bioinfo/MapSplice2UserGuide#CommandLine"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: check_index_consistency not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return

def cat_files(files_tobe_cat, cated_file):
    if os.path.exists(cated_file) and rerun_all == 0:
        return cated_file
    bowtie2sam_cmd = ["cat"]
    for sam_file in files_tobe_cat:
        bowtie2sam_cmd = bowtie2sam_cmd + [sam_file]
    cated_file_fs = open(cated_file, "w")
    if DEBUG == 1:
        print >> sys.stderr, "[%s] > %s combine to file" % (bowtie2sam_cmd, cated_file)
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=cated_file_fs)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: combine to final sam"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cat not found on this system"
        exit(1)
    return cated_file


def cat_multithread_files(files_tobe_cat, start_thread, end_thread, cated_file, remove_original_file):
    combine_file = []    
    for a in range(start_thread, end_thread):
        combine_file = combine_file + [files_tobe_cat + "." + str(a)]
    cat_files(combine_file, cated_file)
    #if remove_original_file > 0:
    #    remove_files(combine_file)
    return cated_file
    
    
def remove_files(files_tobe_remove):
    if DEBUG > 0:
        return
    bowtie2sam_cmd = ["rm"]
    for sam_file in files_tobe_remove:
        bowtie2sam_cmd = bowtie2sam_cmd + [sam_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] remove files" % bowtie2sam_cmd
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: remove files failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: rm not found on this system"
        exit(1)


def copy_file(files_tobe_copy, copied_file):
    if os.path.exists(copied_file) and rerun_all == 0:
        return copied_file
    bowtie2sam_cmd = ["cp",
                      files_tobe_copy,
                      copied_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] > %s copy file" % (bowtie2sam_cmd, copied_file)
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: copy file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cp not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return copied_file
    
def remove_pair_no(input_file, output_file):
    print >> sys.stderr, "[%s] Formatting SAM file" % (right_now())
    if os.path.exists(output_file) and rerun_all == 0:
        return output_file
    bowtie2sam_cmd = [bin_dir + "RemovePairNo",
                      input_file,
                      output_file]
    if DEBUG == 1:
        print >> sys.stderr, "Formatting SAM file [%s]" % (bowtie2sam_cmd)
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Formatting SAM file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: RemovePairNo not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return output_file

def sort_by_name_c(files_tobe_sort, sorted_file):
    print >> sys.stderr, "[%s] Sorting file" % (right_now())
    if os.path.exists(sorted_file) and rerun_all == 0:
        return sorted_file
    env = os.environ.copy()
    env.update({"LC_ALL": "C"})
    bowtie2sam_cmd = ["sort",
                      "-k1,1",
                      "-S", "3500000",
                      "-o", sorted_file, 
                      "-T", temp_dir]
    for sam_file in files_tobe_sort:
        bowtie2sam_cmd = bowtie2sam_cmd + [sam_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % bowtie2sam_cmd
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, env = env)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Sorting file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return sorted_file


def check_reads_format(reads_files, minimum_read_len, pair_end):
    print >> sys.stderr, "[%s] Checking read format" % (right_now())
    mapsplice_log = open(logging_dir + "check_reads_format.log", "a")
    read_files_array = reads_files.split(',')
    dividereads_cmd = [bin_dir + "check_reads_format"]
    for read_file in read_files_array:
        dividereads_cmd = dividereads_cmd + [read_file]
    dividereads_cmd = dividereads_cmd + [str(0)] + [str(pair_end)] + [str(minimum_read_len)] 
    if DEBUG == 1:
        print >> sys.stderr, "[%s] Checking read format" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Checking read format failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: check_reads_format not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)


def read_chromo_sizes(all_sams_path, chromo_size_file, chrom_names, chrom_head, chromo_fai_file):
    print >> sys.stderr, "[%s] Checking reference sequence length" % (right_now())
    if os.path.exists(chromo_size_file) and os.path.exists(chrom_names) and rerun_all == 0:
        return chromo_size_file
    mapsplice_log = open(logging_dir + "read_chromo_sizes.log", "w")
    merge_sam_cmd = [bin_dir + "read_chromo_size"]  
    for sam_file in all_sams_path:
        merge_sam_cmd = merge_sam_cmd + [sam_file]
    merge_sam_cmd = merge_sam_cmd + [chromo_fai_file] + [chrom_head] + [chromo_size_file] + [chrom_names]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] Checking reference sequence length" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Checking reference sequence length failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: read_chromo_size not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return chromo_size_file


def fqfa2sam(samfile, reads_files, flag, log_file):
    print >> sys.stderr, "[%s] Converting unmapped reads to sam" % (right_now())
    if os.path.exists(samfile) and rerun_all == 0:
        return samfile    
    mapsplice_log = open(log_file, "w")
    merge_sam_cmd = [bin_dir + "reads2unmappedsam"]  
    merge_sam_cmd = merge_sam_cmd + [samfile]
    merge_sam_cmd = merge_sam_cmd + [str(flag)]
    for read_file in reads_files:
        merge_sam_cmd = merge_sam_cmd + [read_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert unmapped reads to sam failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: reads2unmappedsam not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return samfile


def collectstats(final_stats, 
                 normal_stats, 
                 fusion_stats, 
                 unmapped_stats,
                 pair_end,
                 do_fusion_flag,
                 do_repeat,
                 do_annotate,
                 candidate_fusion_file,
                 well_annotated_fusion_file,
                 not_well_annotated_fusion_file,
                 circular_RNAS_file,
                 unmapped_sams,                 
                 log_file):
    print >> sys.stderr, "[%s] Collecting stats of read alignments and junctions" % (right_now())
    mapsplice_log = open(log_file, "w")
    merge_sam_cmd = [bin_dir + "collectstats"]  
    merge_sam_cmd = merge_sam_cmd + [final_stats] + [normal_stats] + [fusion_stats] + [unmapped_stats] 
    merge_sam_cmd = merge_sam_cmd + [str(pair_end)] + [str(do_fusion_flag)] + [str(do_repeat)] + [str(do_annotate)]
    merge_sam_cmd = merge_sam_cmd + [candidate_fusion_file] + [well_annotated_fusion_file] + [not_well_annotated_fusion_file] + [circular_RNAS_file]
    for sam_file in unmapped_sams:
        merge_sam_cmd = merge_sam_cmd + [sam_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Collect stats of reads alignments and junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: collectstats not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return final_stats


def SetUnmappedBitFlag(unmappedsamfile, unmappedsetbit, log_file):
    print >> sys.stderr, "[%s] Setting unmapped paired end reads bit flag" % (right_now())
    if os.path.exists(unmappedsetbit) and rerun_all == 0:
        return unmappedsetbit    
    mapsplice_log = open(log_file, "w")
    merge_sam_cmd = [bin_dir + "SetUnmappedBitFlag"]  
    merge_sam_cmd = merge_sam_cmd + [unmappedsamfile]
    merge_sam_cmd = merge_sam_cmd + [unmappedsetbit]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Set unmapped paired end reads bit flag failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SetUnmappedBitFlag not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return unmappedsetbit


def sam2bam(bowtie_mapped_sam, bam_file, allchromos_file_fai):
    print >> sys.stderr, "[%s] Converting bowtie sam file to bam format" % (right_now())
    bowtie_mapped_bam = bam_file
    if os.path.exists(bowtie_mapped_bam) and rerun_all == 0:
        return bowtie_mapped_bam
    sam2import_cmd = [bin_dir + "samtools", "view", "-S", "-b", "-o",
                      bowtie_mapped_bam, bowtie_mapped_sam]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] " % sam2import_cmd
    try:    
        retcode = subprocess.call(sam2import_cmd)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: import sam to bam failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: samtools not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return bowtie_mapped_bam


def FilterFusionSamByFusionJunc(all_sams_path, junc_bed, remained_sam, filtered_sam, log_file):
    print >> sys.stderr, "[%s] Filtering fusion alignments by filtered fusion junctions" % (right_now())
    if os.path.exists(remained_sam) and os.path.exists(filtered_sam) and rerun_all == 0:
        return (remained_sam, filtered_sam)
    mapsplice_log = open(log_file, "w")
    merge_sam_cmd = [bin_dir + "FilterFusionAlignmentsByFilteredFusions"] + [junc_bed] + [remained_sam] + [filtered_sam]
    for sam_file in all_sams_path:
        merge_sam_cmd = merge_sam_cmd + [sam_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s] FilterFusionAlignmentsByFilteredFusions" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filtering fusion alignments By filtered fusion junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterFusionAlignmentsByFilteredFusions not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
            
        exit(1)
    return (remained_sam, filtered_sam)


def junc_db(junc_file, chrom_files_dir, min_anchor_width, max_anchor, max_threshold, synthetic_file, log_file):
    print >> sys.stderr, "[%s] Generating synthetic junction sequences" % (right_now())
    if os.path.exists(synthetic_file) and rerun_all == 0:
        return synthetic_file
    chrom_files_dir = chrom_files_dir + "/"
    syn_junc_cmd = [bin_dir + "junc_db", str(min_anchor_width), str(max_anchor), str(max_threshold),
                    junc_file, chrom_files_dir, synthetic_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % syn_junc_cmd
    mapsplice_log = open(log_file, "w")
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Generating synthetic junction sequences failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: junc_db not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return synthetic_file


def fusion_junc_db(junc_file, fusion_junc_file, chrom_files_dir, min_anchor_width, 
                   max_anchor, max_threshold_each, max_threshold_total, synthetic_file, log_file):
    print >> sys.stderr, "[%s] Synthetic fusion junctions sequence" % (right_now())
    if os.path.exists(synthetic_file) and rerun_all == 0:
        return synthetic_file
    chrom_files_dir = chrom_files_dir + "/"
    syn_junc_cmd = [bin_dir + "junc_db_fusion",
                    str(min_anchor_width),
                    str(max_anchor),
                    str(max_threshold_each),
                    str(max_threshold_total),
                    junc_file,
                    "1",
                    fusion_junc_file,
                    "1",
                    chrom_files_dir,
                    synthetic_file]
    if DEBUG == 1:
        print >> sys.stderr, "[%s]" % syn_junc_cmd
    mapsplice_log = open(log_file, "w")
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=mapsplice_log)
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Synthetic fusion junctions sequence failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: junc_db_fusion not found on this system"
            print >> sys.stderr, "Please re-build MapSplice by running 'make' in the MapSplice directory"
        exit(1)
    return synthetic_file


def print_arguments(argu_file):

    #print >> sys.stderr, "Printing argugment to log"
    
    """
    argu_log = open(argu_file, "w")
            
    print >> argu_log, "output_dir=[%s]" % (output_dir)
        
    print >> argu_log, "splice_mismatches=[%s]" % (splice_mismatches)
        
    print >> argu_log, "FASTA_file_extension=[%s]" % (FASTA_file_extension)
        
    print >> argu_log, "min_intron_length=[%s]" % (min_intron_length)
        
    print >> argu_log, "max_intron_length=[%s]" % (max_intron_length)
            
    print >> argu_log, "read_width=[%s]" % (read_width)
        
    print >> argu_log, "flank_case=[%s]" % (flank_case)
        
    print >> argu_log, "fusion_flank_case=[%s]" % (fusion_flank_case)
        
    print >> argu_log, "chromosome_files_dir=[%s]" % (chromosome_files_dir)
            
    print >> argu_log, "bwt_idx_prefix=[%s]" % (bwt_idx_prefix)
            
    print >> argu_log, "bowtie_threads=[%s]" % (bowtie_threads)
            
    print >> argu_log, "max_hits=[%s]" % (max_hits)
            
    print >> argu_log, "boundary=[%s]" % (boundary)
            
    print >> argu_log, "unmapped_reads=[%s]" % (unmapped_reads)
        
    print >> argu_log, "output_bam =[%s]" % (output_bam)
            
    print >> argu_log, "is_paired=[%s]" % (is_paired)
        
    print >> argu_log, "seg_len=[%s]" % (seg_len)

    print >> argu_log, "format_flag=[%s]" % (format_flag)
            
    print >> argu_log, "append_mismatch=[%s]" % (append_mismatch)
            
    print >> argu_log, "collect_stat=[%s]" % (collect_stat)
            
    print >> argu_log, "rm_temp=[%s]" % (rm_temp)
            
    print >> argu_log, "do_fusion=[%s]" % (do_fusion)
            
    print >> argu_log, "run_mapper=[%s]" % (run_mapper)
            
    print >> argu_log, "max_insert=[%s]" % (max_insert)

    print >> argu_log, "max_delete=[%s]" % (max_delete)
            
    print >> argu_log, "minimum_map_len=[%s]" % (min_map_len)
            
    print >> argu_log, "do_filter_fusion_by_repeat=[%s]" % (do_filter_fusion_by_repeat)
            
    print >> argu_log, "rerun_all=[%s]" % (rerun_all)
            
    print >> argu_log, "DEBUG=[%s]" % (DEBUG)
            
    print >> argu_log, "input_reads_1=[%s]" % (input_reads_1)
            
    print >> argu_log, "input_reads_2=[%s]" % (input_reads_2)
            
    argu_log.close()
    """




def main(argv=None):
    if argv is None: #Allow main function to be called from the interactive Python prompt
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hvo:s:m:x:p:c:i:I:1:2:k:", 
                                        ["version",
                                         "help",  
                                         "splice-mis=",
                                         "segment-mismatches=",
                                         "min-intron=",
                                         "max-intron=",
                                         "non-canonical",
                                         "fusion",
                                         "fusion-non-canonical",
                                         "bowtie-index=",
                                         "threads=",
                                         "chromosome-dir=",
                                         "seglen=",
                                         "output=",
                                         "max-hits=",
                                         "filtering=",
                                         "keep-tmp",
                                         "not-rerun-all",
                                         "DEBUG",
                                         "run-MapPER",
                                         "ins=",
                                         "del=",
                                         "bam",
                                         "min-map-len=",
                                         "min-len=",
                                         "min-entropy=",
                                         "min-fusion-entropy=",
                                         "qual-scale=",
                                         "filter-fusion-by-repeat=",
                                         "gene-gtf=",
                                         "max-append-mis="])
        except getopt.error, msg:
            raise Usage(msg)    
        if len(opts) == 0:
            raise Usage(use_message)
        
        check_executables()
                
        global rerun_all
        global DEBUG
        global output_dir
        global logging_dir
        global temp_dir
        global original_dir
        global filteredbest_dir
        global remap_dir
        global formated_chrom_dir
        global formated_reads_dir
        global fusion_dir
        global cluster_dir
        global cluster_result_dir
        global cluster_data_dir
        global cluster_data_parsedPER_dir
        global fusion_result_fusionRead_dir
        global filter_repeats_dir
        global bin_dir
        global chromosome_files_dir
        global splice_mismatches_double
        global splice_mismatches_single
        global FASTA_file_extension
        global min_intron_length
        global max_intron_length
        global read_width
        global min_read_len
        global flank_case
        global fusion_flank_case
        global qual_scale
        global bwt_idx_prefix
        global bowtie_threads
        global max_insert
        global max_delete
        global is_paired
        global input_reads_1
        global input_reads_2
        global output_bam
        global do_filter_fusion_by_repeat
        global max_hits
        global filtering_level
        global format_flag
        global is_fastq
        global do_fusion
        global rm_temp
        global run_mapper
        global collect_stat
        global seg_len
        global min_map_len
        global gene_gtf_file
        global junction_tobe_syn
        global append_mismatch
        global remap_mismatch
        global boundary
        global unmapped_reads_sams
        global min_entropy_repeats
        global min_entropy
        global build_bowtie_index_dir
        global seg_len_overided
        
        # option processing
        for option, value in opts:
            if option in ("-o", "--output"):
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
                temp_dir = output_dir + "tmp/"
                original_dir = temp_dir + "original/"
                filteredbest_dir = temp_dir + "best/"
                remap_dir = temp_dir + "remap/"
                fusion_dir = temp_dir + "fusion/"
                cluster_dir = temp_dir + "cluster/"
                cluster_result_dir = cluster_dir + "result/"
                cluster_data_dir = cluster_dir + "data/"
                cluster_data_parsedPER_dir = cluster_data_dir + "parsedPER/"              
                formated_chrom_dir = temp_dir + "formated_chrom/"
                formated_reads_dir = temp_dir + "formated_reads/"
                filter_repeats_dir = temp_dir + "filter_repeats/"
            if option in ("-v", "--version"):
                raise Usage(ver_message)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option == "--ins":
                max_insert = int(value)
                if max_insert < 0 or max_insert > 10:
                    print >> sys.stderr, "Error: arg to --ins must be in range [0, 10]"
                    exit(1)
            if option == "--del":
                max_delete = int(value)
                if max_delete < 0:
                    print >> sys.stderr, "Error: arg to --del must be >= 0"
                    exit(1)            
            if option == "--max-append-mis":
                append_mismatch = int(value)
                if append_mismatch < 0:
                    print >> sys.stderr, "Error: arg to --max-append-mis must >= 0"
                    exit(1)
            if option in ("-m", "--splice-mis"):
                splice_mismatches_single = int(value)
                if splice_mismatches_single < 0 or splice_mismatches_single > 2:
                    print >> sys.stderr, "Error: arg to -m/--splice-mis must be in range of [0, 2]"
                    exit(1)
            if option in ("-i", "--min-intron"):
                min_intron_length = int(value)
                if min_intron_length <= 0:
                    print >> sys.stderr, "Error: arg to -i/--min-intron must > 0"
                    exit(1)                
            if option in ("-I", "--max-intron"):
                max_intron_length = int(value)
                if max_intron_length <= 0:
                    print >> sys.stderr, "Error: arg to -I/--max-intron must > 0"
                    exit(1)
            if option in ("-p", "--threads"):
                bowtie_threads = int(value)
                if bowtie_threads <= 0:
                    print >> sys.stderr, "Error: arg to -p/--threads must > 0"
                    exit(1)
            if option in ("-s", "--seglen"):
                seg_len = int(value)
                seg_len_overided = True
                if seg_len <= 0:
                    print >> sys.stderr, "Error: arg to -s/--seglen must > 0"
                    exit(1)
            if option == "--non-canonical":
                flank_case = 0
            if option == "--fusion":
                do_fusion = 1
            if option == "--fusion-non-canonical":
                fusion_flank_case = 0
                do_fusion = 1
            if option in ("-x", "--bowtie-index"):
                bwt_idx_prefix = value
            if option in ("-c", "--chromosome-dir"):
                chromosome_files_dir = value
                chromosome_files_dir = chromosome_files_dir + "/"
            if option == "--bam":
                output_bam = 1                                    
            if option in ("-k", "--max-hits"):
                max_hits = int(value)
                if max_hits <= 0:
                    print >> sys.stderr, "Error: arg to -k/--max-hits must be > 0"
                    exit(1)
            if option == "--filtering":
                filtering_level = int(value)
                if filtering_level < 1 or filtering_level > 2:
                    print >> sys.stderr, "Error: arg to --filtering can only be 1 and 2"
                    exit(1)
            if option == "--min-map-len":
                min_map_len = int(value)
                if min_map_len <= 0:
                    print >> sys.stderr, "Error: arg to --min-map-len must > 0"
                    exit(1)
            if option == "--min-len":
                min_read_len = int(value)
                if min_read_len <= 0:
                    print >> sys.stderr, "Error: arg to --min-len must > 0"
                    exit(1)
            if option == "--gene-gtf":
                gene_gtf_file = value                
            if option == "--min-entropy":
                min_entropy = float(value)
            if option == "--min-fusion-entropy":
                min_entropy_repeats = float(value)
            if option == "--qual-scale":
                if value == "phred64":
                    qual_scale = value
                elif value == "phred33":
                    qual_scale = value
                elif value == "solexa64":
                    qual_scale = value
                else:
                    print >> sys.stderr, "Error: --qual-scale can only be phred64 or phred33 or solexa64"
                    exit(1)                
            if option == "--keep-tmp":
                rm_temp = 0
            if option == "--not-rerun-all":
                rerun_all = 0
            if option == "--run-MapPER":
                run_mapper = 1
            if option == "--DEBUG":
                DEBUG = 1
            if option == "-1":
                input_reads_1 = value
            if option == "-2":
                input_reads_2 = value
        if max_delete >= min_intron_length:
            print >> sys.stderr, "Error: maximum deletion length(--del) must < minimum intron length(-i)"
            exit(1)
        if do_fusion == 1:
            min_map_len = 0
        if bwt_idx_prefix == "":
            build_bowtie_index_dir = output_dir + "bowtie_index/"
            bwt_idx_prefix = build_bowtie_index_dir + "built_bowtie_index"
            #print >> sys.stderr, "Error: -x/--bowtie-index not specified"
            #exit(1)
        if chromosome_files_dir == "":
            print >> sys.stderr, "Error: -c/--chromosome-dir not specified"
            exit(1)    
        if input_reads_1 == "":
            print >> sys.stderr, "Error: -1 not specified"
            exit(1)
        if input_reads_1 != "" and input_reads_2 != "":
            is_paired = 1        
        print_arguments(logging_dir + "argu_log")   
        
        start_time = datetime.now() 
        print >> sys.stderr
        print >> sys.stderr, "-----------------------------------------------" 
        print >> sys.stderr, "[%s] Beginning Mapsplice run (%s)" % (right_now(), get_version())
        print >> sys.stderr, "[%s] Bin directory: %s " % (right_now(), bin_dir)
        
        prepare_output_dir()
        no_original_junction = False
        
        # Validate all the input files, check all prereqs
        splitted_reads1 = input_reads_1.split(',')
        for read_file in splitted_reads1:
            check_file_existence(read_file)
        splitted_reads2 = []
        if input_reads_2 != "":
            splitted_reads2 = input_reads_2.split(',')
            if len(splitted_reads1) != len(splitted_reads2):
                print >> sys.stderr, "Error: number of read files not consistent in -1 and -2"
                exit(0)
            for read_file in splitted_reads2:
                check_file_existence(read_file)
        check_file_existence(chromosome_files_dir)
        check_bowtie_index(bwt_idx_prefix, chromosome_files_dir, 0, FASTA_file_extension, "reference sequence")
        inspect_bowtie_index(bwt_idx_prefix, logging_dir + "bowtie_refnames")
        ####Read chromosome size     
        all_chromos_path = read_dir_by_suffix(chromosome_files_dir, FASTA_file_extension)   
        chromo_fai_file = temp_dir + "chromo.fai"
        chromo_size_file = read_chromo_sizes(all_chromos_path, 
                                             temp_dir + "chrom_sizes", 
                                             fusion_dir + "chrName.txt", 
                                             temp_dir + "chrom_head",
                                             chromo_fai_file)
        check_index_consistency(logging_dir + "bowtie_refnames", temp_dir + "chrom_sizes")
        
        if input_reads_2 != "":
            combined_reads = ""
            for ix in range(len(splitted_reads1)):
                combined_reads = combined_reads + splitted_reads1[ix]
                combined_reads = combined_reads + ','
                combined_reads = combined_reads + splitted_reads2[ix]
                if ix < len(splitted_reads1) - 1:
                    combined_reads = combined_reads + ','
            check_reads_format(combined_reads, min_read_len, is_paired)
        else:
            check_reads_format(input_reads_1, min_read_len, is_paired)
        qual_scale_detected = ""
        (maxlen, format_flag, qual_scale_detected, low_support_thresh) = extract_maxlen(logging_dir + "check_reads_format.log")
        read_width = maxlen
        if maxlen < 50 and seg_len_overided == False:
            seg_len = max(18, maxlen/2) 
        if format_flag == "-q":
            is_fastq = 1
        else:
            is_fastq = 0 
        if is_fastq > 0: 
            if qual_scale == "":
                if qual_scale_detected != "unknown":
                    qual_scale = qual_scale_detected
                else:
                    print >> sys.stderr, "Error: Auto detect quality scale failed, please use --qual-scale to specify manually"
                    exit(0)
            else:
                if qual_scale_detected != "unknown" and qual_scale_detected != qual_scale:
                    print >> sys.stderr, "Warning: Auto detected quality scale [%s] not consitent with user specified [%s], user specified is used" \
                    % (qual_scale_detected, qual_scale)
                else:
                    print >> sys.stderr, "Auto detect quality scale failed, user specified [%s] is used" % qual_scale
                    
        #####Run mapsplice multi-thread version without remapping
        do_fusion_and_pairend = False    
        original_sam = call_mapsplice_multithreads(segment_mismatches,
                                                   splice_mismatches_double,
                                                   splice_mismatches_single, 
                                                   format_flag,
                                                   max_hits * 10,
                                                   bowtie_threads,
                                                   seg_len,
                                                   chromo_size_file,
                                                   chromosome_files_dir,
                                                   flank_case,
                                                   max_insert,
                                                   max_delete,
                                                   min_map_len,
                                                   original_dir + "original.sam",
                                                   bwt_idx_prefix,
                                                   input_reads_1,
                                                   input_reads_2,
                                                   "", #original_dir + "original_unmapped",
                                                   original_dir + "debug_info",
                                                   "",
                                                   no_original_junction,
                                                   False,
                                                   max_intron_length,
                                                   append_mismatch,
                                                   min_read_len,
                                                   min_intron_length,
                                                   qual_scale,
                                                   logging_dir + "mapsplice_original.log")            
        cat_multithread_files(original_dir + "original.sam", 1, bowtie_threads + 1, original_dir + "original.sam", 1)
        
        #####Convert sam alignments to junctions
        sam2juncarray([original_sam], original_dir + "ori.all_junctions.txt", chromosome_files_dir, 
                 maxlen, 1, 350000000, "", 1, logging_dir + "sam2junc_ori.all_junctions.log")
        remove_files([original_dir + "original.sam"])
    
        #####Filtering original junction
        max_seg = 2
        filter_junc_by_min_mis_lpq(original_dir + "ori.all_junctions.txt", 
                                   filteredbest_dir + "ori.all_junctions.filtered_by_min_mis_lpq.remained.txt", 
                                   filteredbest_dir + "ori.all_junctions.filtered_by_min_mis_lpq.filtered.txt", 
                                   remap_mismatch, float(max_seg + 1) / float(10),
                                   logging_dir + "ori.all_junctions_filtered_by_min_mis_lpq.log")
        entropy_weight = 0.097718
        lpq_weight = 0.66478
        ave_mis_weight = -0.21077
        min_score = 0.719
        (remained_junc, filtered_out_junc) = filterjuncbyROCarguNoncanon(filteredbest_dir + "ori.all_junctions.filtered_by_min_mis_lpq.remained.txt", 
                                                                         filteredbest_dir + "best_junction.txt",
                                                                         filteredbest_dir + "best_junction_semi_non_canon_filtered_by_ROCargu.txt", 
                                                                         entropy_weight, 
                                                                         lpq_weight, 
                                                                         ave_mis_weight, 
                                                                         min_score, 
                                                                         logging_dir + "best_junction_semi_non_canon_remained_ROC.log")
        if remained_junc == "" and filtered_out_junc == "":
            no_original_junction = True
        junction_tobe_syn = filteredbest_dir + "best_junction.txt"


        ##### generate synthetic junction sequence
        if no_original_junction == False:
            junc_db(junction_tobe_syn, 
                    chromosome_files_dir, 
                    2, 
                    38, 
                    400, 
                    remap_dir + "synthetic_alljunc_sequenc.txt",
                    logging_dir + "junc_db.log")
            check_bowtie_index(remap_dir + "syn_idx_prefix", 
                               remap_dir + "synthetic_alljunc_sequenc.txt", 
                               rerun_all, 
                               FASTA_file_extension,
                               "junction synthetic sequence")
        syn_idx_prefix = remap_dir + "syn_idx_prefix"
        if is_paired > 0 and do_fusion > 0:
            do_fusion_and_pairend = True
        
        ##### Run mapsplice multi-thread version remapping  
        remapped_sam = call_mapsplice_multithreads(segment_mismatches,
                                                   splice_mismatches_double,
                                                   splice_mismatches_single,
                                                   format_flag,
                                                   max_hits * 10,
                                                   bowtie_threads,
                                                   seg_len,
                                                   chromo_size_file,
                                                   chromosome_files_dir,
                                                   flank_case,
                                                   max_insert,
                                                   max_delete,
                                                   min_map_len,
                                                   remap_dir + "remapped.sam",
                                                   bwt_idx_prefix,
                                                   input_reads_1,
                                                   input_reads_2,
                                                   remap_dir + "remap_unmapped",
                                                   remap_dir + "debug_info",
                                                   syn_idx_prefix,
                                                   no_original_junction,
                                                   do_fusion_and_pairend,
                                                   max_intron_length,
                                                   append_mismatch,
                                                   min_read_len,
                                                   min_intron_length,
                                                   qual_scale,
                                                   logging_dir + "mapsplice_remap.log")
        cat_multithread_files(remap_dir + "remapped.sam", 1, bowtie_threads + 1, remap_dir + "remapped.sam", 1)
        cat_multithread_files(remap_dir + "remap_unmapped.1", 1, bowtie_threads + 1, remap_dir + "remap_unmapped.1", 1)
        if do_fusion == 0:
            fqfa2sam(remap_dir + "remap_unmapped.1.sam", [remap_dir + "remap_unmapped.1"], is_fastq, logging_dir + "remap_unmapped.1_2sam.log")     
            unmapped_reads_sams = unmapped_reads_sams + [remap_dir + "remap_unmapped.1.sam"]    
        if is_paired > 0:
            cat_multithread_files(remap_dir + "remap_unmapped.2", 1, bowtie_threads + 1, remap_dir + "remap_unmapped.2", 1)
            if do_fusion == 0:
                fqfa2sam(remap_dir + "remap_unmapped.2.sam", [remap_dir + "remap_unmapped.2"], is_fastq, logging_dir + "remap_unmapped.2_2sam.log")  
                unmapped_reads_sams = unmapped_reads_sams + [remap_dir + "remap_unmapped.2.sam"]
        
        
        ##### Run alignment handler filtering on remapped sam alignments
        max_mate_dist = 50000
        filtered_alignment_base = remap_dir + "_filtered_normal_alignments"
        filtered_alignment_append = "_filtered_normal_alignments"
        filter_flag = 12 + 32
        if filtering_level == 2:
            filter_flag = filter_flag + 2048
        if is_paired > 0 and do_fusion > 0:
            filter_flag = filter_flag + 256
        min_anchor = 0
        min_junction_anchor = 10
        min_mismatch = 5
        add_soft_clip = 1
        mate_dist_sd = 100
        max_anchor_diff = 50
        intron_dist_sd = 500
        encompassing_fusion_region_extension = 50000
        input_sam_file = remapped_sam
        min_coverage = 0
        fragment_length = 400
        fragment_length_sd = 100
        avearge_fragment_length = 225
        min_isoform_length = read_width / 2
        min_encompass_count = 1
        run_alignment_handler_multi("", 
                              is_paired,
                              max_mate_dist,
                              max_hits * 10,
                              filtered_alignment_base,
                              filtered_alignment_append,
                              read_width,
                              "",
                              filter_flag,
                              max_insert,
                              max_delete,
                              min_anchor,
                              min_junction_anchor,
                              min_mismatch,
                              add_soft_clip,
                              mate_dist_sd,
                              intron_dist_sd,
                              max_anchor_diff,
                              chromo_size_file,
                              encompassing_fusion_region_extension,
                              bowtie_threads,
                              min_coverage,
                              fragment_length,
                              fragment_length_sd,
                              avearge_fragment_length,
                              boundary,
                              min_isoform_length,
                              min_encompass_count,
                              min_entropy,
                              low_support_thresh,
                              input_sam_file,
                              logging_dir + "alignment_handler_remap.log",
                              logging_dir + "alignment_handler_remap.err")
        
        unmapped_reads_sams = unmapped_reads_sams + [remap_dir + "_filtered_normal_alignments.unmapped"]
        fusion_paired_alignments = input_sam_file + "_filtered_normal_alignments" + ".fusion_paired.comb"
        combine_fusion_paired = [input_sam_file + "_filtered_normal_alignments" + ".fusion_paired"] \
                              + [input_sam_file + "_filtered_normal_alignments" + ".bothunspliced.fusion_paired"]              
        cat_files(combine_fusion_paired, fusion_paired_alignments)
        single_alignments = input_sam_file + "_filtered_normal_alignments" + ".single.comb"
        combine_single = [input_sam_file + "_filtered_normal_alignments" + ".single"] \
                              + [input_sam_file + "_filtered_normal_alignments" + ".bothunspliced.single"]                    
        cat_files(combine_single, single_alignments)
        paired_alignments = [input_sam_file + "_filtered_normal_alignments" + ".paired"] \
                          + [input_sam_file + "_filtered_normal_alignments" + ".bothunspliced.paired"]
        remove_files([input_sam_file + "_filtered_normal_alignments"])
        remove_files([input_sam_file])
        
        ##### search fusion alignments
        filtered_fusion_alignment_base = ""
        if do_fusion > 0:
            if is_paired > 0:
                #--------------------cluster--------------------#
                parseCluster(fusion_paired_alignments, cluster_dir, logging_dir + "remap_parseCluster.log")
                #--------------------generate cluster region--------------------#
                cluster(cluster_dir, logging_dir + "remap_cluster.log")
            cluster_region = cluster_result_dir + "cluster.txt"
            ####Generate unmapped paired reads
            unmapped_fq1 = remap_dir + "remap_unmapped.1"
            unmapped_fq2 = ""
            if is_paired > 0:
                unmapped_fq2 = remap_dir + "remap_unmapped.2"
            
            ####Run MapSplice fusion
            no_original_fusion = False
            call_mapsplice_multithreads_fusion(segment_mismatches,
                                               splice_mismatches_double,
                                               splice_mismatches_single,
                                               format_flag,
                                               max_hits * 10,
                                               bowtie_threads,
                                               seg_len,
                                               chromo_size_file,
                                               chromosome_files_dir,
                                               flank_case,
                                               fusion_flank_case,
                                               max_insert,
                                               max_delete,
                                               min_map_len,
                                               fusion_dir + "normal.sam",
                                               bwt_idx_prefix,
                                               unmapped_fq1,
                                               unmapped_fq2,
                                               "", #fusion_dir + "fusion_original_unmapped",
                                               fusion_dir + "debug_info",
                                               syn_idx_prefix,
                                               "",
                                               cluster_region,
                                               fusion_dir + "fusion_alignments_original.sam",
                                               no_original_junction,
                                               no_original_fusion,
                                               max_intron_length,
                                               append_mismatch,
                                               min_read_len,
                                               min_intron_length,
                                               qual_scale,
                                               logging_dir + "mapsplice_fusion_original.log")
            
            ####convert fusion alignments to fusion junctions and do filtering
            fusionsam2junc_filteranchor_newfmt(fusion_dir + "fusion_alignments_original.sam", 
                                               fusion_dir + "original_fusion_junction.txt", maxlen, "", 1, 
                                               chromosome_files_dir, logging_dir + "fusionsam2junc_original_fusion_junction.log")
            fusion_junction_bef_remap = fusion_dir + "original_fusion_junction.txt"
            (remained_junc, filtered_out_junc) = filteroriginalfusion(fusion_dir + "original_fusion_junction.txt", 
                                 fusion_dir + "original_fusion_junction.remained.txt", 
                                 fusion_dir + "original_fusion_junction.filtered.txt", 
                                 2, float(max_seg + 1) / float(10),
                                 logging_dir + "ori.filteroriginalfusion.log")#
            if remained_junc == "" and filtered_out_junc == "":
                no_original_fusion = True

            ####generate fusion junction synthetic sequence
            fusion_syn_idx_prefix = ""
            if no_original_fusion == False:
                fusion_junc_db(junction_tobe_syn, fusion_dir + "original_fusion_junction.remained.txt", chromosome_files_dir, 2, 
                               38, 20, 50, fusion_dir + "fusion_synthetic_sequence.txt", logging_dir + "fusion_junc_db.log")
                check_bowtie_index(fusion_dir + "syn_fusion_idx_prefix", 
                                   fusion_dir + "fusion_synthetic_sequence.txt", 
                                   rerun_all, 
                                   FASTA_file_extension,
                                   "fusion junction synthetic sequence")
                fusion_syn_idx_prefix = fusion_dir + "syn_fusion_idx_prefix"
                
            #### Run MapSplice fusion Remapping 
            call_mapsplice_multithreads_fusion(segment_mismatches,
                                               splice_mismatches_double,
                                               splice_mismatches_single,
                                               format_flag,
                                               max_hits * 10,
                                               bowtie_threads,
                                               seg_len,
                                               chromo_size_file,
                                               chromosome_files_dir,
                                               flank_case,
                                               fusion_flank_case,
                                               max_insert,
                                               max_delete,
                                               min_map_len,
                                               fusion_dir + "normal.sam",
                                               bwt_idx_prefix,
                                               unmapped_fq1,
                                               unmapped_fq2,
                                               fusion_dir + "fusion_unmapped",
                                               fusion_dir + "debug_info",
                                               syn_idx_prefix,
                                               fusion_syn_idx_prefix,
                                               cluster_region,
                                               fusion_dir + "fusion_alignments_remap.sam",
                                               no_original_junction,
                                               no_original_fusion,
                                               max_intron_length,
                                               append_mismatch,
                                               min_read_len,
                                               min_intron_length,
                                               qual_scale,
                                               logging_dir + "mapsplice_fusion.log")
            
            remapped_filtered_fusion = [fusion_dir + "fusion_alignments_remap.sam"]
            sort_by_name_c(remapped_filtered_fusion, fusion_dir + "fusion_alignments.sam")
        
            ####Filter fusion junction and fusion alignments by primary filters
            fusionsam2junc_filteranchor_newfmt(fusion_dir + "fusion_alignments.sam", 
                                               fusion_dir + "remapped_fusion_junction.txt", maxlen, "", 1, 
                                               chromosome_files_dir, logging_dir + "fusionsam2junc_remapped_fusion_junction.log")
            fusion_min_mis = 2
            fusion_min_hits = 1
            fusion_min_lpq = 0
            fusion_min_entropy = 0.001
            filterremappedfusion(fusion_dir + "remapped_fusion_junction.txt", 
                                 fusion_dir + "remapped_fusion_junction.remained.txt", 
                                 fusion_dir + "remapped_fusion_junction.filtered.txt", 
                                 fusion_min_mis, 
                                 fusion_min_hits, 
                                 fusion_min_lpq, 
                                 fusion_min_entropy, 
                                 logging_dir + "filterremappedfusion.log")
            fusion_alignments_tobe_filtered = [fusion_dir + "fusion_alignments.sam"]
            FilterFusionSamByFusionJunc(fusion_alignments_tobe_filtered, 
                                        fusion_dir + "remapped_fusion_junction.remained.txt",
                                        fusion_dir + "fusion_alignments.remained.sam",
                                        fusion_dir + "fusion_alignments.filtered.sam", 
                                        logging_dir + "FilterFusionSamByFusionJunc.log")
            unmapped_reads_sams = unmapped_reads_sams + [fusion_dir + "fusion_alignments.filtered.sam.unmapped"]
            remapped_filtered_fusion = fusion_dir + "fusion_alignments.remained.sam"
    
    
            ####Generate Combined alignments 
            comb_fusion_alignments = []
            if is_paired > 0:
                comb_fusion_alignments = [single_alignments] + [fusion_paired_alignments] + [remapped_filtered_fusion]
            else:
                comb_fusion_alignments = [remapped_filtered_fusion]  
            comb_fusion_alignments_sorted = sort_by_name_c(comb_fusion_alignments, fusion_dir + "combined_fusion_normal.sam")
            
            #combine fusion unmapped reads
            cat_multithread_files(fusion_dir + "fusion_unmapped.1", 1, bowtie_threads + 1, fusion_dir + "fusion_unmapped.1", 1)
            fqfa2sam(fusion_dir + "fusion_unmapped.1.sam", [fusion_dir + "fusion_unmapped.1"], is_fastq, logging_dir + "fusion_unmapped.1_2sam.log")  
            unmapped_reads_sams = unmapped_reads_sams + [fusion_dir + "fusion_unmapped.1.sam"]
            if is_paired > 0:
                cat_multithread_files(fusion_dir + "fusion_unmapped.2", 1, bowtie_threads + 1, fusion_dir + "fusion_unmapped.2", 1)
                fqfa2sam(fusion_dir + "fusion_unmapped.2.sam", [fusion_dir + "fusion_unmapped.2"], is_fastq, logging_dir + "fusion_unmapped.2_2sam.log")
                unmapped_reads_sams = unmapped_reads_sams + [fusion_dir + "fusion_unmapped.2.sam"]
        

            ####Filter alignments with fusion
            filtered_fusion_alignment_base = fusion_dir + "_filtered_fusion_alignments"
            filtered_fusion_alignment_append = "_filtered_fusion_alignments"
            input_sam_file = comb_fusion_alignments_sorted
            filter_flag = 12 + 32 + 128 + 1024
            if filtering_level == 2:
                filter_flag = filter_flag + 2048
            fragment_length = 400
            fragment_length_sd = 100
            avearge_fragment_length = 225
            min_isoform_length = read_width / 2
            min_encompass_count = 1
            if maxlen >= 75:
                min_encompass_count = 0
            run_alignment_handler_multi("", 
                                        is_paired,
                                        max_mate_dist,
                                        max_hits * 10,
                                        filtered_fusion_alignment_base,
                                        filtered_fusion_alignment_append,
                                        read_width,
                                        "",
                                        filter_flag,
                                        max_insert,
                                        max_delete,
                                        min_anchor,
                                        min_junction_anchor,
                                        min_mismatch,
                                        add_soft_clip,
                                        mate_dist_sd,
                                        intron_dist_sd,
                                        max_anchor_diff,
                                        chromo_size_file,
                                        encompassing_fusion_region_extension,
                                        bowtie_threads,
                                        min_coverage,
                                        fragment_length,
                                        fragment_length_sd,
                                        avearge_fragment_length,
                                        boundary,
                                        min_isoform_length,
                                        min_encompass_count,
                                        min_entropy,
                                        low_support_thresh,
                                        input_sam_file,
                                        logging_dir + "alignment_handler_fusion.log",
                                        logging_dir + "alignment_handler_fusion.err")
            fusion_alignments = input_sam_file + "_filtered_fusion_alignments"
            unmapped_reads_sams = unmapped_reads_sams + [fusion_dir + "_filtered_fusion_alignments.unmapped"]
        
        if is_paired > 0:
            comb_unmapped_alignments_sorted = sort_by_name_c(unmapped_reads_sams, temp_dir + "combined_unmapped.sam")
            SetUnmappedBitFlag(temp_dir + "combined_unmapped.sam", temp_dir + "combined_unmapped_setbitflag.sam", logging_dir + "SetUnmappedBitFlag.log")
            unmapped_reads_sams = [temp_dir + "combined_unmapped_setbitflag.sam"]

        final_alignments = []
        if is_paired > 0 and do_fusion > 0:
            if no_original_fusion == False:
                final_alignments = [fusion_alignments] + paired_alignments
            else:
                final_alignments = paired_alignments + [single_alignments] + [fusion_paired_alignments]
        elif is_paired > 0 and do_fusion == 0:
            final_alignments = paired_alignments + [single_alignments] + [fusion_paired_alignments]
        elif is_paired == 0 and do_fusion > 0:
            if no_original_fusion == False:
                final_alignments = [fusion_alignments] + [single_alignments]
            else:
                final_alignments = [single_alignments]
        elif is_paired == 0 and do_fusion == 0:
            final_alignments = [single_alignments]
        
        final_alignments = final_alignments + unmapped_reads_sams
        add_head = [temp_dir + "chrom_head"] + final_alignments
        cat_files(add_head, temp_dir + "final_alignments_headed.sam")
        if output_bam > 0:
            if is_paired > 0:
                remove_pair_no(temp_dir + "final_alignments_headed.sam", temp_dir + "final_alignments_headed.sam.pair_no_removed")
            else:
                copy_file(temp_dir + "final_alignments_headed.sam", temp_dir + "final_alignments_headed.sam.pair_no_removed")
            sam2bam(temp_dir + "final_alignments_headed.sam.pair_no_removed", output_dir + "alignments.bam", chromo_fai_file)
        else:
            if is_paired > 0:
                remove_pair_no(temp_dir + "final_alignments_headed.sam", output_dir + "alignments.sam")
            else:
                copy_file(temp_dir + "final_alignments_headed.sam", output_dir + "alignments.sam")
        copy_file(filtered_alignment_base + ".fil.junc", output_dir + "junctions.txt")
        copy_file(filtered_alignment_base + ".fil.junc.del", output_dir + "deletions.txt")
        copy_file(filtered_alignment_base + ".fil.junc.ins", output_dir + "insertions.txt")
        
        do_repeat = 0
        do_annotate = 0
        candidate_fusion_file = output_dir + "fusions_candidates.txt"
        well_annotated_fusion_file = output_dir + "fusions_well_annotated.txt"
        not_well_annotated_fusion_file = output_dir + "fusions_not_well_annotated.txt"
        circular_RNAS_file = output_dir + "circular_RNAs.txt"
        if do_fusion > 0:
            copy_file(filtered_fusion_alignment_base + ".ori.junc.fusion.encompass.no.compass", output_dir + "fusions_raw.txt")
            min_anchor_repeat = 15
            if maxlen / 4 < min_anchor_repeat:
                min_anchor_repeat = 10  
            max_anchor_diff = maxlen / 2

            ##### filter fusion alignments by repeats
            filterbyrepeats(output_dir + "fusions_raw.txt",
                            output_dir + "junctions.txt",
                            output_dir + "fusions_candidates.txt",
                            chromosome_files_dir,
                            seg_len,
                            100,
                            max_anchor_diff,
                            min_anchor_repeat,
                            min_entropy_repeats,
                            bowtie_threads,
                            bwt_idx_prefix,
                            chromo_size_file,
                            logging_dir)
            do_repeat = 1
            
            ##### filter fusion by gene
            if gene_gtf_file != "":
                filterbyannotation(output_dir + "fusions_candidates.txt",
                                   output_dir + "fusions_well_annotated.txt",
                                   output_dir + "fusions_not_well_annotated.txt",
                                   output_dir + "circular_RNAs.txt",
                                   gene_gtf_file,
                                   logging_dir)
                do_annotate = 1
        
        if collect_stat > 0:
            collectstats(output_dir + "stats.txt", 
                         filtered_alignment_base + ".final_stat.txt", 
                         filtered_fusion_alignment_base + ".final_stat.txt", 
                         logging_dir + "SetUnmappedBitFlag.log",
                         is_paired,
                         do_fusion,
                         do_repeat,
                         do_annotate,
                         candidate_fusion_file,
                         well_annotated_fusion_file,
                         not_well_annotated_fusion_file,
                         circular_RNAS_file,
                         unmapped_reads_sams,
                         logging_dir + "collectstats.log")        
                         
        if rm_temp and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, True)
            
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr
        print >> sys.stderr, "[%s] Mapsplice finished running (time used: %s)" % (right_now(), duration)
        print >> sys.stderr, "-----------------------------------------------" 
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2

if __name__ == "__main__":
    sys.exit(main())
