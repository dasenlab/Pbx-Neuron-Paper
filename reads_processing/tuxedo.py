# RNA-Seq tuxedo pipeline
# updated 09/25/2014
# Lisa Cohen
# Genome Technology Center
# NYU Langone Medical Center
# tophat --> cufflinks --> cuffquant --> cuffdiff

import os

def get_bowtie1_string_single(aligned_files_dir,reference,fastq_file,sam_filename):
    # format:
    # bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
    # manual: http://bowtie-bio.sourceforge.net/manual.shtml#bowtie-options-m
    # -v is the alignment mismatches allowed
    # -m is the suppression of reads beyond the number of matches, so -m 1 will only allow 1 read match
    #sam_filename=aligned_files_dir+filename+".sam"
    bowtie1_string="bowtie -v 2 -m 1 -S "+reference+" "+fastq_file+" "+sam_filename
    return bowtie1_string

def get_bowtie1_string_paired(aligned_files_dir,reference,fastq_file_list,sam_filename):
    read1=fastq_file_list[0]
    read2=fastq_file_list[1]
    bowtie1_string="bowtie -v 2 -m 1 -S "+reference+" -1 "+read1+" -2 "+read2+" "+sam_filename
    return bowtie1_string

def get_bowtie2_string_single(aligned_files_dir,fastq_file,filename):
    sam_filename=aligned_files_dir+filename+".sam"
    #print sam_filename
    bowtie2_string=["bowtie2 -x "+reference+" "+fastq_file+" -S "+sam_filename]
    return bowtie2_string

def get_bowtie2_string_paired(aligned_files_dir,fastq_file1,fastq_file2,filename):
# http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# format:
# bowtie2 -x pmirabHI4320 -1 ChipIP1_CGATGT_L001_R1_001.fastq -2 ChipIP1_CGATGT_L001_R2_001.fastq -S ChipIP1_CGATGT_L001_R1_001.sam
    sam_filename=aligned_files_dir+filename+".sam"
    #print sam_filename
    bowtie2_string="bowtie2 -x "+reference+" -1 "+fastq_file1+" -2 "+fastq_file2+" -S "+sam_filename
    return bowtie2_string

def get_tophat_string_annotation(ref,annotation,filenames,output_dir):
# Format: tophat -p 32 -G annotation -o <outputdir> ref reads_1.fq reads_2.fq
# info on turning off coverage search:
# https://www.biostars.org/p/49224/
# Coverage search feature will only look for novel junctions. It saves memory and time to turn it off.
    tophat_string="tophat -p 8 --GTF "+annotation+" -o "+output_dir+" "+ref+" "+filenames 
    return tophat_string

def get_tophat_string(ref,annotation,filenames,output_dir):
# Manual: http://ccb.jhu.edu/software/tophat/manual.shtml
# tophat -p 32 -o <outputdir> ref reads_1.fq reads_2.fq
# -p is the number of threads
    tophat_string="tophat -p 8 -o "+output_dir+" "+ref+" "+filenames 
    return tophat_string

def get_tophat_output(tophat_outputdir):
    tophat_output=tophat_outputdir+"accepted_hits.bam"
    #print "The tophat output file exists:",tophat_output
    #print os.path.isfile(tophat_output)
    return tophat_output

def get_dupremoved_output(base_filename,tophat_outputdir):
# to-do, move this to samtools module
    dupremoved_output=tophat_outputdir+base_filename+".dedup.bam"
    #print os.path.isfile(dupremoved_output)
    #print dupremoved_output
    return dupremoved_output

def get_cufflinks_string(cufflinks_outputdir,tophat_output_string,annotation):
# cufflinks manual:
# http://cole-trapnell-lab.github.io/cufflinks/manual/
    cufflinks_string="cufflinks -g "+annotation+" -p 8 -o "+cufflinks_outputdir+" "+tophat_output_string
    #print cufflinks_string
    return cufflinks_string

def get_cuffmerge_string(reference_fasta,annotation,assembly_filename,cuffmerge_dir):
# don't need to do cuffmerge unless you actually want to merge gtf from different reference-guided assemblies
# get gene names by giving one standard .gtf to cufflinks
# cuffmerge -o <outputdirpath> -g genes.gtf -s genome.fa -p 8 assemblies.txt
    cuffmerge_string="cuffmerge -o "+cuffmerge_dir+"merge -g "+annotation+" -s "+reference_fasta+" -p 6 "+assembly_filename
    #print cuffmerge_string
    return cuffmerge_string

def get_cuffcompare_string(annotation,cuffcompare_outputdir,cufflinks_outputdir):
# The -R flag does ??
# cuffcompare -o cuffcompare_output -R -r genome.gtf cufflinks_output/transcripts.gtf
    cuffcompare_string="cuffcompare -o "+cuffcompare_outputdir+" -R -r "+annotation+" "+cufflinks_outputdir+"transcripts.gtf"
    #print "This cufflinks output file exists:",cufflinks_outputdir+"transcripts.gtf"
    #print os.path.isfile(cufflinks_outputdir+"transcripts.gtf")
    return cuffcompare_string

def get_cuffquant_string(cuffquant_outputdir,tophat_output,annotation):
    # cuffquant [options]* <annotation.(gtf/gff)> <aligned_reads.(sam/bam)>
    # requires cuffmerge output
    cuffquant_string="cuffquant -p 8 -o "+cuffquant_outputdir+" "+annotation+" "+tophat_output
    #print cuffquant_string
    return cuffquant_string

def get_cuffquant_cxb_output(cuffquant_dir):
    cxb_filename=cuffquant_dir+"abundances.cxb"
    print "This file exists:",cxb_filename
    print os.path.isfile(cxb_filename)
    return cxb_filename

def get_cuffdiff_string(basedir,cuffdiff_outputdir,groupsoffiles,annotation,labels):
    # Run cuffdiff by using the merged transcriptome assembly along with the cxb files from cuffquant for each replicate:
    # needs groups of files, each group of files is separated by comma, groups are separated from each other with a space, from cuffquant
    cuffdiff_string="cuffdiff -o "+cuffdiff_outputdir+" -p $NSLOTS "+" -L "+labels+" "+annotation+" "+groupsoffiles
    return cuffdiff_string

def create_assembly_file(transcripts_filename,assembly_filename):
# input is list of samples
# Creates a file called assemblies.txt that lists the transcripts assembly file for each sample
# This file is required input for cuffmerge, which is required input for cuffquant
# One line per sample
# format:
# /ifs/home/cohenl06/data/sequencing/for_steve/cufflinks/JMS.140529.B1/transcripts.gtf
    with open(assembly_filename,"a") as assembly_file:
        assembly_file.write(transcripts_filename+"\n")

def assembly_file_check(assembly_file):
    with open(assembly_file) as assemblyfile:
        for line in assemblyfile:
            filename=line.rstrip("\n")
            print filename
            print os.path.isfile(filename)

def alignment_summary_single(tophat_outputdir):
    align_sum_file=tophat_outputdir+"align_summary.txt"
    if os.path.isfile(align_sum_file)==True:
        with open(align_sum_file) as datafile:
            for line in datafile:
                line_data=line.split()
                print line_data
                if len(line_data)!=0:
                    first=line_data[0]
                    if first=="Left:":
                        reads_input=line_data[1]
                        #print "Reads input:",reads_input
                    if first=="Mapped:":
                        reads_mapped=line_data[1]
                        #print "Reads mapped:",reads_mapped
                        percent_mapped_input_raw=line_data[2]
                        percent_mapped_input=percent_mapped_input_raw[1:]
                        #print percent_mapped_input
                    if first=="of":
                        mult_align=line_data[2]
                        #print "Number with multiple alignemnts:",mult_align
                        percent_mult_align_raw=line_data[3]
                        percent_mult_align=percent_mult_align_raw[1:-1]
                        #print "Percent with multiple alignments:",percent_mult_align
                        morethan20_align_raw=line_data[7]
                        morethan20_align=morethan20_align_raw[1:]
                        #print "Number with more than 20 align:",morethan20_align
                    if first[-1]=="%":
                        overall_percent_aligned=line_data[0]
                        #print "Overall percent aligned:",overall_percent_aligned
                    else:
                        print "Nothing of importance:",first
                else:
                    print "line_data empty."
            data_list=[reads_input,reads_mapped,percent_mapped_input,mult_align,percent_mult_align,morethan20_align,overall_percent_aligned]
            #print "\t".join(data_list)
            return data_list
    else:
        print "File not ready yet:",align_sum_file
        return "File not ready yet:",align_sum_file


def get_duplicates_removed(tophat_outputdir,sample):
    duplicates_removed_file=tophat_outputdir+"picard_metrics/"+sample+".dups.txt"
    line_data={}
    line_num=0
    if os.path.isfile(duplicates_removed_file)==True:
        with open(duplicates_removed_file) as datafile:
            for line in datafile:
                line_num+=1
                line_data[line_num]=line.split()
        duplicates_removed=line_data[8][8]
        return duplicates_removed

def get_alignment_summary(tophat_outputdir):
    align_sum_file=tophat_outputdir+"align_summary.txt"
    line_data={}
    line_num=0
    if os.path.isfile(align_sum_file)==True:
        with open(align_sum_file) as datafile:
            for line in datafile:
                line_num+=1
                line_data[line_num]=line.split()
        return line_data
    else:
        print "File not ready yet:",align_sum_file
        return "File not ready yet:",align_sum_file

def get_alignment_list(line_data):
    data_list=[]
    left_reads=line_data[2][1]
    left_mapped=line_data[3][1]
    left_percent_mapped=line_data[3][2][1:]
    right_reads=line_data[6][1]
    right_mapped=line_data[7][1]
    right_percent_mapped=line_data[7][2][1:]
    overall_mapped=line_data[9][0]
    pairs_aligned=line_data[11][2]
    data_list=[left_reads,left_mapped,left_percent_mapped,right_reads,right_mapped,right_percent_mapped,overall_mapped,pairs_aligned]
    len(data_list)
    return data_list

def alignment_table_file_tophat(basedir,sample_dictionary):
    alignment_table_filename=basedir+"alignment_stats.txt"
    tophat_dir=basedir+"tuxedo/tophat/"
    header=["Sample","% Duplicates Removed","Left Reads","Left Mapped","Left % Mapped","Right Reads","Right Mapped","Right % Mapped","Overall % Mapped","Pairs Aligned"]
    len(header)
    with open(alignment_table_filename,"w") as datafile:
        datafile.write("\t".join(header))
        datafile.write("\n")
        for sample in sample_dictionary.keys():
            tophat_outputdir=tophat_dir+sample+"/"
            line_data_dictionary=get_alignment_summary(tophat_outputdir)
            data_list=get_alignment_list(line_data_dictionary)
            dups_removed=get_duplicates_removed(tophat_outputdir,sample)
            datafile.write(sample+"\t")
            datafile.write(dups_removed+"\t")
            #print "\t".join(data_list)
            datafile.write("\t".join(data_list))
            datafile.write("\n")
    datafile.close()
    print "Alignment stats written:",alignment_table_filename

def get_bowtie1_alignment_data(filename):
    # parse contents of first four lines from bowtie1 sge file, e.g:
    # # reads processed: 15494922
    # reads with at least one reported alignment: 10529780 (67.96%)
    # reads that failed to align: 1740139 (11.23%)
    # reads with alignments suppressed due to -m: 3225003 (20.81%)
    line_num_list=range(0,4)
    line_data_list=[]
    data_list=[]
    if os.path.isfile(filename)==True:
        with open(filename) as datafile:
            for line_num in line_num_list:
                line_data=next(datafile).split('\t')
                line_data_list.append(line_data[0])
    else:
        print "Filename does not exist:",filename
    reads_data=line_data_list[0].split()
    reads=reads_data[3]
    data_list.append(reads)
    reads_alignment_data=line_data_list[1].split()
    reads_alignment=reads_alignment_data[8]
    data_list.append(reads_alignment)
    percent_reads_alignment=reads_alignment_data[9][1:-2]
    data_list.append(percent_reads_alignment)
    reads_failed_data=line_data_list[2].split()
    reads_failed=reads_failed_data[6]
    data_list.append(reads_failed)
    percent_reads_failed=reads_failed_data[7][1:-2]
    data_list.append(percent_reads_failed)
    suppressed_alignments_data=line_data_list[3].split()
    suppressed_reads=suppressed_alignments_data[8]
    data_list.append(suppressed_reads)
    percent_suppressed_reads=suppressed_alignments_data[9][1:-2]
    data_list.append(percent_suppressed_reads)
    return data_list

#def bowtie1_alignment_file(basedir,sample_dictionary):
