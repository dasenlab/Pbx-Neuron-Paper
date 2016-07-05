import cluster

def get_samtools_bam_string(aligned_files_dir,sam_file):
    bam_file=sam_file[:-4]+".bam"
    sam_bam_string="samtools view -bS "+sam_file+" > "+bam_file
    return sam_bam_string,bam_file

def get_samtools_sort_string(sorted_bam_files_dir,bam_file,sample):
    # format: samtools sort SRR521547.bam SRR521547.sorted
    # output file will automatically add ".bam" to the end of the output file you give
    # specifying ".bam" will cause files to be named as ".bam.bam"
    sorted_bam_filename_in=sorted_bam_files_dir+sample+".sorted"
    samtools_sort_string="samtools sort "+bam_file+" "+sorted_bam_filename_in
    sorted_bam_filename_out=sorted_bam_filename_in+".bam"
    return samtools_sort_string,sorted_bam_filename_out

def get_picard_string(basedir,bam_file,aligned_files_dir,base_filename):
    metrics_dir=aligned_files_dir+"picard_metrics/"
    cluster.check_dir(metrics_dir)
    sge_filename=basedir+"sge_files/picard_"+base_filename+".sge"
    picard_string1='PICARD_SETTINGS="VERBOSITY=WARNING QUIET=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2500000"'
    picard_string2="java -Xmx16G -jar ${PICARD_ROOT}/MarkDuplicates.jar $PICARD_SETTINGS REMOVE_DUPLICATES=true ASSUME_SORTED=false CREATE_INDEX=true INPUT="+bam_file+" OUTPUT="+aligned_files_dir+base_filename+".dedup.bam METRICS_FILE="+metrics_dir+base_filename+".dups.txt"
    return picard_string1,picard_string2

def get_samtools_rmdup_string(basedir,bam_file,aligned_files_dir,base_filename):
    # format: samtools rmdup [-sS] <input.srt.bam> <output.bam>
    samtools_rmdup_string="samtools rmdup" 


def get_RnaSeqMetrics(basedir,bam_file,aligned_files_dir,base_filename):
    # picard RnaSeqMetrics
    # http://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics
    metrics_dir=aligned_files_dir+"picard_metrics/"
    cluster.check_dir(metrics_dir)
    metrics_txt=metrics_dir+base_filename+".txt"
    metrics_pdf=metrics_dir+base_filename+".pdf"
    process_list=[]
    # settings
    process_list.append('REFFLAT="/local/data/iGenomes/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/refFlat.txt.gz"')
    process_list.append('PICARD_SETTINGS="VERBOSITY=WARNING QUIET=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2500000"')
    process_list.append("java -Xmx16G -jar ${PICARD_ROOT}/CollectRnaSeqMetrics.jar \\")
    process_list.append("$PICARD_SETTINGS \\")
    process_list.append("REF_FLAT=${REFFLAT} \\")
    process_list.append("STRAND_SPECIFICITY=NONE \\")
    process_list.append("INPUT="+bam_file+" \\")
    process_list.append("CHART_OUTPUT="+metrics_pdf+" \\")
    process_list.append("OUTPUT="+metrics_txt)	
    return process_list
