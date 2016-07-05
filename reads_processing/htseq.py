import cluster

def get_htseqcount_string(basedir,sorted_bam_filename,gtf_file,sample):
    htseq_counts_dir=basedir+"htseq_counts_third/"
    cluster.check_dir(htseq_counts_dir)
    # use this with gtf file:
    htseq_string="htseq-count --stranded=no --format=bam "+sorted_bam_filename+" "+gtf_file+" > "+htseq_counts_dir+sample+"_counts.txt"
    # use this with gff file:
    #htseq_string="htseq-count --stranded=no --format=bam --idattr=ID --type=gene "+sorted_bam_filename+" "+gtf_file+" > "+htseq_counts_dir+sample+"_counts.txt"
    return htseq_string
