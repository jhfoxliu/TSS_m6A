import sys,os
from sjm_tools import job, check_env

if not sys.argv[1]:
	print("Usage: python %s [directory] [submit(optional)]" % (sys.argv[0]))

#metadata
ref_genome_index = check_env("/home/fox/Database/GLORI_index/Gencode_v45_nochr/",is_path=True)
# transcriptome_A2G_index = check_env()
ref_genome = check_env("/home/fox/Database/Genome/Human/Gencode/v45/GRCh38.primary_assembly.genome.nochr.format.fa")
annot_first_exon_end = check_env("/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.anno.first_exon_end.bed")

#environment
python  = check_env("/home/fox/Software/bin/python3")
python2 = check_env("/home/fox/Software/bin/python2")
samtools = check_env("/home/fox/Software/samtools/1.16/bin/samtools")
hisat2 = check_env("/home/fox/Software/hisat2/2.2.1/hisat2")
hisat2_path = check_env("/home/fox/Software/hisat2/2.2.1/", is_path=True)
cutadapt = check_env("/home/fox/Software/bin/cutadapt")
umi_tools = check_env("/home/fox/Software/bin/umi_tools")
bedtools = check_env("/home/fox/Software/bedtools/2.29.1/bedtools")
bowtie2 = check_env("/home/fox/Software/bowtie2/2.3.4.2/bowtie2")
PATH = "./"

with open(sys.argv[1],'r') as input:
	for line in input.readlines():
		# Read sample list
		if line.strip() == "":
			continue
		name, read1, read2 = line.strip().split("\t")

		# Create workpath and define SJM file
		workpath = PATH + "/" + name + "_gencode/"
		SJM = name + ".CROWN.sjm"

		if os.path.isdir(workpath) == False:
			os.mkdir(workpath)
			
		JOB = job(workpath,SJM) 

		# Sort Bam
		JOB.step_start(step_name="TSS_m6A_annotate_exon",memory="80G")
		JOB.add_process("{python} /home/fox/Scripts/CROWN-seq/CROWN_fetch_first_exon.py -i {name}.TSS_m6A.txt -b {name}.merged.umi.bam -o {name}.TSS_exon_end.pairs.bed".format(python=python, name=name))

		JOB.add_process("{bedtools} sort -i {name}.TSS_exon_end.pairs.bed > {name}.TSS_exon_end.pairs.sorted.bed".format(bedtools=bedtools, name=name))

		JOB.add_process("{bedtools} closest -s -D a -a {name}.TSS_exon_end.pairs.sorted.bed -b {annot_first_exon_end} > {name}.TSS_exon_end_pairs.closest.bed".format(bedtools=bedtools, name=name, annot_first_exon_end=annot_first_exon_end))

		JOB.add_process("{python} /home/fox/Scripts/CROWN-seq/get_TSS_m6A_pair_to_first_exon_end.py -i {name}.TSS_m6A.txt.all.cov20.csv -b {name}.TSS_exon_end_pairs.closest.bed -o {name}.TSS_m6A.txt.all.cov20.annot_exon.csv ".format(python=python, name=name))

		JOB.step_end()

		JOB.job_finish()

		JOB.submit()