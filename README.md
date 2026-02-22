# TSS_m6A
Scripts to quantify internal m6A levels at 5'-isoform-specific manner from CROWN-seq data.

## Dependencies

* Python 3
* Pandas
* Numpy
* Biopython
* Pysam

## Workflow:

(1) Finish CROWN-seq mapping. See https://github.com/jhfoxliu/CROWN-Seq

(2) Run `python CROWN_seq_fetch_all_A_events_multiprocessing_v2.py -p 4 -r {reference_genome} -b {bam_file} -o {name}.txt`

Here: 

{refrence_genome}, the corresponding FASTA file; 

{bam_file}, the CROWN-seq mapping bam file output; 

{name}.txt, a file containing the read number counts for all TSS and A sites. Further filtering is required for this output.

(3) Run `python CROWN_fetch_first_exon.py -i {name}.txt -b {bam_file} -o {name}.bed`

Here: {bam_file} is the same as the one used in (2);

{name}.bed is a file containing the possible positions of the 3' exon edge.

Further correction is needed for the 3' exon edge. This script will also generate a file named **{name}.txt.all.cov20.csv**, which will be used in the final annotation step.

(4) `bedtools sort -i {name}.bed > {name}.sorted.bed`

(5) `bedtools closest -s -D a -a {name}.sorted.bed -b {first_exon_ends}.bed > {name}.sorted.first_exon_closest.bed`

This step annotate the bedfile to the closest 3' exon end based on the provided {first_exon_ends}.bed. (Check CROWN-seq metadata workflow)

You can generate this annotation file by extracting all the first exons from Gencode annotations. 

You can use `sort -u` to remove the duplicates.

(6) `get_TSS_m6A_pair_to_first_exon_end.py -i {name}.txt.all.cov20.csv -b {name}.sorted.first_exon_closest.bed -o {name}.TSS_m6A.txt.all.cov20.annot_exon.csv`

The **{name}.TSS_m6A.txt.all.cov20.annot_exon.csv** file is what you need.

## Contact

Email jil4026@med.cornell.edu or jhfoxliu@gmail.com if you have problems.

## Citation

Pending. Coming soon.

