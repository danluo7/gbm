# gbm

## check quality of reads

    fastqc *.fastq.gz 
    
    python3 -m multiqc . 
    
    mkdir fastqc
    mv *fastqc* fastqc
    
hg38 genome index:    
  
Sources for obtaining gene annotation files formatted for HISAT2: HISAT2 pre-built Genome Index 
http://daehwankimlab.github.io/hisat2/download/

Also available from their FTP site ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/
(this was used by Andrew in Hingtgen's lab)

    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz
    tar -xzvf hg38.tar.gz
  
These 8 files together constitute the index: they are all that is needed to align reads to that reference. 


For mouse mm10 genome hisat pre-built index:

    wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz



Also download the reference genome compatible to Hisat2 from UCSC: 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
for other genomes: https://hgdownload.soe.ucsc.edu/downloads.html

File used by Andrew is:

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz 
    gzip -d hg38.ncbiRefSeq.gtf.gz

For mouse reference UCSC mm10: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/

    wget: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    

note: This is the newest reference but isn't used by hisat2, since it came out after HISAT2 had published their prebuilt indexes: mouse reference UCSC mm39: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/



CAUTION: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz serves as reference only for GENE LEVEL counts and NOT transcript level counts.


Alternatively, UCSC's references are also hosted on S3 bucket: http://daehwankimlab.github.io/hisat2/download/
note: mm10 and GRCm38 are synonymous: https://github.com/kundajelab/atac_dnase_pipelines/issues/143 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/ except: UCSC version will have chr* identifiers in the row names


### alignment using hisat2 (which suceeded hisat and tophat)

setting the Read Group info:

a) Sample and library tags. Can be autmatically detected from current sample naming scheme:
<Sample.ID><Index.Sequence><Lane.ID><Set.number>.fastq

which for this dataset is: 
6931_1_S1_L001_R1_001.fastq.gz

sample.ID: 6931
index.sequence: 1
lane.id: L001
set.number: R1


Also can be identified from the name of a sequence read in the Fastq file:

@(instrument id) @A00201R
:(run number) 402
:(flowcell ID) HCTMMDRXY
:(lane) 1
:(tile) 2102
:(x_pos) 21305
:(y_pos) (read) 28087
:(pair) 2
:(is filtered) N
:(control number) 0
:(UMI) GCTAATAGGA+AGAGCACTAG
:(index sequence)

ID and PU (to enable merging replictes)


ID = Read group identifier = {FLOWCELL BARCODE/ID}.{LANE} HCTMMDRXY.1 or .2
SM = <Sample.ID> 1 through xx (individual sample that was submitted)
LB = <Sample.ID>_<Index.Sequence> (to indentify if same library was used)

PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{library-specific identifier}. This is the most specific definition for a group of reads.


















## Convert SAM files into BAM files

    cd $gbm_data/alignments

    samtools sort -@ 8 -o 1_rep1.bam 1_rep1.sam
    samtools sort -@ 8 -o 1_rep2.bam 1_rep2.sam

    samtools sort -@ 8 -o 2_rep1.bam 2_rep1.sam
    samtools sort -@ 8 -o 2_rep2.bam 2_rep2.sam

    samtools sort -@ 8 -o 3_rep1.bam 3_rep1.sam
    samtools sort -@ 8 -o 3_rep2.bam 3_rep2.sam

    samtools sort -@ 8 -o 4_rep1.bam 4_rep1.sam
    samtools sort -@ 8 -o 4_rep2.bam 4_rep2.sam

    samtools sort -@ 8 -o 5_rep1.bam 5_rep1.sam
    samtools sort -@ 8 -o 5_rep2.bam 5_rep2.sam

    samtools sort -@ 8 -o 6_rep1.bam 6_rep1.sam
    samtools sort -@ 8 -o 6_rep2.bam 6_rep2.sam

    samtools sort -@ 8 -o 7_rep1.bam 7_rep1.sam
    samtools sort -@ 8 -o 7_rep2.bam 7_rep2.sam

    samtools sort -@ 8 -o 8_rep1.bam 8_rep1.sam
    samtools sort -@ 8 -o 8_rep2.bam 8_rep2.sam

    samtools sort -@ 8 -o 9_rep1.bam 9_rep1.sam
    samtools sort -@ 8 -o 9_rep2.bam 9_rep2.sam

    samtools sort -@ 8 -o 10_rep1.bam 10_rep1.sam
    samtools sort -@ 8 -o 10_rep2.bam 10_rep2.sam

    samtools sort -@ 8 -o 11_rep1.bam 11_rep1.sam
    samtools sort -@ 8 -o 11_rep2.bam 11_rep2.sam

    samtools sort -@ 8 -o 12_rep1.bam 12_rep1.sam
    samtools sort -@ 8 -o 12_rep2.bam 12_rep2.sam

    samtools sort -@ 8 -o 13_rep1.bam 13_rep1.sam
    samtools sort -@ 8 -o 13_rep2.bam 13_rep2.sam

    samtools sort -@ 8 -o 14_rep1.bam 14_rep1.sam
    samtools sort -@ 8 -o 14_rep2.bam 14_rep2.sam

    samtools sort -@ 8 -o 15_rep1.bam 15_rep1.sam
    samtools sort -@ 8 -o 15_rep2.bam 15_rep2.sam

    samtools sort -@ 8 -o 16_rep1.bam 16_rep1.sam
    samtools sort -@ 8 -o 16_rep2.bam 16_rep2.sam

    samtools sort -@ 8 -o 17_rep1.bam 17_rep1.sam
    samtools sort -@ 8 -o 17_rep2.bam 17_rep2.sam

    samtools sort -@ 8 -o 18_rep1.bam 18_rep1.sam
    samtools sort -@ 8 -o 18_rep2.bam 18_rep2.sam

    samtools sort -@ 8 -o 19_rep1.bam 19_rep1.sam
    samtools sort -@ 8 -o 19_rep2.bam 19_rep2.sam

    samtools sort -@ 8 -o 20_rep1.bam 20_rep1.sam
    samtools sort -@ 8 -o 20_rep2.bam 20_rep2.sam

    samtools sort -@ 8 -o 21_rep1.bam 21_rep1.sam
    samtools sort -@ 8 -o 21_rep2.bam 21_rep2.sam


Index bam files: make sure that the only files in the directory are the sam and bam/bai files, then:

    find *.bam -exec echo samtools index {} \; | sh

