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


For mouse mm10 genome hisat2 pre-built index:

    mkdir -p $gbm/RNA_REF_FA
    cd $gbm/RNA_REF_FA
    wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
    tar -xzvf mm10_genome.tar.gz

note: GRCm39 = mm39, and GRCm38 = mm10



Also download the reference genome compatible to Hisat2 from UCSC: 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
for other genomes: https://hgdownload.soe.ucsc.edu/downloads.html

File used by Andrew is:

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz 
    gzip -d hg38.ncbiRefSeq.gtf.gz


For mouse reference UCSC mm10: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/

    mkdir -p $gbm/RNA_REF_GTF
    cd $gbm/RNA_REF_GTF
    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
    gzip -d mm10.ncbiRefSeq.gtf.gz

note: mouse reference UCSC mm39 became available in 2020 (9 years newer than 38) but not used by hisat2: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/
reason this matters: mm10 and GRCm38 are synonymous: https://github.com/kundajelab/atac_dnase_pipelines/issues/143 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/ except: UCSC version will have chr* identifiers in the row names. 


CAUTION: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz serves as reference only for GENE LEVEL counts and NOT transcript level counts.


Alternatively, UCSC's references are also hosted on S3 bucket: http://daehwankimlab.github.io/hisat2/download/



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



    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:1 --rg LB:1_TGGTAGAGAT+TGTTGTTCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TGGTAGAGAT+TGTTGTTCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_1_S1_L001_R1_001.fastq.gz -2 $gbm_data/6931_1_S1_L001_R2_001.fastq.gz -S $gbm_data/alignments/1_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:1 --rg LB:1_TGGTAGAGAT+TGTTGTTCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TGGTAGAGAT+TGTTGTTCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_1_S1_L002_R1_001.fastq.gz -2 $gbm_data/6931_1_S1_L002_R2_001.fastq.gz -S $gbm_data/alignments/1_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:2 --rg LB:2_AGTACTCATG+GTAGAGTCAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.AGTACTCATG+GTAGAGTCAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_2_S2_L001_R1_001.fastq.gz -2 $gbm_data/6931_2_S2_L001_R2_001.fastq.gz -S $gbm_data/alignments/2_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:2 --rg LB:2_AGTACTCATG+GTAGAGTCAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.AGTACTCATG+GTAGAGTCAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_2_S2_L002_R1_001.fastq.gz -2 $gbm_data/6931_2_S2_L002_R2_001.fastq.gz -S $gbm_data/alignments/2_rep2.sam

    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:3 --rg LB:3_TACGTGAAGG+GACTGGTTGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TACGTGAAGG+GACTGGTTGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_3_S3_L001_R1_001.fastq.gz -2 $gbm_data/6931_3_S3_L001_R2_001.fastq.gz -S $gbm_data/alignments/3_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:3 --rg LB:3_TACGTGAAGG+GACTGGTTGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TACGTGAAGG+GACTGGTTGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_3_S3_L002_R1_001.fastq.gz -2 $gbm_data/6931_3_S3_L002_R2_001.fastq.gz -S $gbm_data/alignments/3_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:4 --rg LB:4_TGTGGTCCGG+GTTCCGCAGG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TGTGGTCCGG+GTTCCGCAGG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_4_S4_L001_R1_001.fastq.gz -2 $gbm_data/6931_4_S4_L001_R2_001.fastq.gz -S $gbm_data/alignments/4_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:4 --rg LB:4_TGTGGTCCGG+GTTCCGCAGG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TGTGGTCCGG+GTTCCGCAGG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_4_S4_L002_R1_001.fastq.gz -2 $gbm_data/6931_4_S4_L002_R2_001.fastq.gz -S $gbm_data/alignments/4_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:5 --rg LB:5_CCGACAGACT+TACCGAACTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CCGACAGACT+TACCGAACTA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_5_S5_L001_R1_001.fastq.gz -2 $gbm_data/6931_5_S5_L001_R2_001.fastq.gz -S $gbm_data/alignments/5_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:5 --rg LB:5_CCGACAGACT+TACCGAACTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CCGACAGACT+TACCGAACTA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_5_S5_L002_R1_001.fastq.gz -2 $gbm_data/6931_5_S5_L002_R2_001.fastq.gz -S $gbm_data/alignments/5_rep2.sam



    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:6 --rg LB:6_ATCCAGGTAT+ATCAACAGCC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.ATCCAGGTAT+ATCAACAGCC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_6_S6_L001_R1_001.fastq.gz -2 $gbm_data/6931_6_S6_L001_R2_001.fastq.gz -S $gbm_data/alignments/6_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:6 --rg LB:6_ATCCAGGTAT+ATCAACAGCC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.ATCCAGGTAT+ATCAACAGCC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_6_S6_L002_R1_001.fastq.gz -2 $gbm_data/6931_6_S6_L002_R2_001.fastq.gz -S $gbm_data/alignments/6_rep2.sam



    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:7 --rg LB:7_CGATGCGGTT+GGCTCTTGCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CGATGCGGTT+GGCTCTTGCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_7_S7_L001_R1_001.fastq.gz -2 $gbm_data/6931_7_S7_L001_R2_001.fastq.gz -S $gbm_data/alignments/7_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:7 --rg LB:7_CGATGCGGTT+GGCTCTTGCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CGATGCGGTT+GGCTCTTGCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_7_S7_L002_R1_001.fastq.gz -2 $gbm_data/6931_7_S7_L002_R2_001.fastq.gz -S $gbm_data/alignments/7_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:8 --rg LB:8_AATGCGAACA+ATTACTCACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.AATGCGAACA+ATTACTCACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_8_S8_L001_R1_001.fastq.gz -2 $gbm_data/6931_8_S8_L001_R2_001.fastq.gz -S $gbm_data/alignments/8_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:8 --rg LB:8_AATGCGAACA+ATTACTCACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.AATGCGAACA+ATTACTCACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_8_S8_L002_R1_001.fastq.gz -2 $gbm_data/6931_8_S8_L002_R2_001.fastq.gz -S $gbm_data/alignments/8_rep2.sam



    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:9 --rg LB:9_TTCGGTGTGA+AACAAGGCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TTCGGTGTGA+AACAAGGCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_9_S9_L001_R1_001.fastq.gz -2 $gbm_data/6931_9_S9_L001_R2_001.fastq.gz -S $gbm_data/alignments/9_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:9 --rg LB:9_TTCGGTGTGA+AACAAGGCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TTCGGTGTGA+AACAAGGCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_9_S9_L002_R1_001.fastq.gz -2 $gbm_data/6931_9_S9_L002_R2_001.fastq.gz -S $gbm_data/alignments/9_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:10 --rg LB:10_CATTAACTGA+AGAATCTTCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CATTAACTGA+AGAATCTTCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_10_S10_L001_R1_001.fastq.gz -2 $gbm_data/6931_10_S10_L001_R2_001.fastq.gz -S $gbm_data/alignments/10_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:10 --rg LB:10_CATTAACTGA+AGAATCTTCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CATTAACTGA+AGAATCTTCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_10_S10_L002_R1_001.fastq.gz -2 $gbm_data/6931_10_S10_L002_R2_001.fastq.gz -S $gbm_data/alignments/10_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:11 --rg LB:11_CCGGTTCCTA+CGGCAATGGA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CCGGTTCCTA+CGGCAATGGA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_11_S11_L001_R1_001.fastq.gz -2 $gbm_data/6931_11_S11_L001_R2_001.fastq.gz -S $gbm_data/alignments/11_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:11 --rg LB:11_CCGGTTCCTA+CGGCAATGGA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CCGGTTCCTA+CGGCAATGGA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_11_S11_L002_R1_001.fastq.gz -2 $gbm_data/6931_11_S11_L002_R2_001.fastq.gz -S $gbm_data/alignments/11_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:12 --rg LB:12_ATACATCACA+TAACAGTGTT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.ATACATCACA+TAACAGTGTT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_12_S12_L001_R1_001.fastq.gz -2 $gbm_data/6931_12_S12_L001_R2_001.fastq.gz -S $gbm_data/alignments/12_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:12 --rg LB:12_ATACATCACA+TAACAGTGTT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.ATACATCACA+TAACAGTGTT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_12_S12_L002_R1_001.fastq.gz -2 $gbm_data/6931_12_S12_L002_R2_001.fastq.gz -S $gbm_data/alignments/12_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:13 --rg LB:13_TAGAGAATAC+AGAGTGCGGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TAGAGAATAC+AGAGTGCGGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_13_S13_L001_R1_001.fastq.gz -2 $gbm_data/6931_13_S13_L001_R2_001.fastq.gz -S $gbm_data/alignments/13_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:13 --rg LB:13_TAGAGAATAC+AGAGTGCGGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TAGAGAATAC+AGAGTGCGGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_13_S13_L002_R1_001.fastq.gz -2 $gbm_data/6931_13_S13_L002_R2_001.fastq.gz -S $gbm_data/alignments/13_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:14 --rg LB:14_GTCTCGCCAC+CAGCCGATTG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.GTCTCGCCAC+CAGCCGATTG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_14_S14_L001_R1_001.fastq.gz -2 $gbm_data/6931_14_S14_L001_R2_001.fastq.gz -S $gbm_data/alignments/14_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:14 --rg LB:14_GTCTCGCCAC+CAGCCGATTG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.GTCTCGCCAC+CAGCCGATTG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_14_S14_L002_R1_001.fastq.gz -2 $gbm_data/6931_14_S14_L002_R2_001.fastq.gz -S $gbm_data/alignments/14_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:15 --rg LB:15_CGCTGTCTCA+ATGTCGTGGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CGCTGTCTCA+ATGTCGTGGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_15_S15_L001_R1_001.fastq.gz -2 $gbm_data/6931_15_S15_L001_R2_001.fastq.gz -S $gbm_data/alignments/15_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:15 --rg LB:15_CGCTGTCTCA+ATGTCGTGGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CGCTGTCTCA+ATGTCGTGGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_15_S15_L002_R1_001.fastq.gz -2 $gbm_data/6931_15_S15_L002_R2_001.fastq.gz -S $gbm_data/alignments/15_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:16 --rg LB:16_GCTAATAGGA+AGAGCACTAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.GCTAATAGGA+AGAGCACTAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_16_S16_L001_R1_001.fastq.gz -2 $gbm_data/6931_16_S16_L001_R2_001.fastq.gz -S $gbm_data/alignments/16_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:16 --rg LB:16_GCTAATAGGA+AGAGCACTAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.GCTAATAGGA+AGAGCACTAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_16_S16_L002_R1_001.fastq.gz -2 $gbm_data/6931_16_S16_L002_R2_001.fastq.gz -S $gbm_data/alignments/16_rep2.sam



    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:20 --rg LB:20_GCGTGCTGTG+CGGTGACACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.GCGTGCTGTG+CGGTGACACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_20_S20_L001_R1_001.fastq.gz -2 $gbm_data/6931_20_S20_L001_R2_001.fastq.gz -S $gbm_data/alignments/20_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:20 --rg LB:20_GCGTGCTGTG+CGGTGACACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.GCGTGCTGTG+CGGTGACACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_20_S20_L002_R1_001.fastq.gz -2 $gbm_data/6931_20_S20_L002_R2_001.fastq.gz -S $gbm_data/alignments/20_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:21 --rg LB:21_CGAGAGGCGT+GAGACATAAT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CGAGAGGCGT+GAGACATAAT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_21_S21_L001_R1_001.fastq.gz -2 $gbm_data/6931_21_S21_L001_R2_001.fastq.gz -S $gbm_data/alignments/21_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:21 --rg LB:21_CGAGAGGCGT+GAGACATAAT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CGAGAGGCGT+GAGACATAAT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_21_S21_L002_R1_001.fastq.gz -2 $gbm_data/6931_21_S21_L002_R2_001.fastq.gz -S $gbm_data/alignments/21_rep2.sam



the following samples are E0771 mouse breast cancer mets to the brain, so need to align with mouse reference genome

    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:17 --rg LB:17_CCTAACACAG+GGTGGAATAC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CCTAACACAG+GGTGGAATAC -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_17_S17_L001_R1_001.fastq.gz -2 $gbm_data/6931_17_S17_L001_R2_001.fastq.gz -S $gbm_data/alignments/17_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:17 --rg LB:17_CCTAACACAG+GGTGGAATAC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CCTAACACAG+GGTGGAATAC -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_17_S17_L002_R1_001.fastq.gz -2 $gbm_data/6931_17_S17_L002_R2_001.fastq.gz -S $gbm_data/alignments/17_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:18 --rg LB:18_TGCCGGTCAG+GAGGCTCCTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TGCCGGTCAG+GAGGCTCCTA -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_18_S18_L001_R1_001.fastq.gz -2 $gbm_data/6931_18_S18_L001_R2_001.fastq.gz -S $gbm_data/alignments/18_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:18 --rg LB:18_TGCCGGTCAG+GAGGCTCCTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TGCCGGTCAG+GAGGCTCCTA -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_18_S18_L002_R1_001.fastq.gz -2 $gbm_data/6931_18_S18_L002_R2_001.fastq.gz -S $gbm_data/alignments/18_rep2.sam


    hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:19 --rg LB:19_TTAACCTTCG+TAATGGCAAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TTAACCTTCG+TAATGGCAAG -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_19_S19_L001_R1_001.fastq.gz -2 $gbm_data/6931_19_S19_L001_R2_001.fastq.gz -S $gbm_data/alignments/19_rep1.sam
    hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:19 --rg LB:19_TTAACCTTCG+TAATGGCAAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TTAACCTTCG+TAATGGCAAG -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_19_S19_L002_R1_001.fastq.gz -2 $gbm_data/6931_19_S19_L002_R2_001.fastq.gz -S $gbm_data/alignments/19_rep2.sam



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

