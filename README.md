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

note: mouse reference UCSC mm39 became available in 2020 (mm38/mm10 came out in 2011) but not used by hisat2: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/
reason this matters: mm10 and GRCm38 are synonymous: https://github.com/kundajelab/atac_dnase_pipelines/issues/143 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/ except: UCSC version will have chr* identifiers in the row names. 


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

## At this point decide on naming conventions. Order the samples logically (ie. if samples 1 2 3 should be ordered 3 1 2, then in hisat2 put input sample file 1 > 3.bam, and sample file 2 > 1.bam etc. Then make a note of this change.


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

    samtools sort -@ 8 -o 20_rep1.bam 20_rep1.sam
    samtools sort -@ 8 -o 20_rep2.bam 20_rep2.sam

    samtools sort -@ 8 -o 21_rep1.bam 21_rep1.sam
    samtools sort -@ 8 -o 21_rep2.bam 21_rep2.sam
    
mouse ones
    
    samtools sort -@ 8 -o 17_rep1.bam 17_rep1.sam
    samtools sort -@ 8 -o 17_rep2.bam 17_rep2.sam

    samtools sort -@ 8 -o 18_rep1.bam 18_rep1.sam
    samtools sort -@ 8 -o 18_rep2.bam 18_rep2.sam

    samtools sort -@ 8 -o 19_rep1.bam 19_rep1.sam
    samtools sort -@ 8 -o 19_rep2.bam 19_rep2.sam


## Merge the bam files

    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=1.bam INPUT=1_Rep1.bam INPUT=1_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=2.bam INPUT=2_Rep1.bam INPUT=2_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=3.bam INPUT=3_Rep1.bam INPUT=3_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=4.bam INPUT=4_Rep1.bam INPUT=4_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=5.bam INPUT=5_Rep1.bam INPUT=5_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=6.bam INPUT=6_Rep1.bam INPUT=6_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=7.bam INPUT=7_Rep1.bam INPUT=7_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=8.bam INPUT=8_Rep1.bam INPUT=8_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=9.bam INPUT=9_Rep1.bam INPUT=9_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=10.bam INPUT=10_Rep1.bam INPUT=10_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=11.bam INPUT=11_Rep1.bam INPUT=11_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=12.bam INPUT=12_Rep1.bam INPUT=12_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=13.bam INPUT=13_Rep1.bam INPUT=13_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=14.bam INPUT=14_Rep1.bam INPUT=14_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=15.bam INPUT=15_Rep1.bam INPUT=15_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=16.bam INPUT=16_Rep1.bam INPUT=16_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=17.bam INPUT=17_Rep1.bam INPUT=17_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=18.bam INPUT=18_Rep1.bam INPUT=18_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=19.bam INPUT=19_Rep1.bam INPUT=19_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=20.bam INPUT=20_Rep1.bam INPUT=20_Rep2.bam
    java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=21.bam INPUT=21_Rep1.bam INPUT=21_Rep2.bam







## Index all bam files for IGV visualization
Make sure that the only files in the directory are the sam and bam/bai files, then:

can do:
samtools index 1.bam
samtools index 2.bam
or:

    find *.bam -exec echo samtools index {} \; | sh



## use samtoools flagstat to geta basic sumary of an alignment 
for example percent of unmapped reads

	mkdir flagstat
	samtools flagstat 1.bam > flagstat/1.bam.flagstat
	samtools flagstat 2.bam > flagstat/2.bam.flagstat


view the resulting flagstat: 
	
	cat flagstat/1.bam.flagstat
	


## Use Stringtie to generate expression estimates (FPKM and TPM) from the SAM/BAM files generated by HISAT2 and samtools
this takes the bam file, add expression estimates FPKM and TPM using the reference annotation GTF file as a guide, and gives an annotated gtf file for the sample.
    
  
    mkdir -p $gbm_data/expression/stringtie/ref_only/
    cd $gbm_data/expression/stringtie/ref_only/


    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 1/transcripts.gtf -A 1/gene_abundances.tsv $gbm_data/alignments/1.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 2/transcripts.gtf -A 2/gene_abundances.tsv $gbm_data/alignments/2.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 3/transcripts.gtf -A 3/gene_abundances.tsv $gbm_data/alignments/3.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 4/transcripts.gtf -A 4/gene_abundances.tsv $gbm_data/alignments/4.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 5/transcripts.gtf -A 5/gene_abundances.tsv $gbm_data/alignments/5.bam

    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 6/transcripts.gtf -A 6/gene_abundances.tsv $gbm_data/alignments/6.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 7/transcripts.gtf -A 7/gene_abundances.tsv $gbm_data/alignments/7.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 8/transcripts.gtf -A 8/gene_abundances.tsv $gbm_data/alignments/8.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 9/transcripts.gtf -A 9/gene_abundances.tsv $gbm_data/alignments/9.bam

    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 10/transcripts.gtf -A 10/gene_abundances.tsv $gbm_data/alignments/10.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 11/transcripts.gtf -A 11/gene_abundances.tsv $gbm_data/alignments/11.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 12/transcripts.gtf -A 12/gene_abundances.tsv $gbm_data/alignments/12.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 13/transcripts.gtf -A 13/gene_abundances.tsv $gbm_data/alignments/13.bam

    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 14/transcripts.gtf -A 14/gene_abundances.tsv $gbm_data/alignments/14.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 15/transcripts.gtf -A 15/gene_abundances.tsv $gbm_data/alignments/15.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 16/transcripts.gtf -A 16/gene_abundances.tsv $gbm_data/alignments/16.bam

    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 17/transcripts.gtf -A 17/gene_abundances.tsv $gbm_data/alignments/17.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 18/transcripts.gtf -A 18/gene_abundances.tsv $gbm_data/alignments/18.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 19/transcripts.gtf -A 19/gene_abundances.tsv $gbm_data/alignments/19.bam

    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 20/transcripts.gtf -A 20/gene_abundances.tsv $gbm_data/alignments/20.bam
    stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 21/transcripts.gtf -A 21/gene_abundances.tsv $gbm_data/alignments/21.bam




view raw output from stringtie

	less -S 1/transcripts.gtf
	
View transcript records only (column -t makes formatting easier to view). scroll right to see the ### 3 major expression matrices "cov, fpkm, and tpm". ### Cov, or coverage, is average per base covereage over the genomic segment. basically, the coverage of each transcript. it's un-normalized unlike fpkm or tpm.

	grep -v "^#" 1/transcripts.gtf | grep -w "transcript" | column -t | less -S

Limit the view only to transcript records and their expression estimates (FPKM and TPM values). the if statement {if ($3=="transcript") print} means only print the line if the third column is a transcript. then only return field 1, 4, 9 which are the fpkm nad tpm values.

	awk '{if ($3=="transcript") print}' 1/transcripts.gtf | cut -f 1,4,9 | less -S

Could also view gene and trnscript level expression values in the two files generated by stringtie, a transcript.gtf file (transcript level abundance) and a gene_abundances.tsv (a gene level abundance):

	column -t 1/t_data.ctab | less -S
	less -S -x20 1/gene_abundances.tsv

## Parallel to Stringtie, can also run htseq-count on alignments to produce raw counts instead of FPKM/TPM values (what Stringtie outputs) for differential expression analysis


htseq-count basic usage:

	htseq-count [options] <sam_file> <gff_file>

Extra options specified below:

	’–format’ specify the input file format one of BAM or SAM. Since we have BAM format files, select ‘bam’ for this option.
	’–order’ provide the expected sort order of the input file. Previously we generated position sorted BAM files so use ‘pos’.
	’–mode’ determines how to deal with reads that overlap more than one feature. We believe the ‘intersection-strict’ mode is best.
	’–stranded’ specifies whether data is stranded or not. The TruSeq strand-specific RNA libraries suggest the ‘reverse’ option for this parameter.
	’–minaqual’ will skip all reads with alignment quality lower than the given minimum value
	’–type’ specifies the feature type (3rd column in GFF file) to be used. (default, suitable for RNA-Seq and Ensembl GTF files: exon)
	’–idattr’ The feature ID used to identify the counts in the output table. The default, suitable for RNA-SEq and Ensembl GTF files, is gene_id.


	mkdir -p $gbm_data/expression/htseq_counts
	cd $gbm_data/expression/htseq_counts


	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/1.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 1.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/2.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 2.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/3.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 3.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/4.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 4.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/5.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 5.tsv

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/6.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 6.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/7.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 7.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/8.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 8.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/9.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 9.tsv

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/10.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 10.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/11.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 11.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/12.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 12.tsv

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/13.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 13.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/14.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 14.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/15.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 15.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/16.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 16.tsv

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/17.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 17.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/18.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 18.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/19.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 19.tsv

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/20.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 20.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/21.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 21.tsv



Join the results for each replicate together and merge results files into a single matrix for use in edgeR. Since edgeR only does pairwise comparisons, need to pair up the compared samples:

Pairs:

for E0771Br:
(17, 18) vs. 19



	cd $gbm/expression/htseq_counts
	
	join 17.tsv 18.tsv | join - 19.tsv > E0771_gene_read_counts_table_all.tsv


for 011 GBM samples:
- (20, 21) vs. (1, 2)
- (20, 21) vs. 3
- (20, 21) vs. (4, 5)
- (1, 2) vs. 3
- (1, 2) vs. (4, 5)
- 3 vs. (4, 5)
-

    

	join 20.tsv 21.tsv | join - 1.tsv | join - 2.tsv | join - 3.tsv | join - 4.tsv | join - 5.tsv > GBM011_gene_read_counts_table_all.tsv
	
join 20.tsv 21.tsv | join - 1.tsv | join - 2.tsv > GBM011_invitro_vs_slice_gene_read_counts_table.tsv
join 20.tsv 21.tsv | join - 3.tsv > GBM011_invitro_vs_organoids_gene_read_counts_table.tsv




for 024 GBM samples:
- 7 vs. 6
- 7 vs. 8
- 7 vs. 9
- 6 vs. 8
- 6 vs. 9
- 8 vs. 9
-

	join 7.tsv 6.tsv | join - 8.tsv | join - 8.tsv | join - 9.tsv > GBM024_gene_read_counts_table_all.tsv


for UNC lung mets:
- 12 vs. 13
- 12 vs. 10
- 12 vs. 11
- 10 vs. 11
.



	join 12.tsv 13.tsv | join - 10.tsv | join - 11.tsv > UNClung_gene_read_counts_table_all.tsv
	

for UNC GBMs:
- 16 vs. 14
- 16 vs. 15
- 14 vs. 15
.

	join 16.tsv 14.tsv | join - 15.tsv > UNCGBM_gene_read_counts_table_all.tsv
	
	
	
	
Creat a simple text file with just the header that will be used for the table:

	echo "GeneID E0771Br_invitro_1 E0771Br_invitro_2 E0771Br_slices" > E0771_header.txt
	echo "GeneID 011_invitro_1 011_invitro_2 011_slices_1 011_slices_2 011_organoids 011_tissue_1 011_tissue_2" > GBM011_header.txt
	echo "GeneID 024_invitro 024_slices 024_organoids 024_tissue" > GBM024_header.txt
	echo "GeneID UNClung_invitro_p0 UNClung_invitro_p4 UNClung_slice UNClung_tissue" > UNClung_header.txt
	echo "GeneID UNCGBM_invitro UNCGBM_slice UNCGBM_tissue" > UNCGBM_header.txt



Clean up a bit more, add a header, reformat the result as a tab delimited file. note: grep -v "__" is being used to filter out the summary lines at the end of the files that ht-seq count gives to summarize reads that had no feature, were ambiguous, did not align at all, did not align due to poor alignment quality, or the alignment was not unique.

awk -v OFS="\t" '$1=$1' is using awk to replace the single space characters that were in the concatenated version of our header.txt and gene_read_counts_table_all.tsv with a tab character. -v is used to reset the variable OFS, which stands for Output Field Separator. By default, this is a single space. By specifying OFS="\t", we are telling awk to replace the single space with a tab. The '$1=$1' tells awk to reevaluate the input using the new output variable


	cat E0771_header.txt E0771_gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > E0771_gene_read_counts_table_all_final.tsv
	cat GBM011_header.txt GBM011_gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > GBM011_gene_read_counts_table_all_final.tsv
	cat GBM024_header.txt GBM024_gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > GBM024_gene_read_counts_table_all_final.tsv
	cat UNClung_header.txt UNClung_gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > UNClung_gene_read_counts_table_all_final.tsv
	cat UNCGBM_header.txt UNCGBM_gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > UNCGBM_gene_read_count_table_all_final.tsv
	
	head E0771_gene_read_counts_table_all_final.tsv | column -t

remove files no longer needed
	
	rm -f E0771_gene_read_counts_table_all.tsv E0771_header.txt
	
or just delete these files manually 

# Use Ballgown in R for differential expression (DE) analysis using output from Stringtie
## Perform A vs. B comparison, using all replicates, for known (reference only mode) transcripts

(raw counts using htseq output will come later)

	mkdir -p $gbm/de/ballgown/ref_only
	cd $gbm/de/ballgown/ref_only/


Use printf to create/print a table with ids, type (each type of sample is a type), and path to the file, as the header. Then n returns a new line. 
## note: id needs to match the file folder names created by stringtie

Bascially, need a table that needs to look like this to feed into R:

ids type path to file 011_invitro_1 011 $gbm/expression/stringtie/1 011_invitro_2 011 $gbm/expression/stringtie/2 ... ...

this is how the script should look like (without the enters inbetween each line):



pairwise comparisons for 011:

against invitro
printf "\"ids\",\"type\",\"path
\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5
\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20
\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21
\"\n" > 011_tissue_vs_invitro.csv

printf "\"ids\",\"type\",\"path
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20
\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21
\"\n" > 011_slice_vs_invitro.csv

printf "\"ids\",\"type\",\"path
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3
\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20
\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21
\"\n" > 011_organoid_vs_invitro.csv



against slice
printf "\"ids\",\"type\",\"path
\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20
\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n" > 011_invitro_vs_slice.csv


printf "\"ids\",\"type\",\"path
\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n" > 011_tissue_vs_slice.csv


printf "\"ids\",\"type\",\"path
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n" > 011_organoid_vs_slice.csv




against tissue
printf "\"ids\",\"type\",\"path
\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20
\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21
\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5
\"\n" > 011_invitro_vs_tissue.csv


printf "\"ids\",\"type\",\"path
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5
\"\n" > 011_slice_vs_tissue.csv


printf "\"ids\",\"type\",\"path
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3
\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5
\"\n" > 011_organoid_vs_tissue.csv



against organoid
printf "\"ids\",\"type\",\"path
\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20
\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3
\"\n" > 011_invitro_vs_organoid.csv


printf "\"ids\",\"type\",\"path
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3
\"\n" > 011_slice_vs_organoid.csv


printf "\"ids\",\"type\",\"path
\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3
\"\n" > 011_tissue_vs_organoid.csv


printf "\"ids\",\"type\",\"path\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n" > 011_tissue_vs_invitro.csv

printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n" > 011_slice_vs_invitro.csv

printf "\"ids\",\"type\",\"path\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n" > 011_organoid_vs_invitro.csv




	printf "\"ids\",\"type\",\"path\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n" > 011_invitro_vs_slice.csv


	printf "\"ids\",\"type\",\"path\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n" > 011_tissue_vs_slice.csv


	printf "\"ids\",\"type\",\"path\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n" > 011_organoid_vs_slice.csv





	printf "\"ids\",\"type\",\"path\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n" > 011_invitro_vs_tissue.csv


	printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n" > 011_slice_vs_tissue.csv


	printf "\"ids\",\"type\",\"path\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n" > 011_organoid_vs_tissue.csv




	printf "\"ids\",\"type\",\"path\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n" > 011_invitro_vs_organoid.csv


	printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n" > 011_slice_vs_organoid.csv


	printf "\"ids\",\"type\",\"path\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n" > 011_tissue_vs_organoid.csv




pairwise comparisons for 024:

printf "\"ids\",\"type\",\"path
\"\n\"7\",\"024_invitro\",\"$gbm/expression/stringtie/7
\"\n\"6\",\"024_slice\",\"$gbm/expression/stringtie/6
\"\n" > 024_invitro_vs_slice.csv


printf "\"ids\",\"type\",\"path
\"\n\"7\",\"024_invitro\",\"$gbm/expression/stringtie/7
\"\n\"8\",\"024_organoid\",\"$gbm/expression/stringtie/8
\"\n" > 024_invitro_vs_organoid.csv


printf "\"ids\",\"type\",\"path
\"\n\"7\",\"024_invitro\",\"$gbm/expression/stringtie/7
\"\n\"9\",\"024_tissue\",\"$gbm/expression/stringtie/9
\"\n" > 024_invitro_vs_tissue.csv

printf "\"ids\",\"type\",\"path
\"\n\"6\",\"024_slice\",\"$gbm/expression/stringtie/6
\"\n\"8\",\"024_organoid\",\"$gbm/expression/stringtie/8
\"\n" > 024_slice_vs_organoid.csv


printf "\"ids\",\"type\",\"path
\"\n\"8\",\"024_organoid\",\"$gbm/expression/stringtie/8
\"\n\"9\",\"024_tissue\",\"$gbm/expression/stringtie/9
\"\n" > 024_organoid_vs_tissue.csv



pairwise comparison for UNC lung

printf "\"ids\",\"type\",\"path
\"\n\"12\",\"UNClung_invitro_p0\",\"$gbm/expression/stringtie/12
\"\n\"13\",\"UNClung_invitro_p4\",\"$gbm/expression/stringtie/13
\"\n" > UNClung_invitro_vs_slice.csv

printf "\"ids\",\"type\",\"path
\"\n\"12\",\"UNClung_invitro\",\"$gbm/expression/stringtie/12
\"\n\"13\",\"UNClung_invitro\",\"$gbm/expression/stringtie/13
\"\n\"10\",\"UNClung_slice\",\"$gbm/expression/stringtie/10
\"\n" > UNClung_invitro_vs_slice.csv

printf "\"ids\",\"type\",\"path
\"\n\"12\",\"UNClung_invitro\",\"$gbm/expression/stringtie/12
\"\n\"13\",\"UNClung_invitro\",\"$gbm/expression/stringtie/13
\"\n\"11\",\"UNClung_tissue\",\"$gbm/expression/stringtie/11
\"\n" > UNClung_invitro_vs_tissue.csv

printf "\"ids\",\"type\",\"path
\"\n\"10\",\"UNClung_slice\",\"$gbm/expression/stringtie/10
\"\n\"11\",\"UNClung_tissue\",\"$gbm/expression/stringtie/11
\"\n" > UNClung_slice_vs_tissue.csv



pairwise comparison of UNC GBM

printf "\"ids\",\"type\",\"path
\"\n\"16\",\"UNCGBM_invitro\",\"$gbm/expression/stringtie/16
\"\n\"14\",\"UNCGBM_slice\",\"$gbm/expression/stringtie/14
\"\n" > UNCGBM_invitro_vs_slice.csv

printf "\"ids\",\"type\",\"path
\"\n\"16\",\"UNCGBM_invitro\",\"$gbm/expression/stringtie/16
\"\n\"15\",\"UNCGBM_tissue\",\"$gbm/expression/stringtie/15
\"\n" > UNCGBM_invitro_vs_tissue.csv

printf "\"ids\",\"type\",\"path
\"\n\"14\",\"UNCGBM_slice\",\"$gbm/expression/stringtie/14
\"\n\"15\",\"UNCGBM_tissue\",\"$gbm/expression/stringtie/15
\"\n" > UNCGBM_slice_vs_tissue.csv


pairwise comp of E0771

printf "\"ids\",\"type\",\"path
\"\n\"17\",\"E0771Br_invitro\",\"$gbm/expression/stringtie/17
\"\n\"18\",\"E0771Br_invitro\",\"$gbm/expression/stringtie/18
\"\n\"19\",\"E0771Br_slice\",\"$gbm/expression/stringtie/19
\"\n" > E0771Br_invitro_vs_slice.csv



script

pairwise comparisons for 011:

	printf "\"ids\",\"type\",\"path\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n" > 011_invitro_vs_slice.csv

	printf "\"ids\",\"type\",\"path\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n" > 011_invitro_vs_organoid.csv

	printf "\"ids\",\"type\",\"path\"\n\"20\",\"011_invitro\",\"$gbm/expression/stringtie/20\"\n\"21\",\"011_invitro\",\"$gbm/expression/stringtie/21\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n" > 011_invitro_vs_tissue.csv

	printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n" > 011_slice_vs_organoid.csv

	printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n" > 011_slice_vs_tissue.csv

	printf "\"ids\",\"type\",\"path\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/3\"\n\"4\",\"011_tissue\",\"$gbm/expression/stringtie/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/5\"\n" > 011_organoid_vs_tissue.csv




pairwise comparisons for 024:

	printf "\"ids\",\"type\",\"path\"\n\"7\",\"024_invitro\",\"$gbm/expression/stringtie/7\"\n\"6\",\"024_slice\",\"$gbm/expression/stringtie/6\"\n" > 024_invitro_vs_slice.csv

	printf "\"ids\",\"type\",\"path\"\n\"7\",\"024_invitro\",\"$gbm/expression/stringtie/7\"\n\"8\",\"024_organoid\",\"$gbm/expression/stringtie/8\"\n" > 024_invitro_vs_organoid.csv

	printf "\"ids\",\"type\",\"path\"\n\"7\",\"024_invitro\",\"$gbm/expression/stringtie/7\"\n\"9\",\"024_tissue\",\"$gbm/expression/stringtie/9\"\n" > 024_invitro_vs_tissue.csv

	printf "\"ids\",\"type\",\"path\"\n\"6\",\"024_slice\",\"$gbm/expression/stringtie/6\"\n\"8\",\"024_organoid\",\"$gbm/expression/stringtie/8\"\n" > 024_slice_vs_organoid.csv

	printf "\"ids\",\"type\",\"path\"\n\"8\",\"024_organoid\",\"$gbm/expression/stringtie/8\"\n\"9\",\"024_tissue\",\"$gbm/expression/stringtie/9\"\n" > 024_organoid_vs_tissue.csv



pairwise comparison for UNC lung


	printf "\"ids\",\"type\",\"path\"\n\"12\",\"UNClung_invitro_p0\",\"$gbm/expression/stringtie/12\"\n\"13\",\"UNClung_invitro_p4\",\"$gbm/expression/stringtie/13\"\n" > UNClung_invitro_vs_slice.csv

	printf "\"ids\",\"type\",\"path\"\n\"12\",\"UNClung_invitro\",\"$gbm/expression/stringtie/12\"\n\"13\",\"UNClung_invitro\",\"$gbm/expression/stringtie/13\"\n\"10\",\"UNClung_slice\",\"$gbm/expression/stringtie/10\"\n" > UNClung_invitro_vs_slice.csv

	printf "\"ids\",\"type\",\"path\"\n\"12\",\"UNClung_invitro\",\"$gbm/expression/stringtie/12\"\n\"13\",\"UNClung_invitro\",\"$gbm/expression/stringtie/13\"\n\"11\",\"UNClung_tissue\",\"$gbm/expression/stringtie/11\"\n" > UNClung_invitro_vs_tissue.csv

	printf "\"ids\",\"type\",\"path\"\n\"10\",\"UNClung_slice\",\"$gbm/expression/stringtie/10\"\n\"11\",\"UNClung_tissue\",\"$gbm/expression/stringtie/11\"\n" > UNClung_slice_vs_tissue.csv



pairwise comparison of UNC GBM

	printf "\"ids\",\"type\",\"path\"\n\"16\",\"UNCGBM_invitro\",\"$gbm/expression/stringtie/16\"\n\"14\",\"UNCGBM_slice\",\"$gbm/expression/stringtie/14\"\n" > UNCGBM_invitro_vs_slice.csv

	printf "\"ids\",\"type\",\"path\"\n\"16\",\"UNCGBM_invitro\",\"$gbm/expression/stringtie/16\"\n\"15\",\"UNCGBM_tissue\",\"$gbm/expression/stringtie/15\"\n" > UNCGBM_invitro_vs_tissue.csv

	printf "\"ids\",\"type\",\"path\"\n\"14\",\"UNCGBM_slice\",\"$gbm/expression/stringtie/14\"\n\"15\",\"UNCGBM_tissue\",\"$gbm/expression/stringtie/15\"\n" > UNCGBM_slice_vs_tissue.csv


pairwise comp of E0771

	printf "\"ids\",\"type\",\"path\"\n\"17\",\"E0771Br_invitro\",\"$gbm/expression/stringtie/17\"\n\"18\",\"E0771Br_invitro\",\"$gbm/expression/stringtie/18\"\n\"19\",\"E0771Br_slice\",\"$gbm/expression/stringtie/19\"\n" > E0771Br_invitro_vs_slice.csv



start R and load libraries

	R
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)


Load phenotype data from the file just saved in the current working directory (will have to repeat all below steps for each comparison

	pheno_data = read.csv("011_invitro_vs_slice.csv")

Load ballgown data structure and save it to a variable "bg"

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

	bg
ballgown output: ballgown instance with 190734 transcripts and 4 samples


Load all attributes including gene name, put in a table by using a transcript expression function from the ballgown library (texpr), feed it the ballgown library, all of it.

	bg_table = texpr(bg, 'all')

Then pull out just the gene names and gene id's. [, 9:10] means pull out all rows and just column 9 and 10. Keep in mind that the ncbi reference genomes gene id and gene names are the same, so column 9 and 10 will be the same. [, x,x] is essentially how to subset data within R. This table is used in later step for DE after merging the results.

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

Save the ballgown object to a file, R data file, for later use, so later so don't have to load each sample from their directories.

	save(bg, file='bg.rda')
	

# Perform differential expression (DE) analysis with no filtering

creat a results object. Use function stattest, feed the ballgown object, tell it whta we want to look at is transcript data, and the covariate that it's gonna perform the DE on is type, which was created earlier in the table using the printf function. type was H460 vs H460_2G etc. Use getFC to get fold change. Use "meas" to use FPKM as measurement of the foldchange

	results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")

same thing with gene level data

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

	head(results_genes)

if using ncbi's reference genome, this step isn't necessary. "gene name" and "id" are the same. Did this anyway to keep consistent with existing pipelines.

This creates a new table of a merge between original results_gene table, and a specific column in bg_gene_names (which was created earlier using [,9:10]). This is done by matching the column "id" from the results_gene table with "gene_id" column in the bg_gene_names table, and append (via merge) into a new results_genes table.

	results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

save a tab delimited file for both transcript and gene results. separater is a tab, denoted by \t. No quotes around objects to be printed.

	write.table(results_transcripts, "011_invitro_vs_slice_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
	write.table(results_genes, "011_invitro_vs_slice_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)



# Filter low-abundance genes. Removing all TRANSCRIPTS with a variance across the samples of less than one

	bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

Load all attributes including gene name

	bg_filt_table = texpr(bg_filt , 'all')
	bg_filt_gene_names = unique(bg_filt_table[, 9:10])



# Perform DE analysis now using the filtered data

	results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

Output the filtered list of genes and transcripts and save to tab delimited files

	write.table(results_transcripts, "011_invitro_vs_slice_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
	write.table(results_genes, "011_invitro_vs_slice_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

Identify the significant genes with q-value < 0.05. Note that q-value is what most people will use to filter in a large dataset.

	sig_transcripts = subset(results_transcripts, results_transcripts$qval<0.23)
	sig_genes = subset(results_genes, results_genes$qval<0.23)

	head(sig_genes)

	nrow(sig_genes)
	
output: [1] 3856

Output the signifant gene results to a pair of tab delimited files

	write.table(sig_transcripts, "011_invitro_vs_slice_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
	write.table(sig_genes, "011_invitro_vs_slice_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

	quit()
	y


	grep -v feature 011_invitro_vs_slice_gene_results_filtered.tsv | wc -l

	grep -v feature 011_invitro_vs_slice_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t 	#Higher abundance in invitro
	grep -v feature 011_invitro_vs_slice_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t 	#Higher abundance in slices

save the results into a new file with just the names of the genes (column 6)

	grep -v feature 011_invitro_vs_slice_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > 011_invitro_vs_slice_DE_genes.txt
	head 011_invitro_vs_slice_DE_genes.txt



## In vitro vs. tissue R code


R
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)



	pheno_data = read.csv("011_invitro_vs_tissue.csv")



	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

	bg



	bg_table = texpr(bg, 'all')



	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)



	save(bg, file='bg.rda')
	



	results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")



	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

	head(results_genes)




	results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))



	write.table(results_transcripts, "011_invitro_vs_tissue_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
	write.table(results_genes, "011_invitro_vs_tissue_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)





	bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)



	bg_filt_table = texpr(bg_filt , 'all')
	bg_filt_gene_names = unique(bg_filt_table[, 9:10])





	results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))



	write.table(results_transcripts, "011_invitro_vs_tissue_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
	write.table(results_genes, "011_invitro_vs_tissue_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)



	sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
	sig_genes = subset(results_genes, results_genes$pval<0.05)

	head(sig_genes)

	nrow(sig_genes)
	




	write.table(sig_transcripts, "011_invitro_vs_tissue_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
	write.table(sig_genes, "011_invitro_vs_tissue_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

	quit()
	n


	grep -v feature 011_invitro_vs_tissue_gene_results_filtered.tsv | wc -l

	grep -v feature 011_invitro_vs_tissue_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t 	#Higher abundance in invitro
	grep -v feature 011_invitro_vs_tissue_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t 	#Higher abundance in tissue



	grep -v feature 011_invitro_vs_tissue_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > 011_invitro_vs_tissue_DE_genes.txt
	head 011_invitro_vs_tissue_DE_genes.txt




## slices vs. tissue



