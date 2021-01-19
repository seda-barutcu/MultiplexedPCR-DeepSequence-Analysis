# MultiplexedPCR-DeepSequence-Analysis
Identify and count gene-specific and non-specific amplicon reads from high-throughput sequencing of multiplexed-PCR

Systen requirements Python2 or Python3.

C19_SparSeq_amplicon_sequence_finder.py : This script is formated specifically for analyzing multiplexed PCR fastq files from the output of C19_SparSeq pipeline.

C19SS_pilot_R1_nonSpecAmp.py : This code processes each fastq file given as a name-list in "filelist_R1.txt" text file (line75), and outputs the non-specific amplicons binned to each individual forward primer in the multiplexed pcr.
In addition, the total numbers of the non-specific amlicons are summed up by the type of sample such as "water, Twist, LTRI(patient-sample)" for the pilot cohort processed for C19-SparSeq, and reported for each primer per sample category.

C19SS_pilot_R2_nonSpecAmp.py : This code processes each fastq file given as a name-list in "filelist_R2.txt" text file (line75), and outputs the non-specific amplicons binned to each individual forward primer in the multiplexed pcr.
In addition, the total numbers of the non-specific amlicons are summed up by the type of sample such as "water, Twist, LTRI(patient-sample)" for the pilot cohort processed for C19-SparSeq, and reported for each primer per sample category.
