#install.packages("BiocManager")
library(BiocManager)
#install(c("sangerseqR","annotate")) 

setwd("C:/Users/jason/Dropbox/My PC (DESKTOP-IOO274E)/Desktop/Data/Lemna_singleinoculations/16S seqs (CAGEF)")

library(sangerseqR)

ITS<-read.abif("C1_F_D07_029_2022-09-21.ab1")

ITSseq <- sangerseq(ITS)
str(ITSseq)

library(annotate)

SeqX<-makeBaseCalls(ITSseq) # Call 

SeqX@primarySeq

SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
head(SeqXBlastDF)

fileNamesC <- c("C1_F_D07_029_2022-09-21.ab1", "C1_R_A01_008_2022-09-21.ab1","C2_F_E07_028_2022-09-21.ab1","C2_R_B01_007_2022-09-21.ab1","C3_F_F07_027_2022-09-21.ab1","C3_R_C01_006_2022-09-21.ab1",
               "C5_F_G07_026_2022-09-21.ab1","C5_R_D01_005_2022-09-21.ab1","C7_F_H07_025_2022-09-21.ab1",
               "C7_R_E01_004_2022-09-21.ab1","C11_F_A09_040_2022-09-21.ab1","C11_R_F01_003_2022-09-21.ab1",
               "C13_F_B09_039_2022-09-21.ab1","C13_R_G01_002_2022-09-21.ab1")

for(i in 1:6){
  seq<-read.abif(fileNamesC[i])
  assign(paste("seq",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seq",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seq",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blast.DF",i,sep=""),SeqXBlastDF)
}

#OK, so we have an issue with the 3rd sequence, 5 and 6 too...

for(i in 9:14){
  seq<-read.abif(fileNamesC[i])
  assign(paste("seq",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seq",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seq",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blast.DF",i,sep=""),SeqXBlastDF)
}

for(i in 7:8){
  seq<-read.abif(fileNamesC[i])
  assign(paste("seq",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seq",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seq",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blast.DF",i,sep=""),SeqXBlastDF)
}

# We are having issues with some of the csvs... Can we remove some columns?

write.csv(blast.DF1, file="Blast_1.csv")
blast.DF1.2<-blast.DF1[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF1.2,"Blast_1.2.xlsx")

write.csv(blast.DF2, file="Blast_2.csv")
blast.DF2.2<-blast.DF2[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF2.2,"Blast_2.2.xlsx")

write.csv(blast.DF3, file="Blast_3.csv")
#Failed, try again?

write.csv(blast.DF4, file="Blast_4.csv")
blast.DF4.2<-blast.DF4[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF4.2,"Blast_4.2.xlsx")

write.csv(blast.DF5, file="Blast_5.csv")
#Failed, try again?

write.csv(blast.DF6, file="Blast_6.csv")
#Failed, try again?

write.csv(blast.DF7, file="Blast_7.csv")
blast.DF7.2<-blast.DF7[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF7.2,"Blast_7.2.xlsx")

write.csv(blast.DF8, file="Blast_8.csv")
blast.DF8.2<-blast.DF8[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF8.2,"Blast_8.2.xlsx")

write.csv(blast.DF9, file="Blast_9.csv")
#Failed, try again?

write.csv(blast.DF10, file="Blast_10.csv")
blast.DF10.2<-blast.DF10[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF10.2,"Blast_10.2.xlsx")


write.csv(blast.DF11, file="Blast_11.csv")
blast.DF11.2<-blast.DF11[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF11.2,"Blast_11.2.xlsx")

write.csv(blast.DF12, file="Blast_12.csv")
blast.DF12.2<-blast.DF12[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF12.2,"Blast_12.2.xlsx")

write.csv(blast.DF13, file="Blast_13.csv")
#Failed, try again?

write.csv(blast.DF14, file="Blast_14.csv")
blast.DF14.2<-blast.DF14[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blast.DF14.2,"Blast_14.2.xlsx")

############################################################################################

fileNamesW <- c("W4_F_C09_038_2022-09-21.ab1", "W4_R_H01_001_2022-09-21.ab1","W5_F_D09_037_2022-09-21.ab1","W5_R_A03_016_2022-09-21.ab1","W6_F_E09_036_2022-09-21.ab1","W6_R_B03_015_2022-09-21.ab1",
                "W7_F_F09_035_2022-09-21.ab1","W7_R_C03_014_2022-09-21.ab1","W7.2_F_G09_034_2022-09-21.ab1",
                "W7.2_R_D03_013_2022-09-21.ab1","W9B_F_H09_033_2022-09-21.ab1","W9B_R_E03_012_2022-09-21.ab1",
                "W10_F_A11_048_2022-09-21.ab1","W10_R_F03_011_2022-09-21.ab1","W15_F_B11_047_2022-09-21.ab1",
                "W15_R_G03_010_2022-09-21.ab1","W17_F_C11_046_2022-09-21.ab1","W17_R_H03_009_2022-09-21.ab1")
  
for(i in 1:4){
  seq<-read.abif(fileNamesW[i])
  assign(paste("seqW",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqW",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqW",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastW.DF",i,sep=""),SeqXBlastDF)
}

write.csv(blastW.DF1, file="BlastW_1.csv")
blastW.DF1.2<-blastW.DF1[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF1.2,"BlastW_1.2.xlsx")

write.csv(blastW.DF2, file="BlastW_2.csv")
blastW.DF2.2<-blastW.DF2[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF2.2,"BlastW_2.2.xlsx")

write.csv(blastW.DF3, file="BlastW_3.csv")
blastW.DF3.2<-blastW.DF3[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF3.2,"BlastW_3.2.xlsx")

write.csv(blastW.DF4, file="BlastW_4.csv")
blastW.DF4.2<-blastW.DF4[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF4.2,"BlastW_4.2.xlsx")

for(i in 5:8){
  seq<-read.abif(fileNamesW[i])
  assign(paste("seqW",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqW",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqW",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastW.DF",i,sep=""),SeqXBlastDF)
}

write.csv(blastW.DF5, file="BlastW_5.csv")
blastW.DF5.2<-blastW.DF5[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF5.2,"BlastW_5.2.xlsx")

write.csv(blastW.DF6, file="BlastW_6.csv")
blastW.DF6.2<-blastW.DF6[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF6.2,"BlastW_6.2.xlsx")

write.csv(blastW.DF7, file="BlastW_7.csv")
#Didn't work, no BLAST hits.

write.csv(blastW.DF8, file="BlastW_8.csv")
blastW.DF8.2<-blastW.DF8[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF8.2,"BlastW_8.2.xlsx")

for(i in 9:12){
  seq<-read.abif(fileNamesW[i])
  assign(paste("seqW",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqW",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqW",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastW.DF",i,sep=""),SeqXBlastDF)
}

#12 didn't work

write.csv(blastW.DF9, file="BlastW_9.csv")
blastW.DF9.2<-blastW.DF9[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF9.2,"BlastW_9.2.xlsx")

write.csv(blastW.DF10, file="BlastW_10.csv")
blastW.DF10.2<-blastW.DF10[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF10.2,"BlastW_10.2.xlsx")

write.csv(blastW.DF11, file="BlastW_11.csv")
blastW.DF11.2<-blastW.DF11[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF11.2,"BlastW_11.2.xlsx")

for(i in 13:16){
  seq<-read.abif(fileNamesW[i])
  assign(paste("seqW",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqW",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqW",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastW.DF",i,sep=""),SeqXBlastDF)
}

# 16 didn't work

write.csv(blastW.DF13, file="BlastW_13.csv")
blastW.DF13.2<-blastW.DF13[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF13.2,"BlastW_13.2.xlsx")

write.csv(blastW.DF14, file="BlastW_14.csv")
blastW.DF14.2<-blastW.DF14[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF14.2,"BlastW_14.2.xlsx")

write.csv(blastW.DF16, file="BlastW_16.csv")
blastW.DF16.2<-blastW.DF16[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF16.2,"BlastW_16.2.xlsx")

for(i in 17:18){
  seq<-read.abif(fileNamesW[i])
  assign(paste("seqW",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqW",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqW",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastW.DF",i,sep=""),SeqXBlastDF)
}

write.csv(blastW.DF18, file="BlastW_18.csv")
blastW.DF18.2<-blastW.DF18[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastW.DF18.2,"BlastW_18.2.xlsx")

####################################################################

fileNamesR <- c("R1_F_D11_045_2022-09-21.ab1", "R1_R_A05_024_2022-09-21.ab1","R2_F_E11_044_2022-09-21.ab1",
                "R2_R_B05_023_2022-09-21.ab1","R3_F_F11_043_2022-09-21.ab1","R3_R_C05_022_2022-09-21.ab1",
                "R4_F_G11_042_2022-09-21.ab1","R4_R_D05_021_2022-09-21.ab1","R5_F_H11_041_2022-09-21.ab1",
                "R5_R_E05_020_2022-09-21.ab1","R6_F_A02_008_2022-09-21.ab1","R6_R_F05_019_2022-09-21.ab1",
                "R7_F_B02_007_2022-09-21.ab1","R7_R_G05_018_2022-09-21.ab1","R8_F_C02_006_2022-09-21.ab1",
                "R8_R_H05_017_2022-09-21.ab1","R9_F_D02_005_2022-09-21.ab1","R9_R_A07_032_2022-09-21.ab1",
                "R10_F_E02_004_2022-09-21.ab1","R10_R_B07_031_2022-09-21.ab1",
                "R11_F_F02_003_2022-09-21.ab1","R11_R_C07_030_2022-09-21.ab1")

for(i in 1:4){
  seq<-read.abif(fileNamesR[i])
  assign(paste("seqR",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqR",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqR",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastR.DF",i,sep=""),SeqXBlastDF)
}

write.csv(blastR.DF1, file="BlastR_1.csv")
blastR.DF1.2<-blastR.DF1[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF1.2,"BlastR_1.2.xlsx")

write.csv(blastR.DF2, file="BlastR_2.csv")
blastR.DF2.2<-blastR.DF2[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF2.2,"BlastR_2.2.xlsx")

write.csv(blastR.DF4, file="BlastR_4.csv")
blastR.DF4.2<-blastR.DF4[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF4.2,"BlastR_4.2.xlsx")

for(i in 6:8){
  seq<-read.abif(fileNamesR[i])
  assign(paste("seqR",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqR",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqR",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastR.DF",i,sep=""),SeqXBlastDF)
}

# I aborted 5

write.csv(blastR.DF6, file="BlastR_6.csv")
blastR.DF6.2<-blastR.DF6[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF6.2,"BlastR_6.2.xlsx")

write.csv(blastR.DF7, file="BlastR_7.csv")
blastR.DF7.2<-blastR.DF7[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF7.2,"BlastR_7.2.xlsx")

write.csv(blastR.DF8, file="BlastR_8.csv")
blastR.DF8.2<-blastR.DF8[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF8.2,"BlastR_8.2.xlsx")

for(i in 9:12){
  seq<-read.abif(fileNamesR[i])
  assign(paste("seqR",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqR",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqR",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastR.DF",i,sep=""),SeqXBlastDF)
}

#11 failed.

write.csv(blastR.DF9, file="BlastR_9.csv")
blastR.DF9.2<-blastR.DF9[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF9.2,"BlastR_9.2.xlsx")

write.csv(blastR.DF10, file="BlastR_10.csv")
blastR.DF10.2<-blastR.DF10[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF10.2,"BlastR_10.2.xlsx")

write.csv(blastR.DF12, file="BlastR_12.csv")
blastR.DF12.2<-blastR.DF12[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF12.2,"BlastR_12.2.xlsx")

for(i in 13:16){
  seq<-read.abif(fileNamesR[i])
  assign(paste("seqR",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqR",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqR",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastR.DF",i,sep=""),SeqXBlastDF)
}

#15 didn't really work

write.csv(blastR.DF13, file="BlastR_13.csv")
blastR.DF13.2<-blastR.DF13[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF13.2,"BlastR_13.2.xlsx")

write.csv(blastR.DF14, file="BlastR_14.csv")
blastR.DF14.2<-blastR.DF14[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF14.2,"BlastR_14.2.xlsx")

write.csv(blastR.DF15, file="BlastR_15.csv")
blastR.DF15.2<-blastR.DF15[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF15.2,"BlastR_15.2.xlsx")

write.csv(blastR.DF16, file="BlastR_16.csv")
blastR.DF16.2<-blastR.DF16[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF16.2,"BlastR_16.2.xlsx")

for(i in 17:20){
  seq<-read.abif(fileNamesR[i])
  assign(paste("seqR",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqR",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqR",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastR.DF",i,sep=""),SeqXBlastDF)
}

# No hits for 20

write.csv(blastR.DF17, file="BlastR_17.csv")
blastR.DF17.2<-blastR.DF17[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF17.2,"BlastR_17.2.xlsx")

write.csv(blastR.DF18, file="BlastR_18.csv")
blastR.DF18.2<-blastR.DF18[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF18.2,"BlastR_18.2.xlsx")

write.csv(blastR.DF19, file="BlastR_19.csv")
blastR.DF19.2<-blastR.DF19[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF19.2,"BlastR_19.2.xlsx")

write.csv(blastR.DF20, file="BlastR_20.csv")
blastR.DF20.2<-blastR.DF20[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF20.2,"BlastR_20.2.xlsx")

for(i in 21:22){
  seq<-read.abif(fileNamesR[i])
  assign(paste("seqR",i,sep=""),seq)
  
  san.seq <- sangerseq(seq)
  assign(paste("san.seqR",i,sep=""),san.seq)
  
  SeqX<-makeBaseCalls(san.seq) # Call 
  
  assign(paste("prim_seqR",i,sep=""),SeqX@primarySeq)
  
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame')
  assign(paste("blastR.DF",i,sep=""),SeqXBlastDF)
}

write.csv(blastR.DF21, file="BlastR_21.csv")
blastR.DF21.2<-blastR.DF21[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF21.2,"BlastR_21.2.xlsx")

write.csv(blastR.DF22, file="BlastR_22.csv")
blastR.DF22.2<-blastR.DF22[,c(10,11,9,14,15,24,26,27)]
write_xlsx(blastR.DF22.2,"BlastR_22.2.xlsx")

