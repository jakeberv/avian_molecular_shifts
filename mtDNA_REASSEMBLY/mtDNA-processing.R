require(readxl)

#https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    #if the value is NA, skip it
    if(!is.na(as.character(data[rowNum,"seq"])) == TRUE){
    	fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    	fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
    }
    
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


data<-read_excel(path="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/MITOGENOME-Berv-Aug-2021.xlsx", sheet="jake_mod")

#data must be sorted according to mtdna sort order

mtdna.12s <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`seq 1`, 200))
mtdna.16s <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`seq 2`, 200))
mtdna.ND1 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND1_seq`, 200))
mtdna.ND2 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND2_seq`, 200))
mtdna.COX1 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`COX1_seq`, 200))
mtdna.COX2 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`COX2_seq`, 200))
mtdna.ATP8 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ATP8_seq`, 200))
mtdna.ATP6 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ATP6_seq`, 200))
mtdna.COX3 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`COX3_seq`, 200))
mtdna.ND3 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND3_seq`, 200))
mtdna.ND4L <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND4L_seq`, 200))
mtdna.ND4 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND4_seq`, 200))
mtdna.ND5 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND5_seq`, 200))
mtdna.CYTB <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`CYTB_seq`, 200))
mtdna.ND6 <- data.frame(name=tail(data$`NAMES FOR EXPORT`,200), seq=tail(data$`ND6_seq`, 200))

#write the rRNAs
writeFasta(data=mtdna.12s, filename="mtdna.12s.fasta")
writeFasta(data=mtdna.16s, filename="mtdna.16s.fasta")
#write the protein coding genes
writeFasta(data=mtdna.ND1, filename="mtdna.ND1.fasta")
writeFasta(data=mtdna.ND2, filename="mtdna.ND2.fasta")
writeFasta(data=mtdna.COX1, filename="mtdna.COX1.fasta")
writeFasta(data=mtdna.COX2, filename="mtdna.COX2.fasta")
writeFasta(data=mtdna.ATP8, filename="mtdna.ATP8.fasta")
writeFasta(data=mtdna.ATP6, filename="mtdna.ATP6.fasta")
writeFasta(data=mtdna.COX3, filename="mtdna.COX3.fasta")
writeFasta(data=mtdna.ND3, filename="mtdna.ND3.fasta")
writeFasta(data=mtdna.ND4L, filename="mtdna.ND4L.fasta")
writeFasta(data=mtdna.ND4, filename="mtdna.ND4.fasta")
writeFasta(data=mtdna.ND5, filename="mtdna.ND5.fasta")
writeFasta(data=mtdna.CYTB, filename="mtdna.CYTB.fasta")
writeFasta(data=mtdna.ND6, filename="mtdna.ND6.fasta")








