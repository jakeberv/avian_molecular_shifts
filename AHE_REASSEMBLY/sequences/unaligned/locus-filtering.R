

require(seqinr)
require(readxl)
require(stringi)

#path to directory where fasta files are
filepath<-"~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/min2x/qual20_cov2_haplo"

#path to directory where fasta files will be written
outpath<-"~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/min2x/min50bp_min10p"

#this file path goes to the "best only" haplotype directory

#get list of all file paths in the directory
file.dir<-grep("[.]fasta",list.files(filepath,full.names=T),value=T)

#get list of shortened file paths
file.dir.short<-grep("[.]fasta",list.files(filepath,full.names=F),value=T)

#create counter for file.dir
count<-seq(1, length(file.dir))

#read in the excel doc with translation table of codes and names
translation<-read_excel(path="~/jsb439@cornell.edu/AnchoredEnrichment/RY-coding/translation-table.xlsx", sheet="Sheet3")

#function to rename and sort sequences
renamer<-function(fas){
	
	#store names of sequences in order
	tmp<-names(fas)
		
	#create a data frame with rownames equal to sequence names
	tmp<-as.data.frame(tmp)
	rownames(tmp)<-tmp[,1]
	
	#create new data frame of new labels from spreadsheet
	tmp2<-as.data.frame(translation$new.tip.labels.haplo)
	rownames(tmp2)<-translation$codes.haplo
	
	#use rownames to swap codes for new tip labels
	reordered.names<-as.character(as.data.frame(tmp2[rownames(tmp),])[,1])
	
	#return new vector of names in order
	#return(reordered.names)
	
	names(fas)<-reordered.names
	
	#sort sequences by name
	sort<-sort(names(fas))
	fas<-fas[sort]
	
	return(fas)
}


#set up for loop to operate on directory of files
for(i in 1:length(file.dir)){
	
	print(paste("operating on ", file.dir.short[i]))
	print(paste((count[i]/length(file.dir)*100), " % complete"))
	
	#testing code for filtering
	
	fas<-read.fasta(file.dir[i], forceDNAtolower=F, as.string=T, set.attributes=F)

	remove<-names(fas) %in% translation$codes.haplo

	print(paste(sum(as.numeric(names(fas) %in% translation$codes.haplo))," sequences retained", " out of ", length(names(fas)), " sequences"))

	#remove sequences not in spreadsheet master list
	fas.pruned<-fas[remove]

	#output renamed sequences in same order as fas.pruned
	fas.renamed <- renamer(fas.pruned)


	#trim leading/trailing Ns
	for (j in names(fas.renamed)){
		fas.renamed[names(fas.renamed)==j] <- trimws(fas.renamed[names(fas.renamed)==j], whitespace="N")
	}

	#first pass filter to remove sequences that are now zero length
	#filter out sequences within each fasta file that are less than 1 bp
	bp<-1 #set bp filter

	for (j in names(fas.renamed)){
	
		totalcount<-nchar(fas.renamed[[j]])
		
		if ((totalcount) < bp){
		
			fas.renamed[names(fas.renamed)==j]<-NULL
		
		#print(paste(j, " sequence was removed because is is zero length after trimming"))
		}
	}
	
	#second pass filter
	#filter out sequences within each fasta file that have high % of N
	Nfilter<-0.4 #set N filter (above this % threshold, sequence is deleted)

	for (j in names(fas.renamed)){
	
		Ncount<-lengths(regmatches(fas.renamed[names(fas.renamed)==j], gregexpr("N", fas.renamed[names(fas.renamed)==j])))
		totalcount<-nchar(fas.renamed[[j]])
		
		if ((Ncount/totalcount) > Nfilter){
		
			fas.renamed[names(fas.renamed)==j]<-NULL
		
		print(paste(j, " sequence was removed because high %N"))
		}
	}
	
	#third pass filter to remove short sequences which remain
	#filter out sequences within each fasta file that are less than 100 bp
	#can be more stringent later on
	
	bp<-50 #set bp filter (if less than this threshold, removed)

	for (j in names(fas.renamed)){
	
		totalcount<-nchar(fas.renamed[[j]])
		
		if ((totalcount) < bp){
		
			fas.renamed[names(fas.renamed)==j]<-NULL
		
		print(paste(j, " sequence was removed because < ", paste(bp), " bp"))
		}
	}


	#set up output files
	out.string<-paste(outpath,file.dir.short[i],sep="/")
	
	
	#fourth pass filter to require at least 10% of target individuals to have data
	#can be more stringent later on
	individualfilter<-0.10
	
	#length(names(fas))
	
	if (length(names(fas.renamed)) > (individualfilter * 400)){
		
		write.fasta(fas.renamed,file.out=out.string, names=names(fas.renamed), as.string=T)

	}  else {
		print(paste(file.dir.short[i], " did not contain enough data"))
	}
	
}



#as.SeqFastadna #can return GC content