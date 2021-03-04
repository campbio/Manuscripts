#################
# read.maf
#
# Read in a MAF file into a data.frame, data.table, or GenomicRanges object
#
#################
read.maf = function(file, class=c("dt", "df", "gr"), keys=c("Chromosome", "Start_position", "End_position", "Strand", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"), standardize=TRUE, showProgress=FALSE, ...) {
	suppressPackageStartupMessages(require(data.table))
	suppressPackageStartupMessages(require(GenomicRanges))
    
	class = match.arg(class)
	
	if(class == "df") {
	
	  temp.maf = suppressWarnings(fread(file, data.table=FALSE, showProgress=showProgress, ...))
  	  if(standardize == TRUE) { standardize.maf.columns(temp.maf) }
  	  
	} else if (class == "dt" | class == "gr") {
	
	  temp.maf = suppressWarnings(fread(file, data.table=TRUE, showProgress=showProgress, ...))
  	  if(standardize == TRUE) { standardize.maf.columns(temp.maf) } ## Standardize before setting keys
	  setkeyv(temp.maf, keys)	  
	  
    }
    
    if(class == "gr") {
      maf.gr = GRanges(seqnames=Rle(temp.maf$Chromosome), IRanges(names=temp.maf$Chromosome, start=temp.maf$Start_position, end=temp.maf$End_position), strand=temp.maf$Strand, temp.maf)
      return(maf.gr)
    } else {
      return(temp.maf)
    }
}



##################
# standarize.maf.columns
#
# Some programs such as MutSig will change the column names of the maf. This function automatically checks and fixes the columns
#
##################
standardize.maf.columns = function(maf) {
  colnames.fix = c("Start_Position", "End_Position")
  colnames.correct = c("Start_position", "End_position")

  colnames.match = match(colnames.fix, colnames(maf))
  ix = !is.na(colnames.match)
  if(sum(ix) > 0) { colnames(maf)[colnames.match[ix]] = colnames.correct[ix] }
  return(maf)
}




##################
# maf.coding
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a coding mutation
#
##################
maf.coding = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins")	

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}



#################
# maf.exonic
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a coding or a silent mutation
#
##################
maf.exonic = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("Silent", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins")	

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}



#################
# maf.genic
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a coding, silent, UTR, or intronic mutation
#
##################
maf.genic = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("3'UTR", "5'UTR", "Intron", "Silent", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins")	

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}


#################
# maf.truncating
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a truncating mutations (not including splice_site)
#
##################
maf.truncating = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "De_novo_Start_OutOfFrame")

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}


#################
# lof.pval
#
# Determines if a gene in each maf is enriched for loss of function (i.e. truncating) mutations 
# using the Fisher's exact test.
#
##################
lof.pval = function(maf) {
	maf$Hugo_Symbol = as.factor(maf$Hugo_Symbol)
 
	t.ix = maf.truncating(maf, index=TRUE)

	## Calculate totals for numbers of LoF mutations vs. other mutations in each gene
	gene.tot = aggregate(rep(1, nrow(maf)), by=list(maf$Hugo_Symbol), sum)
	gene.tr = aggregate(t.ix, by=list(maf$Hugo_Symbol), sum)
	gene.nt = aggregate(!t.ix, by=list(maf$Hugo_Symbol), sum)
	gene.other = t(sapply(1:nrow(gene.tot), function(i) c(sum(gene.tr$x[-i]), sum(gene.nt$x[-i]))))

	gene.ta = cbind(n1=gene.tr$x, n2=gene.nt$x, n3=gene.other[,1], n4=gene.other[,2])
	rownames(gene.ta) = gene.tot[,1]

	## Apply Fisher's exact test to each row
	fet = function(v, alt="greater") {
	  f = fisher.test(matrix(v, ncol=2), alternative=alt)
	  d = c(OR=as.numeric(f$estimate), Pvalue=as.numeric(f$p.value))
	}
	res = t(apply(gene.ta, 1, fet))

	res = data.frame(Gene=gene.tot[,1], gene.ta, res)
	colnames(res) = c("Gene", "LOF_in_Gene", "Other_Mut_in_Gene", "Total_LOF_minus_Gene", "Total_Other_Mut_minus_Gene", "OR", "P")
	return(res)
}




make.maf.id = function(maf, type=c("site", "patient", "custom"), sep="_", custom=NULL) {

  ## Takes a maf and generates an ID that should be unique for each site or patient
  ## Written by: Josh Campbell
  ## 2-4-2016

  ## Parameters: ##########################################################################################
  # 
  # maf			A data.frame or matrix with Mutation Annotation Format
  # type		Determines the columns that will be used to make the ID. 
  #				site: Chromosome, Start_position, End_position, Reference_Allele, Tumor_Seq_Allele2
  #				patient: Same as site but with the addition of Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode to make the ID scheme unique to each patient
  #				custom: Uses column names from the "custom" variable.
  #	sep			String used to collapse the selcted columns for each row of the maf
  #	custom		Can pass a custom set of columns to use for collapsing the maf rows.
  #
  #########################################################################################################

  ## Values: ###############################################################################################
  # make.maf.id
  #
  # Returns a vector equal to the number of rows in the maf where each entry is a combination
  # of the column scheme chosen with type
  #
  #########################################################################################################

  type = match.arg(type)
  
  if(type == "site") {
    cols = c("Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele2")
  } else if (type == "patient") {
    cols = c("Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode")
  } else if (type == "custom") {
    cols = custom
    if (length(setdiff(cols, colnames(maf))) > 0 | is.null(custom)) {
      stop("Customs columns to generate ids need to match colnames in maf.")
    }
  }
  
  maf.site.id = apply(maf[,cols], 1, paste, collapse=sep)
  maf.site.id = gsub(" ", "", maf.site.id)
  return(maf.site.id)
}  






###################
# mutation.context.snv192
#
# Summarizes the 192 mutation contexts per sample
#
##################

mutation.context.snv192 = function(maf) {

  require(Biostrings)
  
  final.mut.type = rep(NA, nrow(maf))
  final.mut.context = rep(NA, nrow(maf))
  final.mut.strand = rep(NA, nrow(maf))
    
  ## Get mutation type
  initial.maf.type = paste(maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, sep="")
  
  ## Get mutation context info for those on "+" strand
  forward.change = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind = maf$Variant_Type == "SNP" & initial.maf.type %in% forward.change

  final.mut.type[ind] = paste(maf$Reference_Allele[ind], ">", maf$Tumor_Seq_Allele2[ind], sep="")
  final.mut.context[ind] = toupper(substring(as.character(maf$ref_context[ind]), 10,12))
  final.mut.strand[ind] = ifelse(maf$Transcript_Strand[ind] == "+", "U", "T")

  ## Get mutation context info for those on "-" strand
  rev.change = c("A>G","A>T","A>C","G>T","G>C","G>A")
  ind = maf$Variant_Type == "SNP" & initial.maf.type %in% rev.change

  ## Reverse complement the context so only 6 mutation categories instead of 12
  rev.context = reverseComplement(DNAStringSet(maf$ref_context[ind]))
  rev.refbase = reverseComplement(DNAStringSet(maf$Reference_Allele[ind]))
  rev.altbase = reverseComplement(DNAStringSet(maf$Tumor_Seq_Allele2[ind]))

  final.mut.type[ind] = paste(as.character(rev.refbase), ">", as.character(rev.altbase), sep="")
  final.mut.context[ind] = toupper(substring(as.character(rev.context), 10,12))
  final.mut.strand[ind] = ifelse(maf$Transcript_Strand[ind] == "-", "U", "T")
  
  
  maf.mut.id = paste(final.mut.type, final.mut.context, final.mut.strand, sep="_")
  Tumor_ID = as.factor(maf$Tumor_Sample_Barcode)

  ## Define all mutation types for 196 substitution scheme
  b1 = rep(rep(c("A", "C", "G", "T"), each=24), 2)
  b2 = rep(rep(c("C", "T"), each=12), 8)
  b3 = rep(c("A", "C", "G", "T"), 48)
  mut.trinuc = apply(cbind(b1, b2, b3), 1, paste, collapse="")
  mut.type = rep(rep(rep(forward.change, each=4), 4), 2)
  mut.strand = rep(c("T", "U"), each=96)

  mut.id = apply(cbind(mut.type, mut.trinuc, mut.strand), 1, paste, collapse="_")  

  Mutation = factor(maf.mut.id, levels=mut.id)
  
  mut.table = xtabs(~ Mutation + Tumor_ID)
  mut.summary = data.frame(Mutation, Type=final.mut.type, Context=final.mut.context, Strand=final.mut.strand, stringsAsFactors=FALSE)
  
  return(list(mutation_table=mut.table, maf_mutations=mut.summary))  
}


###################
# mutation.context.snv96
#
# Summarizes the 96 mutation contexts per sample
#
##################

mutation.context.snv96 = function(maf) {

  require(Biostrings)
  
  final.mut.type = rep(NA, nrow(maf))
  final.mut.context = rep(NA, nrow(maf))
  final.mut.strand = rep(NA, nrow(maf))
    
  ## Get mutation type
  initial.maf.type = paste(maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, sep="")
  
  ## Get mutation context info for those on "+" strand
  forward.change = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind = maf$Variant_Type == "SNP" & initial.maf.type %in% forward.change

  final.mut.type[ind] = paste(maf$Reference_Allele[ind], ">", maf$Tumor_Seq_Allele2[ind], sep="")
  final.mut.context[ind] = toupper(substring(as.character(maf$ref_context[ind]), 10,12))

  ## Get mutation context info for those on "-" strand
  rev.change = c("A>G","A>T","A>C","G>T","G>C","G>A")
  ind = maf$Variant_Type == "SNP" & initial.maf.type %in% rev.change

  ## Reverse complement the context so only 6 mutation categories instead of 12
  rev.context = reverseComplement(DNAStringSet(maf$ref_context[ind]))
  rev.refbase = reverseComplement(DNAStringSet(maf$Reference_Allele[ind]))
  rev.altbase = reverseComplement(DNAStringSet(maf$Tumor_Seq_Allele2[ind]))

  final.mut.type[ind] = paste(as.character(rev.refbase), ">", as.character(rev.altbase), sep="")
  final.mut.context[ind] = toupper(substring(as.character(rev.context), 10,12))
  
  
  maf.mut.id = paste(final.mut.type, final.mut.context, sep="_")
  Tumor_ID = as.factor(maf$Tumor_Sample_Barcode)

  ## Define all mutation types for 96 substitution scheme
  b1 = rep(c("A", "C", "G", "T"), each=24)
  b2 = rep(rep(c("C", "T"), each=12), 4)
  b3 = rep(c("A", "C", "G", "T"), 24)
  mut.trinuc = apply(cbind(b1, b2, b3), 1, paste, collapse="")
  mut.type = rep(rep(forward.change, each=4), 4)

  mut.id = apply(cbind(mut.type, mut.trinuc), 1, paste, collapse="_")  

  Mutation = factor(maf.mut.id, levels=mut.id)
  
  mut.table = xtabs(~ Mutation + Tumor_ID)
  mut.summary = data.frame(Mutation, Type=final.mut.type, Context=final.mut.context, stringsAsFactors=FALSE)
  
  return(list(mutation_table=mut.table, maf_mutations=mut.summary))  
}

