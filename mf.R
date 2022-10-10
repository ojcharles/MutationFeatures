# Generate for each change at each residue, the mutation features.
# These can then be fed into drugres classification etc.


# -------------------- Setup
### runtime vars
infasta = paste0("/query/HCMV_UL54.fasta")
#v_evals = c("1E-100", "1E-50", "1e-7")
v_eval = 1e-7

# passed arguments
# //todo make a rgument
blast_db_name = "uniref50_virus.fasta"

library(stringr)
library(Peptides)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyr)
library(readr)
library(FSelector)
library(protr)
library(ape)
library(Biostrings)

#source("R/functions.R")


# ------------------------------------------------------------ setup empty df
inseq = as.character(unlist(readFASTA(infasta)))
inseq_vector = as.vector(str_split_fixed(inseq, pattern = "", n = nchar(inseq)))

# generate a table of location, wild type, mutant type residue
locs = 1:length(inseq_vector)
amino_acids = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M","N", "P", "Q", "R", "S", "T", "V", "W", "Y")
wt_mt = expand.grid(inseq_vector, amino_acids)
df = cbind(locs, wt_mt)
colnames(df) = c("loc", "wt", "mt")
head(df)



##### is there a PDB?
pdb_file = gsub(".fasta", ".pdb", infasta )
use_pdb = file.exists(pdb_file)



# ------------------------------------------------------------Evolutionary Features
# -------------------- PSSM

tdir = "/tmp"
pssm_file = paste0(tdir, "/pssm.txt")
temp_blast = paste0(tdir,"/psiblast.fa")
temp_blast_msa = paste0(tdir,"/psiblast_msa.fa")

# search for homology
command = paste0("psiblast -num_threads 6 -query ", infasta," -db db/", blast_db_name, " -num_iterations 2 -out_ascii_pssm ",
                              pssm_file," -save_pssm_after_last_round  -out ",temp_blast," -outfmt '6 qseqid sseqid sseq qstart' -inclusion_ethresh ",v_eval ," -evalue ", v_eval," ;")
system(command)

# convert blastoutput to fasta
out_fasta = c("")
t = read.table(temp_blast)
for(i in 1:nrow(t)){
  seq = paste0(paste0(rep("-", t$V4[i] - 1), collapse = ""), t$V3[i])
  out_fasta = c(out_fasta, paste0(">", t$V2[i]))
  out_fasta = c(out_fasta, seq)
}
writeLines(out_fasta, temp_blast)

# align outputted fasta
command = paste0("mafft --add  ", temp_blast," --keeplength --thread 4 ",infasta, " > ", temp_blast_msa," 2>",tdir,"/err.txt" )
system(command)


# function to take pssm file -> pssm score matrix
pssmfile2df = function(pssm_file){
  x=read.delim(pssm_file,skip = 2,sep = "",header = FALSE)
  cols = x[1,1:20]
  x=x[-1,-c(1,23:44)]
  d=which(x=="Lambda")
  if(length(d)!=0){
    x=x[-c(d:dim(x)[1]),]
  }
  x=x[,-1]
  colnames(x)=cols
  rownames(x)=NULL
  return(x)
}
pssm = pssmfile2df(pssm_file)


# //todo slow, profile.
df$pssm_wt = 0; df$pssm_mt = 0; df$pssm_diff = 0;
cols = colnames(pssm)
for(r in 1:nrow(df)){
  pos = df$loc[r]
  wtAA = df$wt[r]
  mtAA = df$mt[r]
  sus2 = pssm[pos,]
  wtcol = grep(wtAA, cols)
  mtcol = grep(mtAA, cols)
  df$pssm_wt[r] = pssm[pos,wtcol]
  df$pssm_mt[r] = pssm[pos,mtcol]
  df$pssm_diff[r] = abs(as.numeric(df$pssm_mt[r]) - as.numeric(df$pssm_wt[r]))
  df$pssm_mean[r] = mean(as.numeric(pssm[pos,]))
}


# -------------------- sequence conservation
# # https://compbio.cs.princeton.edu/conservation/
command = paste0("python2 /mflibs/conservation_code/score_conservation.py -m /mflibs/conservation_code/matrix/blosum62.bla -p FALSE -g 0.99 ",temp_blast_msa," > ",tdir,"/conservation.txt")
system(command)
conservation = data.frame(read.table(paste0(tdir,"/conservation.txt"),header = F, sep = "\t")[,1:2])
colnames(conservation) = c("loc", "jsdiv")
conservation$jsdiv = as.numeric(conservation$jsdiv)
conservation[conservation$jsdiv < 0,2] = 0 # handle NA
df = merge(df, conservation, by = "loc", all.x = T)


### psi-blast conservation score
# //todo depthnorm could do with a more informative value
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
a = Biostrings::readAAStringSet(temp_blast_msa)
a = as.matrix(a)
psi_depth = data.frame(loc = 1:ncol(a), psi_depth = 0, psi_unique = 0)
for(i in 1:ncol(a)){
  t = as.character(a[,i])
  t = t[t != "-"]
  psi_depth$psi_depth[i] = round(length(t),0)
  psi_depth$psi_depthnorm[i] = round(length(t) / nrow(a),2)
  psi_depth$psi_unique[i] = length(unique(t))
  
}
df = merge(df, psi_depth, by = "loc", all.x = T)

  
 
  
  

  
  
  
  
#--------------------  grantham / blossum
# grantham
grantham = readr::read_tsv("/mflibs/grantham.tsv") %>%
  tidyr:: gather(SECOND,SCORE, -FIRST) %>% dplyr::filter(SCORE > 0)  
df$grantham = 0
for(i in 1:nrow(grantham)){
  wt = grantham$FIRST[i]
  mt = grantham$SECOND[i]
  df[df$wt ==wt & df$mt == mt,]$grantham = grantham$SCORE[i]
  df[df$mt ==wt & df$wt == mt,]$grantham = grantham$SCORE[i]
}

#blossum
data(AABLOSUM62)
blosum = data.frame(AABLOSUM62)
blosum$wt = rownames(blosum)
blosum = reshape2::melt(blosum, id.vars = "wt")
colnames(blosum) = c("wt", "mt", "blosum")
df$blosum62 = 0
for(i in 1:nrow(blosum)){
  wt = blosum$wt[i]
  mt = blosum$mt[i]
  df[df$wt ==wt & df$mt == mt,]$blosum62 = blosum$blosum[i]
  df[df$mt ==wt & df$wt == mt,]$blosum62 = blosum$blosum[i]
}





#--------------------  physical properties
# wdwradius
# //todo vwd from, wvd to?
vdw = read.csv("/mflibs/vdw_radius.csv", skip = 1)
df$vdw_radius_delta = 0
for(r in 1:nrow(df)){
  wtAA = df$wt[r]
  mtAA = df$mt[r]
  df$vdw_radius_delta[r] = round(abs(vdw[which(vdw$AA == wtAA),2] - vdw[which(vdw$AA == mtAA),2]),1)
}


# dydrogen moment
df$from_hydro = hmoment(df$wt)
df$to_hydro = hmoment(df$mt)
df$change_hydro = abs(df$to_hydro - df$from_hydro)
# hydrophobicity
df$from_hydrophobicity = hydrophobicity(df$wt)
df$to_hydrophobicity = hydrophobicity(df$mt)
df$change_hydrophobicity = abs(df$to_hydrophobicity - df$from_hydrophobicity)
# isoelectric point
df$from_PI = pI(df$wt)
df$to_PI = pI(df$mt)
df$change_PI = df$to_hydro - df$from_hydro
# molweight
df$from_mw = mw(df$wt)
df$to_mw = mw(df$mt)
df$change_mw = abs(df$to_mw - df$from_mw)
# charge
df$from_charge = charge(df$wt)
df$to_charge = charge(df$mt)
df$change_charge = abs(df$to_charge - df$from_charge)







# ------------------------------------------------------------ Structural features
# #--------------------  DSSP
if(use_pdb){
  
  parse.dssp <- function(file){
    ## --------------- Reading the dssp file ------------------ ##
    con <- file(file, 'r')
    
    counter <- 0
    resnum <- c()
    respdb <- c()
    chain <- c()
    aa <- c()
    ss <- c()
    sasa <- c()
    phi <- c()
    psi <- c()
    
    while(TRUE){
      line <- readLines(con, n = 1)
      counter <- counter + 1
      
      if (counter == 1){
        l <- strsplit(line, split = "")[[1]]
        l <- paste(l, collapse = "")
        if ("have bz2" %in% l){
          first_valid_line <- 29 # dssp file coming from the API
        } else {
          first_valid_line <- 28 # dssp file coming from the sync
        }
      }
      
      if (counter > first_valid_line & length(line) != 0){
        a <- strsplit(line, split = "")[[1]]
        resnum <- c(resnum, paste(a[1:5], collapse = ""))
        respdb <- c(respdb, paste(a[6:10], collapse = ""))
        chain <- c(chain, paste(a[11:12], collapse = ""))
        aa <- c(aa, paste(a[13:14], collapse = ""))
        ss <- c(ss, paste(a[15:17], collapse = ""))
        sasa <- c(sasa, paste(a[36:38], collapse = ""))
        phi <- c(phi, paste(a[104:109], collapse = ""))
        psi <- c(psi, paste(a[110:115], collapse = ""))
      }
      if (length(line) == 0){
        break
      }
    }
    close(con)
    
    ## ------ Setting the variable types ------------- ##
    resnum <- as.numeric(resnum)
    respdb <- as.numeric(respdb)
    chain <- gsub(" ", "", chain)
    aa <- gsub(" ", "", aa)
    ss <- gsub("   ", "C", ss)
    ss <- gsub(" ", "", ss)
    
    ## -------- Building the dataframe ---------------- ##
    df <- as.data.frame(matrix(c(resnum, respdb, chain, aa,
                                ss, sasa, phi, psi), ncol = 8),
                        stringsAsFactors = FALSE)
    
    colnames(df) <- c('resnum', 'respdb', 'chain', 'aa', 'ss',
                      'asa', 'phi', 'psi')
    
    df$resnum <- as.numeric(df$resnum)
    df$respdb <- as.numeric(df$respdb)
    df$asa <- as.numeric(df$asa)
    df$phi <- as.numeric(df$phi)
    df$psi <- as.numeric(df$psi)
    
    ## --------------- Remove empty lines between chains ------------- ##
    badlines <- c()
    for (i in 1:nrow(df)){
      if (df$aa[i] == '!' | df$aa[i] == 'X'){
        badlines <- c(badlines, i)
      }
    }
    if (length(badlines) != 0){
      df <- df[-badlines,]
      df$resnum <- 1:nrow(df)
    }
    
    
    
    ## --------------- ASA -> RSA ------------- ##
    ##### rsa
    asa_lookup = c('A'=129.0, 'R'=274.0, 'N'=195.0, 'D'=193.0, 'C'=167.0,
                  'E'=223.0, 'Q'=225.0, 'G'=104.0, 'H'=224.0, 'I'=197.0,
                  'L'=201.0, 'K'=236.0, 'M'=224.0, 'F'=240.0, 'P'=159.0,
                  'S'=155.0, 'T'=172.0, 'W'=285.0, 'Y'=263.0, 'V'=174.0)
    
    # for each residue, get asa
    df$rsa = 0.0
    for(i in 1:nrow(df)){
      asa = df$asa[i]
      total_surface_area_for_residue = as.numeric(asa_lookup[which(names(asa_lookup) == df$aa[i])])
      df$rsa[i] = asa / total_surface_area_for_residue
    }
    
    return(df[,c(1,5,6,7,8,9)])
    
  }
  dssp_exec = "/root/.local/bin/mkdssp"
  t = readLines(pdb_file)
  missing_header = sum(grepl("HEADER", t)) == 0
  missing_cryst1 = sum(grepl("CRYST1", t)) == 0
  # push example header if misssing
  if( missing_header || missing_cryst1 ){
    t = c("HEADER    A HEADER                             01-JAN-22   1ABC  ",
        "TITLE     CRYSTAL STRUCTURE OF SOMETHING IMPORTANT",
        "CRYST1  103.917  125.550  220.577  90.00  90.00  90.00 P 21 21 21    8 ",
        t )
    pdb_file = gsub(".fasta", "fixed_.pdb", infasta )
  }
  writeLines(t, pdb_file)
  command = paste0(dssp_exec, " ",pdb_file," --output-format=dssp > ",tdir,"/dssp.txt")
  system(command)
  dssp = parse.dssp( paste0(tdir,"/dssp.txt") )
  names(dssp)[2:6] = paste0("dssp_",names(dssp)[2:6])

  df = merge(df, dssp, by.x = "loc", by.y = "resnum")


}







# ------------------------------------------------------------ output
write.csv(df, paste0(tdir, "/df.csv"),row.names = F)



















## #--------------------   AlphaFold2 pLDDT ~ disorder
#residue_disorder = data.frame(loc = 1:length(inseq_vector), pdb_plddt_order_dssp = 0)
#
## handle that alphafold returns not nicely formatted pdb with lost space between cols
#x = readLines(pdb_file)
#which_to_look_at = grep(pattern = "[A-Z]{1}[0-9]{3,5}" , x)
#x = x[which_to_look_at]
#for(line in which_to_look_at){
#  problem_str = stringr::str_extract(x[line], "[A-Z]{1}[0-9]{3,5}")
#  rest = stringr::str_split(x[line], "[A-Z]{1}[0-9]{3,5}",simplify = T)
#  to_insert = stringr::str_split(problem_str, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",simplify = T)
#  newline = paste(rest[1] , to_insert[1] , to_insert[2], rest[2])
#  x[line] = newline
#}
#x = x[1:(length(x) -2)]
#x = read.delim(text = x ,sep = "",header = F)
#for(i in 1:max(x[,6])){
#  # i is residue
#  residue_disorder$loc[i] = i
#  residue_disorder$pdb_plddt_order_dssp[i] = max(x[x$V6 == i,]$V11)
#}
#dat = merge(dat, residue_disorder, by = "loc")







# #--------------------  suspect
# sus = read.table(paste0("data/", virus,"_", gene,"_suspect.txt"),sep = "\t", stringsAsFactors = F, header = T)
# sus = melt(sus, id.vars = "pos")
# 
# # remove not covered mutations - as we used a structure  -- can ignore this as alphafold now
# dat = dat[dat$loc >= min(sus$pos),]
# dat = dat[dat$loc <= max(sus$pos),]
# 
# # add col for from_Score
# del = c()
# for(r in 1:nrow(dat)){
#   pos = dat$loc[r]
#   if(!pos %in% unique(sus$pos)){;next}
#   wtAA = dat$from[r]
#   mtAA = dat$to[r]
#   sus2 = sus[sus$pos == pos,]
#   sus2$value = sus2$value / 10
#   sus2$value = as.numeric(substr(sus2$value,1,1)) + 1 # bin values
#   dat$sus_wt[r] = sus2[sus2$variable == wtAA,3]
#   dat$sus_mt[r] = sus2[sus2$variable == mtAA,3]
#   dat$sus_ave[r] = mean(sus2$value) # i.e. how sensetive is pos
#   dat$sus_diff[r] = max(sus2$value) - min(sus2$value) # ie. whats the relative change of this?
# }




#--------------------  netsurfp
# conda = "py3"
# "netsurfp2 --csv netsurf.csv hhblits /home/oscar/lib/netsurfp2/models/hhsuite.pb temp.fasta netsurf/"

#net = read.csv(paste0("data/",virus, "_", gene, "_netsurfp.csv"))
#names(net)[4:21] = paste0("net_",names(net)[4:21]) # we will use dssp to generate features that are identical also
#dat = merge(dat, net, by.x = "loc", by.y = "n", all.x = T)








# 
# #--------------------  DDG by way of foldx - SLOW
# # there is something odd uccuring, the ddg values seem to follow an upward trend over time. potentially need to remove intermediate files within loop.
# dat$foldx_ddg = 0
# for(r in 1:nrow(dat)){
#   foldx_mutant = paste0(dat$from[r], #WT
#          "A", # CHAIN - always A for alphafold
#          dat$loc[r], #pos
#          dat$to[r] # to
#          )
# 
#   command = paste0("wsl cp ./data/", virus, "_", gene,".pdb ./temp/temp.pdb ;
#                     cp ./lib/rotabase.txt ./temp/
#                     cp -R ./lib/molecules ./temp/
#                     cd ./temp ;
#                    ~/tools/foldx_20221231 --pdb temp.pdb -c PositionScan --positions ", foldx_mutant)
#   system(command)
# 
#   # parse outputs
#   foldx_ddg = read.table("./temp/PS_temp_scanning_output.txt")[2,2]
#   dat$foldx_ddg[r] = foldx_ddg
# }
  
  
  
  
  
  
  
  
  