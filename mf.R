# Generate for each change at each residue, the mutation features.
# These can then be fed into drugres classification etc.


# -------------------- Setup
### runtime vars
args = commandArgs(trailingOnly=TRUE)
infasta = as.character(args[1]) #"/query/HCMV_UL97.fasta"
blast_db_name = as.character(args[2]) # "uniref50.fasta"
threads = as.numeric(args[3]) # 32
v_eval = as.character(args[4]) # 1e-7 # psiblast e value


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
library(bio3d)

#source("R/functions.R")


# ------------------------------------------------------------ setup empty df
inseq = as.character(unlist(readFASTA(infasta)))
inseq_vector = as.vector(str_split_fixed(inseq, pattern = "", n = nchar(inseq)))

# generate a table of location, wild type, mutant type residue
nlocs = length(inseq_vector)
locs = 1:nlocs
amino_acids = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M","N", "P", "Q", "R", "S", "T", "V", "W", "Y")
wt_mt = expand.grid(inseq_vector, amino_acids)
df = cbind(locs, wt_mt)
colnames(df) = c("loc", "wt", "mt")

##### is there a PDB?
pdb_file = gsub(".fasta", ".pdb", infasta )
use_pdb = file.exists(pdb_file)

# outfile
out_file = gsub(".fasta", "_MF.csv", infasta )



# ------------------------------------------------------------Evolutionary Features
# -------------------- PSSM
tdir = "/tmp"
pssm_file = paste0(tdir, "/seq_evol_pssm.txt")
temp_blast = paste0(tdir,"/psiblast.fa")
temp_blast_msa = paste0(tdir,"/psiblast_msa.fa")

# search for homology
command = paste0("psiblast -num_threads ", threads ," -query ", infasta," -db db/", blast_db_name, " -num_iterations 2 -out_ascii_pssm ",
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

# deduplicate fasta
system( paste0("awk '/^>/{f=!d[$1];d[$1]=1}f' ", temp_blast, " > /tmp/temp.fa") )
system( paste0("cp /tmp/temp.fa ", temp_blast) )

# align outputted fasta
command = paste0("mafft --add  ", temp_blast," --keeplength --thread ",  threads ," ",infasta, " > ", temp_blast_msa," 2>",tdir,"/err.txt" )
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

df$seq_evol_pssm_wt = 0; df$seq_evol_pssm_mt = 0; df$seq_evol_pssm_diff = 0;
cols = colnames(pssm)
for(r in 1:nrow(df)){
  pos = df$loc[r]
  wtAA = df$wt[r]
  mtAA = df$mt[r]
  sus2 = pssm[pos,]
  wtcol = grep(wtAA, cols)
  mtcol = grep(mtAA, cols)
  df$seq_evol_pssm_wt[r] = pssm[pos,wtcol]
  df$seq_evol_pssm_mt[r] = pssm[pos,mtcol]
  df$seq_evol_pssm_diff[r] = abs(as.numeric(df$seq_evol_pssm_mt[r]) - as.numeric(df$seq_evol_pssm_wt[r]))
  df$seq_evol_pssm_mean[r] = mean(as.numeric(pssm[pos,]))
}





# get the move ave along a vector, do this before merge by loc
mav =  function(x, n = 5){
  # only is sensible for odd mav values
  if( n %% 2 == 0){stop("CONSERVATION: cannot run an even moving average")}
  y = apply(embed(x, n), 1, mean)
  # fill NA start and end
  n_to_pad = floor(n / 2) # by taking the mav we lose the first and last n_to_pad values
  pad_start = rep(y[1] , n_to_pad)
  pad_end = rep(y[length(y)] , n_to_pad)
  y = c(pad_start , y , pad_end)
  return(y)
}





# -------------------- princeton conservation
# ref: https://compbio.cs.princeton.edu/conservation/
command = paste0("python2 /mflibs/conservation_code/score_conservation.py -m /mflibs/conservation_code/matrix/blosum62.bla -p FALSE -g 0.99 ",temp_blast_msa," > ",tdir,"/seq_evol_conservation.txt")
system(command)
conservation = data.frame(read.table(paste0(tdir,"/seq_evol_conservation.txt"),header = F, sep = "\t")[,1:2])
colnames(conservation) = c("loc", "seq_evol_conservation")
conservation$seq_evol_conservation = as.numeric(conservation$seq_evol_conservation)
conservation[conservation$seq_evol_conservation < 0,2] = 0 # handle NA
conservation[,1] = conservation[,1] + 1 # reindex to start at 1 not 0
conservation$seq_evol_conservation_ma5 = mav(conservation[,2],5)
conservation$seq_evol_conservation_ma11 = mav(conservation[,2],11)
df = merge(df, conservation, by = "loc", all.x = T)





# -------------------- psi-blast conservation
# //todo depthnorm could do with a more informative value
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
a = Biostrings::readAAStringSet(temp_blast_msa)
a = as.matrix(a)
psi_depth = data.frame(loc = 1:ncol(a), seq_evol_psi_depth = 0, seq_evol_psi_unique = 0)
for(i in 1:ncol(a)){
  t = as.character(a[,i])
  t = t[t != "-"]
  psi_depth$seq_evol_psi_depth[i] = round(length(t),0)
  psi_depth$seq_evol_psi_depthnorm[i] = round(length(t) / nrow(a),2)
  psi_depth$seq_evol_psi_unique[i] = length(unique(t))
  psi_depth$seq_evol_psi_depth_ma5 = mav(psi_depth$seq_evol_psi_depth,5)
  psi_depth$seq_evol_psi_depthnorm_ma5 = mav(psi_depth$seq_evol_psi_depthnorm,5)
  psi_depth$seq_evol_psi_unique_ma5 = mav(psi_depth$seq_evol_psi_unique,5) 
}
df = merge(df, psi_depth, by = "loc", all.x = T)


 


#--------------------  residue-residue coevolution
# get loca-locb-coupling
command = paste0("bash /scripts/msa2coupling.sh -i=", temp_blast_msa ,
  " -o=/tmp/seq_evol_covar_coupling.tab")
system(command)
seq_coevol = read.table("/tmp/seq_evol_covar_coupling.tab")[,c(1,3,6)]
colnames(seq_coevol) = c("loc", "locb", "coupling")
# //todo, coupling to a p2rank loc?
# //todo make image of loca-locb matrix?
library(dplyr)
seq_coevol2 = data.frame()
for(i in 1:nlocs){
  which_locs = c( which(seq_coevol$loc == i) , which(seq_coevol$locb == i) )
  t_coevol = seq_coevol[which_locs,]
  seq_evol_max_residue_coupling = max(t_coevol$coupling)
  seq_evol_residue_coupling_mean = mean(t_coevol$coupling)
  seq_evol_top10th_coupling_mean = t_coevol %>% 
  slice(1 : as.integer( max(locs) / 10))   %>%
  summarise(seq_evol_top10th_coupling_mean = mean(coupling))

  seq_coevol2 = rbind(seq_coevol2, 
    data.frame(loc = i,
    seq_evol_max_residue_coupling,
    seq_evol_residue_coupling_mean,
    seq_evol_top10th_coupling_mean) )
}
df = merge(df, seq_coevol2, by = "loc", all.x = T)





#--------------------  Amino acid substitution scores
# grantham
grantham = readr::read_tsv("/mflibs/grantham.tsv") %>%
  tidyr:: gather(SECOND,SCORE, -FIRST) %>% dplyr::filter(SCORE > 0)  
df$seq_evol_grantham = 0
for(i in 1:nrow(grantham)){
  wt = grantham$FIRST[i]
  mt = grantham$SECOND[i]
  df[df$wt ==wt & df$mt == mt,]$seq_evol_grantham = grantham$SCORE[i]
  df[df$mt ==wt & df$wt == mt,]$seq_evol_grantham = grantham$SCORE[i]
}





#blossum
data(AABLOSUM62)
blosum = data.frame(AABLOSUM62)
blosum$wt = rownames(blosum)
blosum = reshape2::melt(blosum, id.vars = "wt")
colnames(blosum) = c("wt", "mt", "blosum")
df$seq_evol_blosum62 = 0
for(i in 1:nrow(blosum)){
  wt = blosum$wt[i]
  mt = blosum$mt[i]
  df[df$wt ==wt & df$mt == mt,]$seq_evol_blosum62 = blosum$blosum[i]
  df[df$mt ==wt & df$wt == mt,]$seq_evol_blosum62 = blosum$blosum[i]
}





# proximity to N and C terminus
df$seq_struc_proximity50_N_C = 0
df[df$loc < 50,]$seq_struc_proximity50_N_C = 1
df[df$loc > (nlocs - 50),]$seq_struc_proximity50_N_C = 1

df$seq_struc_proximity10_N_C = 0
df[df$loc < 10,]$seq_struc_proximity50_N_C = 1
df[df$loc > (nlocs - 10),]$seq_struc_proximity50_N_C = 1





# sequence diversity measures - from bio3d
t_msa = bio3d::read.fasta(temp_blast_msa)
t = bio3d::entropy(t_msa)
t1 = data.frame(loc = 1:nlocs,
  seq_evol_SHentropy = t$H,
  seq_evol_SHentropy_norm = t$H.norm,
  seq_evol_SHentropy10 = t$H.10,
  seq_evol_SHentropy10_norm = t$H.10.norm,
  seq_evol_conservation_bio3d = bio3d::conserv(t_msa, method = "similarity", sub.matrix = "bio3d"),
  seq_evol_conservation_blosum62 = bio3d::conserv(t_msa, method = "similarity", sub.matrix = "blosum62")
  #seq_evol_conservation_pam30 = bio3d::conserv(t_msa, method = "similarity", sub.matrix = "pam30") # fails
)
df = merge(df, t1, by = "loc", all.x = T)





# generate a HMMprofile for the alignment, and return the " -log(emmission_probability)" using HMMER3
# this simple appends the columns, not by mt or wt
read_hmmprofile <- function(file){
  text = readLines(file)
  
  start = grep("HMM", text)[2]
  end = grep("//", text)
  
  if(length(start) == 0 || length(end) == 0) {stop("malformed hmm profile")}
  
  # parser currently only handles /f formatted .hmm files
  if( grepl("HMMER3/f", text[1]) ){
    text = text[start:end]
    
    
    which_emmission = grep(" [0-9]{1,9} ", text)
    which_transition = which_emmission + 2
    
    if(length(start) == 0 || length(end) == 0) {stop("the parser cannot find residue emmission or transition probabilities")}
    
    df_emmission = read.table(text = text[which_emmission])[,1:21]
    df_transition = read.table(text = text[which_transition])
    
    df = cbind(df_emmission,
               df_transition)
    
    colnames(df) = c("position",
                     "A",
                     "C",
                     "D",
                     "E",
                     "F",
                     "G",
                     "H",
                     "I",
                     "K",
                     "L",
                     "M",
                     "N",
                     "P",
                     "Q",
                     "R",
                     "S",
                     "T",
                     "V",
                     "W",
                     "Y",
                     "m->m",
                     "m->i",
                     "m->d",
                     "i->m",
                     "i->i",
                     "d->m",
                     "d->d")
  }
  
  return(df)
}
system("hmmbuild -n query_hmm  --symfrac 0 /tmp/psiblast_msa.hmm /tmp/psiblast_msa.fa ")
d_hmmer = read_hmmprofile("/tmp/psiblast_msa.hmm")
colnames(d_hmmer) = paste0("seq_evol_HmmEmmProb_", colnames(d_hmmer) )
df = merge(df,d_hmmer[,1:21], by.x = "loc", by.y = "seq_evol_HmmEmmProb_position")





# residues in Pfam domain - binary boolean
command = paste0("/scripts/Seq2PfamResidues.sh -i=", infasta, " -o=/tmp/pfam_domains.out")
system(command)
# read in key line from table
text = readLines("/tmp/pfam_domains.out")
# if no hits
if(length(text) <= 13){
    pfam_from = 0
    pfam_to = 0
}else{
    # else hit found - only care about best
    text = text[c(2,4)]
    t_pfam = read.table(text = text)[,]
    pfam_from = as.numeric(t_pfam[1,20])
    pfam_to = as.numeric(t_pfam[1,21])
    rm(t_pfam)
}
df$seq_evol_pfam_domain = 0
df[df$loc %in% pfam_from:pfam_to,]$seq_evol_pfam_domain = 1




# ------------------------------------------------------------physicochemical Features
# define a maping of AA -> proterty vector, then just apply
# wt, mt, diff
vdw = read.csv("/mflibs/vdw_radius.csv", skip = 1)
# protscale key features
seq_phys_bulkiness = read.table("/mflibs/protscale/bulkiness.tsv", header = T)
seq_phys_recognition_factors = read.table("/mflibs/protscale/recognition_factor.tsv", header = T)
seq_phys_buried_residues = read.table("/mflibs/protscale/fraction_buried.tsv", header = T)
seq_evol_relative_mutability = read.table("/mflibs/protscale/relative_mutability.tsv", header = T)
seq_phys_polarity = read.table("/mflibs/protscale/polarity.tsv", header = T)
protscale = data.frame(
  # key protscale
  AA = seq_phys_bulkiness[,1],
  seq_phys_bulkiness = as.numeric(seq_phys_bulkiness[,2]),
  seq_phys_recognition_factors = as.numeric(seq_phys_recognition_factors[,2]),
  seq_phys_buried_residues = as.numeric(seq_phys_buried_residues[,2]),
  seq_evol_relative_mutability = as.numeric(seq_evol_relative_mutability[,2]),
  seq_phys_polarity = as.numeric(seq_phys_polarity[,2]),

  # R protr
  seq_phys_hydrophobicity = hydrophobicity(seq_phys_bulkiness[,1]),  # hydrophobicity
  seq_phys_hmoment = hmoment(seq_phys_bulkiness[,1]),                # dydrogen moment
  seq_phys_isolectric = pI(seq_phys_bulkiness[,1]),                  # isoelectric point
  seq_phys_molweight = mw(seq_phys_bulkiness[,1]),                   # molecular weight
  seq_phys_charge = charge(seq_phys_bulkiness[,1]),                  # charge

  # other
  seq_phys_vdw_radius = vdw[,2]
)
newcols = c( paste0(names(protscale[1,-1]), "_wt"),
              paste0(names(protscale[1,-1]), "_mt"),
              paste0(names(protscale[1,-1]), "_diff") )
n_newcols = length(newcols)
physdat = data.frame(matrix(ncol = n_newcols , nrow = nrow(df)))
colnames(physdat) = newcols
for(aa in protscale$AA){
  physdat[ which(df$wt == aa) , 1:(n_newcols / 3)] = protscale[protscale$AA == aa,-1]
  physdat[ which(df$mt == aa) , ((n_newcols / 3)+1) : (2*(n_newcols / 3)) ] = protscale[protscale$AA == aa,-1]
}
physdat[ , (2*(n_newcols / 3)) : n_newcols ] = (physdat[ , 1: (n_newcols / 3)]) - (physdat[ , ( (n_newcols / 3)+1) : (2*(n_newcols / 3)) ])
df = cbind(df,physdat)





# ------------------------------------------------------------ Structural features
if(1 == 1){
  struc = list()
  # --------------------  from sequence
  ### disorder
  command = paste0("/scripts/Seq2Disorder.sh -i=", infasta,
    " -o=/tmp/seq2disorder.csv")
  system(command)
  struc$seq2disorder = read.csv("/tmp/seq2disorder.csv", header = F)
  colnames(struc$seq2disorder) = c("seq_struc_disorder")
  struc$seq2disorder$loc = 1:nrow(struc$seq2disorder)
  df = merge(df, struc$seq2disorder, by = "loc", all.x = T)

  ### secondary structure
  command = paste0("/scripts/Seq2SecStruc.sh -i=", infasta,
    " -o=/tmp/seq2ss.csv")
  system(command)
  struc$seq2SecStruc = read.csv("/tmp/seq2ss.csv")[,c(1,3,4,5,6)]
  colnames(struc$seq2SecStruc) = c("loc", "seq_struc_seq2ss_ss",
  "seq_struc_seq2ss_C", "seq_struc_seq2ss_E", "seq_struc_seq2ss_H")
  df = merge(df, struc$seq2SecStruc, by = "loc", all.x = T)
}





# --------------------  PDB -> structural features
### DSSP
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
    d <- as.data.frame(matrix(c(resnum, respdb, chain, aa,
                                ss, sasa, phi, psi), ncol = 8),
                        stringsAsFactors = FALSE)
    
    colnames(d) <- c('loc', 'pdb_struc_dssp_respdb',
                      'pdb_struc_dssp_chain', 'pdb_struc_dssp_aa',
                      'pdb_struc_dssp_ss', 'pdb_struc_dssp_asa',
                      'pdb_struc_dssp_phi', 'pdb_struc_dssp_psi')
    
    d$loc <- as.numeric(d$loc)
    d$pdb_struc_dssp_respdb <- as.numeric(d$pdb_struc_dssp_respdb)
    d$pdb_struc_dssp_asa <- as.numeric(d$pdb_struc_dssp_asa)
    d$pdb_struc_dssp_phi <- as.numeric(d$pdb_struc_dssp_phi)
    d$pdb_struc_dssp_psi <- as.numeric(d$pdb_struc_dssp_psi)
        
    ## --------------- Remove empty lines between chains ------------- ##
    badlines <- c()
    for (i in 1:nrow(d)){
      if (d$pdb_struc_dssp_aa[i] == '!' | d$pdb_struc_dssp_aa[i] == 'X'){
        badlines <- c(badlines, i)
      }
    }
    if (length(badlines) != 0){
      d <- d[-badlines,]
      d$loc <- 1:nrow(d)
    }
    
    ## --------------- ASA -> RSA ------------- ##
    ##### rsa
    asa_lookup = c('A'=129.0, 'R'=274.0, 'N'=195.0, 'D'=193.0, 'C'=167.0,
                  'E'=223.0, 'Q'=225.0, 'G'=104.0, 'H'=224.0, 'I'=197.0,
                  'L'=201.0, 'K'=236.0, 'M'=224.0, 'F'=240.0, 'P'=159.0,
                  'S'=155.0, 'T'=172.0, 'W'=285.0, 'Y'=263.0, 'V'=174.0)
    
    # for each residue, get asa
    d$pdb_struc_dssp_rsa = 0.0
    for(i in 1:nrow(d)){
      asa = d$pdb_struc_dssp_asa[i]
      total_surface_area_for_residue = as.numeric(asa_lookup[which(names(asa_lookup) == d$pdb_struc_dssp_aa[i])])
      d$pdb_struc_dssp_rsa[i] = asa / total_surface_area_for_residue
    }
    
    return(d[,c(1,5,6,7,8,9)])
  }

  dssp_exec = "/usr/local/bin/mkdssp"
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

  df = merge(df, dssp, by.x = "loc", by.y = "loc")
}





# --------------------  Residue clustering, how close are closest [2,5] residues
if(use_pdb){
  command = paste0("python3 /scripts/pdb2ResDistMatrix.py ", pdb_file, " /tmp/pdb_struc_mean_k_closest_residues.csv")
  system(command)
  res_clust = read.csv("/tmp/pdb_struc_mean_k_closest_residues.csv")
  colnames(res_clust) = c("loc","pdb_struc_mean_2_closest_residues","pdb_struc_mean_5_closest_residues")
  df = merge(df, res_clust, by = "loc", all.x = T)
}





# --------------------  protein ligand binding site complex
### p2rank
# predict where the most likely binding pocket is.  which residues are involved?
if(use_pdb){
  command = paste0("/tools/p2rank_2.4/prank predict -f ",pdb_file," -o /tmp/p2rank") 
  system(command)
  tfile = list.files("/tmp/p2rank", "*.pdb_residues.csv", full.names = T)
  tdf = read.csv(tfile)
  #residue part of key ligand site? 0 is not, 1 is yes
  ligand_interracting_locs = tdf[tdf$pocket == 1,2]
  df$ligand_p2rank_best_pocket = 0
  df[df$loc %in% ligand_interracting_locs,]$ligand_p2rank_best_pocket = 1
  # generic zscore over all pockets
  tdf2 = tdf[,c(2,5,6)]
  colnames(tdf2) = c("loc", "pdb_ligand_p2rank_zscore", "pdb_ligand_p2rank_prob")
  df = merge(df, tdf2, by = "loc", all.x = T)
}





# -------------------- Protein structure, normal mode analysis
# ref Skjaerven, L. et al. (2014) BMC Bioinformatics 15, 399. Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
if(use_pdb){
  b3d_pdb <- bio3d::read.pdb( pdb_file)
  b3d_modes <- bio3d::nma(b3d_pdb)
  b3d_nma_fluct = b3d_modes$fluctuations
  # t = bio3d::deformation.nma(b3d_modes) # non-trivial to assign value to residue
  # t = bio3d::gnm(b3d_pdb) # # non-trivial to assign value to residue
  t = bio3d::torsion.pdb(b3d_pdb)
  #t$alpha # handle NA
  #t$omega

  t1 = data.frame(loc = 1:nlocs,
    pdb_md_nma_fluctuations = b3d_nma_fluct,
    pdb_struc_torson_alpha = t$alpha,
    pdb_struc_torson_omega = t$omega
    )
  # first and last residues do not have certain angles, as theres no neighbour
  t1$pdb_struc_torson_alpha[1] = 0
  t1$pdb_struc_torson_alpha[(nlocs - 1):nlocs] = 0
  t1$pdb_struc_torson_omega[nlocs] = 0
  df = merge(df, t1, by = "loc", all.x = T)
}




 
# -------------------- Protein sequence natural language embedding
if(1==2){
  command = paste0("bash /scripts/Seq2ProtLangRep.sh -i=",infasta, " -o=/tmp/natlang/prot5")
  system(command)

  # # append protein vector to each residue row
  # t_seq2protlangrep = as.numeric(unlist(read.csv("/tmp/natlang/prot5_protein.csv", header = F)))
  # for( c in 1:length(t_seq2protlangrep) ){
  #   df = cbind(df, rep(t_seq2protlangrep , nrow(df)) )
  # }

  # append residue vector to each residue row
  t_seq2residuelangrep = read.csv("/tmp/natlang/prot5_residue.csv", header = F)
  colnames( t_seq2residuelangrep) = paste0("seq2residuelangrep_", 1:ncol(t_seq2residuelangrep))
  t_seq2residuelangrep = cbind( loc = 1:nrow(t_seq2residuelangrep),
                            t_seq2residuelangrep)
  df = merge(df, t_seq2residuelangrep, by = "loc")
}





# ------------------------------------------------------------ output
write.csv(df, out_file, row.names = F)









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
  
  
  
  
  
  
  
  
  
