# Generate for each change at each residue, the mutation features.
# These can then be fed into drugres classification etc.


# -------------------- Setup
### runtime vars
# //todo make pass argument
infasta = paste0("/query/HCMV_UL54.fasta")
#v_evals = c("1E-100", "1E-50", "1e-7")
v_eval = 1e-7

# passed arguments
# //todo make pass argument
blast_db_name = "uniref50.fasta"
threads = 32

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





# -------------------- princeton conservation
# ref: https://compbio.cs.princeton.edu/conservation/
command = paste0("python2 /mflibs/conservation_code/score_conservation.py -m /mflibs/conservation_code/matrix/blosum62.bla -p FALSE -g 0.99 ",temp_blast_msa," > ",tdir,"/seq_evol_conservation.txt")
system(command)
conservation = data.frame(read.table(paste0(tdir,"/seq_evol_conservation.txt"),header = F, sep = "\t")[,1:2])
colnames(conservation) = c("loc", "seq_evol_conservation")
conservation$seq_evol_conservation = as.numeric(conservation$seq_evol_conservation)
conservation[conservation$seq_evol_conservation < 0,2] = 0 # handle NA
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
  
}
df = merge(df, psi_depth, by = "loc", all.x = T)


 


#--------------------  residue-residue coevolution
# get loca-locb-coupling
command = paste0("bash /app/msa2coupling.sh -i=", temp_blast_msa ,
  " -o=/tmp/seq_evol_covar_coupling.tab")
system(command)
seq_coevol = read.table("/tmp/seq_evol_covar_coupling.tab")[,c(1,3,6)]
colnames(seq_coevol) = c("loc", "locb", "coupling")
# //todo, coupling to a p2rank loc?
# //todo make image of loca-locb matrix?
library(dplyr)
seq_coevol2 = seq_coevol[,c(1,3)] %>% 
            group_by(loc) %>%
            summarise(seq_evol_max_residue_coupling = max(coupling))
df = merge(df, seq_coevol2, by = "loc", all.x = T)
seq_coevol2 = seq_coevol[,c(1,3)] %>% 
            group_by(loc) %>%
            summarise(seq_evol_coupling_mean = mean(coupling))
df = merge(df, seq_coevol2, by = "loc", all.x = T)
seq_coevol2 = seq_coevol[,c(1,3)] %>% 
            arrange(desc(coupling)) %>%   
            group_by(loc) %>% 
            slice(1 : as.integer( max(locs) / 10))   %>%
            summarise(seq_evol_top10th_coupling_mean = mean(coupling))
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





# ------------------------------------------------------------Physiochemical Features
# wdwradius
vdw = read.csv("/mflibs/vdw_radius.csv", skip = 1)
df$seq_phys_vdw_radius_wt = 0
df$seq_phys_vdw_radius_mt = 0
df$seq_phys_vdw_radius_delta = 0
for(r in 1:nrow(df)){
  wtAA = df$wt[r]
  mtAA = df$mt[r]
  df$seq_phys_vdw_radius_wt = round(vdw[which(vdw$AA == wtAA),2],1)
  df$seq_phys_vdw_radius_mt = round(vdw[which(vdw$AA == mtAA),2],1)
  df$seq_phys_vdw_radius_delta[r] = round(abs(vdw[which(vdw$AA == wtAA),2] - vdw[which(vdw$AA == mtAA),2]),1)
}





# dydrogen moment
df$seq_phys_wt_hydro = hmoment(df$wt)
df$seq_phys_mt_hydro = hmoment(df$mt)
df$change_hydro = abs(df$seq_phys_mt_hydro - df$seq_phys_wt_hydro)
# hydrophobicity
df$seq_phys_wt_hydrophobicity = hydrophobicity(df$wt)
df$seq_phys_mt_hydrophobicity = hydrophobicity(df$mt)
df$change_hydrophobicity = abs(df$seq_phys_mt_hydrophobicity - df$seq_phys_wt_hydrophobicity)
# isoelectric point
df$seq_phys_wt_PI = pI(df$wt)
df$seq_phys_mt_PI = pI(df$mt)
df$change_PI = df$seq_phys_mt_hydro - df$seq_phys_wt_hydro
# molweight
df$seq_phys_wt_mw = mw(df$wt)
df$seq_phys_mt_mw = mw(df$mt)
df$change_mw = abs(df$seq_phys_mt_mw - df$seq_phys_wt_mw)
# charge
df$seq_phys_wt_charge = charge(df$wt)
df$seq_phys_mt_charge = charge(df$mt)
df$change_charge = abs(df$seq_phys_mt_charge - df$seq_phys_wt_charge)





# ------------------------------------------------------------ Structural features
struc = list()
# --------------------  from sequence
### disorder
command = paste0("/app/Seq2Disorder.sh -i=", infasta,
  " -o=/tmp/seq2disorder.csv")
system(command)
struc$seq2disorder = read.csv("/tmp/seq2disorder.csv")
colnames(struc$seq2disorder) = c("seq_struc_disorder")
struc$seq2disorder$loc = 1:nrow(struc$seq2disorder)
df = merge(df, struc$seq2disorder, by = "loc", all.x = T)

### secondary structure
command = paste0("/app/Seq2SecStruc.sh -i=", infasta,
  " -o=/tmp/seq2ss.csv")
system(command)
struc$seq2SecStruc = read.csv("/tmp/seq2ss.csv")[,c(1,3,4,5,6)]
colnames(struc$seq2SecStruc) = c("loc", "seq_struc_seq2ss_ss",
"seq_struc_seq2ss_C", "seq_struc_seq2ss_E", "seq_struc_seq2ss_H")
df = merge(df, struc$seq2SecStruc, by = "loc", all.x = T)





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

  df = merge(df, dssp, by.x = "loc", by.y = "resnum")
}





# --------------------  Residue clustering, how close are closest [2,5] residues
if(use_pdb){
  command = paste0("python3 /app/pdb2ResDistMatrix.py ", pdb_file, " /tmp/pdb_struc_mean_k_closest_residues.csv")
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
  ligand_interracting_locs = tdf[tdf$pocket==1,2]
  df$ligand_p2rank_best_pocket = 0
  df[ligand_interracting_locs,]$ligand_p2rank_best_pocket = 1
  # generic zscore over all pockets
  tdf2 = tdf[,c(2,5,6)]
  colnames(tdf2) = c("loc", "pdb_ligand_p2rank_zscore", "pdb_ligand_p2rank_prob")
  df = merge(df, tdf2, by = "loc", all.x = T)
}





# --------------------  Protein structure, normal mode analysis
# ref Skjaerven, L. et al. (2014) BMC Bioinformatics 15, 399. Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
if(use_pdb){}
b3d_pdb <- bio3d::read.pdb( pdb_file)
b3d_modes <- bio3d::nma(b3d_pdb)
b3d_nma_fluct = b3d_modes$fluctuations
# t = bio3d::deformation.nma(b3d_modes) # non-trivial to assign value to residue
# t = bio3d::gnm(b3d_pdb) # # non-trivial to assign value to residue
t = bio3d::torsion.pdb(b3d_pdb)
t$alpha # handle NA
t$omega

t1 = data.frame(loc = 1:nlocs,
  pdb_md_nma_fluctuations = b3d_nma_fluct,
  pdb_struc_torson_alpha = t$alpha,
  pdb_struc_torson_omega = t$omega
  )
df = merge(df, t1, by = "loc", all.x = T)










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
  
  
  
  
  
  
  
  
  
