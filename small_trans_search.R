## default values
case_in_control = 0
vals_threshold = 6

## read command line variables
args = commandArgs(trailingOnly=TRUE)

sample_name = args[1]
path_to_sample = args[2]
DIR = args[3]
outdir_path = args[4]

write.table("chr1 chr2 pos1 pos2 str_sum_intra str_sum_inter trans_search_value", paste(outdir_path, "/", sample_name,"_small_translocations, sep = ""), row.names =F, col.names = F, append = T)

## cis_trans norm coef

cis_trans_norm <- function(path_to_sample, DIR){
  d1 = matrix(0, nrow = 22, ncol = 22)
  d2 = matrix(0, nrow = 22, ncol = 22)
  d3 = matrix(0, nrow = 22, ncol = 22)
  d4 = matrix(0, nrow = 22, ncol = 22)
  d5 = matrix(0, nrow = 22, ncol = 22)
  d6 = matrix(0, nrow = 22, ncol = 22)
  d7 = matrix(0, nrow = 22, ncol = 22)
  d8 = matrix(0, nrow = 22, ncol = 22)
  d9 = matrix(0, nrow = 22, ncol = 22)
  d10 = matrix(0, nrow = 22, ncol = 22)
  for (chr_intra in 1:22) {
    #print(chr_intra)
    case <- read.table(
      paste(path_to_sample, "/chr", as.character(chr_intra), "_track", sep="")
      , stringsAsFactors = F, head=F)
    control <- read.table(
      paste(DIR, "/control_samples/chr",as.character(chr_intra),"_track_control", sep="")
      , stringsAsFactors = F, head=F)
    
    #inter_case = case[, -c(chr_inter + 3)]
    #inter_case = rowSums(inter_case[,4:24])
    
    inter_control = control[, -c(chr_intra + 1)]
    inter_control = rowSums(inter_control[,2:22])
    #cis_trans_coef
    
    a1 = case[, 4:25]
    a1 = a1[, -c(chr_intra)]
    a2 = case[, 3 + chr_intra]
    b1 = control[, 2:23]
    b1 = b1[, -c(chr_intra)]
    b2 = control[, 1 + chr_intra]
    
    for (chr_inter in 1:22) {
      if(chr_inter!=chr_intra){
        
        d1[chr_intra, chr_inter] = sum(case[, 3 + chr_inter])
        d2[chr_intra, chr_inter] = sum(case[, 3 + chr_intra])
        d3[chr_intra, chr_inter] = sum(control[, 1 + chr_inter])
        d4[chr_intra, chr_inter] = sum(control[, 1 + chr_intra])
        d5[chr_intra, chr_inter] = sum(case[, 4:25])
        d6[chr_intra, chr_inter] = sum(control[, 2:23])
        d7[chr_intra, chr_inter] = sum(a1)
        d8[chr_intra, chr_inter] = sum(a2)
        d9[chr_intra, chr_inter] = sum(b1)
        d10[chr_intra, chr_inter] = sum(b2)
        
      }
    }
  }
  d = (d1/d2)/(d3/d4)
  d1 = (d7/d8)/(d9/d10)
  d1[which(is.na(d1))] = 0
  d[which(d==Inf | d==-Inf | is.na(d))] = 0
  md = colSums(d)/21
  md1 = rowSums(d1)/21
  md11 = md1/mean(md1)
  ctm = cbind(md11)%*%rbind(md)
  return(1/ctm)
}



cis_trans_norm_coef = cis_trans_norm(path_to_sample, DIR)

for (chr_intra in 1:22) {
  ## read data
  case <- read.table(
    paste(path_to_sample, "/chr", as.character(chr_intra), "_track", sep="")
    , stringsAsFactors = F, head=F)
  control <- read.table(
    paste(DIR, "/control_samples/chr",as.character(chr_intra),"_track_control", sep="")
    , stringsAsFactors = F, head=F)
  
  ## control correction
  if(case_in_control == 1){
    control[, 1:23] = control[, 1:23] - case[, 3:25]
  }
  
  for (chr_inter in 1:22) {
    if(chr_inter!=chr_intra){
      control_intra = control[, 1 + chr_intra]
      control_inter = control[, 1 + chr_inter]
      case_intra = case[, 3 + chr_intra]
      case_inter = case[, 3 + chr_inter]
      
      ## down and up threshold for noize removing
      # down threshold
      zero_fraction = length(which(control_intra==0))/length(control_intra)
      th_down = as.numeric(quantile(control_intra, 0.4 + zero_fraction))/2
      # sequencing deph normalization
      s = case_intra/control_intra
      s[which(is.na(s) | s==Inf)] = 0
      th_down = th_down*(sum(case_intra)/sum(control_intra))
      # up threshold
      th_down = 15
      th_up = as.numeric(quantile(control_intra, 0.95))
      
      ## removing noize
      w = which(case_inter<th_down | case_inter>th_up )
      case_inter[w] = 0
      
      ## vals
      vals = (case_inter/case_intra)/(control_inter/control_intra)
      # cis_trans normalization
      vals = vals*cis_trans_norm_coef[chr_intra, chr_inter]
      # NA and Inf removing
      vals[which((vals == Inf) | is.na(vals))] = 0
      which_under_th = which(vals>vals_threshold)
      if(length(which_under_th)!=0){
        vals_res = cbind(chr_intra, chr_inter ,case[, 1:2], case[, 3 + chr_intra], case[, 3 + chr_inter], vals)
        write.table(vals_res[which_under_th,], paste(outdir_path, "/", sample_name,"_small_translocations, sep = ""), row.names =F, col.names = F, append = T)
      }
    }
  }
}

