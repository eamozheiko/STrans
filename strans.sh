#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"



#####################################################################################################
### PATHS TO INPUT FILES
#####################################################################################################
## default
CAPTURE_TRACK_PATH=${DIR}/MedExome_hg19_capture_targets.bed
FS_TRACK_PATH=${DIR}/fs_track_100kb.txt
CHR_SIZES_PATH=${DIR}/chr_sizes_hg19
CONTROL_COVERAGE_TRACK_PATH=${DIR}/control_samples
BINSIZE=10000

## edit here
SAMPLE_NAME= 
VALID_PAIRS_FILE_PATH=
OUTDIR_FOLDER_PATH=



mkdir ${OUTDIR_FOLDER_PATH}
cd ${OUTDIR_FOLDER_PATH}


cat ${VALID_PAIRS_FILE_PATH} | awk '{
if(NR!=1){
  if($1!="chrX" && $3!="chrX" && $1!="chrY" && $3!="chrY"){
    if($1==$3){
      if(($4-$2)>1000000 || ($2-$4)>1000000){
        print $1 " " $2 " " $3 " " $4
      }
    } else {
      print $1 " " $2 " " $3 " " $4
    }

  }
}
}' > allValidPairs${SAMPLE_NAME}
time sed -i -r 's/chr//g' allValidPairs${SAMPLE_NAME}


echo ${BINSIZE} > binsize

#generate bins
for chr in {1..22}
do
echo ${chr} > chr
cat ${CHR_SIZES_PATH} | awk 'BEGIN{getline chr < "chr"; getline binsize < "binsize";}{
  if(NR==chr){
  max = $2/binsize
  max = max - max%1
  max = max + 1
  for (i = 1; i <= max; ++i){
    bin_st=(i-1)*binsize
    bin_en=bin_st+binsize
    print bin_st " " bin_en
    }
  }

}' > chr${chr}_bins
done

# fist step
for chr1 in {1..22}
do
echo ${chr1}
echo ${chr1} > chr1
echo ${chr1} > chr2

## bin st
cat chr${chr1}_bins | awk '{
  print $1
}' > bin_st

cat chr${chr1}_bins | awk '{
  print $2
}' > bin_en

cat bin_en | tail -n1 > bin_en_max


cat allValidPairs${SAMPLE_NAME} | awk 'BEGIN{
    k = 0
    getline chr1 < "chr1"
    getline bin_st < "bin_st"
    getline bin_en < "bin_en"
    getline bin_en_max < "bin_en_max"
    val = 1
  }{
    if($1==$3 && $1==chr1 && ($4-$2)>500){
      coord = $2
      while (coord>bin_en){
        print bin_st " " bin_en " " 0
        getline bin_st < "bin_st"
        getline bin_en < "bin_en"
        val = 0
      }
      val = val + 1
    }
  
  }END{
      while(bin_en<=bin_en_max && bin_en>0 && k!=1){
        print bin_st " " bin_en " " 0
        if(bin_en==bin_en_max){k=1}
        getline bin_st < "bin_st"
        getline bin_en < "bin_en"
        
      }
  }' > chr${chr1}_track
done




###############################################################################################
for chr1 in {1..22}
do
for chr2 in {1..22}
do
echo ${chr1}
echo ${chr2}
echo ${chr1} > chr1
echo ${chr2} > chr2


## bin st
cat chr${chr1}_bins | awk '{
  print $1
}' > bin_st

cat chr${chr1}_bins | awk '{
  print $2
}' > bin_en

cat bin_en | tail -n1 > bin_en_max


# inter tracks
  # case1


cat allValidPairs${SAMPLE_NAME} | awk 'BEGIN{
  getline chr1 < "chr1"
  getline chr2 < "chr2"
}{
  if($1==$3 && chr1==$1 && chr1==chr2){
      print $4
      print $2
  }
  if($1==chr2 && $3==chr1 && chr1!=chr2){
    print $4
  }
  if($1==chr1 && $3==chr2 && chr1!=chr2){
    print $2
  }
}' > inter_track_pre1
cat inter_track_pre1 | sort -k1,1n > inter_track_pre2
cat inter_track_pre2 | awk 'BEGIN{
  k = 0
  val = 1
  getline bin_st < "bin_st"
  getline bin_en < "bin_en"
  getline bin_en_max < "bin_en_max"
}{
  coord = $1
  while (coord>bin_en){
    print val
    getline bin_st < "bin_st"
    getline bin_en < "bin_en"
    val = 0
  }
  val = val + 1
}END{
  while(bin_en<=bin_en_max && bin_en>0 && k!=1){
    print 0
    if(bin_en==bin_en_max){k=1}
    getline bin_st < "bin_st"
    getline bin_en < "bin_en"
      
  }
}' > inter_track
rm inter_track_pre1
rm inter_track_pre2




cat chr${chr1}_track > chr${chr1}_track_pre
cat chr${chr1}_track_pre | awk '{
  getline inter_track < "inter_track"
  print $0 " " inter_track
}' > chr${chr1}_track
rm chr${chr1}_track_pre


done
done

#rm
for chr in {1..22}
do
rm chr${chr}_bins
done
rm bin_st
rm bin_en
rm chr
rm inter_track

# trans search
Rscript ${DIR}/small_trans_search.R ${SAMPLE_NAME} ${OUTDIR_FOLDER_PATH} ${DIR} ${OUTDIR_FOLDER_PATH}

exit































