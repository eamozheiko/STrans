# STrans

This script useful for finding small translocations at 10kb resolution from Hi-C data. It receives as input valid_pairs_file.  valid_pairs_file format is :


```
    chromosome1 position1 chromosome2 position2
```

example:
```
1 10000 X 20000
```

## running strans

Enter paths to input files such as sample_name, valid_pairs_file_path and outdir_folder_path at the biginning of the strans.sh file.

to run program from command line do:
```
./strans.sh
```
output file will be named as:

${SAMPLE_NAME}_small_translocations

and will be in output folder

## requirements 

R, Shell


