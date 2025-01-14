#!/bin/bash

Func_Usage()
{
cat << __EOF
========================================================================================
Usage: `basename $0` 
    1. <Script.R>	This Rscipt will be used. Mostly sample place as the shell script
    2. <CountsFile>     Counts file with first col as features
    3. <SAMPLE_FILE>    Sample information file. Should be a csv like below:
                        #SampleName,GroupName,ReadConfig,Strand,PI
                        Mut_123,MUT,PE,fr-stranded,PIKA
                        Mut_265,MUT,PE,fr-stranded,PIKA
    4. <OUT_DIR>        /FullPath/ to OutDir
    5. <EDGER_BCV>      edgeR BCV values can be:
                        edgeR_BCV == 2    for Biological (any species)
                        edgeR_BCV == 0.01 for Technical (any species)
                        edgeR_BCV == 0.1  for None (inbreeding organism i.e. mouse, fly etc..)
                        edgeR_BCV == 0.4  for None (outbreeding organism i.e. human)
    6. <CONTRAST_FILE>  /FullPath/ to contrast file. Should be a csv like below:
                        condition_A,condition_B
                        WT,MUT
    7. <TPM_file>       /FullPath/ to TPM file.
    8. <AnnotationFile> File with annotations. Relevant column should be named as FEATURE_NAME.
                        Can be none.
========================================================================================
__EOF
}

script=$1
COUNTS_FILE=$2
SAMPLE_FILE=$3
OUT_DIR=$4
EDGER_BCV=$5
CONTRAST_FILE=$6
TPM_file=$7
AnnotationFile=$8

######################################################################################
NofExpectedArguments=8
if [ $# -ne $NofExpectedArguments ]
then
	Func_Usage
	echo -e "Number of arguments specifed:\t $# \n"
	echo -e "Number of arguments required:\t $NofExpectedArguments \n"
	echo -e "\nArguments provided:\n $*"
	exit 1
fi
for((i=1;i<=$#;i++)); do echo -e "$i.\t${!i}";done;echo
Start=`date +'%x %R'`
start=`date +%s`
######################################################################################

ModulesToUnload=""
ModulesToLoad="  gcc/8.2.0 R/4.3.2"
for i in `echo $ModulesToLoad`; do module unload `dirname $i` ;done
module unload $ModulesToUnload
module load $ModulesToLoad \
|| { (>&2 echo -e "\nERROR: $ModulesToLoad modules could not be loaded\n"); exit 1 ;}
echo -e "Modules loaded:\n$ModulesToLoad\n"

if [ "$1" == "-h" ]
then
	Func_Usage
	exit
fi

######################################################################################
ulimit -n 10000

mkdir -p $OUT_DIR

Rscript $script \
	$COUNTS_FILE \
	$SAMPLE_FILE \
	$OUT_DIR \
	$EDGER_BCV \
	$CONTRAST_FILE \
	$TPM_file \
	$AnnotationFile

######################################################################################
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )

echo -e "\n\nStarted $0 at:\n$Start\nEnded at:\n`date +'%x %R'`"
echo -e "\nTotal runtime was (in secs) :\n$runtime"
