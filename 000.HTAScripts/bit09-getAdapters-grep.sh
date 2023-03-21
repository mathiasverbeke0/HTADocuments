#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: ./bit09-getAdapters.sh /home/user/fastqc/ 4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this script requires 2 parameters:\n
		1. Path of the folder with the FastQC output files.\n
		2. Maximum number of overrepresented sequences to find.\n\n
        Does a grep -A [max_num] on fastqc_data.txt files.\n";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 2 ]; then
	usage
fi
################################################################################
# INPUT FOLDER CONTAINING .fastq.gz files
################################################################################
inputFolder=$1;
# Remove trailing slash if this is last char
len=${#inputFolder};
lastPos=$(expr $len - 1);
lastChar=${inputFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	inputFolder=${inputFolder:0:$lastPos};
fi
################################################################################
# RUN AND GET OVERREPRESENTED SEQUENCES
################################################################################
numberLines=$(expr $2 + 1);
for i in $( find ${inputFolder} -name 'fastqc_data.txt'); do
	echo "### $i ###";
	grepCommand="grep -A ${numberLines} 'Overrepresented' $i";
    grepCommand="$grepCommand | grep -B ${numberLines} 'END'";
	# echo -e "$myCommand\n";
	output=$(eval $grepCommand);
	echo -e "$output\n";
done
################################################################################