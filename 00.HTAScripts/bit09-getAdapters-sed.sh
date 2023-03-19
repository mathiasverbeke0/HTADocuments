#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: ./bit09-getAdapters.sh /home/user/fastqc/
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this script requires 1 parameters:\n
		1. Path of the folder with the FastQC output files.\n\n
        Does a sed on fastqc_data.txt files.\n";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 1 ]; then
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
for i in $( find ${inputFolder} -name 'fastqc_data.txt'); do
	echo "### $i ###";
    # Get lines from 'Overrepresented' to 'END_MODULE' with sed 
	sedCommand="sed -n -e '/Overrepresented/,/END_MODULE/ p' $i";
    # Remove first line with 'Overrepresented'
    sedCommand="$sedCommand | sed '1,1d'";
    # Remove last line with 'END_MODULE'
    sedCommand="$sedCommand | sed '\$d'";
	#echo -e "$sedCommand\n";
	output=$(eval $sedCommand);
	echo -e "$output\n";
done
################################################################################