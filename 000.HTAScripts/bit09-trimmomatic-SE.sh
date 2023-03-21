#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: bash bit09-trimmomatic-SE.sh /home/user/sra/ /home/user/trimmomatic/ 
#             /opt/Trimmomatic-0.39/adapters/adapters.fa 4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this Trimmomatic script requires 4 parameters:\n
		1. Path of the folder with fastq.gz input files.\n
		2. Path of the output folder.\n
        3. The adapter file (including full path)
		4. Number of threads to use.\n\n
        Run on Single End (SE) read files (*.fastq.gz).\n
		Uses Sanger fastq -phred33 scores.";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 4 ]; then
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
# CREATE FOLDER (IF NOT EXISTS)
################################################################################
outFolder=$2;
# Remove trailing slash if this is last char
len=${#outFolder};
lastPos=$(expr $len - 1);
lastChar=${outFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	outFolder=${outFolder:0:$lastPos};
fi
mkdir -p ${outFolder}
################################################################################
# RUN TRIMMOMATIC
################################################################################
# Path to Trimmomatic
pathTrimmomatic="/opt/Trimmomatic-0.39/trimmomatic-0.39.jar";
# Part of input/file name found with "find $1 -name '*.fastq.gz'"
tmpPart='"{}"';
tmpPart2='${IN}';
tmpPart3='$(basename ${IN})';
# Compose command: find files, pipe to threads and use IN for Trommomatic
trimCommand="find $inputFolder -name '*.fastq.gz' | xargs --max-procs=$4";
trimCommand="$trimCommand -I {} sh -c 'IN=$tmpPart;";
# Add java call, first part of Trimmomatic command
trimCommand="$trimCommand java -jar ${pathTrimmomatic} SE -threads 2 -phred33";
# Add log file
trimCommand="$trimCommand -trimlog ${outFolder}/${tmpPart3}_log.txt";
# Add sample to process
trimCommand="$trimCommand ${tmpPart2} ${outFolder}/${tmpPart3}";
# Add other Trimmomatic parameters
trimCommand="$trimCommand ILLUMINACLIP:${3}:2:30:10 LEADING:3 TRAILING:3";
trimCommand="$trimCommand SLIDINGWINDOW:4:15 MINLEN:35';";
# Show command
echo -e "$trimCommand";
# Execute
output=$(eval $trimCommand);
# Show output
echo -e "$output\n";
################################################################################
