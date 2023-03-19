#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: bash bit09-filtering-samtools.sh /home/user/hisat2/ 20
#             /home/user/hisat2_filtered/ 4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="This filtering script requires 4 parameters:\n
		1. Path of the folder with mapping files to run samtools on.\n
        2. MAPping Quality value to filter (use default 20).\n
        3. Path of the output folder.\n
		4. Number of threads to use (max. 8 on BIT server!).\n\n
        Run samtools (on HISAT2 output).";
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
# OUTPUT FOLDER (CREATE IF NOT EXISTS)
################################################################################
outFolder=$3;
# Remove trailing slash if this is last char
len=${#outFolder};
lastPos=$(expr $len - 1);
lastChar=${outFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	outFolder=${outFolder:0:$lastPos};
fi
mkdir -p ${outFolder}
################################################################################
# RUN SAMTOOLS
################################################################################
# MAPping Quality and number of threads to use
mapq=$2;
threads=$4;
################################################################################
for i in $(ls ${inputFolder}/*.bam); do
    inputFile="${i}";
    # Remove .bam
    posKeep=$(expr ${#i} - 4);
    baseNameTmp=${i:0:$posKeep};
    # Remove path to get sample name for output folder
	baseName=${baseNameTmp/"${inputFolder}/"/""};
    # Show inputfile
    echo "### Inputfile: $inputFile (basename: $baseName) ###\n";
    # SAMTOOLS VIEW command
    viewCommand="samtools view -bq ${mapq} $inputFile";
    viewCommand="$viewCommand > ${outFolder}/${baseName}_filtered.bam;";
    echo -e "$viewCommand";
    outputViewCommand=$(eval $viewCommand);
    echo -e "$outputViewCommand";
    # SAMTOOLS SORT command
    sortCommand="samtools sort ${outFolder}/${baseName}_filtered.bam";
    sortCommand="$sortCommand -o ${outFolder}/${baseName}_filtered_sorted.bam;";
    echo -e "$sortCommand";
    outputSortCommand=$(eval $sortCommand);
    echo -e "$outputSortCommand";
    # SAMTOOLS INDEX command
    indexCommand="samtools index ${outFolder}/${baseName}_filtered_sorted.bam;";
    echo -e "$indexCommand";
    outputIndexCommand=$(eval $indexCommand);
    echo -e "$outputIndexCommand\n";
done
################################################################################