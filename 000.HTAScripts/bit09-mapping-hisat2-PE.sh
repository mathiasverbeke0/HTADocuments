#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: bash bit09-mapping-hisat2-PE.sh 
#             /home/user/trimmomatic
#             /data/igenomes/hisat2-index-mm10/genome
#             /home/user/hisat2/
#             4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="This mapping script requires 4 parameters:\n
		1. Path of the folder with input files to run HISAT2 on.\n
        2. Path of the genome or transcriptome index (and prefix) of the files.\n
        3. Path of the output folder (with trailing slash).\n
		4. Number of threads to use (max. 8 on BIT server!).\n\n
        Run HISAT2 to map paired-end (PE) reads to the reference transcriptome.\n
        Files expected to be *_1P.fastq.gz and *_2P.fastq.gz.\n
		Typically you run this on trimmed reads e.g. after Trimmomatic.";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 4 ]; then
	usage
fi
################################################################################
# 1. INPUT FOLDER CONTAINING .fastq.gz files
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
# 2. LOCATION OF INDEX FILES
################################################################################
pathIndex=$2;
################################################################################
# 5. CREATE FOLDER (IF NOT EXISTS)
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
# 6. RUN TOPHAT2 USING NUMBER OF THREADS
################################################################################
# Concatenate _1 filenames
for i in $(ls ${inputFolder}/*_1P.fastq.gz); do
	inputFile1="${i}";
    # Remove _1P.fastq.gz and make filename 2
    posKeep=$(expr ${#i} - 12);
	baseNameTmp=${i:0:$posKeep};
    inputFile2="${baseNameTmp}_2P.fastq.gz";
    # Remove path to get sample name for output folder
	baseName=${baseNameTmp/"$inputFolder"/""};
    # Show inputfiles
    echo "### Inputfiles: $inputFile1 $inputFile2 ###";
    echo "### Running HISAT2 ###";
    # Compose command
    hisat2Command="hisat2 -x ${pathIndex} --min-intronlen 50 -p $4";
    hisat2Command="$hisat2Command -1 ${inputFile1} -2 ${inputFile2}";
    hisat2Command="$hisat2Command -S ${outFolder}${baseName}.sam";
    hisat2Command="$hisat2Command --un-conc ${outFolder}${baseName}_unmapped.fastq";
    # Show command
    echo -e "$hisat2Command";
    # Execute
    output=$(eval $hisat2Command);
    # Show output
    echo -e "$output";
    # Convert sam to bam
    sam2bamCommand="samtools view -bS ${outFolder}${baseName}.sam > ${outFolder}${baseName}.bam";
    echo -e "$sam2bamCommand";
    output2=$(eval $sam2bamCommand);
    echo -e "$output2";
    # Delete large .sam files (only keep .bam)
    samrmCommand="rm ${outFolder}${baseName}.sam";
    echo -e "$samrmCommand";
    output3=$(eval $samrmCommand);
    echo -e "$output3";
done
################################################################################