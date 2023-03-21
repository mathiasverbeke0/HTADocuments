#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: bash bit09-htseqcount.sh /home/pacoh/hisat2_filtered/ 
#             /data/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf
#             exon gene_id /home/pacoh/htseqcount/ 4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this HTSeqCount script requires 6 parameters:\n
		1. Path of the folder with the mapping files to run HTSeqCount on.\n
        2. Path of the GTF file.\n
		3. Path of the output folder.\n
		4. Feature type to use from GTF file (e.g. CDS, exon).\n
		5. Attribute to use as label for counting (e.g. gene_id, gene_name, transcript_id).\n
		6. Number of threads to use.\n";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 6 ]; then
	usage
fi
################################################################################
# INPUT FOLDER
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
# RUN HTSEQCOUNT
################################################################################
# GTF annotation file
pathGTF=$2;
feature=$4;
attribute=$5;
threads=$6;
# Part of input/file name
tmpPart='"{}"';
tmpPart2='${IN}';
tmpPart3='$(basename ${IN%.*})';
# Compose command
htseqCommand="find ${inputFolder} -name '*_sorted.bam'";
htseqCommand="$htseqCommand | xargs --max-procs=${threads} -I {} sh -c 'IN=$tmpPart;";
htseqCommand="$htseqCommand htseq-count -f bam -m intersection-strict -s no ";
htseqCommand="$htseqCommand -a 10 -t ${feature} -i ${attribute} ${tmpPart2} ${pathGTF}";
htseqCommand="$htseqCommand > ${outFolder}/counts_${tmpPart3}.txt';";
# Show command
echo -e "$htseqCommand";
# Execute
output=$(eval $htseqCommand);
# Show output
echo -e "$output\n";
################################################################################