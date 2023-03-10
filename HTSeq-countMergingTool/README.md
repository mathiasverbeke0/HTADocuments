# HTSeq-count Merging Tool
HTSeq-count Merging Tool is a command-line application that takes two HTSeq-count files as input and merges them based on a shared gene ID field. The output is a new file that contains only the gene ID and the summed read counts.

The tool first constructs two data frames from the input files and then merges them into a new data frame using the gene ID field as the key. The read counts for each gene ID in both files are summed, and a new data frame is then written to an output file.

The application is written in R and can be run from the command line. The input files should be tab-separated text files with the gene ID in the first column and the read count in the second column. The output file will be a tab-separated text file with the gene ID in the first column and the summed read count in the second column.

## Use case
This application was developed to merge two HTSeq-count files generated from RNA-seq analysis, to obtain a merged output file that can be used for downstream differential expression analysis. In RNA-seq analysis, raw reads from a sample are first trimmed, filtered, and mapped to a reference genome or transcriptome, and then the mapped reads are counted per gene using a tool such as HTSeq. However, if there are multiple fastq files per sample, HTSeq-count files will be generated separately for each fastq file. Thus, merging these HTSeq-count files is a necessary step to combine the gene counts for each sample, which can then be used for differential expression analysis.

## Requirements
This script requires R to be installed on your machine.

## Usage
To use this script, run the following command in the command line:
```bash
Rscript HTSeq-countMergingTool.R file1 file2
```
Where file1 and file2 are the paths to the HTSeq-count files you want to merge. If the file names do not end with .txt, the script will raise an error and terminate. 

## Input
This script takes two HTSeq-count files with .txt extension as input. The files should be tab-delimited text files, where the first column is the gene_id and the second column is the read count for that gene_id. There should be no header row in either file.

## Output
The output of this script is a tab-delimited text file containing the gene_id and summed read counts. The file name will be the same as the first input file, but with "_merged" appended to the end of the file name.

## Contributions
Contributions to this tool are more than welcome. If you find a bug or have a feature request, please feel free to open an issue on the GitHub repository or to submit a pull request. All contributions are greatly appreciated!