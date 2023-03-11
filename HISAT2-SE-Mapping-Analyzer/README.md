# HISAT2 Single-End Read Alignment Summary Analyzer
The HISAT2 Single-End Read Alignment Summary Analyzer is a Python command-line application that takes a HISAT2 single-end read alignment summary file as input and extracts the ...
* percentage of unpaired reads
* percentage of reads that were aligned 0 times
* percentage of reads that were aligned exactly 1 time
* percentage of reads that were aligned more than 1 time
* overall alignment rate
The application then outputs these percentages to a CSV file and, optionally, an Excel file.

## Requirements
This application requires the following packages to be installed:
* Python (version 3.6 or later)
* pandas (version 1.3.3 or later)
* argparse (version 1.4.0 or later)

You can install pandas and argparse using pip, by running the following command:
```bash
pip install pandas argparse
```

Note that this application has been tested on Linux, and should work on other Unix-based systems. It has not been tested on Windows.

## Usage
```bash
usage: python HISAT2-SE-Mapping-Analyzer.py [-h] [-o OUTDIR] [-e] filename

positional arguments:
  filename              Provide the name of the HISAT2 single-end read alignment summary file.

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Provide the name of the output directory.
  -e, --excel           Generate Excel files.
```

## Output
The application generates a CSV file with the following columns:
* SampleNo
* % unpaired
* % aligned 0 times
* % aligned exactly 1 time
* % aligned >1 times
* % overall alignment rate

In addition to the CSV file, the application will also generate an Excel file with the same columns if the -e or --excel option is used.

## Examples
```bash
$ python HISAT2-SE-Mapping-Analyzer.py sample_summary.txt -o results/
```
In this example, the HISAT2-SE-Mapping-Analyzer.py script is run on the file sample_summary.txt, and the output is saved to a directory called results (which must exist prior to running the script). The script will generate a CSV file with the results in the results directory.

```bash
$ python HISAT2-SE-Mapping-Analyzer.py -e -o results/ sample_summary.txt
```
In this example, the user wants to analyze the HISAT2 single-end read alignment summary file sample_summary.txt. The results will be saved to an output directory called results/. The -e flag is included to generate an Excel file along with the CSV file.

## Future improvements
* Add support for paired-end read alignment summary files
* Implement support for additional output file formats, such as JSON and YAML

## Contributing
Contributions are always welcome! If you find any issues or have suggestions for how to improve this tool, please feel free to open an issue or submit a pull request. To contribute, please follow these steps:

1. Fork the repository
2. Create a new branch for your changes
3. Make your changes and test them
4. Submit a pull request