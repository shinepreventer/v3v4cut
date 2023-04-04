~~~markdown
# Sequencing Data Processing Script

This is a python script for processing sequencing data files in fastq format. It can trim the sequences based on the V4 start codes, and calculate and output the sequencing quality and base pair counts before and after trimming.

## Usage

To use this script, you need to provide three arguments: input_dir, output_dir and log_file. The input_dir is the directory that contains the fastq files to be processed. The output_dir is the directory that will store the trimmed fastq files. The log_file is the file that will record the statistics of the processing. You can also use the default values for these arguments, which are:

- input_dir: H:/v4test
- output_dir: H:/v4test/v4new
- log_file: H:/v4test/log.txt

To run the script, you can use the following command:

```bash
python trim_v4.py --input_dir H:/v4test --output_dir H:/v4test/v4new --log_file H:/v4test/log.txt
~~~

You can also omit the arguments if you want to use the default values:

```bash
python trim_v4.py
```

The script will process all the fastq files in the input directory, and write the trimmed sequences to the output directory. It will also write the statistics of each file to the log file, such as the number of removed and kept base pairs before and after trimming.



## trim.py   This file does not contain a comparison of sequencing quality information before and after cut
