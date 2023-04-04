# 16sRDNA v4 region trimming and filtering script

This is a python script that uses cutadapt and multiprocessing modules to trim and filter 16sRDNA raw data (.fastq format) v4 region.

## Requirements

- Python 3.6 or higher
- cutadapt 4.3 or higher

## Usage

```
python trim_v4.py -i INPUT -o OUTPUT -p PRIMER [-q QUAL] [-l LEN]
```

- INPUT: Input fastq file or prefix of paired-end fastq files. For example, sample.fastq for single-end data or sample for paired-end data (sample_R1.fastq and sample_R2.fastq).
- OUTPUT: Output fastq file or prefix of paired-end fastq files. For example, sample_trimmed.fastq.gz for single-end data or sample_trimmed for paired-end data (sample_trimmed_R1.fastq.gz and sample_trimmed_R2.fastq.gz).
- PRIMER: Primer sequence or pair of primer sequences separated by comma. For example, GTGCCAGCMGCCGCGGTAA for single-end data or GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT for paired-end data.
- QUAL: Minimum average quality score for filtering (default: 20.0).
- LEN: Minimum sequence length for filtering (default: 100).

## Example

```
python trim_v4.py -i sample -o sample_trimmed -p GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT -q 25 -l 150
```

This command will trim and filter the paired-end data sample_R1.fastq and sample_R2.fastq using the primer sequences GTGCCAGCMGCCGCGGTAA and GGACTACHVGGGTWTCTAAT, and write the results to sample_trimmed_R1.fastq.gz and sample_trimmed_R2.fastq.gz. It will also filter out the records that have an average quality score lower than 25 or a sequence length shorter than 150.
