# Import modules
import argparse
import logging
import gzip
import multiprocessing
import cutadapt

# Create a logger object with a name
logger = logging.getLogger("cut")
# Set the logging level to INFO
logger.setLevel(logging.INFO)
# Create a file handler object to write logs to log.txt
file_handler = logging.FileHandler("log.txt")
# Set the format of the log messages
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# Add the formatter to the file handler
file_handler.setFormatter(formatter)
# Add the file handler to the logger
logger.addHandler(file_handler)

# A simple quality filtering function based on Bio.SeqIO module
# Input: a SeqRecord object, a minimum average quality score, a minimum sequence length
# Output: a filtered SeqRecord object or None
def quality_filter(record, min_qual, min_len):
    # Calculate the average quality score of the record
    avg_qual = sum(record.letter_annotations["phred_quality"]) / len(record)
    # Check if the average quality score and the sequence length meet the criteria
    if avg_qual >= min_qual and len(record) >= min_len:
        # Return the record if passed
        return record
    else:
        # Return None if failed
        return None

# A function that takes a file name and a primer sequence as arguments and performs the primer trimming and quality filtering for single-end data
def process_file_se(file_name, primer):
    # Get the input file name from the file name
    input_file = file_name + ".fastq"
    # Get the output file name from the file name
    output_file = file_name + "_trimmed.fastq.gz"
    # Create a cutadapt adapter object with the primer sequence and the --anywhere option
    adapter = cutadapt.adapters.Adapter(primer, where="anywhere")
    # Create a cutadapt statistics object to record the trimming results
    stats = cutadapt.statistics.Statistics()
    # Try to read the input file and write the output file with gzip compression
    try:
        # Open the input file with SeqIO.parse
        with SeqIO.parse(input_file, "fastq") as input_handle:
            # Open the output file with gzip.open and SeqIO.write
            with gzip.open(output_file, "wt") as output_handle:
                # Loop through each record in the input file
                for record in input_handle:
                    # Trim the primer from the record using the adapter object
                    trimmed_record = adapter.match_and_trim(record)
                    # Check if the record was trimmed
                    if trimmed_record:
                        # Filter the record by quality and length
                        filtered_record = quality_filter(trimmed_record, min_qual, min_len)
                        # Check if the record passed the filtering
                        if filtered_record:
                            # Write the filtered record to the output file
                            SeqIO.write(filtered_record, output_handle, "fastq")
                            # Update the statistics object with the trimming result
                            stats.update(filtered_record)
        # Log a message indicating the completion of the trimming and filtering for this file
        logger.info(f"Finished trimming and filtering {input_file} and writing to {output_file}")
        # Log a message indicating the statistics of the trimming and filtering for this file
        logger.info(f"Trimming and filtering statistics: {stats}")
    # Handle possible exceptions
    except FileNotFoundError:
        # Log an error message if the input file is not found
        logger.error(f"Input file {input_file} not found")
    except ValueError:
        # Log an error message if the input file is not in fastq format
        logger.error(f"Input file {input_file} not in fastq format")
    except Exception as e:
        # Log an error message for any other exception
        logger.error(f"An unexpected error occurred: {e}")

# A function that takes a pair of file names and a pair of primer sequences as arguments and performs the primer trimming and quality filtering for paired-end data
def process_file_pe(file_name_1, file_name_2, primer_1, primer_2):
    # Get the input file names from the file names
    input_file_1 = file_name_1 + ".fastq"
    input_file_2 = file_name_2 + ".fastq"
    # Get the output file names from the file names
    output_file_1 = file_name_1 + "_trimmed.fastq.gz"
    output_file_2 = file_name_2 + "_trimmed.fastq.gz"
    # Create a cutadapt adapter object with the primer sequence and the --anywhere option for each end
    adapter_1 = cutadapt.adapters.Adapter(primer_1, where="anywhere")
    adapter_2 = cutadapt.adapters.Adapter(primer_2, where="anywhere")
    # Create a cutadapt statistics object to record the trimming results for each end
    stats_1 = cutadapt.statistics.Statistics()
    stats_2 = cutadapt.statistics.Statistics()
    # Try to read the input files and write the output files with gzip compression
    try:
        # Open the input files with SeqIO.parse
        with SeqIO.parse(input_file_1, "fastq") as input_handle_1, SeqIO.parse(input_file_2, "fastq") as input_handle_2:
            # Open the output files with gzip.open and SeqIO.write
            with gzip.open(output_file_1, "wt") as output_handle_1, gzip.open(output_file_2, "wt") as output_handle_2:
                # Loop through each pair of records in the input files
                for record_1, record_2 in zip(input_handle_1, input_handle_2):
                    # Trim the primer from each record using the adapter object
                    trimmed_record_1 = adapter_1.match_and_trim(record_1)
                    trimmed_record_2 = adapter_2.match_and_trim(record_2)
                    # Check if both records were trimmed
                    if trimmed_record_1 and trimmed_record_2:
                        # Filter each record by quality and length
                        filtered_record_1 = quality_filter(trimmed_record_1, min_qual, min_len)
                        filtered_record_2 = quality_filter(trimmed_record_2, min_qual, min_len)
                        # Check if both records passed the filtering
                        if filtered_record_1 and filtered_record_2:
                            # Write the filtered records to the output files
                            SeqIO.write(filtered_record_1, output_handle_1, "fastq")
                            SeqIO.write(filtered_record_2, output_handle_2, "fastq")
                            # Update the statistics objects with the trimming results
                            stats_1.update(filtered_record_1)
                            stats_2.update(filtered_record_2)
        # Log a message indicating the completion of the trimming and filtering for these files
        logger.info(f"Finished trimming and filtering {input_file_1} and {input_file_2} and writing to {output_file_1} and {output_file_2}")
        # Log a message indicating the statistics of the trimming and filtering for these files
        logger.info(f"Trimming and filtering statistics: {stats}")
    # Handle possible exceptions
    except FileNotFoundError:
        # Log an error message if any of the input files is not found
        logger.error(f"Input file {input_file} not found")
    except ValueError:
        # Log an error message if any of the input files is not in fastq format
        logger.error(f"Input file {input_file} not in fastq format")
    except Exception as e:
        # Log an error message for any other exception
        logger.error(f"An unexpected error occurred: {e}")
# Create an argument parser object to get command line arguments
parser = argparse.ArgumentParser(description="A script to trim 16sRDNA raw data (.fastq format) v4 region using cutadapt")
# Add arguments for input file, output file, primer sequence, minimum average quality score and minimum sequence length
parser.add_argument("-i", "--input", type=str, required=True, help="Input fastq file or prefix of paired-end fastq files")
parser.add_argument("-o", "--output", type=str, required=True, help="Output fastq file or prefix of paired-end fastq files")
parser.add_argument("-p", "--primer", type=str, required=True, help="Primer sequence or pair of primer sequences separated by comma")
parser.add_argument("-q", "--qual", type=float, default=20.0, help="Minimum average quality score for filtering (default: 20.0)")
parser.add_argument("-l", "--len", type=int, default=100, help="Minimum sequence length for filtering (default: 100)")
# Parse the arguments and get their values
args = parser.parse_args()
input_file = args.input
output_file = args.output
primer = args.primer
min_qual = args.qual
min_len = args.len
# Check if the input file contains a dot
if "." in input_file:
    # Assume single-end data
    # Check if the output file contains a dot
    if "." in output_file:
        # Use the output file name as is
        output_file_name = output_file
    else:
        # Add .fastq.gz to the output file name
        output_file_name = output_file + ".fastq.gz"
    # Process the input file using the primer sequence
    process_file_se(input_file, primer)
else:
    # Assume paired-end data
    # Check if the primer contains a comma
    if "," in primer:
        # Split the primer into two sequences
        primer_1, primer_2 = primer.split(",")
    else:
        # Use the same primer sequence for both ends
        primer_1 = primer_2 = primer
    # Check if the output file contains a dot
    if "." in output_file:
        # Use the output file name as prefix
        output_file_prefix = output_file
    else:
        # Use the output file name as is
        output_file_prefix = output_file
    # Add _R1 and _R2 to the input and output file names
    input_file_1 = input_file + "_R1.fastq"
    input_file_2 = input_file + "_R2.fastq"
    output_file_1 = output_file_prefix + "_R1.fastq.gz"
    output_file_2 = output_file_prefix + "_R2.fastq.gz"
    # Process the input files using the primer sequences
    process_file_pe(input_file_1, input_file_2, primer_1, primer_2)
