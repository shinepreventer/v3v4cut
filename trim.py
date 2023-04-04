# 这是一个用于处理测序数据的脚本
# 作者：wzy
# 日期：2023-04-04
import os
import threading
import queue

input_dir = 'H:/v4test'
output_dir = 'H:/v4test/v4new'
log_file = 'H:/v4test/log.txt'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

v4_start_codes = [
    'GTGACAGCAGCCGCGGTAA',
    'GTGACAGCCGCCGCGGTAA',
    'GTGACAGCGGCCGCGGTAA',
    'GTGACAGCTGCCGCGGTAA',
    'GTGCCAGCAGCCGCGGTAA',
    'GTGCCAGCCGCCGCGGTAA',
    'GTGCCAGCGGCCGCGGTAA',
    'GTGCCAGCTGCCGCGGTAA'
]

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# Create a queue to store the files to be processed
file_queue = queue.Queue()

# Create a lock object to synchronize the access to the log file
log_lock = threading.Lock()

# Define a function to process each file
def process_file():
    # Get a file from the queue
    file = file_queue.get()
    # Open input and output files
    with open(file.path, 'r') as input_file, open(os.path.join(output_dir, file.name), 'w') as output_file:
        # Initialize counters for removed and kept bases
        removed_bp = 0
        kept_bp = 0
        # Loop through each line in input file
        for line in input_file:
            # If line starts with @ or +, it contains header or quality score information
            if line.startswith('@') or line.startswith('+'):
                # Write line to output file without any modification
                output_file.write(line)
            # Otherwise, line contains sequence data
            else:
                # Remove any whitespace characters from line
                seq = line.strip()
                # If file name ends with .2.fastq, it is a reverse sequence
                if file.name.endswith('.2.fastq'):
                    # Reverse and complement the sequence using dictionary
                    rev_seq = ''.join(complement.get(base, base) for base in reversed(seq))
                    # Find the index of the first occurrence of any pattern in v4_start_codes
                    v4_start = next((rev_seq.find(code) for code in v4_start_codes if rev_seq.find(code) != -1), -1)
                    # If index is not -1, pattern is found
                    if v4_start != -1:
                        # Calculate how many bases are removed from end of sequence
                        removed_bp += len(seq) - (v4_start + len(v4_start_codes[0]))
                        # Calculate how many bases are kept from start of sequence
                        kept_bp += v4_start + len(v4_start_codes[0])
                        # Write only the kept part of sequence to output file, reversing and complementing it again
                        output_file.write(''.join(complement.get(base, base) for base in reversed(rev_seq[v4_start:v4_start + len(v4_start_codes[0])])) + '\n')
                    # If index is -1, pattern is not found
                    else:
                        # All bases are removed from sequence
                        removed_bp += len(seq)
                # If file name does not end with .2.fastq, it is a forward sequence
                else:
                    # Find the index of the first occurrence of any pattern in v4_start_codes
                    v4_start = next((seq.find(code) for code in v4_start_codes if seq.find(code) != -1), -1)
                    # If index is not -1, pattern is found
                    if v4_start != -1:
                        # Calculate how many bases are removed from start of sequence
                        removed_bp += v4_start
                        # Calculate how many bases are kept from end of sequence
                        kept_bp += len(seq) - v4_start
                        # Write only the kept part of sequence to output file
                        output_file.write(seq[v4_start:] + '\n')
                    # If index is -1, pattern is not found
                    else:
                        # All bases are removed from sequence
                        removed_bp += len(seq)
        # Acquire the lock before writing to log file
        log_lock.acquire()
 # Open the log file and write the file name and the counters
        with open(log_file, 'a') as log:
            log.write(f'{file.name}\t{removed_bp}\t{kept_bp}\n')
        # Release the lock after writing to log file
        log_lock.release()
    # Indicate that the task is done
    file_queue.task_done()

# Create a list to store the threads
threads = []

# Define the number of threads to use
num_threads = 8

# Loop through the number of threads
for i in range(num_threads):
    # Create a thread that runs the process_file function
    thread = threading.Thread(target=process_file)
    # Add the thread to the list
    threads.append(thread)
    # Start the thread
    thread.start()

# Loop through all files in input directory
for file in os.scandir(input_dir):
    # Only process files that end with .fastq
    if file.name.endswith('.fastq'):
        # Put the file in the queue
        file_queue.put(file)

# Wait for the queue to be empty
file_queue.join()

# Print a message to indicate the code is done
print('Code finished successfully.')
