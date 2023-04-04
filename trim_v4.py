
import os

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

with open(log_file, 'w') as log:
    log.write('filename\tremoved_bp\tkept_bp\n')
    for filename in os.listdir(input_dir):
        if filename.endswith('.fastq'):
            with open(os.path.join(input_dir, filename), 'r') as input_file:
                with open(os.path.join(output_dir, filename), 'w') as output_file:
                    removed_bp = 0
                    kept_bp = 0
                    for line in input_file:
                        if line.startswith('@'):
                            output_file.write(line)
                        elif line.startswith('+'):
                            output_file.write(line)
                        else:
                            seq = line.strip()
                            if filename.endswith('.2.fastq'):
                                rev_seq = ''.join(complement.get(base, base) for base in reversed(seq))
                                v4_start = -1
                                for code in v4_start_codes:
                                    v4_start = rev_seq.find(code)
                                    if v4_start != -1:
                                        break
                                if v4_start != -1:
                                    removed_bp += len(seq) - (v4_start + len(code))
                                    kept_bp += v4_start + len(code)
                                    output_file.write(''.join(complement.get(base, base) for base in reversed(rev_seq[v4_start:])) + '\n')
                                else:
                                    removed_bp += len(seq)
                            else:
                                v4_start = -1
                                for code in v4_start_codes:
                                    v4_start = seq.find(code)
                                    if v4_start != -1:
                                        break
                                if v4_start != -1:
                                    removed_bp += v4_start
                                    kept_bp += len(seq) - v4_start
                                    output_file.write(seq[v4_start:] + '\n')
                                else:
                                    removed_bp += len(seq)
                    log.write(f'{filename}\t{removed_bp}\t{kept_bp}\n')
