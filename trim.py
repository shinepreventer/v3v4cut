# 这是一个用于处理测序数据的脚本
# 作者：wzy
# 日期：2023-04-04

# 导入os模块，用于操作文件和目录
import os
# 导入logging模块，用于记录日志信息
import logging
# 导入biopython模块，用于处理测序数据
from Bio import SeqIO

# 定义常量，表示输入目录，输出目录和日志文件的路径
INPUT_DIR = 'H:/v4test'
OUTPUT_DIR = 'H:/v4test/v4new'
LOG_FILE = 'H:/v4test/log.txt'

# 定义常量，表示8个v4起始序列
V4_START_CODES = [
    'GTGACAGCAGCCGCGGTAA',
    'GTGACAGCCGCCGCGGTAA',
    'GTGACAGCGGCCGCGGTAA',
    'GTGACAGCTGCCGCGGTAA',
    'GTGCCAGCAGCCGCGGTAA',
    'GTGCCAGCCGCCGCGGTAA',
    'GTGCCAGCGGCCGCGGTAA',
    'GTGCCAGCTGCCGCGGTAA'
]
# 定义一个函数，用于处理一个测序数据文件
def process_file(filename):
    # 使用os.path.join函数拼接输入文件和输出文件的完整路径
    input_file = os.path.join(INPUT_DIR, filename)
    output_file = os.path.join(OUTPUT_DIR, filename)
    # 初始化移除的碱基数和保留的碱基数为0
    removed_bp = 0
    kept_bp = 0
    # 使用biopython的SeqIO模块来读取和写入测序数据文件
    with open(input_file, 'r') as input_handle:
        with open(output_file, 'w') as output_handle:
            # 遍历输入文件中的每一条记录，每条记录包含标识符，序列和质量分数
            for record in SeqIO.parse(input_handle, 'fastq'):
                # 获取序列字符串
                seq = str(record.seq)
                # 如果文件名以.2.fastq结尾，说明是反向测序数据，需要反转并互补序列
                if filename.endswith('.2.fastq'):
                    seq = str(record.reverse_complement().seq)
                # 初始化v4起始位置为-1，表示没有找到
                v4_start = -1
                # 遍历8个v4起始序列，查找是否在序列中出现
                for code in V4_START_CODES:
                    v4_start = seq.find(code)
                    # 如果找到了，就跳出循环
                    if v4_start != -1:
                        break
                # 如果找到了v4起始位置，就计算移除的碱基数和保留的碱基数，并将保留的部分写入输出文件
                if v4_start != -1:
                    removed_bp += v4_start
                    kept_bp += len(seq) - v4_start
                    output_handle.write(seq[v4_start:] + '\n')
                # 否则，就将整个序列都算作移除的碱基数
                else:
                    removed_bp += len(seq)
    # 返回文件名，移除的碱基数和保留的碱基数
    return filename, removed_bp, kept_bp
    # 使用logging模块来配置日志信息的格式和级别，并指定日志文件的路径
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO, filename=LOG_FILE)
# 在日志文件中写入表头，包含文件名，移除的碱基数和保留的碱基数
logging.info('filename\tremoved_bp\tkept_bp')
# 如果输出目录不存在，就创建它，并使用try-except语句来处理可能出现的异常，并在日志文件中记录异常信息
try:
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    # 遍历输入目录下的所有文件
    for filename in os.listdir(INPUT_DIR):
        # 如果文件名以.fastq结尾，说明是测序数据文件
        if filename.endswith('.fastq'):
            # 使用process_file函数来处理文件，并获取返回的结果
            filename, removed_bp, kept_bp = process_file(filename)
            # 在日志文件中写入结果
            logging.info(f'{filename}\t{removed_bp}\t{kept_bp}')
except Exception as e:
    # 如果出现异常，就在日志文件中记录异常信息，并打印到屏幕上
    logging.error(e)
    print(e)
