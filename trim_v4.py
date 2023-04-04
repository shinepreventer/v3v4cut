# 这是一个用于处理测序数据的脚本
import os
import logging
from Bio import SeqIO
import argparse
import multiprocessing

# 定义常量
INPUT_DIR = 'H:/v4test'
OUTPUT_DIR = 'H:/v4test/v4new'
LOG_FILE = 'H:/v4test/log.txt'

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
    # 拼接文件路径
    input_file = os.path.join(args.input_dir, filename)
    output_file = os.path.join(args.output_dir, filename)
    # 初始化计数器
    removed_bp = 0
    kept_bp = 0
    # 读写文件
    with open(input_file, 'r') as input_handle:
        with open(output_file, 'w') as output_handle:
            # 处理每条记录
            for record in SeqIO.parse(input_handle, 'fastq'):
                # 获取序列
                seq = str(record.seq)
                # 判断反向测序
                if filename.endswith('.2.fastq'):
                    seq = str(record.reverse_complement().seq)
                # 查找v4起始位置
                v4_start = -1
                for code in V4_START_CODES:
                    v4_start = seq.find(code)
                    if v4_start != -1:
                        break
                # 计算移除和保留的碱基数，并写入输出文件
                if v4_start != -1:
                    removed_bp += v4_start
                    kept_bp += len(seq) - v4_start
                    output_handle.write(seq[v4_start:] + '\n')
                else:
                    removed_bp += len(seq)
    # 返回结果
    return filename, removed_bp, kept_bp

# 创建解析器对象，并添加参数
parser = argparse.ArgumentParser(description='Process sequencing data files.')
parser.add_argument('--input_dir', type=str, default='H:/v4test', help='The input directory.')
parser.add_argument('--output_dir', type=str, default='H:/v4test/v4new', help='The output directory.')
parser.add_argument('--log_file', type=str, default='H:/v4test/log.txt', help='The log file.')
# 解析参数，并赋值给args变量
args = parser.parse_args()

# 配置日志信息，并写入表头
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO, filename=args.log_file)
logging.info('filename\tremoved_bp\tkept_bp')
# 处理异常情况，并创建输出目录
try:
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # 创建进程池，并指定进程数为CPU核心数
    pool = multiprocessing.Pool(multiprocessing.cpu_count()) #pool = multiprocessing.Pool(4) # 创建一个包含4个进程的进程池
    # 定义一个函数，用于记录日志信息
    def log_result(result):
        filename, removed_bp, kept_bp = result
        logging.info(f'{filename}\t{removed_bp}\t{kept_bp}')
    # 遍历输入目录下的所有文件，并使用process_file函数来处理文件，并将结果传递给log_result函数
    for filename in os.listdir(args.input_dir):
        if filename.endswith('.fastq'):
            pool.apply_async(process_file, args=(filename,), callback=log_result)
    # 关闭进程池并等待所有进程结束
    pool.close()
    pool.join()
except Exception as e:
    # 记录并打印异常信息
    logging.error(e)
    print(e)
GTGACAGCAGCCGCGGTAA,GTGACAGCCGCCGCGGTAA,GTGACAGCGGCCGCGGTAA,GTGACAGCTGCCGCGGTAA,GTGCCAGCAGCCGCGGTAA,GTGCCAGCCGCCGCGGTAA,GTGCCAGCGGCCGCGGTAA,GTGCCAGCTGCCGCGGTAA
