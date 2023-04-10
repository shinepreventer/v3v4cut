# 导入os模块
import os

# 指定输入文件夹和输出文件夹
input_dir = "/media/wzy/work/v4test"
output_dir = "/media/wzy/work/v4test/v4new"

# 如果输出文件夹不存在，创建它
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 指定正向和反向引物序列
forward_primer = "GTGCCAGCMGCCGCGGTAA"
reverse_primer = "GGACTACHVGGGTWTCTAAT"

# 遍历输入文件夹中的所有文件
for file in os.listdir(input_dir):
    # 如果文件是fastq格式，进行裁剪
    if file.endswith(".fastq"):
        # 指定输入文件和输出文件的完整路径
        input_file = os.path.join(input_dir, file)
        output_file = os.path.join(output_dir, file)
        # 构造cutadapt命令
        cutadapt_command = f"cutadapt -g {forward_primer} -a {reverse_primer} -o {output_file} {input_file} --report minimal"
        # 执行cutadapt命令，并将输出重定向到log文件中
        os.system(cutadapt_command + f" >> {output_dir}/log.txt")
