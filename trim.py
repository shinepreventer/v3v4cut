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

# 定义一个函数来判断文件类型
def check_file_type(file):
    # 打开文件
    with open(file) as f:
        # 创建一个空集合来存储序列ID
        ids = set()
        # 遍历每一行
        for line in f:
            # 如果是序列ID行
            if line.startswith("@"):
                # 去掉换行符和最后一个数字
                id = line.strip().rsplit("/", 1)[0]
                # 如果id已经在集合中，说明是双端
                if id in ids:
                    return "paired"
                # 否则，将id加入集合
                else:
                    ids.add(id)
        # 如果遍历完没有发现重复的id，说明是单端
        return "single"

# 遍历输入文件夹中的所有文件
for file in os.listdir(input_dir):
    # 如果文件是fastq格式，进行裁剪
    if file.endswith(".fastq"):
        # 指定输入文件和输出文件的完整路径
        input_file = os.path.join(input_dir, file)
        output_file = os.path.join(output_dir, file)
        # 如果文件名包含_1，说明是双端文件的第一个文件
        if "_1" in file:
            # 构造另一个输入文件和输出文件的完整路径，将_1替换为_2
            input_file2 = input_file.replace("_1", "_2")
            output_file2 = output_file.replace("_1", "_2")
            # 构造cutadapt命令，指定正向和反向引物序列，以及另一个输出文件  -m 10为去除切除后小于10bp的序列
            cutadapt_command = f"cutadapt -g {forward_primer} -o {output_file} {input_file} -m 10 --report minimal"
        # 如果文件名包含_2，说明是双端文件的第二个文件
        elif "_2" in file:
            # 构造cutadapt命令，只指定反向引物序列
            cutadapt_command = f"cutadapt -a {reverse_primer} -o {output_file} {input_file} -m 10 --report minimal"
        # 否则，调用函数判断文件类型
        else:
            file_type = check_file_type(input_file)
            # 如果是双端合并后的文件
            if file_type == "paired":
                # 构造cutadapt命令，指定正向和反向引物序列
                cutadapt_command = f"cutadapt -g {forward_primer} -a {reverse_primer} -o {output_file} {input_file} -m 10 --report minimal"
            # 如果是单独的左端文件
            elif file_type == "single":
                # 构造cutadapt命令，只指定正向引物序列
                cutadapt_command = f"cutadapt -g {forward_primer} -o {output_file} {input_file} -m 10 --report minimal"
        # 执行cutadapt命令，并将输出重定向到log文件中
        os.system(cutadapt_command + f" >> {output_dir}/log.txt")
