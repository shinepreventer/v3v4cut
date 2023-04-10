# 16SrDNA测序v4区裁剪脚本

这是一个用python写的脚本，用于从16SrDNA测序（v3v4区）的下机文件原始序列中提取出v4区部分，并保留原始信息和正常的格式，使其可以导入qiime2进行后续的分析。

## 依赖

- python 3.7或更高版本
- cutadapt 4.3或更高版本
- os模块

## 使用方法

1. 将您的原始文件放在/media/wzy/work/v4test目录下，确保它们是fastq格式。
2. 修改脚本中的引物序列，以匹配您的实验设计。
3. 运行脚本：trim.py
4. 检查输出文件和log文件，它们将被保存在/media/wzy/work/v4test/v4new目录下。
5. 将输出文件导入qiime2进行后续的分析。

## 注意事项

- 本脚本使用cutadapt软件来实现序列裁剪功能，如果您没有安装cutadapt，请先安装它。
- 本脚本没有进行任何质量控制或过滤操作，如果您需要对您的数据进行进一步的处理，请在运行本脚本之前或之后自行完成。
- 本脚本没有进行任何错误处理或异常处理，如果您遇到任何问题，请检查您的输入文件和参数是否正确，或者联系作者。
- cutadapt安装请参考https://cutadapt.readthedocs.io/en/stable/installation.html

## 联系方式

如果您有任何问题、建议或反馈，请发送邮件至wzyttkx@gmail.com，谢谢！
