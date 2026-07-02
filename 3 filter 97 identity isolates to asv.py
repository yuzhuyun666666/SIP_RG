import os
import pandas as pd

# 文件夹路径，存放所有的txt文件
folder_path = 'derep_isolates_alignment'
output_csv = 'blast_tsb_97_derep.csv'

# 初始化一个空的DataFrame来存储所有满足条件的行
all_filtered_data = pd.DataFrame()

# 遍历文件夹中的所有txt文件
for filename in os.listdir(folder_path):
    if filename.endswith('.txt'):
        file_path = os.path.join(folder_path, filename)
        # 读取txt文件，假设文件以空格或制表符分隔
        data = pd.read_csv(file_path, sep='\t', header=None)

        # 筛选第三列数据大于97的行
        filtered_data = data[data[2] > 97]

        # 将筛选后的数据追加到总的DataFrame中
        all_filtered_data = pd.concat([all_filtered_data, filtered_data])

# 将合并后的数据保存到一个CSV文件中
all_filtered_data.to_csv(output_csv, index=False, header=False)
