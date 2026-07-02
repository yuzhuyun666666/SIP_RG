import os
import subprocess

def run_blastn(query_file, db, output_file):
    blastn_cmd = [
        'blastn', '-query', query_file, '-db', db,
        '-out', output_file, '-outfmt', '6 qseqid sseqid pident length qlen qstart qend sstart send qcovs'
    ]
    subprocess.run(blastn_cmd, check=True)

def calculate_coverage(length, qlen):
    return (length / qlen) * 100

def process_blast_output(output_file):
    identity_list = []
    coverage_list = []

    with open(output_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            qseqid, sseqid, pident, length, qlen, slen, qstart, qend, sstart, send = fields
            identity = float(pident)
            coverage = calculate_coverage(int(length), int(qlen))
            identity_list.append(identity)
            coverage_list.append(coverage)

    return identity_list, coverage_list

def process_files(input_dir, db, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            query_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename.rsplit('.', 1)[0] + "_blast_output.txt")
            
            # Run BLAST
            run_blastn(query_file, db, output_file)
            
            # Process BLAST output
            identity_list, coverage_list = process_blast_output(output_file)
            print(f"Results for {filename}:")
            print(f"Identity: {identity_list}")
            print(f"Coverage: {coverage_list}")

# 输入和输出目录路径
input_dir = "isolate_seperated_asv_fasta_files"
output_dir = "derep_isolates_alignment"
db = "sequences_db/sequences"

# 处理所有文件
process_files(input_dir, db, output_dir)

print("所有文件已处理完成。")
