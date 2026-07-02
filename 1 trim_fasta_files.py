import os
import regex as re  
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def iupac_to_regex(seq):
    iupac_dict = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    return ''.join([iupac_dict[base] for base in seq])

def find_v4b_region(seq, forward_primer_regex, reverse_primer_regex):
    forward_match = re.search(forward_primer_regex, str(seq), flags=re.BESTMATCH)
    reverse_match = re.search(reverse_primer_regex, str(seq), flags=re.BESTMATCH)
    
    if forward_match and reverse_match:
        return seq[forward_match.start():reverse_match.end()]
    else:
        return None

def process_files(trimmomatic_path, input_dir, output_dir, primers_file, forward_primer, reverse_primer):
    forward_primer_regex = f"({iupac_to_regex(forward_primer)}){{e<=1}}"
    reverse_primer_regex = f"({iupac_to_regex(str(Seq(reverse_primer).reverse_complement()))}){{e<=1}}"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    failed_files = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            input_fasta = os.path.join(input_dir, filename)
            output_fasta = os.path.join(output_dir, filename.rsplit('.', 1)[0] + "_v4b_extracted.fasta")
            
            records = []
            for record in SeqIO.parse(input_fasta, "fasta"):
                v4b_seq = find_v4b_region(record.seq, forward_primer_regex, reverse_primer_regex)
                if v4b_seq:
                    new_record = SeqRecord(v4b_seq, id=record.id, description="V4b region")
                    records.append(new_record)
            
            if records:
                SeqIO.write(records, output_fasta, "fasta")
                print(f"V4b区域提取完成: {output_fasta}")
            else:
                print(f"未找到V4b区域: {filename}")
                failed_files.append(filename)
    
    if failed_files:
        failed_output_tsv = os.path.join(output_dir, "failed_files.tsv")
        with open(failed_output_tsv, "w") as tsvfile:
            tsvfile.write("Failed Files\n")
            for failed_file in failed_files:
                tsvfile.write(f"{failed_file}\n")
        print(f"未成功提取V4b区域的文件已保存到 {failed_output_tsv}")

# 运行代码
input_dir = "syncom"
output_dir = "genomic_16sfull"
primers_file = "adapters.fa"
trimmomatic_path = "C:\\Users\\YZY\\Trimmomatic-0.39\\trimmomatic-0.39.jar"
forward_primer = ""
reverse_primer = ""

process_files(trimmomatic_path, input_dir, output_dir, primers_file, forward_primer, reverse_primer)
print("所有文件已处理完成。")
