
from io import StringIO, BytesIO
from Bio import SeqIO
def read_fasta_file(file):
    sequences = []
    stringio = StringIO(file.getvalue().decode("utf-8"))
    for record in SeqIO.parse(stringio, "fasta"):
        sequences.append(record)
    return sequences

def save_to_fasta(targets, CDSs):
    fasta_content = ""
    for target, cds in zip(targets, CDSs):
        fasta_content += f">{target}\n"
        for i in range(0, len(cds), 50):
            fasta_content += str(cds[i:i + 50]) + "\n"
    return fasta_content


# DNA配列をアミノ酸配列に変換する関数
def dna_to_protein(dna_sequence,codon_dict):
    amino_acids = []
    for i in range(0, len(dna_sequence) - 2, 3):  # 3つずつ処理
        codon = dna_sequence[i:i+3]
        if codon in codon_dict:
            amino_acids.append(codon_dict[codon])  # コドンからアミノ酸を取得
        else:
            amino_acids.append('?')  # 無効なコドンの場合
    return ''.join(amino_acids)

# コドン使用率を計算する関数
def calculate_codon_usage(cds_sequences,codon_dict):
    codon_usage = {codon: 0 for codon in codon_dict.keys()}
    total_codons = 0

    for sequence in cds_sequences:
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in codon_usage:
                codon_usage[codon] += 1
                total_codons += 1

    amino_acid_usage = {aa: 0 for aa in set(codon_dict.values())}
    for codon, aa in codon_dict.items():
        amino_acid_usage[aa] += codon_usage[codon]

    codon_freq = {codon: (count / amino_acid_usage[codon_dict[codon]]) if amino_acid_usage[codon_dict[codon]] else 0 for codon, count in codon_usage.items()}
    return codon_usage, codon_freq
 