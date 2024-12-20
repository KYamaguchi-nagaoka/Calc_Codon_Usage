import re
import streamlit as st
from Bio import SeqIO
from io import StringIO, BytesIO
import pandas as pd

# コドン辞書
codon_dict = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# UをTに置換した新しい辞書を作成
codon_dict_with_t = {key.replace('U', 'T'): value for key, value in codon_dict.items()}
codon_dict = codon_dict_with_t

# コドン使用率を計算する関数
def calculate_codon_usage(cds_sequences):
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

def main():
    st.title("CDS Codon Usage Analyzer")
    
    # ファイルのアップロード
    CDS_file = st.file_uploader("Drop CDS file", type=["fasta"])

    if CDS_file:
        # 配列情報を格納するリスト
        CDS = []

        # Read CDS file
        cds_stringio = StringIO(CDS_file.getvalue().decode("utf-8"))
        for record in SeqIO.parse(cds_stringio, "fasta"):
            CDS.append(str(record.seq))

        # コドン使用率を計算
        codon_usage, codon_freq = calculate_codon_usage(CDS)
        
        # 結果をDataFrameに保存
        data = []
        for codon in codon_usage:
            amino_acid = codon_dict[codon]
            count = codon_usage[codon]
            freq = codon_freq[codon]
            data.append([amino_acid, codon, count, freq])
        
        df = pd.DataFrame(data, columns=["amino_acid", "codon", "count", "freq"])
        df['freq'] = df['freq'].round(9)  # 小数点以下の桁数を指定

        # CSVに変換
        csv = df.to_csv(index=False)
        b_io = BytesIO(csv.encode())

        # ダウンロードボタンの追加
        st.download_button(label="Download Codon Usage CSV", data=b_io, file_name="codon_usage.csv", mime="text/csv")
        st.success("Codon usage analysis completed successfully.")
        st.dataframe(df)  # 解析結果を表示

if __name__ == '__main__':
    main()
