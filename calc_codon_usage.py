import re
import streamlit as st
from Bio import SeqIO
from io import StringIO, BytesIO
import pandas as pd
from sub_folder.function import *

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
        codon_usage, codon_freq = calculate_codon_usage(CDS,codon_dict)
        
        # 結果をDataFrameに保存
        data = []
        for codon in codon_usage:
            amino_acid = codon_dict[codon]
            count = codon_usage[codon]
            freq = codon_freq[codon]
            data.append([amino_acid, codon, count, freq])
        
        df = pd.DataFrame(data, columns=["amino_acid", "codon", "count", "freq"])
        df['freq'] = df['freq'].round(3)  # 小数点以下の桁数を指定

        # CSVに変換
        csv = df.to_csv(index=False)

# ファイル名の入力
        filename_input = st.text_input("Enter filename for the CSV file:")
        
        # ダウンロードボタンの追加
        if filename_input:
            b_io = BytesIO(csv.encode())
            st.download_button(label="Download Codon Usage CSV", data=b_io, file_name=f"{filename_input}.csv", mime="text/csv")
            st.success("Codon usage analysis completed successfully.")
        
        st.dataframe(df)  # 解析結果を表示
if __name__ == '__main__':
    main()
