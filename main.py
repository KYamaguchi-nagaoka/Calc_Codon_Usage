from Bio import SeqIO
import pandas as pd
import polars as pl
from gtfparse import read_gtf
import streamlit as st
from io import StringIO, BytesIO
import re
import pyarrow
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


# タブを作成
tab_titles = ['Create CDS Sequences', 'Verify the CDS Sequences', 'Calculate Codon Usage']
tab1, tab2, tab3 = st.tabs(tab_titles)
 
# 各タブにコンテンツを追加
with tab1:
    st.header("Make CDS from All_DNA_seq and gtf_file")
    # ファイルのアップロード
    DNA_file = st.file_uploader("drop ALL_DNA_sequence file", type=["fasta"])
    gtf_file = st.file_uploader("drop gtf file", type=["gtf"])
    if DNA_file and gtf_file:
        sequences = read_fasta_file(DNA_file)
        df = read_gtf(gtf_file)
        df = df[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'protein_id', 'exon_id']]
        
        # Clean and convert protein_id to integer
        df = df.filter(pl.col("protein_id").is_not_null() & (pl.col("protein_id") != ""))
        df = df.with_columns(pl.col("protein_id").cast(pl.Int64).alias("protein_id"))

        df = df.filter(pl.col("exon_id").is_not_null() & (pl.col("exon_id") != ""))
        df = df.with_columns(pl.col("exon_id").cast(pl.Float64).alias("exon_id"))

        CDSs = []
        targets = []

        # 全ユニークなprotein_idに対してループ
        for target in df['protein_id'].unique():
            # 現在のターゲットに一致する行を抽出
            protein_cds = df.filter(pl.col("protein_id") == target)

            # exon_id の ".以降の数字" でソート
            protein_cds = protein_cds.sort("start")

            # DNA配列から対応するCDSを抽出し連結
            complete_cds = ""
            for row in protein_cds.iter_rows(named=True):
                start = row["start"]
                end = row["end"]
                strand = row["strand"]
                frame = row["frame"]
                gene_id = row["seqname"]

                # `sequences` から該当するSeqRecordを検索
                matching_seq = [seq for seq in sequences if seq.id == gene_id]
                if not matching_seq:
                    raise ValueError(f"Gene ID '{gene_id}' not found in sequences.")
                
                # 配列を取得
                sub_seq = matching_seq[0].seq[start - 1 : end]  # 1-based indexに対応

                # 連結
                complete_cds += sub_seq
            
            # strand を考慮
            if strand == "-":
                complete_cds = complete_cds.reverse_complement()

            # 結果を保存
            targets.append(target)
            CDSs.append(complete_cds)

        # Save the CDSs to a FASTA file
        # ファイル名の入力
        filename_input = st.text_input("Enter filename the FASTA file:")
        fasta_content = save_to_fasta(targets, CDSs)
        # ダウンロードボタンの追加
        if filename_input:
            b_io = BytesIO(fasta_content.encode('utf-8'))
            st.download_button(label="Download FASTA", data=b_io, file_name=f"{filename_input}.fasta", mime="text/fasta")
            st.success("CDS extraction and file generation completed successfully.")
 
with tab2:
    st.header('Verify the CDS Sequences')
    # ファイルのアップロード
    CDS_file = st.file_uploader("Drop CDS file", type=["fasta"])
    protein_file = st.file_uploader("Drop protein file", type=["fasta"])

    if CDS_file and protein_file:
        # 配列情報を格納する辞書
        CDS = {}
        proteins_seq = {}

        # Read CDS file
        cds_stringio = StringIO(CDS_file.getvalue().decode("utf-8"))
        for record in SeqIO.parse(cds_stringio, "fasta"):
            protein_sequence = dna_to_protein(str(record.seq),codon_dict)
            CDS[record.id] = protein_sequence

        # Read protein file
        protein_stringio = StringIO(protein_file.getvalue().decode("utf-8"))
        for record in SeqIO.parse(protein_stringio, "fasta"):
            proteins_seq[record.id] = str(record.seq)
        
        # 一致しない配列を確認
        errors = []
        for key in CDS:
            if key in proteins_seq and CDS[key] != proteins_seq[key]:
                errors.append(key)

        st.write(f"Protein と一致しなかった CDS の数: {len(errors)}")
        st.write(f"一致しなかった CDS: {errors}")
    
    
 
with tab3:
    st.header('Calculate Codon Usage')
        # ファイルのアップロード
    CDS_file2 = st.file_uploader("Drop checked CDS file", type=["fasta"])

    if CDS_file2:
        # 配列情報を格納するリスト
        CDS = []

        # Read CDS file
        cds_stringio = StringIO(CDS_file2.getvalue().decode("utf-8"))
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
    

