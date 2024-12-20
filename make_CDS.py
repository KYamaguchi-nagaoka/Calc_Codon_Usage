from Bio import SeqIO
import pandas as pd
import polars as pl
from gtfparse import read_gtf
import streamlit as st
from io import StringIO, BytesIO

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

def main():
    st.title("Make CDS from All_DNA_seq and gtf_file")
    
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
        filename_input = st.text_input("Enter filename for the FASTA file:")
        fasta_content = save_to_fasta(targets, CDSs)
        b_io = BytesIO(fasta_content.encode('utf-8'))
        st.download_button(label="Download FASTA", data=b_io, file_name=filename_input+".fasta", mime="text/fasta")
        st.success("CDS extraction and file generation completed successfully.")

if __name__ == '__main__':
    main()
