import re
import streamlit as st
from Bio import SeqIO
from io import StringIO
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
    st.title("CDS Checker")
    
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

if __name__ == '__main__':
    main()


