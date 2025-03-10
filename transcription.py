import os
import tkinter as tk
from tkinter import filedialog

def find_genes(dna_sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    genes = []
    
    i = 0
    while i < len(dna_sequence) - 2:
        codon = dna_sequence[i:i+3]
        if codon == start_codon:
            for j in range(i + 3, len(dna_sequence) - 2, 3):
                stop_codon = dna_sequence[j:j+3]
                if stop_codon in stop_codons:
                    gene_seq = dna_sequence[i:j+3]
                    genes.append((i, j + 3, gene_seq))
                    i = j + 3
                    break
        i += 1
    
    return genes

# Таблица кодонов для трансляции
codon_table = {
   "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*", 
    "TGA": "*", "TGT": "C", "TGC": "C", "TGG": "W"
}

def translate_gene(dna_gene):
    protein = "".join(codon_table.get(dna_gene[i:i+3], "?") for i in range(0, len(dna_gene), 3))
    return protein

def read_dna_from_file():
    root = tk.Tk()
    root.withdraw()
    filename = filedialog.askopenfilename(title="Выберите файл с ДНК", filetypes=[("Text Files", "*.txt")])
    if not filename:
        return None
    with open(filename, "r") as file:
        dna_sequence = file.read().strip().replace(" ", "").replace("\n", "")
    return dna_sequence

def save_proteins_to_file(proteins):
    unique_proteins = {protein for protein in proteins if len(protein) >= 50}  # Фильтр по длине
    filename = "proteins_filtered.txt"
    with open(filename, "w") as file:
        for i, protein in enumerate(unique_proteins):
            file.write(f">Protein_{i+1}\n{protein}\n")
    print(f"Найдено {len(unique_proteins)} уникальных белков длиной >= 50 аминокислот. Они сохранены в {filename}.")

dna_sequence = read_dna_from_file()
if dna_sequence:
    genes = find_genes(dna_sequence)
    proteins = [translate_gene(gene_seq) for _, _, gene_seq in genes]
    save_proteins_to_file(proteins)
else:
    print("Файл не был выбран.")