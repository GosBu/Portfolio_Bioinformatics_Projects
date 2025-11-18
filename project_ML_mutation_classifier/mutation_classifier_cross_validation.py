#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 19:06:34 2025

@author: mago
"""

"""
Projekt 2: Przewidywanie mutacji w sekwencjach DNA związanych z chorobami genetycznymi

Cel projektu: Celem jest stworzenie modelu klasyfikacyjnego, który określi, czy dana mutacja w sekwencji DNA jest patogenna.\
    
"""

import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# Pobranie danych (wczytanie pliku z danymi)
df = pd.read_csv("/Users/mago/Downloads/dane_projekt2.csv")


# Eksploracja danych
# Wyświetlenie podstawowych statystyk
print(df.describe())
print(df.info())

# Ile braków w każdej kolumnie
print(df.isna().sum())

# Długości sekwencji DNA
lengths = df["DNA_Sequence"].str.len()
print(lengths.head())       # Podejrzenie kilku pierwszych
print(lengths.describe())       # Statystyki

# Ile braków w każdej kolumnie
print(df.isna().sum())

# Najczęstsze pozycje mutacji
position_counts = df["Mutation_Position"].value_counts().sort_values(ascending=False)
print("\nNajczęstsze pozycje mutacji:\n", position_counts.head(15))

# Najrzadsze pozycje mutacji
rare_positions = df["Mutation_Position"].value_counts()
print("\nNajrzadsze pozycje mutacji:\n", rare_positions.tail(15))

# Całkowita liczba mutacji w każdym genie
gene_counts = df["Gene"].value_counts()
print("\nLiczba mutacji w każdym genie:\n", gene_counts.head(10))

# Liczba mutacji patogennych w każdym genie
pathogenic_counts = df[df["Pathogenicity"] == "pathogenic"]["Gene"].value_counts()
print("\nLiczba mutacji patogennych w każdym genie: \n", pathogenic_counts.head(10))

# Liczba mutacji patogennych w każdym genie
pathogenic = df[df["Pathogenicity"] == "pathogenic"]        # Filtrujemy mutacje patogenne

pathogenic_changes = pathogenic.groupby(["Gene", "Reference_Base", "Alternate_Base"]).size().sort_values(ascending=False)       # Grupujemy po genie i zmianie nukleotydu, liczymy wystąpienia

# Wyświetlamy top 10
print("\nNajczęstsze mutacje patogenne (Gen + zmiana nukleotydu):\n", pathogenic_changes.head(10))

# Procent mutacji patogennych w każdym genie
percent_pathogenic = (pathogenic_counts / gene_counts * 100).round(2)
print("\nProcent mutacji patogennych w każdym genie:\n", percent_pathogenic.head(10))

# Liczba wystąpień każdej mutacji (Ref>Alt) oraz patogenność
mutation_counts = df.groupby(["Reference_Base", "Alternate_Base", "Pathogenicity"]).size().sort_values(ascending=False)
print("\nNajczęstsze mutacje (Ref>Alt : Pathogenicity):\n", mutation_counts.head(20))


# Ekstrakcja cech - przetwarzanie danych

# Funkcja do one-hot encoding dla DNA_Sequence (bianrne kodowanie nukleotydów)
def one_hot_sequence(seq):
    mapping = {"A": [1,0,0,0],
               "C": [0,1,0,0],
               "G": [0,0,1,0],
               "T": [0,0,0,1]}
    return [bit for nucleotide in seq for bit in mapping[nucleotide]]

# One-hot encoding dla Reference_Base i Alternate_Base
def one_hot_base(base):
    mapping = {"A": [1,0,0,0],
               "C": [0,1,0,0],
               "G": [0,0,1,0],
               "T": [0,0,0,1]}
    return mapping[base]

# Kodowanie sekwencji DNA
ecoding_seq = df["DNA_Sequence"].apply(one_hot_sequence)
seq_cols = [f"Pos{i}_{nt}" for i in range(50) for nt in ["A","C","G","T"]]
ecoding_seq_df = pd.DataFrame(ecoding_seq.tolist(), columns=seq_cols)

# Kodowanie Reference_Base
ecoding_ref = df["Reference_Base"].apply(one_hot_base)
ref_cols = [f"Ref_{nt}" for nt in ["A","C","G","T"]]
ecoding_ref_df = pd.DataFrame(ecoding_ref.tolist(), columns=ref_cols)

# Kodowanie Alternate_Base
ecoding_alt = df["Alternate_Base"].apply(one_hot_base)
alt_cols = [f"Alt_{nt}" for nt in ["A","C","G","T"]]
ecoding_alt_df = pd.DataFrame(ecoding_alt.tolist(), columns=alt_cols)

# Kodowanie Mutation_Position - zostawiamy jako liczbę
ecoding_pos_df = df[["Mutation_Position"]]

# Łączymy wszystkie cechy 
X_final = pd.concat([ecoding_seq_df, ecoding_ref_df, ecoding_alt_df, ecoding_pos_df], axis=1)
print(X_final.head())

# Kodowanie Pathogenicity
ecoding_pathogenicity = pd.get_dummies(df["Pathogenicity"], prefix="Patho")

# Przypisujemy cechy X i etykiety y
X = X_final
y = df["Pathogenicity"].map({"benign":0, "pathogenic":1})

# Skalowanie Mutation_Position
scaler = MinMaxScaler()
X_final["Mutation_Position"] = scaler.fit_transform(X_final[["Mutation_Position"]])


# Podział na zbiory treningowy i testowy (70% treningowe + walidacyjne, 30% testowe)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=62)

# Wypisanie informacji
total = len(X)
print(f"\nTrain: {len(X_train)} próbek ({len(X_train)/total:.1%})")
print(f"\nTest:  {len(X_test)} próbek ({len(X_test)/total:.1%})")


# Budowa modelu ML (z ważeniem klas)
model = LogisticRegression(solver='liblinear', max_iter=800, class_weight='balanced', random_state=62)
kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=62)         # podobny rozkład klas dla każdej części 


# Ewaluacja modelu
scores = cross_val_score(model, X, y, cv=kf, scoring='roc_auc')         # Miara zdolności modelu do rozróżniania klas
print("ROC AUC dla każdej iteracji:", scores)
print("Średni ROC AUC:", np.mean(scores))


# Interpretacja wyników
# Rysowanie wykresu "Ilość mutacji z podziałem patogenne i nie"
# Dane do wykresu
total_mutations = df["Gene"].value_counts()                  # wszystkie mutacje
pathogenic_mutations = df[df["Pathogenicity"] == "pathogenic"]["Gene"].value_counts()  # mutacje patogenne
genes = total_mutations.index 

# Tworzenie wykresu słupkowego
plt.figure(figsize=(12,6))
bar_width = 0.4
x = range(len(genes))
plt.bar(x, total_mutations, width=bar_width, label="Całkowite mutacje", color="purple")
plt.bar([i + bar_width for i in x], pathogenic_mutations, width=bar_width, label="Mutacje patogenne", color="green")
plt.xticks([i + bar_width/2 for i in x], genes, rotation=45)
plt.xlabel("Gen")
plt.ylabel("Liczba mutacji")
plt.title("Ilość mutacji z podziałem patogenne i nie")
plt.legend()
plt.tight_layout()
plt.show()
