# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 23:43:41 2025

@author: magor
"""


# Tworzenie modułu o nazwie biologia.py
def opis_komorki():
    return "Komórka to podstawowa jednostka życia"

def  licz_nukleotydy(sekwencja):
    from collections import Counter
    licznik_nukleotydow = Counter(sekwencja)
    return f"Zliczone nukleotydy: {licznik_nukleotydow}"
    

# Utworzenie nowego katalogu przy użyciu modułu os
# Utworzenie pliku, do którego zpisujemy wyniki z funkcji licz_nukleotydy

import os

nowy_katalog = "dane_bio"

if not os.path.exists(nowy_katalog):
    os.mkdir(nowy_katalog)
    print(f"Katalog '{nowy_katalog}' został utworzony.")

sciezka_do_pliku_1 = os.path.join(nowy_katalog, "nukleotydy.txt")
sekwencja = "AGCTTAGCTAAGGCT"

with open(sciezka_do_pliku_1, 'w') as plik:
    plik.write(licz_nukleotydy(sekwencja))
    print(f"Plik '{sciezka_do_pliku_1}' został utworzony i zapisany.")


# Zapisanie do pliku daty i czasu utworzenia pliku

import datetime

data_czas_utworzenia_pliku = datetime.datetime.now()
sciezka_do_pliku_2 = os.path.join(nowy_katalog, "nukleotydy.txt")

with open(sciezka_do_pliku_2, 'a') as plik:
    plik.write(f"Data i czas utworzenia pliku: {data_czas_utworzenia_pliku}")
    print(f"Data i czas utworzenia pliku: {data_czas_utworzenia_pliku}")
