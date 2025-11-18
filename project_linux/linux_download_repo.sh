#!/bin/bash

# Utworzenie katalogu o zadanej nazwie

# Komunikat dla uzytkownika
echo "Musimy utworzyc nowy katalog. Jak chcesz nazwac nowy katalog, w ktorym beda przechowywane pobrane pliki?"

# Odczyt podanej przez uzytkownika nazwy
read NAZWA_KATALOGU

# Sprawdzamy czy nazwa zostala wpisana
if [ -z "$NAZWA_KATALOGU" ]; then
	echo "Blad - aby uruchomic skrypt musisz dodac nawze tworzonego katalogu. Nie podano nazwy, wiec skrypt zakonczy dzialanie."
	exit 1
fi

# Tworzenie katalogu o podanej nazwie
mkdir -p "$NAZWA_KATALOGU"



# Pobieranie danych bezposrednio z repozytorium GitHub za pomoca adresu https (bez kolonowania rezpoytorium)

# Komunikat dla uzytkownika
echo "Pobrane pliki beda zapisane jako archiwum zip. Podaj nazwe nowo utworzonego archium zip, w ktorym przechowamy pobrane pliki "

# Odczyt podanej nazwy
read PLIK_ZIP

# Sprawdzamy czy nazwa zostala wspiana
if [ -z "$PLIK_ZIP" ]; then
	echo "Blad, nie podales nazwy pliku. Skrypt zakonczy dzialanie."
	exit 1
fi

# Komunikat dla uzytkownika
echo "Wklej skopiowany link (https://github.com/USER/REPO) repozytorium, ktore chcesz pobrac: "

# Oczyt wprowadzonego linku
read LINK

# Sprawdamy, czy link zostal podany
if [ -z "$LINK" ]; then
	echo "Blad, link nie zostal wprowadzony. Skrypt zakonczy dzialanie."
	exit 1
fi

if [[ "$LINK" == *.git ]]; then
	echo "Blad w linku - link nie moze konczyc sie na .git. Poprawny link ma strukture: https://github.com/USER/REPO). Skrypt zakonczy dzialanie."
	exit 1
fi

# Zapisanie linku w odpowiednim formacie, aby pobrac repozytorium jako archiwum zip
LINK_ZIP="${LINK}/archive/refs/heads/main.zip"

# Pobranie plikow z podanego linku i zapisanie ich od razu w formacie archiwum zip
wget "$LINK_ZIP" -O "${NAZWA_KATALOGU}/${PLIK_ZIP}.zip"

# Sprawdzamy, czy pobranie pliku sie powiodlo
if [ $? -eq 0 ]; then
	echo "Pobrano repozyorium jako archiwum zip."
# Komunikat dla uzytkowania pojawiajacy sie po zakonczeniu dzialania skryptu

	echo "Rozmiar pobranego pliku: "
	du -h "${NAZWA_KATALOGU}/${PLIK_ZIP}.zip"
else
	echo "Blad pobierania repozytorium."
	exit 1
fi
