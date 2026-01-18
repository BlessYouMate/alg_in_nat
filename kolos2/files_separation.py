import pandas as pd
import glob
import os

def create_summary_excel():
    output_filename = "Wyniki_Zbiorcze.xlsx"
    csv_files = glob.glob("*_results.csv")
    
    if not csv_files:
        print("Nie znaleziono plików *_results.csv w katalogu!")
        return

    print(f"Znaleziono {len(csv_files)} plików. Rozpoczynam tworzenie Excela...")

    # Sortujemy pliki, żeby były w ładnej kolejności w Excelu
    csv_files.sort()

    with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
        for file in csv_files:
            try:
                # Wczytanie danych
                df = pd.read_csv(file)
                
                # Nazwa arkusza z nazwy pliku (usuwamy _results.csv)
                sheet_name = file.replace("_results.csv", "")
                # Excel ma limit 31 znaków na nazwę arkusza, przycinamy w razie czego
                sheet_name = sheet_name[:31]
                
                print(f"Przetwarzanie: {sheet_name}...")

                # Przygotowanie struktury pod Excela (układ kolumnowy)
                final_df = pd.DataFrame()

                generations = [20, 50, 100, 500]
                
                for gen in generations:
                    # Wybieramy dane tylko dla konkretnej generacji
                    subset = df[df['generation'] == gen][['f1', 'f2']].reset_index(drop=True)
                    
                    # Nazywamy kolumny jasn, np. "Gen 20 f1"
                    subset.columns = [f'Gen {gen} - f1', f'Gen {gen} - f2']
                    
                    # Dodajemy do głównej tabeli (concat łączy "obok siebie", axis=1)
                    final_df = pd.concat([final_df, subset], axis=1)
                    
                    # Opcjonalnie: dodajemy pustą kolumnę dla odstępu (wizualne oddzielenie)
                    final_df[f'sep_{gen}'] = "" 

                # Usuwamy ostatnią pustą kolumnę separatora
                if not final_df.empty:
                    final_df = final_df.iloc[:, :-1]

                # Zapis do arkusza
                final_df.to_excel(writer, sheet_name=sheet_name, index=False)
                
            except Exception as e:
                print(f"Błąd przy przetwarzaniu pliku {file}: {e}")

    print(f"\nGotowe! Plik zapisano jako: {output_filename}")

if __name__ == "__main__":
    # Upewnij się, że masz zainstalowane biblioteki:
    # pip install pandas openpyxl
    create_summary_excel()