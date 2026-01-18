import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

def plot_benchmark_results():
    # Znajdź wszystkie pliki z wynikami
    files = glob.glob("*_results.csv")

    if not files:
        print("Nie znaleziono plików CSV. Uruchom program C++.")
        return

    # Kolory dla poszczególnych checkpointów
    colors = {
        20: 'lightgray',
        50: 'limegreen',
        100: 'orange',
        500: 'blue'
    }

    alphas = {
        20: 0.5,
        50: 0.7,
        100: 0.8,
        500: 1.0
    }

    labels = {
        20: 'Gen 20',
        50: 'Gen 50',
        100: 'Gen 100',
        500: 'Gen 500 (Final)'
    }

    for file in files:
        print(f"Generowanie wykresu dla: {file}")
        try:
            data = pd.read_csv(file)

            # Zwiększyłem nieco szerokość figury, żeby legenda ładnie wyglądała obok
            plt.figure(figsize=(11, 8))

            # Grupujemy po generacjach
            for gen in [20, 50, 100, 500]:
                subset = data[data['generation'] == gen]
                if not subset.empty:
                    plt.scatter(
                        subset['f1'],
                        subset['f2'],
                        c=colors.get(gen, 'black'),
                        label=labels.get(gen, str(gen)),
                        s=15 if gen == 500 else 10,
                        alpha=alphas.get(gen, 1.0)
                    )

            # Tytuł i opisy na podstawie nazwy pliku
            title_name = file.replace("_results.csv", "")
            plt.title(f'Pareto Front Evolution: {title_name}')
            plt.xlabel('$f_1(x)$')
            plt.ylabel('$f_2(x)$')
            plt.grid(True, linestyle='--', alpha=0.5)
            
            # --- ZMIANA: LEGENDA POZA WYKRESEM ---
            # bbox_to_anchor=(1.05, 1) -> przesunięcie w prawo o 5% szerokości wykresu
            # loc='upper left' -> lewy górny róg legendy styka się z tym punktem
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            # -------------------------------------

            # --- WAŻNE: DOPASOWANIE UKŁADU ---
            # Bez tego legenda zostałaby ucięta przy zapisie do pliku!
            plt.tight_layout()
            # ---------------------------------

            # Zapisz jako PNG
            output_file = title_name + "_plot.png"
            plt.savefig(output_file, dpi=300)
            plt.close()

        except Exception as e:
            print(f"Błąd przy pliku {file}: {e}")

if __name__ == "__main__":
    plot_benchmark_results()