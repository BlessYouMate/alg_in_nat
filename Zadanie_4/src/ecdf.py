import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt

# Ustalenie ścieżek
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.path.abspath('.') # Fallback dla interaktywnych notatników

# Ustaw folder 'results' jako podfolder miejsca, gdzie jest skrypt
base_results_path = os.path.join(script_dir, "results")

print(f"Rozpoczęcie analizy")
print(f"Katalog bazowy skryptu: {script_dir}")
print(f"Szukanie wyników w: {base_results_path}")


# Konfiguracja

# 1. Definicja stałych
F1_OPT = -5.5562  # Optimum dla F1
F2_OPT = 0.0     # Optimum dla F2

# 2. Definicja 51 Prógów Jakości (Delta)
THRESHOLDS = np.logspace(1, -8, 51)

# 3. Definicja wzorców plików (względnych do folderu 'results')
PATHS = {
    'binary': [
        (os.path.join(base_results_path, 'f1', 'bin', '*.txt'), F1_OPT),
        # (os.path.join(base_results_path, 'f1', 'gray', '*.txt'), F1_OPT),
        (os.path.join(base_results_path, 'f2', 'bin', '*.txt'), F2_OPT),
        # (os.path.join(base_results_path, 'f2', 'gray', '*.txt'), F2_OPT)
    ],
    'real': [
        (os.path.join(base_results_path, 'f1', 'real', '*.txt'), F1_OPT),
        (os.path.join(base_results_path, 'f2', 'real', '*.txt'), F2_OPT)
    ]
}

# 4. Definicja progów DO WYKRESU (różne dla różnych zestawów)
PLOT_THRESHOLDS = {
    'binary': [
        1.25892541179417, 0.831763771102671, 0.549540873857625, 
        0.363078054770101, 0.239883291901949, 0.158489319246111, 
        0.10471285480509, 0.0691830970918936, 0.0457088189614875,
        0.0301995172040202, 0.0199526231496888, 0.0131825673855641,
        0.0087096358995608, 0.00575439937337157, 0.00380189396320561,
        0.00251188643150958, 0.00165958690743756
    ],
    'real': [
        0.0301995172040202,
        0.0199526231496888,
        0.0131825673855641,
        0.0087096358995608
    ]
}

def generate_ecdf_data(file_groups, thresholds):
    # Wczytanie i konwersja na "Dane Surowe - Delta"
    all_delta_series = []
    print(" Wczytuję dane surowe")
    
    total_files_found = 0
    min_file_length = float('inf') # Zmienna do śledzenia długości plików

    for pattern, f_opt in file_groups:
        file_paths = glob.glob(pattern)
        print(f" znalaziono {len(file_paths)} plików dla wzorca {pattern}")
        total_files_found += len(file_paths)
        
        for file_path in file_paths:
            # Wczytaj CAŁY plik .txt
            raw_fitness = pd.read_csv(file_path, header=None).squeeze("columns")
            
            # Śledź najkrótszą długość pliku
            if len(raw_fitness) < min_file_length:
                min_file_length = len(raw_fitness)
            
            # Przelicz fitness na błąd (delta)
            delta_series = raw_fitness - f_opt
            
            # Upewnij się, że delta nigdy nie jest ujemna (błędy numeryczne blisko 0)
            delta_series[delta_series < 0] = 0
            
            all_delta_series.append(delta_series)

    # Sprawdzenie, czy cokolwiek znaleziono
    if total_files_found == 0 or len(all_delta_series) == 0:
        print(f"KRYTYCZNY BŁĄD: Nie znaleziono żadnych plików .txt w oczekiwanych lokalizacjach")
        print(f"Sprawdź ścieżki w sekcji 'PATHS' skryptu")
        return None, 0

    # DYNAMICZNE USTALENIE MAX_ITERATIONS
    # Użyj długości najkrótszego pliku jako limitu
    max_iters = min_file_length 
    print(f"  automatycznie wykryto długość serii: {max_iters} iteracji (wierszy)")
    
    # Przytnij wszystkie serie do tej samej, minimalnej długości
    all_delta_series_truncated = [s.iloc[:max_iters] for s in all_delta_series]

    # Połącz wszystkie serie w jedną dużą tabelę (DataFrame)
    raw_delta_df = pd.concat(all_delta_series_truncated, axis=1)
    raw_delta_df.columns = [f"run_{i+1}" for i in range(len(raw_delta_df.columns))]
    raw_delta_df.index.name = "Iteracja"
    
    total_runs = len(raw_delta_df.columns)
    print(f" połączono w tabelę 'Dane Surowe - Delta' (Wymiary: {raw_delta_df.shape})")

    # Stworzenie Tabeli "Czasów Trafień"
    print(" Obliczam czasy trafień")
    
    hit_time_df = pd.DataFrame(index=thresholds, columns=raw_delta_df.columns, dtype=float)
    hit_time_df.index.name = "Próg Jakości (Delta)"
    
    for run_name in raw_delta_df.columns:
        run_data = raw_delta_df[run_name] 
        for th in thresholds:
            hit_indices = np.where(run_data <= th)[0]
            
            if len(hit_indices) > 0:
                first_hit_time = hit_indices[0]
                hit_time_df.loc[th, run_name] = first_hit_time
            else:
                hit_time_df.loc[th, run_name] = max_iters + 1 # Użyj dynamicznego max_iters

    # Stworzenie Tabeli do Wykresu ECDF
    print("  Obliczam finalne dane ECDF")
    
    # Użyj dynamicznego max_iters
    ecdf_plot_df = pd.DataFrame(index=range(max_iters), columns=thresholds)
    ecdf_plot_df.index.name = "Iteracja"
    
    for it in range(max_iters):
        for th in thresholds:
            success_count = (hit_time_df.loc[th] <= it).sum()
            proportion = success_count / total_runs
            ecdf_plot_df.loc[it, th] = proportion

    print(" zakończono przetwarzanie")
    return (raw_delta_df, pd.DataFrame(thresholds, columns=["Progi Jakości"]), hit_time_df, ecdf_plot_df), max_iters

def plot_ecdf(ecdf_data, max_iters, title, filename, thresholds_to_plot):
    print(f"  Rysuję wykres: {title}")
    
    plt.figure(figsize=(12, 8))
    
    for th in thresholds_to_plot:
        # Znajdź najbliższy próg w danych
        actual_th = ecdf_data.columns[np.abs(ecdf_data.columns - th).argmin()]
    
        plt.plot(ecdf_data.index, ecdf_data[actual_th], label=rf"Próg $\Delta f \leq$ {actual_th:.1e}")

    plt.title(f"Wykres ECDF - {title}", fontsize=16)
    plt.xlabel("Iteracje (Generacje)", fontsize=12)
    plt.ylabel("Proporcja rozwiązanych uruchomień", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legenda z boku, bo dla Binary jest dużo linii
    plt.tight_layout() # Dopasowanie marginesów
    plt.ylim(-0.05, 1.05)
    plt.xlim(0, max_iters) # Użyj dynamicznego max_iters
    
    plt.savefig(filename, dpi=150)
    print(f" zapisano wykres jako: {filename}")
    plt.close()


# GŁÓWNA LOGIKA SKRYPTU

for report_type in ['binary', 'real']:
    print(f"\nPrzetwarzanie zestawu: {report_type.upper()}")
    
    # 1. Wygeneruj wszystkie 4 tabele
    processed_data, detected_max_iters = generate_ecdf_data(
        PATHS[report_type],
        THRESHOLDS
    )
    
    # 2. Sprawdź, czy dane zostały wczytane poprawnie
    if processed_data is None:
        print(f"Przerwano przetwarzanie zestawu {report_type.upper()} z powodu braku plików")
        continue 

    raw_delta, thresholds_df, hit_time, ecdf_data = processed_data
    
    # 3. Zapisz raport Excel z 4 zakładkami
    output_excel_file = os.path.join(script_dir, f"ECDF_{report_type.upper()}_Report.xlsx")
    print(f" Zapisuję raport Excel: {output_excel_file}")
    
    with pd.ExcelWriter(output_excel_file) as writer:
        raw_delta.to_excel(writer, sheet_name="1. Dane Surowe (Delta)")
        thresholds_df.to_excel(writer, sheet_name="2. Progi Jakosci")
        hit_time.to_excel(writer, sheet_name="3. Czasy Trafien (Iteracje)")
        ecdf_data.to_excel(writer, sheet_name="4. Dane ECDF (Wykres)")

    print(f" pomyślnie zapisano {output_excel_file}")
    
    # 4. Wygeneruj i zapisz wykres
    output_plot_file = os.path.join(script_dir, f"ECDF_{report_type.upper()}_Plot.png")
    
    # Pobierz odpowiednie progi ze słownika PLOT_THRESHOLDS
    my_thresholds = PLOT_THRESHOLDS[report_type]
    
    plot_ecdf(ecdf_data, detected_max_iters, f"Zestaw {report_type.upper()} (F1+F2)", output_plot_file, my_thresholds)


print("\nAnaliza zakończona pomyślnie")