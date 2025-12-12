import os
import numpy as np
import pandas as pd
import glob

# ==========================================
# KONFIGURACJA ŚCIEŻEK (WERSJA KULOODPORNA)
# ==========================================
# Pobieramy katalog, w którym fizycznie leży ten plik skryptu
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 1. Sprawdź, czy results jest obok skryptu
path_option_1 = os.path.join(SCRIPT_DIR, "results")
# 2. Sprawdź, czy results jest piętro wyżej (jeśli skrypt jest w podfolderze)
path_option_2 = os.path.join(os.path.dirname(SCRIPT_DIR), "results")

if os.path.exists(path_option_1):
    BASE_DIR = path_option_1
    OUTPUT_FILE = os.path.join(SCRIPT_DIR, "Wyniki.xlsx")
    print(f"[OK] Znaleziono dane w: {BASE_DIR}")
elif os.path.exists(path_option_2):
    BASE_DIR = path_option_2
    # Zapiszemy wynik piętro wyżej, w głównym folderze projektu
    OUTPUT_FILE = os.path.join(os.path.dirname(SCRIPT_DIR), "Wyniki.xlsx")
    print(f"[OK] Znaleziono dane (piętro wyżej) w: {BASE_DIR}")
else:
    print("\n[BŁĄD KRYTYCZNY] Nie znaleziono folderu 'results' ani tu, ani piętro wyżej!")
    print(f"Szukano w:\n 1. {path_option_1}\n 2. {path_option_2}")
    exit()

print(f"Plik wynikowy trafi do: {OUTPUT_FILE}")
print("---------------------------\n")

# ==========================================
# POZOSTAŁA KONFIGURACJA
# ==========================================
DIMENSIONS = [5, 15, 30]
FUNCTIONS = ["rosenbrock", "salomon", "whitley"]
OPTIMUM = 0.0

# Progi jakości (51 progów, logarytmicznie co 0.2)
exponents = np.arange(-8.0, 2.01, 0.2)
THRESHOLDS = np.power(10.0, exponents)

# Wymóg: "kolumny po 41 liczb"
NUM_POINTS = 41 

# Twoje wybrane progi do wykresów
CHOSEN_THRESHOLDS = {
    # n=5: Algorytm jest świetny, bierzemy bardzo wyśrubowane progi
    5:  [1.58489319246114,
        0.158489319246114,
        0.0630957344480202,
        0.0158489319246113,
        0.00000001],

    # n=15: Widać walkę, bierzemy szeroki zakres od łatwego do trudnego
    15: [6.31e+01, 6.31e+00, 6.31e-01, 6.31e-02, 1.58e-02],

    # n=30: Jest ciężko, bierzemy łatwe progi, żeby pokazać początkowy postęp
    30: [
        100.0, 
        6.30957344480205, 
        0.630957344480203, 
        0.0630957344480202, 
        0.0251188643150961
    ]
}

# ==========================================
# FUNKCJE
# ==========================================

def find_closest_threshold(target, thresholds):
    idx = (np.abs(thresholds - target)).argmin()
    return thresholds[idx]

def load_and_aggregate_dim(dim):
    all_runs = []
    dim_str = f"n{dim}"
    
    print(f" Szukam danych dla n={dim} w podfolderach...")
    
    for func in FUNCTIONS:
        folder = os.path.join(BASE_DIR, func, dim_str)
        if not os.path.exists(folder):
            print(f"   [INFO] Brak folderu: {folder}")
            continue
            
        files = glob.glob(os.path.join(folder, "*.txt"))
        
        if len(files) > 0:
            # print(f"   -> {func}: znaleziono {len(files)} plików")
            for fpath in files:
                try:
                    df = pd.read_csv(fpath, header=None)
                    vals = df[0].values - OPTIMUM
                    vals[vals < 0] = 0.0
                    all_runs.append(vals)
                except: pass
        else:
            print(f"   [WARN] Folder pusty: {folder}")
    
    if not all_runs: return None, 0

    max_len = max(len(r) for r in all_runs)
    n_runs = len(all_runs)
    data_matrix = np.zeros((n_runs, max_len))
    
    for i, run in enumerate(all_runs):
        l = len(run)
        data_matrix[i, :l] = run
        if l < max_len:
            data_matrix[i, l:] = run[-1]
            
    return data_matrix, max_len

def calculate_resampled_ecdf(data_matrix, max_len):
    n_runs = data_matrix.shape[0]
    indices = np.linspace(0, max_len - 1, NUM_POINTS).astype(int)
    budget_axis = indices 
    ecdf_matrix = np.zeros((NUM_POINTS, len(THRESHOLDS)))
    data_T = data_matrix.T 
    
    for i, time_idx in enumerate(indices):
        current_state = data_T[time_idx, :]
        for j, th in enumerate(THRESHOLDS):
            success_count = np.sum(current_state <= th)
            ecdf_matrix[i, j] = success_count / n_runs
            
    return ecdf_matrix, budget_axis

# ==========================================
# MAIN
# ==========================================

def main():
    print(f"Generowanie pliku '{OUTPUT_FILE}'...")
    
    try:
        writer = pd.ExcelWriter(OUTPUT_FILE, engine='xlsxwriter')
        workbook = writer.book
    except ImportError:
        print("[BŁĄD] Brak biblioteki xlsxwriter. Zainstaluj: pip install xlsxwriter")
        return

    for dim in DIMENSIONS:
        print(f"\nPrzetwarzanie n={dim}...")
        
        data_matrix, max_len = load_and_aggregate_dim(dim)
        if data_matrix is None:
            print(f" [SKIP] Nie udało się wczytać żadnych danych dla n={dim}")
            continue
            
        ecdf_vals, budget_axis = calculate_resampled_ecdf(data_matrix, max_len)
        
        df = pd.DataFrame(ecdf_vals, columns=THRESHOLDS)
        df.insert(0, "Budżet", budget_axis)
        
        sheet_name = str(dim)
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # --- WYKRESY W EXCELU ---
        worksheet = writer.sheets[sheet_name]
        chart = workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
        
        chosen = CHOSEN_THRESHOLDS[dim]
        print(f"  Dodawanie wykresu do Excela dla progów: {chosen}")
        
        for target in chosen:
            real_th = find_closest_threshold(target, THRESHOLDS)
            try:
                col_idx = df.columns.get_loc(real_th)
                chart.add_series({
                    'name':   f'Próg: {real_th:.2e}',
                    'categories': [sheet_name, 1, 0, NUM_POINTS, 0],
                    'values':     [sheet_name, 1, col_idx, NUM_POINTS, col_idx],
                })
            except KeyError:
                pass

        # --- TYTUŁY I OPISY OSI (ZMIENIONE WG ŻYCZENIA) ---
        chart.set_title ({'name': f'Krzywe ECDF (n={dim})'})
        
        chart.set_x_axis({
            'name': 'Liczba Generacji',  # Zamiast "Liczba wywołań"
            'min': 0, 
            'max': max_len
        })
        
        chart.set_y_axis({
            'name': 'Frakcja udanych przebiegów', # Zamiast "Prawdopodobieństwo"
            'min': 0, 
            'max': 1.0,
            'major_unit': 0.1  # Podziałka co 0.1 (10%)
        })
        
        chart.set_size({'width': 720, 'height': 480})
        
        # Wstaw wykres obok tabeli
        worksheet.insert_chart('BA2', chart)

    writer.close()
    print(f"\nGotowe! Otwórz plik 'Wyniki.xlsx'.")

if __name__ == "__main__":
    main()