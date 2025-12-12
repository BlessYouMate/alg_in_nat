import os
import numpy as np
import pandas as pd
import glob

# ==========================================
# KONFIGURACJA ŚCIEŻEK (WERSJA KULOODPORNA)
# ==========================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 1. Sprawdź, czy results jest obok skryptu
path_option_1 = os.path.join(SCRIPT_DIR, "results")
# 2. Sprawdź, czy results jest piętro wyżej
path_option_2 = os.path.join(os.path.dirname(SCRIPT_DIR), "results")

if os.path.exists(path_option_1):
    BASE_DIR = path_option_1
    OUTPUT_FILE = os.path.join(SCRIPT_DIR, "Wyniki.xlsx")
    print(f"[OK] Znaleziono dane w: {BASE_DIR}")
elif os.path.exists(path_option_2):
    BASE_DIR = path_option_2
    OUTPUT_FILE = os.path.join(os.path.dirname(SCRIPT_DIR), "Wyniki.xlsx")
    print(f"[OK] Znaleziono dane (piętro wyżej) w: {BASE_DIR}")
else:
    print("\n[BŁĄD KRYTYCZNY] Nie znaleziono folderu 'results' ani tu, ani piętro wyżej!")
    exit()

print(f"Plik wynikowy trafi do: {OUTPUT_FILE}")
print("---------------------------\n")

# ==========================================
# POZOSTAŁA KONFIGURACJA
# ==========================================
DIMENSIONS = [5, 15, 30]
FUNCTIONS = ["rosenbrock", "salomon", "whitley"]
OPTIMUM = 0.0

# Progi jakości
exponents = np.arange(-8.0, 2.01, 0.2)
THRESHOLDS = np.power(10.0, exponents)

# Wymóg: "kolumny po 41 liczb"
NUM_POINTS = 41 

# Twoje wybrane progi do wykresów
CHOSEN_THRESHOLDS = {
    # n=5: Algorytm jest świetny, bierzemy bardzo wyśrubowane progi
    5:  [1.58489319246114,
        0.100000000000001,
        0.0630957344480202,
        0.0158489319246113,
        0.00000001],

    # n=15: Widać walkę, bierzemy szeroki zakres
    15: [15.8489319246114, 
         0.630957344480203, 
         0.0630957344480202, 
         0.0251188643150961, 
         0.0158489319246113],

    # n=30: Jest ciężko, bierzemy łatwe progi
    30: [
        63.0957344480205,
        0.630957344480203, 
        0.398107170553503, 
        0.0398107170553503, 
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
    """
    Zwraca:
    1. macierz numpy (do obliczeń ECDF)
    2. max_len (długość osi czasu)
    3. DataFrame (dane surowe do zapisu)
    """
    all_runs = []
    all_names = [] # Lista nazw kolumn dla danych surowych
    dim_str = f"n{dim}"
    
    print(f" Szukam danych dla n={dim} w podfolderach...")
    
    for func in FUNCTIONS:
        folder = os.path.join(BASE_DIR, func, dim_str)
        if not os.path.exists(folder):
            print(f"   [INFO] Brak folderu: {folder}")
            continue
            
        files = glob.glob(os.path.join(folder, "*.txt"))
        
        if len(files) > 0:
            for fpath in files:
                try:
                    df = pd.read_csv(fpath, header=None)
                    vals = df[0].values - OPTIMUM
                    vals[vals < 0] = 0.0
                    all_runs.append(vals)
                    
                    # Tworzymy nazwę kolumny: np. "rosenbrock_run_1"
                    run_name = os.path.basename(fpath).replace(".txt", "")
                    all_names.append(f"{func}_{run_name}")
                except: pass
        else:
            print(f"   [WARN] Folder pusty: {folder}")
    
    if not all_runs: return None, 0, None

    max_len = max(len(r) for r in all_runs)
    n_runs = len(all_runs)
    
    # Macierz do obliczeń ECDF
    data_matrix = np.zeros((n_runs, max_len))
    
    # Słownik do stworzenia DataFrame z danymi surowymi
    raw_data_dict = {}
    
    for i, run in enumerate(all_runs):
        l = len(run)
        
        # 1. Wypełnianie macierzy numpy (z paddingiem ostatnią wartością)
        data_matrix[i, :l] = run
        if l < max_len:
            data_matrix[i, l:] = run[-1]
            
        # 2. Wypełnianie słownika dla danych surowych (też z paddingiem, żeby Excel był równy)
        # Jeśli wolisz puste komórki w Excelu zamiast powtórzeń, można tu użyć np. np.nan
        padded_run = data_matrix[i, :] 
        raw_data_dict[all_names[i]] = padded_run
            
    # Tworzenie DataFrame z surowymi danymi
    raw_df = pd.DataFrame(raw_data_dict)
    
    return data_matrix, max_len, raw_df

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
        
        # Pobieramy teraz także raw_df
        data_matrix, max_len, raw_df = load_and_aggregate_dim(dim)
        
        if data_matrix is None:
            print(f" [SKIP] Nie udało się wczytać żadnych danych dla n={dim}")
            continue
            
        ecdf_vals, budget_axis = calculate_resampled_ecdf(data_matrix, max_len)
        
        df = pd.DataFrame(ecdf_vals, columns=THRESHOLDS)
        df.insert(0, "Budżet", budget_axis)
        
        # --- ZAPIS ARKUSZA ECDF (wymagany) ---
        sheet_name = str(dim)
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # --- ZAPIS ARKUSZA RAW (dodatkowy, z danymi surowymi) ---
        sheet_raw = f"{dim}_RAW"
        print(f"  Zapisywanie danych surowych do arkusza '{sheet_raw}'...")
        raw_df.to_excel(writer, sheet_name=sheet_raw, index_label="Generacja")
        
        # --- WYKRESY W EXCELU (w arkuszu ECDF) ---
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

        # Tytuły i opisy osi
        chart.set_title ({'name': f'Krzywe ECDF (n={dim})'})
        
        chart.set_x_axis({
            'name': 'Liczba Generacji', 
            'min': 0, 
            'max': max_len
        })
        
        chart.set_y_axis({
            'name': 'Frakcja udanych przebiegów', 
            'min': 0, 
            'max': 1.0,
            'major_unit': 0.1 
        })
        
        chart.set_size({'width': 720, 'height': 480})
        
        # Wstaw wykres obok tabeli
        worksheet.insert_chart('BA2', chart)

    writer.close()
    print(f"\nGotowe! Otwórz plik 'Wyniki.xlsx'.")

if __name__ == "__main__":
    main()