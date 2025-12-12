import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# --- KONFIGURACJA ---
BASE_DIR = "results"      
OUTPUT_DIR = "plots"      

FUNCTIONS = ["rosenbrock", "salomon", "whitley"]
DIMENSIONS = ["n5", "n15", "n30"]

def load_data_from_folder(folder_path):
    all_runs = []
    if not os.path.exists(folder_path):
        return None
    
    files = [f for f in os.listdir(folder_path) if f.endswith(".txt")]
    if not files:
        return None
    
    # Sortowanie numeryczne
    try:
        files.sort(key=lambda x: int(x.split('_')[1].split('.')[0]) if '_' in x else x)
    except:
        pass

    for file_name in files:
        file_path = os.path.join(folder_path, file_name)
        try:
            with open(file_path, "r") as f:
                # Wczytujemy dane
                values = [float(line.strip()) for line in f if line.strip()]
                if values:
                    all_runs.append(values)
        except Exception:
            pass
    return all_runs

def pad_data(all_runs):
    """Wyrównuje długość przebiegów do najdłuższego."""
    if not all_runs:
        return []
    max_len = max(len(run) for run in all_runs)
    padded_runs = []
    for run in all_runs:
        if len(run) < max_len:
            last_val = run[-1]
            padding = [last_val] * (max_len - len(run))
            padded_runs.append(run + padding)
        else:
            padded_runs.append(run)
    return np.array(padded_runs)

def process_and_plot():
    # Tworzenie głównego folderu wyjściowego
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Generowanie wykresów...")

    for func in FUNCTIONS:
        for dim in DIMENSIONS:
            input_path = os.path.join(BASE_DIR, func, dim)
            all_runs_list = load_data_from_folder(input_path)
            
            if not all_runs_list:
                # print(f"Brak danych dla: {func} {dim}")
                continue

            # Przygotowanie danych
            data_matrix = pad_data(all_runs_list)
            avg_values = np.mean(data_matrix, axis=0)

            # --- RYSOWANIE ---
            fig, ax = plt.subplots(figsize=(10, 6))
            
            x_axis = range(len(avg_values))
            
            # Rysujemy samą linię średnią
            ax.plot(x_axis, avg_values, label=f"Średnia ({len(all_runs_list)} przebiegów)", linewidth=2, color='#0044cc')

            # --- KONFIGURACJA WYGLĄDU ---
            
            # 1. Skala logarytmiczna na Y
            ax.set_yscale('log')
            
            # 2. Gęsta siatka (Grid)
            # 'which="both"' włącza linie dla głównych (10^1, 10^2) i pośrednich wartości
            ax.grid(True, which="major", linestyle="-", linewidth=0.8, color='gray', alpha=0.6)
            ax.grid(True, which="minor", linestyle=":", linewidth=0.5, color='gray', alpha=0.4)

            # 3. Tytuły i osie
            ax.set_title(f"Zbieżność: {func.capitalize()} - {dim}", fontsize=14)
            ax.set_xlabel("Generacja", fontsize=12)
            ax.set_ylabel("Wartość funkcji (Fitness) [log]", fontsize=12)
            
            # 4. Legenda
            ax.legend(loc="upper right", fontsize=10, frameon=True, framealpha=0.9)

            # 5. Dopasowanie zakresów (usuwanie pustej przestrzeni)
            # Ustawiamy dolny limit ciut niżej niż najniższa osiągnięta wartość,
            # żeby wykres nie "szorował" po dolnej krawędzi ramki.
            min_val = np.min(avg_values)
            max_val = np.max(avg_values)
            
            # Margines 10% w dół i w górę na skali log
            if min_val > 0:
                ax.set_ylim(bottom=min_val * 0.5, top=max_val * 2.0)
            
            # Opcjonalnie: Jeśli góra jest absurdalnie wysoka (np. 10^14), przytnij
            # (Odkomentuj jeśli nadal wykresy są "pionową kreską")
            # if max_val > 1e6: ax.set_ylim(top=1e6)

            # --- ZAPIS 1: W podfolderze ---
            sub_folder = os.path.join(OUTPUT_DIR, func, dim)
            os.makedirs(sub_folder, exist_ok=True)
            path_1 = os.path.join(sub_folder, "plot.png")
            plt.savefig(path_1, dpi=300, bbox_inches='tight')

            # --- ZAPIS 2: W głównym folderze ---
            filename_root = f"{func}_{dim}.png"
            path_2 = os.path.join(OUTPUT_DIR, filename_root)
            plt.savefig(path_2, dpi=300, bbox_inches='tight')

            plt.close() 
            print(f"  -> Zapisano: {filename_root}")

    print("Gotowe.")

if __name__ == "__main__":
    process_and_plot()