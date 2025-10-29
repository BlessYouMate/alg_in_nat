import os
import numpy as np
import matplotlib.pyplot as plt

# Ustawienia
base_dir = "results"
functions = ["f1", "f2"]
representations = ["real", "bin", "gray"]
output_dir = "plots"

os.makedirs(output_dir, exist_ok=True)

for func in functions:
    for rep in representations:
        folder_path = os.path.join(base_dir, func, rep)
        
        # Wczytanie wszystkich plików
        all_runs = []
        for file_name in sorted(os.listdir(folder_path)):
            if file_name.endswith(".txt"):
                file_path = os.path.join(folder_path, file_name)
                with open(file_path, "r") as f:
                    values = [float(line.strip()) for line in f if line.strip()]
                    all_runs.append(values)

        if not all_runs:
            print(f"Brak danych w: {folder_path}")
            continue

        min_len = min(len(run) for run in all_runs)
        all_runs = [run[:min_len] for run in all_runs]
        avg_values = np.mean(all_runs, axis=0)

        # Tworzenie wykresu
        plt.figure(figsize=(8, 5))
        plt.plot(range(min_len), avg_values, label=f"{func.upper()} - {rep.capitalize()}", linewidth=2)
        plt.xlabel("Iteracja")
        plt.ylabel("Wartość funkcji F(x)")
        plt.title(f"{func.upper()} ({rep.capitalize()} representation)")

        # Skala i zakres w zależności od funkcji
        if func == "f1":
            plt.ylim([-7, max(avg_values)*1.1])  # zakres od -7 do max
            plt.grid(True, linestyle="--", alpha=0.6)
        elif func == "f2":
            plt.axhline(0, color='red', linestyle='--', linewidth=1.5, label="Minimum = 0")
            plt.grid(True, linestyle="--", alpha=0.6)
        
        plt.legend()
        output_path = os.path.join(output_dir, f"{func}_{rep}.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Wykres zapisany: {output_path}")
