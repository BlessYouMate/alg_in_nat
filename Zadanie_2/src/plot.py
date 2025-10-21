import os
import re
import numpy as np
import matplotlib.pyplot as plt

# Bazowa ścieżka do wyników
base_path = "results"

# Definicja konfiguracji
configs = {
    "Bin": os.path.join(base_path, "bin"),
    "Gray": os.path.join(base_path, "gray")
}

# Oczekiwane wartości n
dimensions = [2, 5, 10]

# Funkcja do wczytywania wyników (dla konkretnego n)
def load_results(path, n):
    runs = []
    if not os.path.exists(path):
        print(f"⚠️ Folder nie istnieje: {path}")
        return runs

    pattern = re.compile(rf"results_n{n}_run\d+\.txt$")
    for file in os.listdir(path):
        if pattern.match(file):
            with open(os.path.join(path, file)) as f:
                values = [float(line.strip()) for line in f if line.strip()]
                runs.append(values)
    return runs

# Uśrednianie długości przebiegów
def pad_and_average(runs):
    if not runs:
        return []
    max_len = max(len(run) for run in runs)
    padded = np.array([np.pad(run, (0, max_len - len(run)), mode='edge') for run in runs])
    return np.mean(padded, axis=0)

# Wczytaj i uśrednij wyniki
mean_results = {}

for method, folder in configs.items():
    for n in dimensions:
        runs = load_results(folder, n)
        if runs:
            mean_results[f"{method}, n={n}"] = pad_and_average(runs)
        else:
            print(f"⚠️ Brak danych dla {method}, n={n}")

plt.figure(figsize=(10,6))
for label, y in mean_results.items():
    plt.plot(y, label=label)

plt.yscale("log")
plt.xlabel("Iteracja")
plt.ylabel("Wartość funkcji F(x)")
plt.title("Porównanie konwergencji (Bin vs Gray dla różnych wymiarów n)")
plt.legend()
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.show()
