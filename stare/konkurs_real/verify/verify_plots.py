import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os

# Lista plików do przetworzenia
FILES = [
    ("grid_rosenbrock.csv", "Rosenbrock (n=2)"),
    ("grid_salomon.csv", "Salomon (n=2)"),
    ("grid_whitley.csv", "Whitley (n=2)")
]

print("Rozpoczynam generowanie wykresów 3D...")

for filename, title in FILES:
    # Sprawdzenie czy plik istnieje
    if not os.path.exists(filename):
        print(f"[BŁĄD] Nie znaleziono pliku: {filename}. Uruchom najpierw kod C++.")
        continue

    try:
        # Wczytanie danych
        data = pd.read_csv(filename)
        
        # Weryfikacja czy mamy dane
        if data.empty:
            print(f"[BŁĄD] Plik {filename} jest pusty.")
            continue

        # Przekształcenie z listy punktów na siatkę (matrix) dla wykresu 3D
        # Oczekujemy 10000 punktów dla siatki 100x100
        size = int(np.sqrt(len(data))) 
        
        X = data['x'].values.reshape(size, size)
        Y = data['y'].values.reshape(size, size)
        Z = data['z'].values.reshape(size, size)

        # Rysowanie
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        # Powierzchnia z mapą kolorów viridis
        surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.9)
        
        ax.set_title(title, fontsize=15)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('f(x, y)')
        
        # Pasek kolorów
        fig.colorbar(surf, shrink=0.5, aspect=10)
        
        # Ustawienie widoku (kąt kamery)
        ax.view_init(elev=30, azim=-45)

        # Zapis do pliku PNG
        output_name = filename.replace("grid_", "verify_").replace(".csv", ".png")
        plt.savefig(output_name, dpi=300)
        print(f" -> Wygenerowano: {output_name}")
        
        plt.close() # Zamknięcie wykresu, żeby zwolnić pamięć
        
    except Exception as e:
        print(f"[BŁĄD] Problem z plikiem {filename}: {e}")

print("Gotowe.")