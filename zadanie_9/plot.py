import matplotlib.pyplot as plt
import pandas as pd

# Wczytanie danych
try:
    data = pd.read_csv("zdt1_results.csv")
    
    # Rysowanie wykresu
    plt.figure(figsize=(8, 6))
    plt.scatter(data['f1'], data['f2'], c='blue', s=10, label='NSGA-II Solutions')
    
    # Ustawienia osi zgodnie z benchmarkiem ZDT1
    plt.xlabel('$f_1(x)$')
    plt.ylabel('$f_2(x)$')
    plt.title('ZDT1 Pareto Front using NSGA-II (20k Evaluations)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Zapisanie wykresu
    plt.savefig("zdt1_plot.png", dpi=300)
    plt.show()
    print("Wykres wygenerowany jako zdt1_plot.png")
    
except FileNotFoundError:
    print("Nie znaleziono pliku zdt1_results.csv. Uruchom najpierw program C++.")