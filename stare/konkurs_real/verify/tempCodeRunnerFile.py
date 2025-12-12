
        # Przekształcenie z listy punktów na siatkę (matrix) dla wykresu 3D
        # Zakładamy, że dane są posortowane lub regularne (pętle w C++ to gwarantują)
        size = int(np.sqrt(len(data))) # Powinno być 100
        X = data['x'].values.reshape(size, size)
        Y = data['y'].values.reshape(size, size)
        Z = data['z'].values.reshape(size, size)

        # Rysowanie
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        # Powierzchnia z mapą kolorów
        surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.9)
        
        ax.set_title(title, fontsize=15)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('f(x, y)')
        
        # Pasek kolorów
        fig.colorbar(surf, shrink=0.5, aspect=10)
        
        # Lepszy kąt widzenia (możesz zmieniać te wartości)
        ax.view_init(elev=30, azim=-45)

        # Zapis
        output_name = filename.replace("grid_", "verify_").replace(".csv", ".png")
        plt.savefig(output_name, dpi=300)
        print(f"Wygenerowano: {output_name}")
        plt.close()
        
    except FileNotFoundError:
        print(f"Brak pliku: {filename}. Uruchom najpierw kod C++!")