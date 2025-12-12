#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- TWOJE IMPLEMENTACJE FUNKCJI (Skopiowane z Twojego kodu) ---

double f_rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
    // Dla N=2 pętla wykona się raz (i=0)
    for (size_t i = 0; i < x.size() - 1; ++i) {
        double t1 = (x[i + 1] - x[i] * x[i]);
        double t2 = (x[i] - 1.0);
        sum += 100.0 * t1 * t1 + t2 * t2;
    }
    return sum;
}

double f_salomon(const std::vector<double>& x) {
    double sum_sq = 0.0;
    for (double val : x) sum_sq += val * val;
    double norm = std::sqrt(sum_sq);
    return 1.0 - std::cos(2.0 * M_PI * norm) + 0.1 * norm;
}

double f_whitley(const std::vector<double>& x) {
    int D = x.size();
    double sum = 0.0;
    for (int j = 0; j < D; ++j) {
        for (int k = 0; k < D; ++k) {
            double temp = 100.0 * std::pow((x[k] - x[j] * x[j]), 2) + std::pow(1.0 - x[j], 2);
            sum += (temp * temp / 4000.0) - std::cos(temp) + 1.0;
        }
    }
    return sum;
}

// --- GENERATOR SIATKI ---

void generate_grid_file(std::string filename, double (*func)(const std::vector<double>&), double low, double high) {
    std::ofstream out(filename);
    // Nagłówek CSV: x, y, z
    out << "x,y,z\n";
    
    int steps = 100; // Wymóg: siatka 100x100
    double step_size = (high - low) / (steps - 1);

    std::cout << "Generowanie: " << filename << " [" << low << ", " << high << "]...\n";

    for (int i = 0; i < steps; ++i) {
        double x = low + i * step_size;
        for (int j = 0; j < steps; ++j) {
            double y = low + j * step_size;
            
            // Obliczamy wartość funkcji dla punktu (x, y)
            std::vector<double> point = {x, y};
            double z = func(point);
            
            out << x << "," << y << "," << z << "\n";
        }
    }
    out.close();
    std::cout << " -> Zapisano.\n";
}

int main() {
    // 1. Rosenbrock
    // Do wykresu bierzemy zakres [-2, 2] (zgodnie z Rys. A.3 w PDF)
    // Zamiast [-30, 30]
    generate_grid_file("grid_rosenbrock.csv", f_rosenbrock, -2.0, 2.0);

    // 2. Salomon
    // Do wykresu bierzemy zakres [-2, 2] (zgodnie z Rys. A.10 w PDF)
    // Zamiast [-100, 100]
    generate_grid_file("grid_salomon.csv", f_salomon, -2.0, 2.0);

    // 3. Whitley
    // Tutaj akurat PDF (Rys A.11) pokazuje zakres ok [-10, 10], 
    // wiec zostawiamy oryginalny, bo wygladal dobrze.
    generate_grid_file("grid_whitley.csv", f_whitley, -10.24, 10.24);

    return 0;
}