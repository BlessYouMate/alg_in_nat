#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using Individual = std::vector<double>;

// ENUM żeby wygodnie wybierać funkcję
enum class Benchmark {
    Rastrigin,
    Schwefel,
    Griewank,
    Levy,
    Eggholder,   // 2D
    DropWave,    // 2D
    Rosenbrock,
    Ackley
};

// --- Implementacje funkcji benchmarkowych ---

double rastrigin(const Individual& x) {
    int n = (int)x.size();
    double A = 10.0;
    double sum = A * n;
    for (double xi : x) sum += (xi * xi - A * std::cos(2.0 * M_PI * xi));
    return sum;
}

double schwefel_fun(const Individual& x) {
    int n = (int)x.size();
    double sum = 0.0;
    for (double xi : x) sum += xi * std::sin(std::sqrt(std::fabs(xi)));
    return 418.9829 * n - sum;
}

double griewank(const Individual& x) {
    int n = (int)x.size();
    double sum = 0.0;
    double prod = 1.0;
    for (int i = 0; i < n; ++i) {
        sum += (x[i] * x[i]) / 4000.0;
        prod *= std::cos(x[i] / std::sqrt((double)i + 1.0));
    }
    return sum - prod + 1.0;
}

double levy(const Individual& x) {
    int n = (int)x.size();
    auto w = [&](int i) { return 1.0 + (x[i] - 1.0) / 4.0; };
    double term1 = std::pow(std::sin(M_PI * w(0)), 2.0);
    double term3 = std::pow(w(n - 1) - 1.0, 2.0) * (1.0 + std::pow(std::sin(2.0 * M_PI * w(n - 1)), 2.0));
    double sum = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double wi = w(i);
        double wi1 = w(i + 1);
        sum += std::pow(wi - 1.0, 2.0) * (1.0 + 10.0 * std::pow(std::sin(M_PI * wi1), 2.0));
    }
    return term1 + sum + term3;
}

// Eggholder (specjalnie 2D)
double eggholder(const Individual& x) {
    if (x.size() < 2) return 1e9;
    double xi = x[0], yi = x[1];
    double a = -(yi + 47.0) * std::sin(std::sqrt(std::fabs(xi / 2.0 + yi + 47.0)));
    double b = -xi * std::sin(std::sqrt(std::fabs(xi - (yi + 47.0))));
    return a + b;
}

// Drop-Wave (2D)
double dropwave(const Individual& x) {
    if (x.size() < 2) return 1e9;
    double xi = x[0], yi = x[1];
    double r = std::sqrt(xi * xi + yi * yi);
    double num = 1.0 + std::cos(12.0 * r);
    double den = 0.5 * (r * r) + 2.0;
    return -num / den;
}

double rosenbrock(const Individual& x) {
    int n = (int)x.size();
    double sum = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double t1 = (x[i + 1] - x[i] * x[i]);
        double t2 = (x[i] - 1.0);
        sum += 100.0 * t1 * t1 + t2 * t2;
    }
    return sum;
}

double ackley_bench(const Individual& x) {
    int d = (int)x.size();
    double a = 20.0, b = 0.2, c = 2.0 * M_PI;
    double sumsq = 0.0, sumcos = 0.0;
    for (double xi : x) { sumsq += xi * xi; sumcos += std::cos(c * xi); }
    return -a * std::exp(-b * std::sqrt(sumsq / d)) - std::exp(sumcos / d) + a + std::exp(1.0);
}

// --- helper: zwraca wskaźnik do funkcji wg enumu ---
double (*getBenchmarkFunction(Benchmark b))(const Individual&) {
    switch (b) {
    case Benchmark::Rastrigin: return rastrigin;
    case Benchmark::Schwefel:  return schwefel_fun;
    case Benchmark::Griewank:  return griewank;
    case Benchmark::Levy:      return levy;
    case Benchmark::Eggholder: return eggholder;
    case Benchmark::DropWave:  return dropwave;
    case Benchmark::Rosenbrock:return rosenbrock;
    case Benchmark::Ackley:    return ackley_bench;
    default: return rastrigin;
    }
}

// --- helper: sugeruje domenę i wymiar dla danego benchmarku ---
void getBenchmarkDomain(Benchmark b, int& dim, double& low, double& high) {
    switch (b) {
    case Benchmark::Rastrigin:
        dim = 10; low = -5.12; high = 5.12; break;
    case Benchmark::Schwefel:
        dim = 10; low = -500; high = 500; break;
    case Benchmark::Griewank:
        dim = 10; low = -600; high = 600; break;
    case Benchmark::Levy:
        dim = 10; low = -10; high = 10; break;
    case Benchmark::Eggholder:
        dim = 2;  low = -512; high = 512; break;
    case Benchmark::DropWave:
        dim = 2;  low = -5.12; high = 5.12; break;
    case Benchmark::Rosenbrock:
        dim = 10; low = -5.0; high = 10.0; break;
    case Benchmark::Ackley:
        dim = 10; low = -32.768; high = 32.768; break;
    default:
        dim = 10; low = -5.12; high = 5.12; break;
    }
}
