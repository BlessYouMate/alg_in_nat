#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>
#include <limits>
#include <iomanip>
#include <filesystem> // Wymaga C++17
#include <fstream>

/*
Strategia Ewolucyjna 
zaimplementowane rozwiazania z wykladu:
    - selekcja (lambda, mu) - wymienamy cale stare pokolenie na nowe
    - mutacja przy uzyciu wektora sigm dla kazdego wymiaru
    - rekombinacja wazona - centroid
    - przycinanie wartosci wychodzacych poza dziedzine do min/max
    - tau i tau_prime - Wsplczynniki uczenia dla samo-adaptacji sigmy. Wyznaczane wg wzorow Schwefela.
dodatkowe rozwiazania spoza wykladu:
    - Mechanizm pedu - Obliczamy wektor momentum jako wykladnicza srednia kroczaca przesuniec centroidy
      i nastepnie dodajemy go do mutacji. Dzieki temu kumulujemy informacje o kierunku,
      w ktorym porusza sie populacja, co pozwala na szybsze przebywanie plaskich obszarow i wydluzanie kroku w dobrym kierunku.
      Jest to uproszczona implementacja koncepcji Evolution Path z algorytmu CMA-ES (Covariance Matrix Adaptation Evolution Strategy).
      Zrodlo: Hansen, N., & Ostermeier, A. (2001). Completely Derandomized Self-Adaptation in Evolution Strategies. Evolutionary Computation, 9(2), 159-195.

    - Strategia Restartow IPOP (Increasing Population) - Gdy algorytm utknie, nastepuje restart ze zwiekszona populacja. 
      Zakladamy, ze jesli mala populacja zawiodla, to nalezy zrestartowac algorytm z wieksza populacja, aby zwiekszyc szanse na globalne zbieznosc.
      Zrodlo: Auger, A., & Hansen, N. (2005). A restart CMA evolution strategy with increasing population size. 2005 IEEE Congress on Evolutionary Computation.

    - Zatrzymanie wzrostu lambdy gdy stanie sie zbyt duza - pomysl autorksi
*/

namespace fs = std::filesystem;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ================================================================
// KONFIGURACJA
// ================================================================
struct Config {
    //wlasnosci minimalizowanej funckji:
    int dim = 30; //liczba wymiarow
    double low = -10.24; //dolne ograniczenie dziedziny
    double high = 10.24; //gorne ograniczenie dziedziny
    long long T_max = 0; //liczba kredytow - wartosc zmieniana w funckji main
    //wlasnosci algorytmu:
    int lambda; //liczba potomkow
    int max_lambda = 1000;
    int mu; //liczba rodzicow
    double startSigma; //startowa wartosc sigma dla wszystkich wymiarow
    double sigma_min = 1e-12;
    double sigma_max = 5.0;
    double momentum_beta = 0.6; // Sila pedu
    int stagnationLimit = 100; //limit iteracji bez poprawy po ktorej wywolujemy restart
    double tau;
    double tau_prime;
};

struct Individual {
    std::vector<double> x;
    std::vector<double> sigma;
    double fitness;
};

std::mt19937 rng(std::random_device{}());

double randUniform(double a, double b) {
    if (a >= b) return a;
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}
double randNormal(double mean, double stddev) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(rng);
}

//obcinanie wartosci spoza dziedziny do min/max
void enforceConstraints(Individual& ind, double low, double high) {
    for (size_t i = 0; i < ind.x.size(); ++i) {
        if (!std::isfinite(ind.x[i])) { ind.x[i] = randUniform(low, high); continue; }
        if (ind.x[i] < low) ind.x[i] = low;
        if (ind.x[i] > high) ind.x[i] = high;
    }
}

//ewaluacja i dekerementowanie kredytow
double evaluate(Individual& ind, double (*func)(const std::vector<double>&), long long& credits) {
    if (credits > 0) {
        ind.fitness = func(ind.x);
        credits--;
    }
    else {
        ind.fitness = std::numeric_limits<double>::infinity();
    }
    return ind.fitness;
}

// ================================================================
// OPERATORY
// ================================================================

Individual calculateCentroid(const std::vector<Individual>& parents, const std::vector<double>& weights, int dim, int mu) {
    Individual c;
    c.x.assign(dim, 0.0);
    c.sigma.assign(dim, 0.0);

    for (int i = 0; i < mu; ++i) {
        for (int d = 0; d < dim; ++d) {
            c.x[d] += weights[i] * parents[i].x[d];
            c.sigma[d] += weights[i] * parents[i].sigma[d];
        }
    }
    return c;
}

void mutate_with_momentum(Individual& ind, const std::vector<double>& momentum, const Config& cfg) {
    double globalNoise = randNormal(0.0, 1.0);

    for (int d = 0; d < cfg.dim; ++d) {
        double localNoise = randNormal(0.0, 1.0);

        // Adaptacja Sigmy z learning rate - inspirowane wykladem
        ind.sigma[d] = ind.sigma[d] * std::exp(cfg.tau_prime * globalNoise + cfg.tau * localNoise);
        if (ind.sigma[d] < cfg.sigma_min) ind.sigma[d] = cfg.sigma_min;
        if (ind.sigma[d] > cfg.sigma_max) ind.sigma[d] = cfg.sigma_max;

        double mean_shift = cfg.momentum_beta * momentum[d];
        ind.x[d] += mean_shift + ind.sigma[d] * randNormal(0.0, 1.0);
    }
    enforceConstraints(ind, cfg.low, cfg.high);
}

//funckja wykorzystywana jedynie do testow
void analyzeResult(const std::vector<double>& bestX) {
    std::cout << "\n=== FINAL ANALYSIS ===\n";
    double distToOpt = 0.0;
    int nearCount = 0;
    std::cout << "Coords (5): [ ";
    for (size_t i = 0; i < bestX.size(); ++i) {
        if (i < 5) std::cout << std::fixed << std::setprecision(2) << bestX[i] << " ";
        double diff = bestX[i] - 1.0;
        distToOpt += diff * diff;
        if (std::abs(diff) < 0.5) nearCount++;
    }
    std::cout << "... ]\n";
    std::cout << "Dist to [1...]: " << std::sqrt(distToOpt) << "\n";
    std::cout << "Near Target: " << nearCount << "/" << bestX.size() << "\n";
}

// ================================================================
// ALGORYTM
// ================================================================
std::vector<double> run_Momentum_ES(double (*objective)(const std::vector<double>&), Config cfg) {
    std::vector<double> history;
    long long credits = cfg.T_max;

    double bestFitnessEver = std::numeric_limits<double>::infinity();
    std::vector<double> bestXEver(cfg.dim);

    int restartCount = 0;

    while (credits > 0) {
        // Restart IPOP - w tej czesci kodu jestesmy tylko na poczatku petli i przy restartach
        if (cfg.lambda >= cfg.max_lambda) {
            cfg.lambda = cfg.lambda + int(4 * std::log((double)cfg.dim));
        }
        else {
            cfg.lambda = (int)(cfg.lambda * 1.5);
        }

        // Selekcja (mu, lambda)
        cfg.mu = std::max(2, cfg.lambda / 2);

        std::vector<double> weights(cfg.mu);
        double sum_w = 0.0;
        for (int i = 0; i < cfg.mu; ++i) { weights[i] = std::log(cfg.mu + 0.5) - std::log(i + 1.0); sum_w += weights[i]; }
        for (int i = 0; i < cfg.mu; ++i) weights[i] /= sum_w;

        //Inicjalizacja
        std::vector<Individual> population(cfg.mu);

        for (int i = 0; i < cfg.mu; ++i) {
            population[i].x.resize(cfg.dim);
            population[i].sigma.resize(cfg.dim);
            for (int d = 0; d < cfg.dim; ++d) {
                population[i].x[d] = randUniform(cfg.low, cfg.high);

                population[i].sigma[d] = cfg.startSigma;
            }

            evaluate(population[i], objective, credits);

            if (population[i].fitness < bestFitnessEver) {
                bestFitnessEver = population[i].fitness;
                bestXEver = population[i].x;
            }
            if (credits <= 0) break;
        }
        if (credits <= 0) break;

        //logging do testow
        //std::cout << ">> RESTART " << restartCount << " (Lam: " << cfg.lambda << ") Edge-Start. Cr: " << credits << "\n";

        // Inicjalizacja Pedu
        std::vector<double> momentum(cfg.dim, 0.0);
        Individual prevCentroid = calculateCentroid(population, weights, cfg.dim, cfg.mu);

        int gen = 0;
        int stagnationCounter = 0;
        double lastIterBest = population[0].fitness;

        while (credits > 0) {
            gen++;

            Individual currentCentroid = calculateCentroid(population, weights, cfg.dim, cfg.mu);

            // Aktualizacja wektora pedu - zgodnie z kierunkiem ruchu populacji
            for (int d = 0; d < cfg.dim; ++d) {
                double move = currentCentroid.x[d] - prevCentroid.x[d];
                momentum[d] = (1.0 - cfg.momentum_beta) * momentum[d] + cfg.momentum_beta * move;
            }
            prevCentroid = currentCentroid;

            //GENEROWANIE POTOMSTWA
            std::vector<Individual> offspring;
            offspring.reserve(cfg.lambda);

            for (int k = 0; k < cfg.lambda; ++k) {
                Individual child = currentCentroid;
                mutate_with_momentum(child, momentum, cfg);
                evaluate(child, objective, credits);
                offspring.push_back(child);
                if (credits <= 0) break;
            }
            if (credits <= 0) break;

            //SELEKCJA (lambda, mu) - zastepujemy cale pokolenie
            std::sort(offspring.begin(), offspring.end(), [](const auto& a, const auto& b) { return a.fitness < b.fitness; });

            for (int i = 0; i < cfg.mu; ++i) population[i] = offspring[i];

            // Aktualizacja rekordu
            if (population[0].fitness < bestFitnessEver) {
                bestFitnessEver = population[0].fitness;
                bestXEver = population[0].x;
                stagnationCounter = 0;
            }

            history.push_back(bestFitnessEver);

            // Statystyki do przerwania
            double avgSigma = 0.0;
            for (const auto& ind : population) for (double s : ind.sigma) avgSigma += s;
            avgSigma /= (cfg.dim * cfg.mu);
            
            if (avgSigma < 1e-4) break;// zbyt mala srednia sigma- utkenlismy

            //Zliczanie Stagnacji
            if (std::abs(lastIterBest - population[0].fitness) < 1e-9) stagnationCounter++;
            else stagnationCounter = 0;
            lastIterBest = population[0].fitness;

            if (stagnationCounter > cfg.stagnationLimit) break;
        }
        restartCount++;
    }

    //loging do testow
    //std::cout << ">>> FINAL: " << std::scientific << bestFitnessEver << "\n";
    //analyzeResult(bestXEver);
    return history;
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

double f_rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
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

int main() {
    std::vector<int> dims = { 5, 15, 30 };
    int runs = 100;
    for (int n : dims) {
        Config cfg;
        cfg.dim = n;
        cfg.T_max = 10000 * n;
        cfg.tau = 1.0 / std::sqrt(2.0 * std::sqrt((double)cfg.dim));
        cfg.tau_prime = 1.0 / std::sqrt(2.0 * (double)cfg.dim);
        
        cfg.startSigma = 2.5;
        
        if (n == 30) {
            cfg.lambda = 60;
            cfg.momentum_beta = 0.9;
        }
        else if (n == 5) {
            cfg.lambda = 8;
            cfg.momentum_beta = 0.6;
        }
        else {
            cfg.lambda = 20;
            cfg.momentum_beta = 0.75;
        }
        std::cout << ">>> Processing Whitley " << n << "D...\n";

        std::string dir = "results/whitley/n" + std::to_string(n);
        fs::create_directories(dir);
        int counter = 0;
        for (int r = 1; r <= runs; ++r) {
            auto history = run_Momentum_ES(f_whitley, cfg);
            std::ofstream out(dir + "/run_" + std::to_string(r) + ".txt");
            for (double val : history) out << val << "\n";
            if (r % 10 == 0) std::cout << ".";
            if (history.back() > -1 && history.back() < 1) counter++;
        }
        std::cout << "number of good scores: "<<counter<<"\n";
        std::cout << " Done.\n";

        std::cout << ">>> Processing Rosenbrock " << n << "D...\n";
        counter = 0;
        dir = "results/rosenbrock/n" + std::to_string(n);

        cfg.high = 30;
        cfg.low = -30;
        cfg.dim = n;
        cfg.T_max = 10000 * n;
        cfg.tau = 1.0 / std::sqrt(2.0 * std::sqrt((double)cfg.dim));
        cfg.tau_prime = 1.0 / std::sqrt(2.0 * (double)cfg.dim);

        cfg.sigma_max = 20;

        cfg.startSigma = 10;

        if (n == 30) {
            cfg.lambda = 80;
            cfg.momentum_beta = 0.9;
            cfg.stagnationLimit = 400;
            cfg.max_lambda = 500;
        }
        else if (n == 5) {
            cfg.lambda = 20;
            cfg.momentum_beta = 0.6;
            cfg.stagnationLimit = 200;
            cfg.max_lambda = 100;

        }
        else {
            cfg.lambda = 40;
            cfg.momentum_beta = 0.8;
            cfg.stagnationLimit = 200;
            cfg.max_lambda = 300;
        }

        
        fs::create_directories(dir);
        for (int r = 1; r <= runs; ++r) {
            auto history = run_Momentum_ES(f_rosenbrock, cfg);
            std::ofstream out(dir + "/run_" + std::to_string(r) + ".txt");
            for (double val : history) out << val << "\n";
            if (r % 10 == 0) std::cout << ".";
            if (history.back() > -1 && history.back() < 1) counter++;
        }
        std::cout << "number of good scores: " << counter << "\n";
        std::cout << " Done.\n";
        
        std::cout << ">>> Processing Salomon " << n << "D...\n";

        cfg.high = 100;
        cfg.low = -100;
        cfg.dim = n;
        cfg.T_max = 10000 * n;
        cfg.tau = 1.0 / std::sqrt(2.0 * std::sqrt((double)cfg.dim));
        cfg.tau_prime = 1.0 / std::sqrt(2.0 * (double)cfg.dim);

        cfg.sigma_max = 50;

        cfg.startSigma = 20;

        cfg.stagnationLimit = 100;

        if (n == 30) {
            cfg.lambda = 120;
            cfg.momentum_beta = 0.5;
            cfg.max_lambda = 500;
            cfg.stagnationLimit = 200;
        }
        else if (n == 5) {
            cfg.lambda = 20;
            cfg.momentum_beta = 0.6;
            cfg.max_lambda = 100;
            cfg.stagnationLimit = 100;
        }
        else {
            cfg.lambda = 40;
            cfg.momentum_beta = 0.7;
            cfg.max_lambda = 200;
            cfg.stagnationLimit = 100;
        }

        counter = 0;
        dir = "results/salomon/n" + std::to_string(n);
        fs::create_directories(dir);
        for (int r = 1; r <= runs; ++r) {
            auto history = run_Momentum_ES(f_salomon, cfg);
            std::ofstream out(dir + "/run_" + std::to_string(r) + ".txt");
            for (double val : history) out << val << "\n";
            if (r % 10 == 0) std::cout << ".";
            if (history.back() > -1 && history.back() < 1) counter++;
        }
        std::cout << "number of good scores: " << counter << "\n";
        std::cout << " Done.\n";
    }
    return 0;
}