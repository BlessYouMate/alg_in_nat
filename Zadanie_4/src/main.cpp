#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>

// ================================================================
// KONFIGURACJA I STRUKTURY
// ================================================================
struct Config {
    int dim = 10;               // wymiar przestrzeni
    int popSize = 50;          // liczba osobników
    double low = -3.0;         // dolna granica
    double high = 3.0;         // górna granica
    double sigma0 = 0.1;       // początkowa wartosc sigma
    double sigma_min = 1e-4;   // minimalna wartosc sigma
    double sigma_max = high / 3;    // maksymalna wartosc sigma
    double p_jump = 0.1;       // szansa na duży skok (Cauchy)
    double crossoverProb = 0.8;// prawdopodobieństwo krzyżowania
    double tournamentP = 0.8;  // prawdopodobieństwo zwycięstwa lepszego osobnika
    int eliteCount = 2;        // ilu najlepszych zachować
    int T_max = 10000;         // maksymalna liczba pokoleń
    int stagnationLimit = 50;  // liczba generacji bez poprawy przed restartem
};

using Individual = std::vector<double>;
std::mt19937 rng(std::random_device{}());

// ================================================================
// FUNKCJE LOSOWE
// ================================================================
double randUniform(double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

double randNormal(double mean = 0.0, double stddev = 1.0) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(rng);
}

double randCauchy(double x0 = 0.0, double gamma = 1.0) {
    std::cauchy_distribution<double> dist(x0, gamma);
    return dist(rng);
}

// ================================================================
// NARZĘDZIA POMOCNICZE
// ================================================================
double reflect(double val, double low, double high) {
    if (val < low) return low + (low - val);
    if (val > high) return high - (val - high);
    return val;
}

// ================================================================
// INICJALIZACJA POPULACJI
// ================================================================
std::vector<Individual> initializePopulation(const Config& cfg) {
    std::vector<Individual> population(cfg.popSize, Individual(cfg.dim));
    for (auto& indiv : population) {
        for (int i = 0; i < cfg.dim; ++i) {
            indiv[i] = randUniform(cfg.low, cfg.high);
        }
    }
    return population;
}

// ================================================================
// FUNKCJA OCENY (fitness) - np. Sphere function
// ================================================================
double evaluate(const Individual& x, double (*objective)(const Individual&)) {
    return objective(x);
}

// ================================================================
// SELEKCJA (turniejowa z prawdopodobieństwem wygranej p)
// ================================================================
Individual tournamentSelection(
    const std::vector<Individual>& population,
    const std::vector<double>& fitness,
    int tournamentSize,
    double p)
{
    std::uniform_int_distribution<int> indexDist(0, (int)population.size() - 1);
    int bestIdx = indexDist(rng);
    for (int i = 1; i < tournamentSize; ++i) {
        int challenger = indexDist(rng);
        // z prawdopodobieństwem p wybieramy lepszego (mniejszy fitness = lepszy)
        bool challengerWins = randUniform(0, 1) < p && fitness[challenger] < fitness[bestIdx];
        if (challengerWins) bestIdx = challenger;
    }
    return population[bestIdx];
}

// ================================================================
// KRZYŻOWANIE (proste arytmetyczne crossover dwóch rodziców)
// ================================================================
Individual crossover(const Individual& parent1, const Individual& parent2, const Config& cfg) {
    Individual child(cfg.dim);
    for (int i = 0; i < cfg.dim; ++i) {
        double alpha = randUniform(0.0, 1.0);
        child[i] = alpha * parent1[i] + (1 - alpha) * parent2[i];
    }
    return child;
}

// ================================================================
// MUTACJA (czyli funkcja sąsiedztwa)
// ================================================================
Individual mutate(const Individual& x, const Config& cfg, double sigma) {
    Individual y = x;
    bool jump = (randUniform(0.0, 1.0) < cfg.p_jump);

    for (int i = 0; i < cfg.dim; ++i) {
        double step;
        if (jump) {
            step = sigma * 5.0 * randCauchy();
        }
        else {
            step = sigma * randNormal();
        }

        y[i] = reflect(x[i] + step, cfg.low, cfg.high);
    }

    return y;
}


// ================================================================
// ADAPTACJA SIGMA wg reguły 1/5 sukcesu
// ================================================================
double adaptSigma(double sigma, double successRate, const Config& cfg) {
    double factorUp = 1.5;
    double factorDown = 0.82;

    if (successRate > 0.2)
        sigma *= factorUp;
    else
        sigma *= factorDown;

    // ograniczenie do dopuszczalnego zakresu
    if (sigma < cfg.sigma_min) sigma = cfg.sigma_min;
    if (sigma > cfg.sigma_max) sigma = cfg.sigma_max;

    return sigma;
}

// ================================================================
// FUNKCJE TESTOWE (cele optymalizacji)
// ================================================================

// Funkcja 1, Dziedzina: [-3, 3]
double f1(const Individual& x) {
    double sum_sq = 0.0;
    for (double xi : x)
        sum_sq += xi * xi;

    double term1 = -5.0 / (1.0 + sum_sq);
    double inner = exp(-5.0 / (1.0 + sum_sq));
    double tan_inner = tan(inner);

    if (fabs(tan_inner) < 1e-100) {
        return term1;
    }

    double term2 = sin(1.0 / tan_inner);
    return term1 + term2;
}

// Funkcja 2 (Ackley), Dziedzina: [-32.768, 32.768]
double f2(const Individual& x) {
    int d = (int)x.size();
    double a = 20.0, b = 0.2, c = 2 * acos(-1.0);

    double sum_sq = 0.0, sum_cos = 0.0;
    for (double xi : x) {
        sum_sq += xi * xi;
        sum_cos += cos(c * xi);
    }

    double term1 = -a * exp(-b * sqrt(sum_sq / d));
    double term2 = -exp(sum_cos / d);

    return term1 + term2 + a + exp(1.0);
}

// ================================================================
// GŁÓWNA PĘTLA ALGORYTMU
// ================================================================

std::vector<double> run_GA_real(double (*objective)(const Individual&), const Config& cfg) {
    std::vector<double> history;
    //Inicjalizacja populacji
    auto population = initializePopulation(cfg);
    std::vector<double> fitness(cfg.popSize);
    int credits = cfg.T_max;
    for (int i = 0; i < cfg.popSize; ++i) {
        fitness[i] = evaluate(population[i], objective);
        credits--;
    }

    double sigma = cfg.sigma0;
    double bestFitness = *std::min_element(fitness.begin(), fitness.end());
    Individual bestIndividual = population[std::min_element(fitness.begin(), fitness.end()) - fitness.begin()];
    int stagnationCounter = 0;

    //Główna pętla
    for (int t = 0; t < credits;) {

        //Selekcja + tworzenie nowej populacji
        std::vector<Individual> offspring;
        int tournamentSize = 2 + (t / 100); // presja rosnie w czasie - coraz wiecej uczestnikow turnieju
        if (tournamentSize > 5) tournamentSize = 5;

        while ((int)offspring.size() < cfg.popSize) {
            Individual parent1 = tournamentSelection(population, fitness, tournamentSize, cfg.tournamentP);
            Individual parent2 = tournamentSelection(population, fitness, tournamentSize, cfg.tournamentP);

            Individual child = parent1;
            if (randUniform(0.0, 1.0) < cfg.crossoverProb) {
                child = crossover(parent1, parent2, cfg);
            }

            child = mutate(child, cfg, sigma);
            offspring.push_back(child);
        }

        //Ewaluacja nowej populacji
        std::vector<double> newFitness(cfg.popSize);
        int successCount = 0;
        for (int i = 0; i < cfg.popSize; ++i) {
            if (t >= credits) break;
            newFitness[i] = evaluate(offspring[i], objective);
            t++;
            if (newFitness[i] < fitness[i % fitness.size()]) {
                successCount++;
            }
        }

        //Adaptacja sigma (Rechenberg 1/5)
        double successRate = (double)successCount / cfg.popSize;
        sigma = adaptSigma(sigma, successRate, cfg);

        //Elityzm – zachowaj najlepszych z poprzedniego i nowowygenerowanego pokolenia
        std::vector<std::pair<double, Individual>> all;
        for (int i = 0; i < cfg.popSize; ++i)
            all.push_back({ newFitness[i], offspring[i] });
        for (int i = 0; i < cfg.popSize; ++i)
            all.push_back({ fitness[i], population[i] });
        std::sort(all.begin(), all.end(), [](auto& a, auto& b) { return a.first < b.first; });

        population.clear();
        fitness.clear();
        for (int i = 0; i < cfg.popSize; ++i) {
            population.push_back(all[i].second);
            fitness.push_back(all[i].first);
        }

        //Aktualizacja najlepszego
        if (fitness[0] < bestFitness) {
            bestFitness = fitness[0];
            bestIndividual = population[0];
            stagnationCounter = 0;
        }
        else {
            stagnationCounter++;
        }

        //Restart po stagnacji
        if (stagnationCounter > cfg.stagnationLimit) {
            std::cout << "[INFO] Restart populacji po stagnacji.\n";
            for (int i = cfg.eliteCount; i < cfg.popSize; ++i) {
                population[i] = mutate(bestIndividual, cfg, sigma); // restart wokół najlepszego
            }
            stagnationCounter = 0;
            sigma = cfg.sigma_max;
        }
        history.push_back(bestFitness);
    }

    //Wynik końcowy
    std::cout << "\nNajlepszy wynik: " << bestFitness << "\n";
    std::cout << "Najlepszy osobnik: ";
    for (double v : bestIndividual) std::cout << v << " ";
    std::cout << std::endl;
    return history;
}
int main() {
    Config cfg;

    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_real(f1, cfg); //pierwszy argument - funkcja celu (f1 lub f2)

        std::string filename_real = "results\\f1\\real\\results_n" + std::to_string(i) + "_run_real.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_real(f2, cfg); //pierwszy argument - funkcja celu (f1 lub f2)

        std::string filename_real = "results\\f2\\real\\results_n" + std::to_string(i) + "_run_real.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }

    return 0;
}
