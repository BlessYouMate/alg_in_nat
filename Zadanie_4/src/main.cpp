#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>
#include "benchmarks.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    double sigma_max = high / 3; // maksymalna wartosc sigma
    double p_jump = 0.1;       // szansa na duży skok (Cauchy)
    double crossoverProb = 0.8;// prawdopodobieństwo krzyżowania
    double tournamentP = 0.8;  // prawdopodobieństwo zwycięstwa lepszego osobnika
    int eliteCount = 2;        // ilu najlepszych zachować
    int T_max = 10000;         // maksymalna liczba pokoleń
    int stagnationLimit = 50;  // liczba generacji bez poprawy przed restartem
    double sigma_restart_multiplier = 3;
    double restartBestFraction = 0.5;

    double diversityThreshold = 0.05;     // próg różnorodności (np. 5% zakresu przestrzeni)
    double maxReplaceFraction = 0.3;      // maksymalnie 30% populacji (poza elitą) podmieniamy
    double minReplaceFraction = 0.0;      // pod koniec już nic nie podmieniamy


    double maxStepFraction = 0.5; //ograniczenie rozmiaru pojedynczego skoku
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
    // Bezpieczne odbicia aż do momentu wejścia do [low, high]
    // Jeśli odbicia nie pomagają (bardzo duże wartości), stosujemy clamp.
    int safety = 0;
    while ((val < low || val > high) && safety < 10) {
        if (val < low) val = low + (low - val);
        else if (val > high) val = high - (val - high);
        ++safety;
    }
    // Jeśli nadal poza obszarem, po prostu przytnij do granicy
    if (val < low) val = low;
    if (val > high) val = high;
    return val;
}

// ================================================================
// MIERNIK RÓŻNORODNOŚCI POPULACJI
// ================================================================
double computePopulationDiversity(const std::vector<Individual>& population, const Config& cfg) {
    int n = cfg.popSize;
    int d = cfg.dim;
    if (n == 0) return 0.0;

    // oblicz środek populacji
    Individual mean(d, 0.0);
    for (const auto& indiv : population) {
        for (int j = 0; j < d; ++j) {
            mean[j] += indiv[j];
        }
    }
    for (int j = 0; j < d; ++j) {
        mean[j] /= n;
    }

    // oblicz średnią wariancję odległości od środka
    double varSum = 0.0;
    for (const auto& indiv : population) {
        double distSq = 0.0;
        for (int j = 0; j < d; ++j) {
            double diff = indiv[j] - mean[j];
            distSq += diff * diff;
        }
        varSum += sqrt(distSq);
    }

    return varSum / n; // średnia odległość osobników od środka populacji
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
            step = sigma * 2.0 * randCauchy();
        }
        else {
            step = sigma * randNormal();
        }
        double range = cfg.high - cfg.low;
        double maxStep = cfg.maxStepFraction * range;
        if (step > maxStep) step = maxStep;
        if (step < -maxStep) step = -maxStep;
        y[i] = reflect(x[i] + step, cfg.low, cfg.high);
    }

    return y;
}


// ================================================================
// ROZRZUCENIE WARTOSCI GDY POPULACJA ZBYT MALO ROZNORODNA
// ================================================================
void adaptiveDiversifyPopulation(std::vector<Individual>& population,
    const Config& cfg,
    double diversity,
    double progress)
{
    // progress ∈ [0,1] – postęp optymalizacji
    if (diversity > cfg.diversityThreshold)
        return; // populacja nadal zróżnicowana – nie ruszamy

    // siła dywersyfikacji maleje z czasem
    double replaceFraction = (1.0 - progress) * cfg.maxReplaceFraction;
    if (replaceFraction < cfg.minReplaceFraction)
        return;

    int elite = cfg.eliteCount;
    int replaceCount = static_cast<int>((cfg.popSize - elite) * replaceFraction);
    if (replaceCount <= 0) return;

    // Wymień ostatnie replaceCount osobników na losowe punkty
    int startIdx = cfg.popSize - replaceCount;
    for (int i = startIdx; i < cfg.popSize; ++i) {
        for (int j = 0; j < cfg.dim; ++j)
            population[i][j] = randUniform(cfg.low, cfg.high);
    }

    /*std::cout << "[DIVERSIFY] generation progress=" << progress
        << " diversity=" << diversity
        << " replaced=" << replaceCount
        << " (" << replaceFraction * 100 << "%)\n";*/
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

// Schwefel
double schwefel(const Individual& ind) {
    double sum = 0.0;
    for (double xi : ind)
        sum += xi * sin(sqrt(fabs(xi))); // bez minusa tutaj!
    return 418.9829 * (double)ind.size() - sum; // minus przed sumą
}

// Michalewicz
double michalewicz(const std::vector<double>& x) {
    double m = 10.0; // parametr “trudności”
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 2 * m);
    }
    return -sum;
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
        // --- Dynamiczna presja selekcyjna ---
        // Presja turnieju rośnie w miarę postępu optymalizacji (progress)
        // i dodatkowo chwilowo zwiększa się przy stagnacji, aby pomóc "wydostać się" z lokalnego minimum.
        double progress = static_cast<double>(t) / cfg.T_max; // [0, 1]
        int baseSize = 2 + static_cast<int>(progress * 3.0); // rośnie płynnie z 2 → 5
        int bonus = std::min(2, stagnationCounter / std::max(1, cfg.stagnationLimit / 4)); // +0..2 w zależności od stagnacji
        int tournamentSize = std::min(6, baseSize + bonus); // ogranicz do 6

        // (opcjonalne logowanie do debugowania)
        /*if (t % 200 == 0) {
            std::cout << "[DEBUG] progress=" << progress
                << " stagnation=" << stagnationCounter
                << " tournamentSize=" << tournamentSize << "\n";
        }*/

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
            //std::cout << "[INFO] Restart populacji po stagnacji.\n";

            // liczba osobników tworzonych przez mutację najlepszego
            int restartBestCount = static_cast<int>((cfg.popSize - cfg.eliteCount) * cfg.restartBestFraction);
            if (restartBestCount < 0) restartBestCount = 0;
            if (restartBestCount > cfg.popSize - cfg.eliteCount)
                restartBestCount = cfg.popSize - cfg.eliteCount;

            // część wokół najlepszego osobnika
            for (int i = cfg.eliteCount; i < cfg.eliteCount + restartBestCount; ++i) {
                population[i] = mutate(bestIndividual, cfg, sigma);
            }

            //część losowa
            for (int i = cfg.eliteCount + restartBestCount; i < cfg.popSize; ++i) {
                for (int j = 0; j < cfg.dim; ++j) {
                    population[i][j] = randUniform(cfg.low, cfg.high);
                }
            }

            // reset stagnacji i parametru sigma
            stagnationCounter = 0;
            sigma *= cfg.sigma_restart_multiplier;
            if (sigma > cfg.sigma_max) sigma = cfg.sigma_max;
        }

        //rozproszeie wartości jesli populacja zbyt malo roznorodna
        double diversity = computePopulationDiversity(population, cfg);
        adaptiveDiversifyPopulation(population, cfg, diversity, progress);

        history.push_back(bestFitness);
        //if (t % 100 == 0) { // co 100 pokoleń
        //    double diversity = computePopulationDiversity(population, cfg);
        //    std::cout << "[GEN " << t << "] Fitness: " << bestFitness
        //        << " | Różnorodność: " << diversity << std::endl;
        //}

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
    cfg.high = 32.768;
    cfg.low = -32.768;
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_real(f2, cfg); //pierwszy argument - funkcja celu (f1 lub f2)

        std::string filename_real = "results\\f2\\real\\results_n" + std::to_string(i) + "_run_real.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    cfg.high = M_PI;
    cfg.low = 0;
    for (int i = 1; i < 2; i++) {
        std::vector<double> history = run_GA_real(michalewicz, cfg); //pierwszy argument - funkcja celu (f1 lub f2)

        std::string filename_real = "results\\f2\\real\\results_n" + std::to_string(i) + "_run_real.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    return 0;
}

//int main() {
//    Config cfg;
//
//    std::vector<Benchmark> tests = {
//        Benchmark::Rastrigin,
//        Benchmark::Schwefel,
//        Benchmark::Griewank,
//        Benchmark::Levy,
//        Benchmark::Eggholder,
//        Benchmark::DropWave,
//        Benchmark::Rosenbrock,
//        Benchmark::Ackley
//    };
//
//    for (auto b : tests) {
//        int dim; double low, high;
//        getBenchmarkDomain(b, dim, low, high);
//        cfg.dim = dim;
//        cfg.low = low;
//        cfg.high = high;
//
//        auto fn = getBenchmarkFunction(b);
//        std::cout << "=== Test: " << (int)b << " dim=" << dim << " range=[" << low << "," << high << "] ===\n";
//
//        // uruchom kilka powtórzeń dla stabilnej statystyki
//        for (int run = 0; run < 5; ++run) {
//            std::vector<double> hist = run_GA_real(fn, cfg);
//            // (opcjonalnie) zapisz hist do pliku results/<nazwa>_run.txt
//            std::string filename = "results/test_" + std::to_string((int)b) + "_run" + std::to_string(run) + ".txt";
//            std::ofstream out(filename);
//            for (double v : hist) out << v << "\n";
//            out.close();
//        }
//    }
//    return 0;
//}
