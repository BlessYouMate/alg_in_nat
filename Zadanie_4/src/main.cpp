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
    int popSize = 100;          // liczba osobników
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


    double maxStepFraction = 0.05; // 5% ograniczenie rozmiaru pojedynczego skoku

    // --- Parametry dla wersji binarnej ---
    double p_cross_binary = 0.9; // Prawdopodobieństwo krzyżowania 
    double p_mut_binary = 2.0 / (10.0 * 16.0); // Prawd. mutacji JEDNEGO bitu 
    bool useGrayCoding = false; // Globalny przełącznik 

    double p_mut_gene = 0.1;
};

// --- DEFINICJE TYPÓW ---
// Genotyp/Fenotyp dla algorytmu rzeczywistego 
using Individual_Real = std::vector<double>;

// Genotyp dla algorytmu binarnego
using Individual_Bin = std::vector<uint16_t>;


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
// KONWERSJA BINARNA / GRAYA (Twoje funkcje)
// ================================================================
uint16_t bin2gray(uint16_t num) {
    return num ^ (num >> 1);
}

uint16_t gray2bin(uint16_t gray) {
    uint16_t num = gray;
    for (uint16_t shift = 1; shift < 16; shift <<= 1)
        num ^= (gray >> shift);
    return num;
}

// Dekoduje Genotyp (Binarny) na Fenotyp (Rzeczywisty)
std::vector<double> decode(const Individual_Bin& v, bool use_gray, double xmin, double xmax) {
    std::vector<double> decoded;
    decoded.reserve(v.size()); // Dobra praktyka
    double domain_width = xmax - xmin;
    for (uint16_t vi : v) {
        uint16_t bin_val = use_gray ? gray2bin(vi) : vi;
        double dec = xmin + (bin_val / 65535.0) * domain_width;
        decoded.push_back(dec);
    }
    return decoded;
}

// ================================================================
// NARZĘDZIA POMOCNICZE (DLA WERSJI RZECZYWISTEJ)
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
// MIERNIK RÓŻNORODNOŚCI POPULACJI (DLA WERSJI RZECZYWISTEJ)
// ================================================================
double computePopulationDiversity(const std::vector<Individual_Real>& population, const Config& cfg) {
    int n = cfg.popSize;
    int d = cfg.dim;
    if (n == 0) return 0.0;

    // oblicz środek populacji
    Individual_Real mean(d, 0.0);
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
// INICJALIZACJA POPULACJI (RZECZYWISTA)
// ================================================================
std::vector<Individual_Real> initializePopulation_real(const Config& cfg) {
    std::vector<Individual_Real> population(cfg.popSize, Individual_Real(cfg.dim));
    for (auto& indiv : population) {
        for (int i = 0; i < cfg.dim; ++i) {
            indiv[i] = randUniform(cfg.low, cfg.high);
        }
    }
    return population;
}

// ================================================================
// FUNKCJA OCENY (fitness) - wrapper dla wersji rzeczywistej
// ================================================================
double evaluate_real(const Individual_Real& x, double (*objective)(const Individual_Real&)) {
    return objective(x);
}

// ================================================================
// SELEKCJA (turniejowa dla wersji RZECZYWISTEJ)
// ================================================================
Individual_Real tournamentSelection_real(
    const std::vector<Individual_Real>& population,
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
// KRZYŻOWANIE (RZECZYWISTE - arytmetyczne)
// ================================================================
Individual_Real crossover_real(const Individual_Real& parent1, const Individual_Real& parent2, const Config& cfg) {
    Individual_Real child(cfg.dim);
    for (int i = 0; i < cfg.dim; ++i) {
        double alpha = randUniform(0.0, 1.0);
        child[i] = alpha * parent1[i] + (1 - alpha) * parent2[i];
    }
    return child;
}

// ================================================================
// MUTACJA (RZECZYWISTA - Gauss/Cauchy)
// ================================================================
Individual_Real mutate_real(const Individual_Real& x, const Config& cfg, double sigma,
    double maxStepFraction, double p_jump)
{
    Individual_Real y = x;
    bool jump = (randUniform(0.0, 1.0) < p_jump);

    for (int i = 0; i < cfg.dim; ++i) {
        double step;
        if (jump) {
            step = sigma * randCauchy(0, 0.2);
        }
        else {
            step = sigma * randNormal();
        }
        double range = cfg.high - cfg.low;
        double maxStep = maxStepFraction * range;
        if (step > maxStep) step = maxStep;
        if (step < -maxStep) step = -maxStep;
        y[i] = reflect(x[i] + step, cfg.low, cfg.high);
    }

    return y;
}


// ================================================================
// ROZRZUCENIE WARTOSCI (DLA WERSJI RZECZYWISTEJ)
// ================================================================
void adaptiveDiversifyPopulation(std::vector<Individual_Real>& population,
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
}


// ================================================================
// ADAPTACJA SIGMA (DLA WERSJI RZECZYWISTEJ)
// ================================================================
double adaptSigma(double sigma, double successRate, const Config& cfg) {
    double factorUp = 1.2;
    double factorDown = 0.9;

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
// ================================================================
// OPERATORY DLA WERSJI BINARNEJ
// ================================================================
// ================================================================

// ================================================================
// INICJALIZACJA POPULACJI (BINARNA)
// ================================================================
std::vector<Individual_Bin> initializePopulation_binary(const Config& cfg) {
    std::vector<Individual_Bin> population(cfg.popSize, Individual_Bin(cfg.dim));
    // Dystrybucja do losowania dowolnej 16-bitowej liczby
    std::uniform_int_distribution<uint16_t> dist(0, 65535);

    for (auto& indiv : population) {
        for (int i = 0; i < cfg.dim; ++i) {
            indiv[i] = dist(rng); // Wypełnij losowymi bitami
        }
    }
    return population;
}

// ================================================================
// FUNKCJA OCENY (BINARNA - z dekodowaniem)
// ================================================================
double evaluate_binary(const Individual_Bin& x_bin, double (*objective)(const Individual_Real&), const Config& cfg) {
    // 1. Dekoduj Genotyp (Binarny) na Fenotyp (Rzeczywisty)
    Individual_Real x_real = decode(x_bin, cfg.useGrayCoding, cfg.low, cfg.high);
    // 2. Oceń Fenotyp
    return objective(x_real);
}

// ================================================================
// SELEKCJA (turniejowa dla wersji BINARNEJ)
// ================================================================
Individual_Bin tournamentSelection_binary(
    const std::vector<Individual_Bin>& population,
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
// KRZYŻOWANIE (BINARNE - Uniform)
// ================================================================
Individual_Bin crossover_binary(const Individual_Bin& parent1, const Individual_Bin& parent2, const Config& cfg) {
    Individual_Bin child(cfg.dim);
    for (int i = 0; i < cfg.dim; ++i) {
        // Dla każdej z 10 liczb (genów) losuj, od którego rodzica ją wziąć
        if (randUniform(0.0, 1.0) < 0.5) {
            child[i] = parent1[i];
        }
        else {
            child[i] = parent2[i];
        }
    }
    return child;
}

// ================================================================
// KRZYŻOWANIE (BINARNE - JEDNOPUNKTOWE "na genach")
// ================================================================
Individual_Bin crossover_binary_one_point(const Individual_Bin& parent1, const Individual_Bin& parent2, const Config& cfg) {
    Individual_Bin child(cfg.dim);
    std::uniform_int_distribution<int> dist(1, cfg.dim - 1);
    int crossover_point = dist(rng); // Losuj punkt cięcia (od 1 do 9)

    // Weź geny od rodzica 1 aż do punktu cięcia
    for (int i = 0; i < crossover_point; ++i) {
        child[i] = parent1[i];
    }
    // Weź resztę genów od rodzica 2
    for (int i = crossover_point; i < cfg.dim; ++i) {
        child[i] = parent2[i];
    }
    return child;
}

// ================================================================
// KRZYŻOWANIE (BINARNE - DWUPUNKTOWE "na genach")
// ================================================================
Individual_Bin crossover_binary_two_point(const Individual_Bin& parent1, const Individual_Bin& parent2, const Config& cfg) {
    Individual_Bin child(cfg.dim);

    std::uniform_int_distribution<int> dist(0, cfg.dim - 1);
    int p1 = dist(rng);
    int p2 = dist(rng);
    if (p1 > p2) std::swap(p1, p2); // Upewnij się, że p1 <= p2

    // Geny od rodzica 1
    for (int i = 0; i < p1; ++i)
        child[i] = parent1[i];

    // Geny od rodzica 2 (pomiędzy punktami)
    for (int i = p1; i < p2; ++i)
        child[i] = parent2[i];

    // Geny znowu od rodzica 1
    for (int i = p2; i < cfg.dim; ++i)
        child[i] = parent1[i];

    return child;
}

// ================================================================
// MUTACJA (BINARNA - Bit-Flip)
// ================================================================
Individual_Bin mutate_binary(const Individual_Bin& x, const Config& cfg) {
    Individual_Bin y = x;
    double p_mut = cfg.p_mut_binary;

    for (int i = 0; i < cfg.dim; ++i) { // Dla każdego z 10 genów 
        for (int k = 0; k < 16; ++k) { // Dla każdego z 16 bitów w genie
            if (randUniform(0.0, 1.0) < p_mut) {
                // Odwróć k-ty bit używając operatora XOR
                y[i] ^= (1 << k);
            }
        }
    }
    return y;
}

// ================================================================
// MUTACJA (BINARNA - Zastąpienie Całego Genu)
// ================================================================
Individual_Bin mutate_binary_gene_reset(const Individual_Bin& x, const Config& cfg) {
    Individual_Bin y = x;
    double p_gene_mut = cfg.p_mut_gene;

    std::uniform_int_distribution<uint16_t> dist(0, 65535);

    for (int i = 0; i < cfg.dim; ++i) { // Dla każdego z 10 genów (uint16_t)
        if (randUniform(0.0, 1.0) < p_gene_mut) {
            // Zastąp CAŁY gen (uint16_t) nową losową wartością
            y[i] = dist(rng);
        }
    }
    return y;
}


// ================================================================
// FUNKCJE TESTOWE
// ================================================================

// Funkcja 1, Dziedzina: [-3, 3]
double f1(const Individual_Real& x) {
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
double f2(const Individual_Real& x) {
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
// GŁÓWNA PĘTLA ALGORYTMU (RZECZYWISTA)
// ================================================================

std::vector<double> run_GA_real(double (*objective)(const Individual_Real&), Config cfg) {
    // Aktualizacja sigmy na wypadek zmiany dziedziny (low/high)
    cfg.sigma_max = (cfg.high - cfg.low) / 3.0;

    std::vector<double> history;
    //Inicjalizacja populacji
    auto population = initializePopulation_real(cfg);
    std::vector<double> fitness(cfg.popSize);
    int credits = cfg.T_max;
    for (int i = 0; i < cfg.popSize; ++i) {
        fitness[i] = evaluate_real(population[i], objective);
        credits--;
    }

    double sigma = cfg.sigma0 * sqrt(cfg.high - cfg.low); // Skaluj sigmę do zakresu
    double bestFitness = *std::min_element(fitness.begin(), fitness.end());
    Individual_Real bestIndividual = population[std::min_element(fitness.begin(), fitness.end()) - fitness.begin()];
    int stagnationCounter = 0;

    //Główna pętla
    for (int t = 0; t < credits;) {
        //Selekcja + tworzenie nowej populacji
        std::vector<Individual_Real> offspring;
        // --- Dynamiczna presja selekcyjna ---
        // Presja turnieju rośnie w miarę postępu optymalizacji (progress)
        // i dodatkowo chwilowo zwiększa się przy stagnacji, aby pomóc "wydostać się" z lokalnego minimum.
        double progress = static_cast<double>(t) / cfg.T_max; // [0, 1]

        double maxStepFraction_dyn = cfg.maxStepFraction * (1.0 - 0.8 * progress);
        double p_jump_dyn = cfg.p_jump * (1.0 - progress);

        int baseSize = 2 + static_cast<int>(progress * 3.0); // rośnie płynnie z 2 → 5
        int bonus = std::min(2, stagnationCounter / std::max(1, cfg.stagnationLimit / 4)); // +0..2 w zależności od stagnacji
        double range = cfg.high - cfg.low;
        int tournamentSize = 2 + (int)(progress * 4);

        // (opcjonalne logowanie do debugowania)
        /*if (t % 200 == 0) {
            std::cout << "[DEBUG] progress=" << progress
                << " stagnation=" << stagnationCounter
                << " tournamentSize=" << tournamentSize << "\n";
        }*/

        while ((int)offspring.size() < cfg.popSize) {
            Individual_Real parent1 = tournamentSelection_real(population, fitness, tournamentSize, cfg.tournamentP);
            Individual_Real parent2 = tournamentSelection_real(population, fitness, tournamentSize, cfg.tournamentP);

            Individual_Real child = parent1;
            if (randUniform(0.0, 1.0) < cfg.crossoverProb) {
                child = crossover_real(parent1, parent2, cfg);
            }

            child = mutate_real(child, cfg, sigma, maxStepFraction_dyn, p_jump_dyn);
            offspring.push_back(child);
        }

        //Ewaluacja nowej populacji
        std::vector<double> newFitness(cfg.popSize);
        int successCount = 0;
        for (int i = 0; i < cfg.popSize; ++i) {
            if (t >= credits) break;
            newFitness[i] = evaluate_real(offspring[i], objective);
            t++;
            if (newFitness[i] < fitness[i % fitness.size()]) {
                successCount++;
            }
        }

        //Adaptacja sigma (Rechenberg 1/5)
        double successRate = (double)successCount / cfg.popSize;
        sigma = adaptSigma(sigma, successRate, cfg);

        //Elityzm – zachowaj najlepszych z poprzedniego i nowowygenerowanego pokolenia
        std::vector<std::pair<double, Individual_Real>> all;
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
                population[i] = mutate_real(population[i], cfg, sigma, maxStepFraction_dyn, p_jump_dyn);

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
    // std::cout << "\nNajlepszy wynik (Real): " << bestFitness << "\n";
    //std::cout << "Najlepszy osobnik: ";
    //for (double v : bestIndividual) std::cout << v << " ";
    //std::cout << std::endl;
    return history;
}

// ================================================================
// GŁÓWNA PĘTLA ALGORYTMU (BINARNA)
// ================================================================
std::vector<double> run_GA_binary(double (*objective)(const Individual_Real&), Config cfg) {
    // Aktualizacja parametrów binarnych na wypadek zmiany dziedziny
    cfg.p_mut_binary = 1.0 / (double)(cfg.dim * 16);

    std::vector<double> history;
    //Inicjalizacja populacji
    auto population = initializePopulation_binary(cfg);
    std::vector<double> fitness(cfg.popSize);
    int credits = cfg.T_max;
    for (int i = 0; i < cfg.popSize; ++i) {
        fitness[i] = evaluate_binary(population[i], objective, cfg);
        credits--;
    }

    double bestFitness = *std::min_element(fitness.begin(), fitness.end());
    int stagnationCounter = 0; // Prosta stagnacja tylko do debugowania

    //Główna pętla
    for (int t = 0; t < credits;) {

        //Selekcja + tworzenie nowej populacji
        std::vector<Individual_Bin> offspring;
        double progress = static_cast<double>(t) / cfg.T_max;
        int tournamentSize = 2 + static_cast<int>(progress * 3.0); // Rośnie od 2 do 5

        while ((int)offspring.size() < cfg.popSize) {
            Individual_Bin parent1 = tournamentSelection_binary(population, fitness, tournamentSize, cfg.tournamentP);
            Individual_Bin parent2 = tournamentSelection_binary(population, fitness, tournamentSize, cfg.tournamentP);

            Individual_Bin child = parent1;
            if (randUniform(0.0, 1.0) < cfg.p_cross_binary) {
                child = crossover_binary(parent1, parent2, cfg);
            }

            child = mutate_binary(child, cfg);
            offspring.push_back(child);
        }

        //Ewaluacja nowej populacji
        std::vector<double> newFitness(cfg.popSize);
        for (int i = 0; i < cfg.popSize; ++i) {
            if (t >= credits) break;
            newFitness[i] = evaluate_binary(offspring[i], objective, cfg);
            t++;
        }

        // W GA nie ma adaptacji sigmy

        //Elityzm – zachowaj najlepszych
        std::vector<std::pair<double, Individual_Bin>> all;
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
            stagnationCounter = 0;
        }
        else {
            stagnationCounter++;
        }

        // W prostym GA nie ma restartu ani dywersyfikacji (można by dodać, ale to uproszczenie)

        history.push_back(bestFitness);
    }

    //Wynik końcowy
    // std::string type = cfg.useGrayCoding ? "Gray" : "Binary";
    // std::cout << "\nNajlepszy wynik (" << type << "): " << bestFitness << "\n";
    return history;
}


// ================================================================
// FUNKCJA MAIN 
// ================================================================
int main() {
    Config cfg;

    // === EKSPERYMENT 1: F1 (Real) ===
    std::cout << "--- Uruchamianie F1 (Real-Valued) ---" << std::endl;
    cfg.low = -3.0;
    cfg.high = 3.0;
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_real(f1, cfg);
        std::string filename_real = "results\\f1\\real\\results_n" + std::to_string(i) + "_run_real.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    std::cout << "\nZakonczono dla f1 real" << "\n";

    // === EKSPERYMENT 2: F1 (Binary) ===
    std::cout << "--- Uruchamianie F1 (Binary) ---" << std::endl;
    cfg.low = -3.0;
    cfg.high = 3.0;
    cfg.useGrayCoding = false; // Kodowanie binarne
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_binary(f1, cfg);
        std::string filename_real = "results\\f1\\bin\\results_n" + std::to_string(i) + "_run_binary.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    std::cout << "\nZakonczono dla f1 bin" << "\n";

    // === EKSPERYMENT 3: F1 (Gray) ===
    std::cout << "--- Uruchamianie F1 (Gray) ---" << std::endl;
    cfg.low = -3.0;
    cfg.high = 3.0;
    cfg.useGrayCoding = true; // Kodowanie Graya
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_binary(f1, cfg);
        std::string filename_real = "results\\f1\\gray\\results_n" + std::to_string(i) + "_run_gray.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    std::cout << "\nZakonczono dla f1 gray" << "\n";


    // === EKSPERYMENT 4: F2 (Real) ===
    std::cout << "--- Uruchamianie F2 (Real-Valued) ---" << std::endl;
    cfg.low = -32.768;
    cfg.high = 32.768;
    // cfg.popSize = 50;
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_real(f2, cfg);
        std::string filename_real = "results\\f2\\real\\results_n" + std::to_string(i) + "_run_real.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    std::cout << "\nZakonczono dla f2 real" << "\n";

    // cfg.popSize = 100;
    // === EKSPERYMENT 5: F2 (Binary) ===
    std::cout << "--- Uruchamianie F2 (Binary) ---" << std::endl;
    cfg.low = -32.768;
    cfg.high = 32.768;
    cfg.useGrayCoding = false; // Kodowanie binarne
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_binary(f2, cfg);
        std::string filename_real = "results\\f2\\bin\\results_n" + std::to_string(i) + "_run_binary.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    std::cout << "\nZakonczono dla f2 bin" << "\n";

    // === EKSPERYMENT 6: F2 (Gray) ===
    std::cout << "--- Uruchamianie F2 (Gray) ---" << std::endl;
    cfg.low = -32.768;
    cfg.high = 32.768;
    cfg.useGrayCoding = true; // Kodowanie Graya
    for (int i = 1; i < 101; i++) {
        std::vector<double> history = run_GA_binary(f2, cfg);
        std::string filename_real = "results\\f2\\gray\\results_n" + std::to_string(i) + "_run_gray.txt";
        std::ofstream out_real(filename_real);
        for (double val : history) out_real << val << "\n";
        out_real.close();
    }
    std::cout << "\nZakonczono dla f2 gray" << "\n";

    return 0;
}