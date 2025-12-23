#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <limits>
#include <iomanip>

// Parametry problemu i algorytmu
const int NUM_VARIABLES = 30;   // Zgodnie z definicją ZDT1 [cite: 4058]
const int NUM_OBJECTIVES = 2;
const int POPULATION_SIZE = 100;
const int MAX_EVALUATIONS = 20000;
const double CROSSOVER_PROB = 0.9;
const double MUTATION_PROB = 1.0 / NUM_VARIABLES;
const double ETA_C = 20.0; // Indeks dystrybucji dla SBX 
const double ETA_M = 20.0; // Indeks dystrybucji dla mutacji 

// Generatory liczb losowych
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

struct Individual {
    std::vector<double> variables;
    std::vector<double> objectives;
    int rank;
    double crowding_distance;
    std::vector<int> dominated_by_indices; // S_p
    int domination_count;                  // n_p

    Individual() : variables(NUM_VARIABLES), objectives(NUM_OBJECTIVES), rank(0), crowding_distance(0.0), domination_count(0) {}
};

// Funkcja ZDT1 [cite: 4058]
void evaluate_zdt1(Individual& ind) {
    double f1 = ind.variables[0];
    
    double g = 0.0;
    for (int i = 1; i < NUM_VARIABLES; ++i) {
        g += ind.variables[i];
    }
    g = 1.0 + 9.0 * (g / (NUM_VARIABLES - 1));
    
    double h = 1.0 - std::sqrt(f1 / g);
    double f2 = g * h;

    ind.objectives[0] = f1;
    ind.objectives[1] = f2;
}

// Inicjalizacja populacji
std::vector<Individual> initialize_population() {
    std::vector<Individual> pop(POPULATION_SIZE);
    for (auto& ind : pop) {
        for (int i = 0; i < NUM_VARIABLES; ++i) {
            ind.variables[i] = dis(gen);
        }
        evaluate_zdt1(ind);
    }
    return pop;
}

// Sprawdzenie dominacji: czy ind1 dominuje ind2?
bool dominates(const Individual& ind1, const Individual& ind2) {
    bool at_least_one_better = false;
    for (int i = 0; i < NUM_OBJECTIVES; ++i) {
        if (ind1.objectives[i] > ind2.objectives[i]) {
            return false; // ind1 jest gorszy w tym kryterium (minimalizacja)
        }
        if (ind1.objectives[i] < ind2.objectives[i]) {
            at_least_one_better = true;
        }
    }
    return at_least_one_better;
}

// Szybkie sortowanie niedominowane [cite: 3237]
std::vector<std::vector<int>> fast_non_dominated_sort(std::vector<Individual>& population) {
    std::vector<std::vector<int>> fronts;
    std::vector<int> first_front;

    for (int i = 0; i < population.size(); ++i) {
        population[i].dominated_by_indices.clear();
        population[i].domination_count = 0;

        for (int j = 0; j < population.size(); ++j) {
            if (i == j) continue;
            if (dominates(population[i], population[j])) {
                population[i].dominated_by_indices.push_back(j);
            } else if (dominates(population[j], population[i])) {
                population[i].domination_count++;
            }
        }

        if (population[i].domination_count == 0) {
            population[i].rank = 1;
            first_front.push_back(i);
        }
    }
    fronts.push_back(first_front);

    int i = 0;
    while (i < fronts.size()) {
        std::vector<int> next_front;
        for (int p_idx : fronts[i]) {
            for (int q_idx : population[p_idx].dominated_by_indices) {
                population[q_idx].domination_count--;
                if (population[q_idx].domination_count == 0) {
                    population[q_idx].rank = i + 2;
                    next_front.push_back(q_idx);
                }
            }
        }
        if (next_front.empty()) break;
        fronts.push_back(next_front);
        i++;
    }
    return fronts;
}

// Obliczanie odległości tłoku (Crowding Distance) [cite: 3317]
void calculate_crowding_distance(std::vector<Individual>& population, const std::vector<int>& front) {
    int l = front.size();
    if (l == 0) return;

    for (int idx : front) {
        population[idx].crowding_distance = 0.0;
    }

    for (int m = 0; m < NUM_OBJECTIVES; ++m) {
        // Sortowanie frontu względem m-tego kryterium
        std::vector<std::pair<double, int>> obj_idx;
        for (int idx : front) {
            obj_idx.push_back({population[idx].objectives[m], idx});
        }
        std::sort(obj_idx.begin(), obj_idx.end());

        population[obj_idx.front().second].crowding_distance = std::numeric_limits<double>::infinity();
        population[obj_idx.back().second].crowding_distance = std::numeric_limits<double>::infinity();

        double min_obj = obj_idx.front().first;
        double max_obj = obj_idx.back().first;
        double range = max_obj - min_obj;

        if (range == 0) continue;

        for (int i = 1; i < l - 1; ++i) {
            population[obj_idx[i].second].crowding_distance += 
                (obj_idx[i + 1].first - obj_idx[i - 1].first) / range;
        }
    }
}

// Operator selekcji turniejowej (Crowded Comparison Operator) 
int tournament_selection(const std::vector<Individual>& pop) {
    int idx1 = std::uniform_int_distribution<>(0, pop.size() - 1)(gen);
    int idx2 = std::uniform_int_distribution<>(0, pop.size() - 1)(gen);

    const Individual& ind1 = pop[idx1];
    const Individual& ind2 = pop[idx2];

    if (ind1.rank < ind2.rank) return idx1;
    if (ind2.rank < ind1.rank) return idx2;
    if (ind1.crowding_distance > ind2.crowding_distance) return idx1;
    return idx2;
}

// Krzyżowanie SBX 
void sbx_crossover(const Individual& p1, const Individual& p2, Individual& c1, Individual& c2) {
    if (dis(gen) <= CROSSOVER_PROB) {
        for (int i = 0; i < NUM_VARIABLES; ++i) {
            if (dis(gen) <= 0.5) {
                if (std::abs(p1.variables[i] - p2.variables[i]) > 1e-14) {
                    double y1 = std::min(p1.variables[i], p2.variables[i]);
                    double y2 = std::max(p1.variables[i], p2.variables[i]);
                    double rand = dis(gen);
                    double beta = 1.0 + (2.0 * (y1 - 0.0) / (y2 - y1));
                    double alpha = 2.0 - std::pow(beta, -(ETA_C + 1.0));
                    
                    double betaq;
                    if (rand <= (1.0 / alpha)) {
                        betaq = std::pow(rand * alpha, (1.0 / (ETA_C + 1.0)));
                    } else {
                        betaq = std::pow(1.0 / (2.0 - rand * alpha), (1.0 / (ETA_C + 1.0)));
                    }
                    
                    c1.variables[i] = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    c2.variables[i] = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
                } else {
                    c1.variables[i] = p1.variables[i];
                    c2.variables[i] = p2.variables[i];
                }
            } else {
                c1.variables[i] = p1.variables[i];
                c2.variables[i] = p2.variables[i];
            }
            // Ograniczenia
            c1.variables[i] = std::max(0.0, std::min(1.0, c1.variables[i]));
            c2.variables[i] = std::max(0.0, std::min(1.0, c2.variables[i]));
        }
    } else {
        c1 = p1;
        c2 = p2;
    }
}

// Mutacja wielomianowa 
void polynomial_mutation(Individual& ind) {
    for (int i = 0; i < NUM_VARIABLES; ++i) {
        if (dis(gen) <= MUTATION_PROB) {
            double y = ind.variables[i];
            double yl = 0.0;
            double yu = 1.0;
            double delta1 = (y - yl) / (yu - yl);
            double delta2 = (yu - y) / (yu - yl);
            double rnd = dis(gen);
            double mut_pow = 1.0 / (ETA_M + 1.0);
            double deltaq;

            if (rnd <= 0.5) {
                double xy = 1.0 - delta1;
                double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (std::pow(xy, (ETA_M + 1.0)));
                deltaq = std::pow(val, mut_pow) - 1.0;
            } else {
                double xy = 1.0 - delta2;
                double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (std::pow(xy, (ETA_M + 1.0)));
                deltaq = 1.0 - std::pow(val, mut_pow);
            }
            ind.variables[i] = y + deltaq * (yu - yl);
            ind.variables[i] = std::max(yl, std::min(yu, ind.variables[i]));
        }
    }
}

int main() {
    std::cout << "Rozpoczynam algorytm NSGA-II dla ZDT1..." << std::endl;
    
    // Krok 1: Inicjalizacja
    auto population = initialize_population();
    int evaluations = POPULATION_SIZE;

    // Główna pętla
    while (evaluations < MAX_EVALUATIONS) {
        // Krok 2: Ewaluacja i tworzenie potomstwa
        std::vector<Individual> offspring;
        while (offspring.size() < POPULATION_SIZE) {
            int p1_idx = tournament_selection(population);
            int p2_idx = tournament_selection(population);
            
            Individual c1, c2;
            sbx_crossover(population[p1_idx], population[p2_idx], c1, c2);
            polynomial_mutation(c1);
            polynomial_mutation(c2);
            
            evaluate_zdt1(c1);
            evaluate_zdt1(c2);
            
            offspring.push_back(c1);
            offspring.push_back(c2);
            evaluations += 2;
        }

        // Krok 3: Łączenie populacji (Rt = Pt + Qt) [cite: 3333]
        std::vector<Individual> R_t = population;
        R_t.insert(R_t.end(), offspring.begin(), offspring.end());

        // Krok 4: Sortowanie niedominowane połączonej populacji [cite: 3334]
        auto fronts = fast_non_dominated_sort(R_t);

        // Krok 5: Tworzenie nowej populacji Pt+1 [cite: 3337]
        std::vector<Individual> next_population;
        int i = 0;
        while (next_population.size() + fronts[i].size() <= POPULATION_SIZE) {
            calculate_crowding_distance(R_t, fronts[i]);
            for (int idx : fronts[i]) {
                next_population.push_back(R_t[idx]);
            }
            i++;
            if (i >= fronts.size()) break;
        }

        // Dopełnianie populacji z ostatniego frontu (sortowanie po crowding distance) [cite: 3343]
        if (next_population.size() < POPULATION_SIZE && i < fronts.size()) {
            calculate_crowding_distance(R_t, fronts[i]);
            
            // Sortowanie malejące po crowding distance
            std::sort(fronts[i].begin(), fronts[i].end(), [&](int a, int b) {
                return R_t[a].crowding_distance > R_t[b].crowding_distance;
            });

            int needed = POPULATION_SIZE - next_population.size();
            for (int k = 0; k < needed; ++k) {
                next_population.push_back(R_t[fronts[i][k]]);
            }
        }

        population = next_population;
    }

    // Zapis wyników do pliku (tylko pierwszy front - rozwiązania Pareto optymalne)
    auto fronts = fast_non_dominated_sort(population);
    std::ofstream outfile("zdt1_results.csv");
    outfile << "f1,f2\n";
    for (int idx : fronts[0]) {
        outfile << population[idx].objectives[0] << "," << population[idx].objectives[1] << "\n";
    }
    outfile.close();

    std::cout << "Zakonczono. Wyniki zapisano w 'zdt1_results.csv'. Uzuj skryptu Python do wizualizacji." << std::endl;

    return 0;
}