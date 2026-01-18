#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <string>
#include <iomanip>

// ================================================================
// KONFIGURACJA EKSPERYMENTU
// ================================================================
enum ProblemType { ZDT1, ZDT2, ZDT3, ZDT4, ZDT6 };

struct Config {
    ProblemType type = ZDT1;
    int num_variables = 30; // Zmieniane na 10, 30, 50
    std::string output_filename = "results.csv";
};

// Globalna konfiguracja (ustawiana w main)
Config current_config;

const int NUM_OBJECTIVES = 2;
const int POPULATION_SIZE = 100;
const int MAX_GENERATIONS = 500;
const double M_PI = 3.14159265358979323846;

// Parametry jDE
const double TAU_1 = 0.1;
const double TAU_2 = 0.1;
const double F_LOWER = 0.1;
const double F_UPPER = 0.9;

// Generatory
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);
std::uniform_int_distribution<> int_dis(0, POPULATION_SIZE - 1);

// ================================================================
// STRUKTURY
// ================================================================
struct Individual {
    std::vector<double> variables;
    std::vector<double> objectives;
    double F;
    double CR;
    int rank;
    double crowding_distance;
    std::vector<int> dominated_by_indices;
    int domination_count;

    Individual() :
        variables(current_config.num_variables),
        objectives(NUM_OBJECTIVES),
        F(0.5), CR(0.9),
        rank(0), crowding_distance(0.0), domination_count(0) {}
};

// ================================================================
// FUNKCJE ZDT
// ================================================================

void evaluate_zdt(Individual& ind) {
    int m = current_config.num_variables;
    double f1, g, h;

    // ZDT1
    if (current_config.type == ZDT1) {
        f1 = ind.variables[0];
        double sum = 0.0;
        for (int i = 1; i < m; ++i) sum += ind.variables[i];
        g = 1.0 + 9.0 * (sum / (m - 1));
        h = 1.0 - std::sqrt(f1 / g);
        ind.objectives[0] = f1;
        ind.objectives[1] = g * h;
    }
    // ZDT2
    else if (current_config.type == ZDT2) {
        f1 = ind.variables[0];
        double sum = 0.0;
        for (int i = 1; i < m; ++i) sum += ind.variables[i];
        g = 1.0 + 9.0 * (sum / (m - 1));
        h = 1.0 - std::pow(f1 / g, 2.0);
        ind.objectives[0] = f1;
        ind.objectives[1] = g * h;
    }
    // ZDT3
    else if (current_config.type == ZDT3) {
        f1 = ind.variables[0];
        double sum = 0.0;
        for (int i = 1; i < m; ++i) sum += ind.variables[i];
        g = 1.0 + 9.0 * (sum / (m - 1));
        h = 1.0 - std::sqrt(f1 / g) - (f1 / g) * std::sin(10.0 * M_PI * f1);
        ind.objectives[0] = f1;
        ind.objectives[1] = g * h;
    }
    // ZDT4
    else if (current_config.type == ZDT4) {
        f1 = ind.variables[0];
        double sum = 0.0;
        for (int i = 1; i < m; ++i) {
            sum += std::pow(ind.variables[i], 2.0) - 10.0 * std::cos(4.0 * M_PI * ind.variables[i]);
        }
        g = 1.0 + 10.0 * (m - 1) + sum;
        h = 1.0 - std::sqrt(f1 / g);
        ind.objectives[0] = f1;
        ind.objectives[1] = g * h;
    }
    // ZDT6
    else if (current_config.type == ZDT6) {
        f1 = 1.0 - std::exp(-4.0 * ind.variables[0]) * std::pow(std::sin(6.0 * M_PI * ind.variables[0]), 6.0);
        double sum = 0.0;
        for (int i = 1; i < m; ++i) sum += ind.variables[i];
        g = 1.0 + 9.0 * std::pow(sum / (m - 1), 0.25);
        h = 1.0 - std::pow(f1 / g, 2.0);
        ind.objectives[0] = f1;
        ind.objectives[1] = g * h;
    }
}

//ograniczenia dziedziny
void enforce_constraints(Individual& ind) {
    // ograniczeia ZDT4
    double lower_bound = 0.0;
    double upper_bound = 1.0;

    // Dla zmiennej x1 zawsze [0, 1]
    if (ind.variables[0] < 0.0) ind.variables[0] = 0.0;
    if (ind.variables[0] > 1.0) ind.variables[0] = 1.0;

    // Dla pozostalych zmiennych
    if (current_config.type == ZDT4) {
        lower_bound = -5.0;
        upper_bound = 5.0;
    }

    for (size_t i = 1; i < ind.variables.size(); ++i) {
        if (ind.variables[i] < lower_bound) ind.variables[i] = lower_bound;
        if (ind.variables[i] > upper_bound) ind.variables[i] = upper_bound;
    }
}

// Inicjalizacja
std::vector<Individual> initialize_population() {
    std::vector<Individual> pop(POPULATION_SIZE);
    double lower = 0.0, upper = 1.0;

    // dla zdt4
    double zdt4_lower = -5.0, zdt4_upper = 5.0;

    for (auto& ind : pop) {
        // x1 zawsze [0, 1]
        ind.variables[0] = dis(gen);

        for (int i = 1; i < current_config.num_variables; ++i) {
            if (current_config.type == ZDT4) {
                std::uniform_real_distribution<> dist(zdt4_lower, zdt4_upper);
                ind.variables[i] = dist(gen);
            }
            else {
                ind.variables[i] = dis(gen);
            }
        }
        ind.F = 0.5;
        ind.CR = 0.9;
        evaluate_zdt(ind);
    }
    return pop;
}

// ================================================================
// ALGORYTM
// ================================================================

// Selekcja Turniejowa z Tlokiem
int crowded_tournament_selection(const std::vector<Individual>& pop) {
    int idx1 = int_dis(gen);
    int idx2 = int_dis(gen);
    const Individual& p1 = pop[idx1];
    const Individual& p2 = pop[idx2];
    if (p1.rank < p2.rank) return idx1;
    if (p2.rank < p1.rank) return idx2;
    if (p1.crowding_distance > p2.crowding_distance) return idx1;
    if (p2.crowding_distance > p1.crowding_distance) return idx2;
    return (dis(gen) < 0.5) ? idx1 : idx2;
}

std::vector<Individual> generate_offspring(const std::vector<Individual>& population) {
    std::vector<Individual> offspring;
    offspring.reserve(POPULATION_SIZE);

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        // Parametry jDE
        double new_F = population[i].F;
        double new_CR = population[i].CR;
        if (dis(gen) < TAU_1) new_F = F_LOWER + dis(gen) * F_UPPER;
        if (dis(gen) < TAU_2) new_CR = dis(gen);

        int r1, r2, r3;

        // HYBRYDOWA STRATEGIA
        // 90% czasu: Hole Filler (Elita + Tlok) -> szybka zbieznosc
        // 10% czasu: Random DE -> Ucieczka z lokalnych optimow (szczegolnie wazne dla ZDT4)
        if (dis(gen) < 0.9) {
            // Hole Filler Strategy
            r1 = crowded_tournament_selection(population);
            do { r2 = crowded_tournament_selection(population); } while (r2 == r1);
            do { r3 = crowded_tournament_selection(population); } while (r3 == r1 || r3 == r2);
        }
        else {
            r1 = int_dis(gen);
            do { r2 = int_dis(gen); } while (r2 == r1);
            do { r3 = int_dis(gen); } while (r3 == r1 || r3 == r2);
        }

        Individual child = population[i];
        child.F = new_F;
        child.CR = new_CR;

        int j_rand = std::uniform_int_distribution<>(0, current_config.num_variables - 1)(gen);

        for (int j = 0; j < current_config.num_variables; ++j) {
            if (dis(gen) < new_CR || j == j_rand) {
                double mutant_val = population[r1].variables[j] +
                    new_F * (population[r2].variables[j] - population[r3].variables[j]);
                child.variables[j] = mutant_val;
            }
            else {
                child.variables[j] = population[i].variables[j];
            }
        }
        enforce_constraints(child);
        evaluate_zdt(child);
        offspring.push_back(child);
    }
    return offspring;
}
bool dominates(const Individual& ind1, const Individual& ind2) {
    bool better = false;
    for (int i = 0; i < NUM_OBJECTIVES; ++i) {
        if (ind1.objectives[i] > ind2.objectives[i]) return false;
        if (ind1.objectives[i] < ind2.objectives[i]) better = true;
    }
    return better;
}

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
            }
            else if (dominates(population[j], population[i])) {
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

void calculate_crowding_distance(std::vector<Individual>& population, const std::vector<int>& front) {
    int l = front.size();
    if (l == 0) return;
    for (int idx : front) population[idx].crowding_distance = 0.0;
    for (int m = 0; m < NUM_OBJECTIVES; ++m) {
        std::vector<std::pair<double, int>> obj_idx;
        for (int idx : front) obj_idx.push_back({ population[idx].objectives[m], idx });
        std::sort(obj_idx.begin(), obj_idx.end());
        population[obj_idx.front().second].crowding_distance = 1e9;
        population[obj_idx.back().second].crowding_distance = 1e9;
        double range = obj_idx.back().first - obj_idx.front().first;
        if (range <= 1e-9) continue;
        for (int i = 1; i < l - 1; ++i) {
            population[obj_idx[i].second].crowding_distance +=
                (obj_idx[i + 1].first - obj_idx[i - 1].first) / range;
        }
    }
}

// Funkcja pomocnicza do zapisu checkpointu
void save_checkpoint(const std::vector<Individual>& population, int generation, std::ofstream& file) {
    std::vector<Individual> temp_pop = population;
    auto fronts = fast_non_dominated_sort(temp_pop);

    for (int idx : fronts[0]) {
        file << temp_pop[idx].objectives[0] << ","
            << temp_pop[idx].objectives[1] << ","
            << generation << "\n";
    }
}

// ================================================================
// MAIN
// ================================================================

void run_experiment(ProblemType type, int dimensions, std::string filename) {
    current_config.type = type;
    current_config.num_variables = dimensions;
    current_config.output_filename = filename;

    std::cout << ">>> Start Experiment: Type=" << type << " Dim=" << dimensions << std::endl;

    std::ofstream outfile(filename);
    outfile << "f1,f2,generation\n";

    auto population = initialize_population();
    auto initial_fronts = fast_non_dominated_sort(population);
    for (const auto& front : initial_fronts) calculate_crowding_distance(population, front);

    for (int gen = 1; gen <= MAX_GENERATIONS; ++gen) {
        auto offspring = generate_offspring(population);
        std::vector<Individual> R_t = population;
        R_t.insert(R_t.end(), offspring.begin(), offspring.end());

        auto fronts = fast_non_dominated_sort(R_t);
        std::vector<Individual> next_population;
        int i = 0;
        while (next_population.size() + fronts[i].size() <= POPULATION_SIZE) {
            calculate_crowding_distance(R_t, fronts[i]);
            for (int idx : fronts[i]) next_population.push_back(R_t[idx]);
            i++;
            if (i >= fronts.size()) break;
        }
        if (next_population.size() < POPULATION_SIZE && i < fronts.size()) {
            calculate_crowding_distance(R_t, fronts[i]);
            std::sort(fronts[i].begin(), fronts[i].end(), [&](int a, int b) {
                return R_t[a].crowding_distance > R_t[b].crowding_distance;
                });
            int needed = POPULATION_SIZE - next_population.size();
            for (int k = 0; k < needed; ++k) next_population.push_back(R_t[fronts[i][k]]);
        }
        population = next_population;

        // Checkpointy: 20, 50, 100, 500
        if (gen == 20 || gen == 50 || gen == 100 || gen == 500) {
            save_checkpoint(population, gen, outfile);
            std::cout << "  Checkpoint saved: Gen " << gen << std::endl;
        }
    }
    outfile.close();
    std::cout << "Saved: " << filename << "\n" << std::endl;
}

int main() {
    std::vector<ProblemType> problems = { ZDT1, ZDT2, ZDT3, ZDT4, ZDT6 };
    std::vector<std::string> prob_names = { "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6" };
    std::vector<int> dims = { 10, 30, 50 };

    for (size_t p = 0; p < problems.size(); ++p) {
        for (int d : dims) {
            std::string fname = prob_names[p] + "_D" + std::to_string(d) + "_results.csv";
            run_experiment(problems[p], d, fname);
        }
    }

    return 0;
}