#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <functional> 
#include <tuple>      
#include <numeric>    
#include <algorithm> 

// UWAGA: PRZED URUCHOMIENIEM UPEWNIJ SIĘ, ŻE ISTNIEJĄ FOLDERY:
// \results\f1\bin
// \results\f1\gray
// \results\f1\real
// \results\f2\bin
// \results\f2\gray
// \results\f2\real

using namespace std;

// --- FUNKCJE CELU ---

// Funkcja 1, Dziedzina: [-3, 3]
double f1(const vector<double>& x) {
    double sum_sq = 0.0;
    for (double xi : x)
        sum_sq += xi * xi;
    double term1 = -5.0 / (1.0 + sum_sq);
    double inner = exp(-5.0 / (1.0 + sum_sq));
    double tan_inner = tan(inner);
    if (abs(tan_inner) < 1e-100) {
        return term1; 
    }
    double term2 = sin(1.0 / tan_inner);
    return term1 + term2;
}

// Funkcja 2 (Ackley), Dziedzina: [-32.768, 32.768]
double f2(const vector<double>& x, double a = 20.0, double b = 0.2, double c = 2 * acos(-1.0)) {
    int d = x.size();
    double sum_sq = 0.0, sum_cos = 0.0;
    for (double xi : x) {
        sum_sq += xi * xi;
        sum_cos += cos(c * xi);
    }
    double term1 = -a * exp(-b * sqrt(sum_sq / d));
    double term2 = -exp(sum_cos / d);
    return term1 + term2 + a + exp(1.0);
}

// --- FUNKCJE POMOCNICZE (START, KODOWANIE, TEMPERATURA) ---

// Punkt startowy to narożnik dziedziny
void create_start(vector<double>& x, int n, double xmax){
  x.clear();
  for(int i = 0; i<n; i++){
    x.push_back(xmax); 
  }
}

uint16_t bin2gray(uint16_t num) {
    return num ^ (num >> 1);
}

uint16_t gray2bin(uint16_t gray) {
    uint16_t num = gray;
    for (uint16_t shift = 1; shift < 16; shift <<= 1)
        num ^= (gray >> shift);
    return num;
}

vector<uint16_t> encode(const vector<double>& x, bool use_gray, double xmin, double xmax) {
    vector<uint16_t> encoded;
    encoded.reserve(x.size());
    double domain_width = xmax - xmin; 
    for (double xi : x) {
        double scaled = (xi - xmin) / domain_width;
        uint16_t code = (uint16_t)round(scaled * 65535.0);
        if (use_gray) code = bin2gray(code);
        encoded.push_back(code);
    }
    return encoded;
}

vector<double> decode(const vector<uint16_t>& v, bool use_gray, double xmin, double xmax) {
  vector<double> decoded;
  double domain_width = xmax - xmin; 
    for (uint16_t vi : v) {
        uint16_t bin_val = use_gray ? gray2bin(vi) : vi;
        double dec = xmin + (bin_val / 65535.0) * domain_width;
        decoded.push_back(dec);
    }
  return decoded;
}

// Parametry dla harmonogramu wyżarzania
const double T_START = 100.0; 
const double T_FINAL = 0.1;

// Funkcja temperatury (harmonogram schładzania)
double temp_func(int k, int k_max) {
    if (k_max <= 0) return T_FINAL;
    double A = (T_START - T_FINAL) * (double)(k_max + 1) / (double)k_max;
    double T_k = A / (double)(k + 1) + T_START - A;
    return max(T_k, T_FINAL);
}

// --- OPERATORY SĄSIEDZTWA ---

// Sąsiedztwo dla reprezentacji BINARNEJ
vector<uint16_t> get_neighbour_binary(const vector<uint16_t>& x_bin, mt19937 &gen) {
    vector<uint16_t> x_new = x_bin; 
    uniform_int_distribution<int> dim_dist(0, x_bin.size() - 1);
    int dim = dim_dist(gen);
    uniform_int_distribution<int> bit_dist(0, 15);
    int bit = bit_dist(gen);
    // flip wybranego bitu
    x_new[dim] ^= (1 << bit); 
    return x_new;
}

// Sąsiedztwo dla reprezentacji RZECZYWISTEJ
vector<double> get_neighbour_real(const vector<double>& x, mt19937 &gen, double xmin, double xmax) {
    
    // Rozkład normalny (Gaussa) N(0, 1)
    normal_distribution<double> noise(0.0, 1.0);

    vector<double> x_prime = x;

    // Losujemy jeden wymiar do zmiany
    uniform_int_distribution<int> dim_dist(0, x.size() - 1);
    int dim_to_change = dim_dist(gen);

    // Zmieniamy wylosowaną współrzędną: x' = x + N(0, 1)
    x_prime[dim_to_change] += noise(gen);
        
    // Pilnujemy granic dziedziny
    if (x_prime[dim_to_change] < xmin) x_prime[dim_to_change] = xmin;
    if (x_prime[dim_to_change] > xmax) x_prime[dim_to_change] = xmax;
    
    return x_prime;
}


// --- GŁÓWNE FUNKCJE ALGORYTMU ---

// Wersja dla reprezentacji BINARNEJ (Bit-flip)
vector<double> run_sa_binary(int n, int credits, bool use_gray, mt19937 &gen,
                                       std::function<double(const vector<double>&)> objective_func,
                                       double xmin, double xmax) 
{
  vector<double> x; 
  create_start(x, n, xmax);
  double Fx = objective_func(x);
  double Fx_best = Fx; 
  vector<double> history; 
  history.push_back(Fx_best); 
  credits -= 1; // 1 kredyt zużyty na start

  vector<uint16_t> x_bin = encode(x, use_gray, xmin, xmax);
  int k_max_iterations = credits; // Reszta kredytów na pętlę

  for(int k = 0; k < k_max_iterations; k++) {
      double T = temp_func(k, k_max_iterations - 1);

      // Używamy sąsiedztwa binarnego
      vector<uint16_t> x_bin_prime = get_neighbour_binary(x_bin, gen);
      vector<double> x_prime = decode(x_bin_prime, use_gray, xmin, xmax);

      double Fx_prime = objective_func(x_prime); // 1 kredyt zużyty
      
      if (Fx_prime < Fx_best) {
          Fx_best = Fx_prime;
      }

      // Kryterium Metropolisa
      if(Fx_prime < Fx){
        x_bin = x_bin_prime;
        x = x_prime;
        Fx = Fx_prime;
      }
      else{
        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(gen);
        double acceptance_prob = exp((Fx - Fx_prime) / T);
        if(u < acceptance_prob){
            x_bin = x_bin_prime;
            x = x_prime;
            Fx = Fx_prime;
        }
      }
      // Zapisujemy NAJLEPSZY wynik
      history.push_back(Fx_best); 
  }
  return history;
}

// Wersja dla reprezentacji RZECZYWISTEJ (Gaussian noise)
vector<double> run_sa_real(int n, int credits, mt19937 &gen,
                                  std::function<double(const vector<double>&)> objective_func,
                                  double xmin, double xmax) 
{
  vector<double> x; 
  create_start(x, n, xmax);
  double Fx = objective_func(x);
  double Fx_best = Fx; 
  vector<double> history; 
  history.push_back(Fx_best); 
  credits -= 1; // 1 kredyt zużyty na start

  int k_max_iterations = credits; // Reszta kredytów na pętlę

  for(int k = 0; k < k_max_iterations; k++) {
      double T = temp_func(k, k_max_iterations - 1);

      // Używamy sąsiedztwa rzeczywistego (bez encode/decode)
      vector<double> x_prime = get_neighbour_real(x, gen, xmin, xmax);

      double Fx_prime = objective_func(x_prime); // 1 kredyt zużyty
      
      if (Fx_prime < Fx_best) {
          Fx_best = Fx_prime;
      }

      // Kryterium Metropolisa
      if(Fx_prime < Fx){
        x = x_prime;
        Fx = Fx_prime;
      }
      else{
        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(gen);
        double acceptance_prob = exp((Fx - Fx_prime) / T);
        if(u < acceptance_prob){
            x = x_prime;
            Fx = Fx_prime;
        }
      }
      // Zapisujemy NAJLEPSZY wynik
      history.push_back(Fx_best);
  }
  return history;
}


// --- FUNKCJA main ---

int main(){
  mt19937 gen(static_cast<unsigned long>(time(nullptr)));

  vector<int> dimensions = {10}; 

  int credits = 10000;
  int repetitions = 100;             

  // Definicja funkcji, ich nazw i domen
  vector<tuple<string, function<double(const vector<double>&)>, double, double>> functions_to_run = {
      make_tuple("f1", f1, -3.0, 3.0),
      // Dopasowanie f2 do wektora który oczekuje 1 argumentu
      make_tuple("f2", [](const vector<double>& v){ return f2(v); }, -32.768, 32.768)
  };

  for (const auto& func_info : functions_to_run) {
      string func_name = get<0>(func_info);
      auto objective_func = get<1>(func_info);
      double xmin = get<2>(func_info);
      double xmax = get<3>(func_info);

      for (int n : dimensions) {
        for (int rep = 1; rep <= repetitions; rep++) {
          
          // --- Wersja binarna ---
          vector<double> history_bin = run_sa_binary(
                n, credits, false, gen, objective_func, xmin, xmax
            );
          string filename_bin = "results\\" + func_name + "\\bin\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
          ofstream out_bin(filename_bin);
          for (double val : history_bin) out_bin << val << "\n";
          out_bin.close();

          // --- Wersja Gray ---
          vector<double> history_gray = run_sa_binary(
                n, credits, true, gen, objective_func, xmin, xmax
            );
          string filename_gray = "results\\" + func_name + "\\gray\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
          ofstream out_gray(filename_gray);
          for (double val : history_gray) out_gray << val << "\n";
          out_gray.close();

          // --- Wersja rzeczywistoliczbowa ---
          vector<double> history_real = run_sa_real(
              n, credits, gen, objective_func, xmin, xmax
          );
          string filename_real = "results\\" + func_name + "\\real\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
          ofstream out_real(filename_real);
          for (double val : history_real) out_real << val << "\n";
          out_real.close();

        }
        cout << "Zakonczono: " << func_name << " dla n = " << n << endl;
      }
  }
  cout << "Zakonczono wszystkie eksperymenty." << endl;
  return 0;
}