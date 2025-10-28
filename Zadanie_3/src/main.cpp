#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include<unordered_set>
#include <numeric>
#include <algorithm>

//Stałe
//MIN i MAX do dziedziny
const double X_MIN = -3;
const double X_MAX = 3;
//Temperatura
const double T_start = 100;
const double T_n = 1;
const int N = 10000; //liczba zmian temperatury



// JEŚLI NIE MA UTWORZYĆ FOLDER \results\bin i \results\gray w \src

using namespace std;

double F(const vector<double>& x) {
    double sum = 0.0;
    for (double xi : x)
        sum += xi * xi;  // xi^2
    return sum;
}


void create_start(vector<double>& x, int n){
  x.clear();
  for(int i = 0; i<n; i++){
    x.push_back(X_MIN);
  }
}

// konwersja bin -> gray
uint16_t bin2gray(uint16_t num) {
    return num ^ (num >> 1);
}

// konwersja gray -> bin
uint16_t gray2bin(uint16_t gray) {
    uint16_t num = gray;
    for (uint16_t shift = 1; shift < 16; shift <<= 1)
        num ^= (gray >> shift);
    return num;
}

vector<uint16_t> encode(const vector<double>& x, bool use_gray) {
    vector<uint16_t> encoded;
    encoded.reserve(x.size());

    for (double xi : x) {
        // przeskalowanie do [0,1]
        double scaled = (xi - X_MIN) / (X_MAX - X_MIN);
        scaled = min(max(scaled, 0.0), 1.0); // zabezpieczenie
        uint16_t code = static_cast<uint16_t>(round(scaled * 65535.0)); // [0,1] -> [0,65535]
        if (use_gray) code = bin2gray(code);
        encoded.push_back(code);
    }

    return encoded;
}

vector<double> decode(const vector<uint16_t>& v, bool use_gray) {
  vector<double> decoded;
    for (uint16_t vi : v) {
        uint16_t bin_val = use_gray ? gray2bin(vi) : vi;
        double dec = X_MIN + (bin_val / 65535.0) * (X_MAX - X_MIN); // odwrotne przeskalowanie
        decoded.push_back(dec);
    }
  return decoded;
}


vector<double> get_neighbour(const vector<double>& x, mt19937 &gen, bool use_gray) {
     //kodowanie x do binarnej reprezentacji
    vector<uint16_t> x_bin = encode(x, use_gray);
    vector<uint16_t> x_new = x_bin; // kopia

    // losuj wymiar, który zmienimy
    uniform_int_distribution<int> dim_dist(0, x_bin.size() - 1);
    int dim = dim_dist(gen);

    // losuj bit do odwrócenia (0–15)
    uniform_int_distribution<int> bit_dist(0, 15);
    int bit = bit_dist(gen);

    // flip wybranego bitu
    x_new[dim] ^= (1 << bit); // maska 000..01 przesunięta o 'bit' (wylosowaną pozycję)

    vector<double> x_prime = decode(x_new, use_gray);
    return x_prime;
}

vector<double> get_neighbour_v2(const vector<double>& x, mt19937 &gen, int n, double T, bool use_gray) {
    //kodowanie x do binarnej reprezentacji
    vector<uint16_t> x_bin = encode(x, use_gray);
    vector<uint16_t> x_new = x_bin;

    uniform_int_distribution<int> bit_count_dist(1, 8); // 1–n bitow na wymiar
    uniform_int_distribution<int> bit_dist(0, 15);      // nr bitu do zmiany

    // zmieniamy każdy wymiar
    for (int dim = 0; dim < n; dim++) {
        int bits_to_flip = bit_count_dist(gen);

        for (int b = 0; b < bits_to_flip; b++) {
            int bit = bit_dist(gen);
            x_new[dim] ^= (1 << bit); // flip bitu w danym wymiarze
        }
    }
    vector<double> x_prime = decode(x_new, use_gray);
    return x_prime;
}

vector<double> get_neighbour_real(const vector<double>& x, mt19937 &gen, int n, double T, bool use_gray) {
    const double RANGE = X_MAX - X_MIN;

    // Zależność: im wyższa temperatura, tym większy krok
    double sigma = (RANGE * 0.2) * (1.0 + T / T_start) * n;
    normal_distribution<double> noise(0.0, sigma);

    vector<double> x_prime = x;

    // Liczba współrzędnych zależna od T
    int num_changed = max(1, int(n * (T / T_start)));
    num_changed = min(num_changed, n);

    vector<int> indices(n);
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), gen);
    indices.resize(num_changed);

    for (int idx : indices) {
        x_prime[idx] += noise(gen);
        if (x_prime[idx] < X_MIN) x_prime[idx] = X_MIN + (X_MIN - x_prime[idx]);
        if (x_prime[idx] > X_MAX) x_prime[idx] = X_MAX - (x_prime[idx] - X_MAX);
    }

    return x_prime;
}

//zad 3
//--------------------------------

// -----------------------------------------------------------
// Funkcja 1:
// f(x) = -5 / (1 + Σ x_i^2) + sin(cot(exp(5 / (1 + Σ x_i^2))))
// gdzie cot(u) = 1 / tan(u)
// -----------------------------------------------------------
double f1(const vector<double>& x) {
    double sum_sq = 0.0;
    for (double xi : x)
        sum_sq += xi * xi;

    double term1 = -5.0 / (1.0 + sum_sq);
    double inner = exp(5.0 / (1.0 + sum_sq));
    double term2 = sin(1.0 / tan(inner)); // cot(u) = 1/tan(u)

    return term1 + term2;
}

// -----------------------------------------------------------
// Funkcja 2 (Ackley):
// f(x) = -a * exp(-b * sqrt(1/d * Σ x_i^2))
//        - exp(1/d * Σ cos(c * x_i))
//        + a + exp(1)
// -----------------------------------------------------------
double f2(const vector<double>& x, double a = 20.0, double b = 0.2, double c = 2 * M_PI) {
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

//funckja do temperatury z wykładu pkt d)
double temp_func(double k){
    double A = log(T_start - T_n)/log(N);
    return T_start - pow(k,A);
}

//funckja do temperatury z wykładu pkt e)
double temp_func_2(double k) {
    return (T_start - T_n) / (1.0 + exp(0.3 * (k - N / 2.0))) + T_n;
}

vector<double> run_local_search(int n, int credits, bool use_gray, mt19937 &gen, vector<double> (*func)(const vector<double>&, mt19937&, int, double, bool)) {
  // punkt startowy
  vector<double> x; 
  create_start(x,n);
  
  double Fx = f1(x);
  vector<double> history; // historia punktów po których skaczemy
  history.push_back(Fx);
  credits -= 1;

 

  double T = T_start;
  int k=1;
  //główna pętla
  uniform_real_distribution<double> dist(0.0, 1.0);
  while(credits > 0){
      //wygeneruj losowego sąsiada (np. flip jednego bitu)
      vector<double> x_prime = func(x, gen,n, T, use_gray);

      //oblicz F(x')
      double Fx_prime = f1(x_prime);
      credits -= 1;

      //jeśli lepszy -> zaakceptuj
      if(Fx_prime < Fx){
        x = x_prime;
        Fx = Fx_prime;
      }
      else{
        double u = dist(gen);
        double acceptance_prob = exp((Fx - Fx_prime) / T); // kryterium metropolis
        if(u < acceptance_prob){
            x = x_prime;
            Fx = Fx_prime;
        }
      }
      history.push_back(Fx);
      T = temp_func_2(k);
      k++;
    }

    return history;
}


int main(){
  // generator losowy z seedem zależnym od czasu
  mt19937 gen(static_cast<unsigned long>(time(nullptr)));
  bool use_gray = false;

  // liczba wymiarów n: 2, 5, 10
  vector<int> dimensions = {2, 5, 10}; 

  // liczba kredytów
  int credits = 10000;

  // liczba powtórzeń algorytmu
  int repetitions = 100;                

  for (int n : dimensions) {
    for (int rep = 1; rep <= repetitions; rep++) {
      // Wersja binarna
      vector<double> history_bin = run_local_search(n, credits, false, gen, get_neighbour_v2);
      string filename_bin = "results\\bin\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
      ofstream out_bin(filename_bin);
      for (double val : history_bin) out_bin << val << "\n";
      out_bin.close();

      // Wersja Gray
      vector<double> history_gray = run_local_search(n, credits, true, gen, get_neighbour_v2);
      string filename_gray = "results\\gray\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
      ofstream out_gray(filename_gray);
      for (double val : history_gray) out_gray << val << "\n";
      out_gray.close();

      //Wersja rzeczywistoliczbowa
      vector<double> history_real = run_local_search(n, credits,false, gen, get_neighbour_real);
      string filename_real = "results\\real\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
      ofstream out_real(filename_real);
      for (double val : history_real) out_real << val << "\n";
      out_real.close();
    }
    cout << "Zakonczono wszystkie powtorzenia dla n = " << n << endl;
  }

  return 0;
}