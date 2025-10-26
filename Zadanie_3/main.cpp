#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <time.h>
#include <string>



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
    x.push_back(10.0);
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
        double scaled = (xi + 10.0) / 20.0; // xi - (xmin) / xmax - (xmin)        [-10,10] -> [0,1]
        // xi - xmin | przesunięcie na początek przedziału do 0
        // / (xmax - xmin) | normalizacja do [0,1]
        uint16_t code = (uint16_t)round(scaled * 65535.0); // [0,1] -> [0,65535]
        // * 65535 | przeskalowanie do docelowego przedziału
        if (use_gray) code = bin2gray(code);
        encoded.push_back(code);
    }

    return encoded;
}

vector<double> decode(const vector<uint16_t>& v, bool use_gray) {
  vector<double> decoded;
    for (uint16_t vi : v) {
        uint16_t bin_val = use_gray ? gray2bin(vi) : vi;
        double dec = -10.0 + (bin_val / 65535.0) * 20.0; // xmin + (code / 65535(2^B-1)) * (xmax-xmin)
        // x / 65535 | przeskalowanie do [0,1]
        // * (xmax - xmin) | rozciągnięcie do przedziału 
        // + xmin przesunięcie do przedziału
        decoded.push_back(dec);
    }
  return decoded;
}


vector<uint16_t> get_neighbour(const vector<uint16_t>& x_bin, mt19937 &gen) {
    vector<uint16_t> x_new = x_bin; // kopia

    // losuj wymiar, który zmienimy
    uniform_int_distribution<int> dim_dist(0, x_bin.size() - 1);
    int dim = dim_dist(gen);

    // losuj bit do odwrócenia (0–15)
    uniform_int_distribution<int> bit_dist(0, 15);
    int bit = bit_dist(gen);

    // flip wybranego bitu
    x_new[dim] ^= (1 << bit); // maska 000..01 przesunięta o 'bit' (wylosowaną pozycję)

    return x_new;
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
const double T_start = 100;

double temp_func(double k){
    double T_n = 1;
    int N = 50; //liczba zmian temperatury
    double A = log(T_start - T_n)/log(N);
    return T_start - pow(k,A);
}

vector<double> run_local_search(int n, int credits, bool use_gray, mt19937 &gen) {
  // punkt startowy
  vector<double> x; 
  create_start(x,n); // x = [10, 10, ..., 10]
  
  double Fx = F(x);
  vector<double> history; // historia punktów po których skaczemy
  history.push_back(Fx);
  credits -= 1;

  //kodowanie x do binarnej reprezentacji
  vector<uint16_t> x_bin = encode(x, use_gray);

  double T = T_start;
  int k=1;
  //główna pętla
  while(credits > 0){
      //wygeneruj losowego sąsiada (np. flip jednego bitu)
      vector<uint16_t> x_bin_prime = get_neighbour(x_bin, gen);
      vector<double> x_prime = decode(x_bin_prime, use_gray);

      //oblicz F(x')
      double Fx_prime = F(x_prime);
      credits -= 1;

      //jeśli lepszy -> zaakceptuj
      if(Fx_prime < Fx){
        x_bin = x_bin_prime;
        x = x_prime;
        Fx = Fx_prime;
      }
      else{
        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(gen);
        double acceptance_prob = exp((Fx - Fx_prime) / T); // kryterium metropolis
        if(u < acceptance_prob){
            x_bin = x_bin_prime;
            x = x_prime;
            Fx = Fx_prime;
        }
      }
      history.push_back(Fx);
      T = temp_func(k);
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
      vector<double> history_bin = run_local_search(n, credits, false, gen);
      string filename_bin = "results\\bin\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
      ofstream out_bin(filename_bin);
      for (double val : history_bin) out_bin << val << "\n";
      out_bin.close();

      // Wersja Gray
      vector<double> history_gray = run_local_search(n, credits, true, gen);
      string filename_gray = "results\\gray\\results_n" + to_string(n) + "_run" + to_string(rep) + ".txt";
      ofstream out_gray(filename_gray);
      for (double val : history_gray) out_gray << val << "\n";
      out_gray.close();
    }
    cout << "Zakonczono wszystkie powtorzenia dla n = " << n << endl;
  }

  return 0;
}