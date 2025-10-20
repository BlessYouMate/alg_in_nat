#include<iostream>
#include<vector>
#include <random>
#include <fstream>
#include <cmath>
using namespace std;

double F(const vector<double>& x) {
    double sum = 0.0;
    for (double xi : x)
        sum += xi * xi;  // xi²
    return sum;
}

void create_start(vector<double>& x, int n){
  x.clear();
  for(int i = 0; i<n; i++){
    x.push_back(10.0);
  }
}

vector<uint16_t> encode(const vector<double>& x) {
    vector<uint16_t> encoded;
    encoded.reserve(x.size());

    for (double xi : x) {
        double scaled = (xi + 10.0) / 20.0;        // [-10,10] → [0,1]
        uint16_t code = (uint16_t)round(scaled * 65535.0); // [0,65535]
        encoded.push_back(code);
    }

    return encoded;
}

vector<double> decode(const vector<uint16_t>& v) {
  vector<double> decoded;
    for (uint16_t vi : v) {
        double dec = -10.0 + (vi / 65535.0) * 20.0;
        decoded.push_back(dec);
    }
  return decoded;
}


vector<uint16_t> get_neighbour(const vector<uint16_t>& x_bin) {
    vector<uint16_t> x_new = x_bin; // kopia

    static random_device rd;
    static mt19937 gen(rd());

    // losuj wymiar, który zmienimy
    uniform_int_distribution<int> dim_dist(0, x_bin.size() - 1);
    int dim = dim_dist(gen);

    // losuj bit do odwrócenia (0–15)
    uniform_int_distribution<int> bit_dist(0, 15);
    int bit = bit_dist(gen);

    // flip wybranego bitu
    x_new[dim] ^= (1 << bit);

    return x_new;
}

int main(){
  // liczba wymiarów n: 2, 5, 10
  int n = 2;
  // liczba kredytów (np. 10_000)
  int credits = 10000;

  // punkt startowy
  vector<double> x;
  create_start(x,n);
  // x = [10, 10, ..., 10]
  double Fx = F(x);
  vector<double> history;
  history.push_back(Fx);
  credits -= 1;

  //kodowanie x do binarnej reprezentacji
  vector<uint16_t> x_bin = encode(x);

  //główna pętla
  while(credits > 0){
    bool found_better = false;

  //pętla szukania pierwszego lepszego sąsiada
    while(!found_better && credits > 0){
      //wygeneruj losowego sąsiada (np. flip jednego bitu)
      vector<uint16_t> x_bin_prime = get_neighbour(x_bin);
      vector<double> x_prime = decode(x_bin_prime);

      //oblicz F(x')
      double Fx_prime = F(x_prime);
      credits -= 1;

      //jeśli lepszy -> zaakceptuj
      if(Fx_prime < Fx){
        x_bin = x_bin_prime;
        x = x_prime;
        Fx = Fx_prime;
        found_better = true;
      }
      history.push_back(Fx);
    }
      //jeśli nie znaleziono lepszego sąsiada -> koniec
      if(!found_better){
        break;
      }
  }

  ofstream out("results_n2.txt");
  for (double val : history) out << val << "\n";
  out.close();

  cout << "Najlepsze znalezione F(x) = " << Fx << endl;
  cout << "Najlepszy punkt: ";
  for (double xi : x) cout << xi << " ";
  cout << endl;
  cout << "Liczba krokopw (wywolan F): " << history.size() << endl;

  return 0;
}