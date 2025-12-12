#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include <set>
#include <ctime>
#include <string>

using namespace std;

// Struktura reprezentująca punkt (rozwiązanie)
struct Point {
    int id;
    vector<double> objectives; 
};

// SPRAWDZENIE DOMINACJI 

// Sprawdza, czy punkt 'b' dominuje punkt 'a' (Maksymalizacja)
// b >= a we wszystkich kryteriach ORAZ b > a w przynajmniej jednym
bool dominates(const Point& b, const Point& a) {
    bool betterInAny = false;
    for (size_t i = 0; i < b.objectives.size(); ++i) {
        if (b.objectives[i] < a.objectives[i]) return false; // Gorszy w czymkolwiek -> Odpada
        if (b.objectives[i] > a.objectives[i]) betterInAny = true;
    }
    return betterInAny;
}

// ALGORYTMY

// 1. Naive Algorithm (Algorytm Naiwny)
vector<Point> naiveAlgorithm(const vector<Point>& P) {
    vector<Point> result;
    // Zewnętrzna pętla po wszystkich punktach
    for (size_t i = 0; i < P.size(); ++i) {
        bool isDominated = false;
        // Wewnętrzna pętla sprawdzająca każdego z każdym
        for (size_t j = 0; j < P.size(); ++j) {
            if (i == j) continue;
            if (dominates(P[j], P[i])) { 
                isDominated = true;
                break;
            }
        }
        if (!isDominated) result.push_back(P[i]);
    }
    return result;
}

// 2. Kung's Algorithm (Algorytm Kunga)
vector<Point> frontRecursive(vector<Point> P) {
    // Warunek stopu
    if (P.size() <= 1) return P;

    size_t mid = P.size() / 2;
    vector<Point> T, B;
    for(size_t i=0; i<mid; ++i) T.push_back(P[i]);
    for(size_t i=mid; i<P.size(); ++i) B.push_back(P[i]);

    // Wywołanie rekurencyjne
    vector<Point> T_prime = frontRecursive(T);
    vector<Point> B_prime = frontRecursive(B);

    // Scalanie
    // T_prime wchodzi zawsze, bo ma lepsze pierwsze kryterium (f1)
    vector<Point> M = T_prime; 
    
    // Sprawdzamy, czy punkty z B_prime są zdominowane przez kogoś z T_prime
    for (const auto& b : B_prime) {
        bool isDominated = false;
        for (const auto& t : T_prime) {
            if (dominates(t, b)) {
                isDominated = true;
                break;
            }
        }
        if (!isDominated) M.push_back(b);
    }
    return M;
}

vector<Point> kungAlgorithm(vector<Point> P) {
    // Sortowanie malejąco po 1. kryterium
    sort(P.begin(), P.end(), [](const Point& a, const Point& b) {
        if (a.objectives[0] != b.objectives[0]) return a.objectives[0] > b.objectives[0];
        if (a.objectives.size() > 1) return a.objectives[1] > b.objectives[1];
        return a.id < b.id;
    });
    return frontRecursive(P);
}

// FUNKCJE POMOCNICZE

// Generuje losowe dane
vector<Point> generateData(int count, int dims) {
    vector<Point> data;
    static mt19937 gen(time(0)); 
    uniform_real_distribution<> dist(0.0, 100.0);
    for (int i = 0; i < count; ++i) {
        Point p;
        p.id = i + 1;
        for (int d = 0; d < dims; ++d) p.objectives.push_back(dist(gen));
        data.push_back(p);
    }
    return data;
}

// Sprawdza identyczność wyników (porównuje zbiory ID)
bool checkEquality(const vector<Point>& A, const vector<Point>& B) {
    if (A.size() != B.size()) return false;
    set<int> idsA, idsB;
    for (const auto& p : A) idsA.insert(p.id);
    for (const auto& p : B) idsB.insert(p.id);
    return idsA == idsB;
}

void runExperiment(ofstream& csvFile, int numPoints, int numDims, string datasetName) {
    cout << "Processing: " << datasetName << " (" << numPoints << " pts, " << numDims << "D)... ";
    
    // 1. Generowanie danych
    auto data = generateData(numPoints, numDims);
    
    // 2. Uruchomienie algorytmów
    auto resultNaive = naiveAlgorithm(data);
    auto resultKung = kungAlgorithm(data);
    
    // 3. Weryfikacja poprawności
    if (checkEquality(resultNaive, resultKung)) {
        cout << "SUCCESS! Algorithms produced the same results (Found " << resultKung.size() << " non-dominated solutions)" << endl;
    } else {
        cout << "ERROR! Algorithms produced different results." << endl;
    }

    // 4. Zapis do CSV
    set<int> paretoIDs;
    for(auto p : resultKung) paretoIDs.insert(p.id);

    for(const auto& p : data) {
        csvFile << datasetName << "," << p.id << "," << (paretoIDs.count(p.id) ? 1 : 0) << ",";
        
        // Pętla zapisująca zawsze 5 kolumn (wypełnia zerami brakujące wymiary)
        for(int k = 0; k < 5; ++k) {
            if (k < numDims) csvFile << p.objectives[k];
            else csvFile << "0";
            
            if (k < 4) csvFile << ",";
        }
        csvFile << "\n";
    }
}

int main() {
    cout << "TASK 1" << endl;
    
    ofstream file("task1.csv");
    // Nagłówek pliku CSV
    file << "Dataset,ID,IsNonDominated,Obj1,Obj2,Obj3,Obj4,Obj5\n";

    // Uruchomienie eksperymentów
    runExperiment(file, 100, 2, "2D_100");
    runExperiment(file, 1000, 5, "5D_1000");

    file.close();
    cout << "\nDone. Results saved to 'task1.csv'." << endl;
    return 0;
}