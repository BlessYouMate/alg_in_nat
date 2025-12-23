#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include <set>
#include <ctime>
#include <string>

using namespace std;

struct Point {
    int id;
    vector<double> objectives;
    int frontId = 0;
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
    for (size_t i = 0; i < P.size(); ++i) {
        bool isDominated = false;
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
    if (P.size() <= 1) return P;

    size_t mid = P.size() / 2;
    vector<Point> T, B;
    for(size_t i=0; i<mid; ++i) T.push_back(P[i]);
    for(size_t i=mid; i<P.size(); ++i) B.push_back(P[i]);

    vector<Point> T_prime = frontRecursive(T);
    vector<Point> B_prime = frontRecursive(B);

    // Scalanie
    // T_prime wchodzi zawsze, bo ma lepsze pierwsze kryterium (f1) dzięki sortowaniu
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

// GENEROWANIE I WCZYTYWANIE DANYCH 

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

vector<Point> loadDataFromFile(string filename) {
    vector<Point> data;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: File " << filename << " not found! Using random data for test.\n";
        return generateData(200, 2);
    }

    double v1, v2;
    int id = 1;
    while (file >> v1 >> v2) {
        Point p;
        p.id = id++;
        p.objectives = { v1, v2 };
        data.push_back(p);
    }
    cout << "Loaded " << data.size() << " points from file.\n";
    return data;
}

// Sprawdza identyczność wyników
bool checkEquality(const vector<Point>& A, const vector<Point>& B) {
    if (A.size() != B.size()) return false;
    set<int> idsA, idsB;
    for (const auto& p : A) idsA.insert(p.id);
    for (const auto& p : B) idsB.insert(p.id);
    return idsA == idsB;
}

// ZADANIA

// Zadanie 1
void runTask1(ofstream& csvFile, int numPoints, int numDims, string datasetName) {
    cout << "Task 1: " << datasetName << "... ";
    auto data = generateData(numPoints, numDims);
    auto resultNaive = naiveAlgorithm(data);
    auto resultKung = kungAlgorithm(data);

    // Weryfikacja poprawności
    if (checkEquality(resultNaive, resultKung)) {
        cout << "OK (Found " << resultKung.size() << " non-dominated solutions)\n";
    } else {
        cout << "ERROR! Algorithms differ.\n";
    }

    set<int> paretoIDs;
    for (auto p : resultKung) paretoIDs.insert(p.id);

    for (const auto& p : data) {
        csvFile << datasetName << "," << p.id << "," << (paretoIDs.count(p.id) ? 1 : 0) << ",0,";

        // Pętla zapisująca zawsze 5 kolumn (wypełnia zerami brakujące wymiary)
        for (int k = 0; k < 5; ++k) {
            if (k < numDims) csvFile << p.objectives[k];
            else csvFile << "0";

            if (k < 4) csvFile << ",";
        }
        csvFile << "\n";
    }
}

// Zadanie 2
void runTask2(ofstream& csvFile, string filename) {
    cout << "Task 2: " << endl;
    vector<Point> currentSet = loadDataFromFile(filename);
    
    int frontCounter = 1;
    
    while (!currentSet.empty()) {
        vector<Point> currentFront = kungAlgorithm(currentSet);
        
        set<int> frontIDs;
        for (auto& p : currentFront) {
            p.frontId = frontCounter;
            frontIDs.insert(p.id);
            
            csvFile << "Task2" << "," << p.id << ",0," << p.frontId << ",";
            csvFile << p.objectives[0] << "," << p.objectives[1] << ",0,0,0\n";
        }
        
        // 2. Filtrowanie znalezionego frontu
        vector<Point> nextSet;
        for (const auto& p : currentSet) {
            // jesli nie ma w znalezionym froncie zapisz do nowej iteracji
            if (frontIDs.find(p.id) == frontIDs.end()) {
                nextSet.push_back(p);
            }
        }

        currentSet = nextSet;
        frontCounter++;
    }

    cout << "Total fronts found: " << (frontCounter - 1) << endl;
}

int main() {
    ofstream file("wyniki.csv");

    file << "Dataset,ID,IsNonDominated,FrontID,Obj1,Obj2,Obj3,Obj4,Obj5\n";

    runTask1(file, 100, 2, "Task1_2D_100");
    runTask1(file, 1000, 5, "Task1_5D_1000");

    runTask2(file, "dane.txt");

    file.close();
    cout << "\nDone! Open 'wyniki.csv' in Excel." << endl;
    return 0;
}