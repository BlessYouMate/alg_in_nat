#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include <set>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;

struct Point {
    int id;
    vector<double> objectives;
    int frontId = 0;
};

// --- FUNKCJE POMOCNICZE DO FORMATOWANIA ---

// Funkcja zamieniająca liczbę na string z PRZECINKIEM (dla Excela w PL)
string plNum(double val) {
    string s = to_string(val);
    replace(s.begin(), s.end(), '.', ',');
    return s;
}

// --- LOGIKA DOMINACJI I ALGORYTMY ---
bool dominates(const Point& b, const Point& a) {
    bool betterInAny = false;
    for (size_t i = 0; i < b.objectives.size(); ++i) {
        if (b.objectives[i] < a.objectives[i]) return false;
        if (b.objectives[i] > a.objectives[i]) betterInAny = true;
    }
    return betterInAny;
}

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

vector<Point> frontRecursive(vector<Point> P) {
    if (P.size() <= 1) return P;
    size_t mid = P.size() / 2;
    vector<Point> T, B;
    for (size_t i = 0; i < mid; ++i) T.push_back(P[i]);
    for (size_t i = mid; i < P.size(); ++i) B.push_back(P[i]);

    vector<Point> T_prime = frontRecursive(T);
    vector<Point> B_prime = frontRecursive(B);

    vector<Point> M = T_prime;
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
    sort(P.begin(), P.end(), [](const Point& a, const Point& b) {
        if (a.objectives[0] != b.objectives[0]) return a.objectives[0] > b.objectives[0];
        if (a.objectives.size() > 1) return a.objectives[1] > b.objectives[1];
        return a.id < b.id;
        });
    return frontRecursive(P);
}

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

bool checkEquality(const vector<Point>& A, const vector<Point>& B) {
    if (A.size() != B.size()) return false;
    set<int> idsA, idsB;
    for (const auto& p : A) idsA.insert(p.id);
    for (const auto& p : B) idsB.insert(p.id);
    return idsA == idsB;
}

// --- FUNKCJE EKSPERYMENTALNE ---
//w zapisywaniu do pliku korzystamy z "polskich" separatowrow - ";" jako separator i "," jako oddzielenie czesi ulamkowej od calkowitej - na potrzeby excel
void runTask1Experiment(ofstream& csvFile, int numPoints, int numDims, string datasetName) {
    cout << "Zadanie 1: " << datasetName << "... ";
    auto data = generateData(numPoints, numDims);
    auto resultNaive = naiveAlgorithm(data);
    auto resultKung = kungAlgorithm(data);

    if (checkEquality(resultNaive, resultKung)) {
        cout << "OK (Znaleziono " << resultKung.size() << ")\n";
    }
    else {
        cout << "BLAD!\n";
    }

    set<int> paretoIDs;
    for (auto p : resultKung) paretoIDs.insert(p.id);

    for (const auto& p : data) {
        csvFile << datasetName << ";" << p.id << ";" << (paretoIDs.count(p.id) ? 1 : 0) << ";0;";

        for (int k = 0; k < 5; ++k) {
            if (k < numDims) csvFile << plNum(p.objectives[k]);
            else csvFile << "0";

            if (k < 4) csvFile << ";";
        }
        csvFile << "\n";
    }
}

vector<Point> loadDataFromFile(string filename) {
    vector<Point> data;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Brak pliku: " << filename << " (generuje losowe)\n";
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
    cout << "Wczytano " << data.size() << " pkt z pliku.\n";
    return data;
}

void runTask2Experiment(ofstream& csvFile, string filename) {
    cout << "Zadanie 2: Obieranie warstw..." << endl;
    vector<Point> currentSet = loadDataFromFile(filename);
    vector<Point> allSortedPoints;
    int frontCounter = 1;

    while (!currentSet.empty()) {
        vector<Point> currentFront = kungAlgorithm(currentSet);
        set<int> frontIDs;
        for (auto& p : currentFront) {
            p.frontId = frontCounter;
            allSortedPoints.push_back(p);
            frontIDs.insert(p.id);
        }
        vector<Point> nextSet;
        for (const auto& p : currentSet) {
            if (frontIDs.find(p.id) == frontIDs.end()) {
                nextSet.push_back(p);
            }
        }
        currentSet = nextSet;
        frontCounter++;
    }

    cout << "Liczba frontow: " << (frontCounter - 1) << endl;

    for (const auto& p : allSortedPoints) {
        csvFile << "Task2_Peeling" << ";" << p.id << ";0;" << p.frontId << ";";
        csvFile << plNum(p.objectives[0]) << ";" << plNum(p.objectives[1]) << ";0;0;0\n";
    }
}

int main() {
    ofstream file("wyniki.csv");

    file << "Dataset;ID;IsNonDominated;FrontID;Obj1;Obj2;Obj3;Obj4;Obj5\n";

    runTask1Experiment(file, 100, 2, "Task1_2D_100");
    runTask1Experiment(file, 1000, 5, "Task1_5D_1000");

    runTask2Experiment(file, "dane.txt");

    file.close();
    cout << "\nGotowe! Otworz 'wyniki.csv' w Excelu." << endl;
    return 0;
}