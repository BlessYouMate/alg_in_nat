#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <cmath>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void runExperiment(int N) {
    const double minCoord = -1.5;
    const double maxCoord = 1.5;
    const double r = 1.0;  // radius of the circle

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(minCoord, maxCoord);

    std::string pointsFilename = "points_" + std::to_string(N) + ".txt";
    std::string resultsFilename = "results_" + std::to_string(N) + ".txt";

    std::ofstream pointsFile(pointsFilename);
    std::ofstream resultsFile(resultsFilename);

    if (!pointsFile || !resultsFile) {
        std::cerr << "Error: unable to open output files for N = " << N << "!" << std::endl;
        return;
    }

    int insideCircle = 0;
    pointsFile << std::fixed << std::setprecision(6);

    // Generate random points
    for (int i = 0; i < N; ++i) {
        double x = dist(gen);
        double y = dist(gen);
        bool inside = (x * x + y * y <= r * r);

        if (inside) insideCircle++;

        pointsFile << "x: " << x << ", y: " << y 
                   << ", inside_circle: " << (inside ? "yes" : "no") << "\n";
    }

    // Estimate the area of the circle
    double areaSquare = (maxCoord - minCoord) * (maxCoord - minCoord);
    double estimatedArea = areaSquare * static_cast<double>(insideCircle) / N;

    resultsFile << "Number of points: " << N << "\n";
    resultsFile << "Points inside the circle: " << insideCircle << "\n";
    resultsFile << "Estimated area of the circle: " << estimatedArea << "\n";
    resultsFile << "Exact area (pi*r^2): " << M_PI * r * r << "\n";
    resultsFile << "Relative error: "
                << std::abs(estimatedArea - M_PI * r * r) / (M_PI * r * r) * 100
                << " %" << std::endl;

    std::cout << "Experiment with " << N << " points completed.\n"
              << "Data saved to '" << pointsFilename << "' and '" << resultsFilename << "'.\n";
}

int main() {
    runExperiment(1000);
    runExperiment(10000);
    runExperiment(100000);

    std::cout << "\nAll experiments completed successfully.\n";
    return 0;
}
