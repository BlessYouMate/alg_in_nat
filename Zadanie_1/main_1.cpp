#include <iostream>
#include <fstream>
#include <random>
#include <vector>

void generate_dataset(int N, int bins, const std::string& name_prefix) {
    double minVal = 0.0;
    double maxVal = 1.0;

    std::vector<double> numbers;
    numbers.reserve(N);

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minVal, maxVal);

    // Generate numbers and save to data file
    std::ofstream data_file(name_prefix + "_data.txt");
    if (!data_file.is_open()) {
        std::cout << "Error: cannot open " << name_prefix << "_data.txt\n";
        return;
    }

    for (int i = 0; i < N; i++) {
        double x = dist(gen);
        numbers.push_back(x);
        data_file << x << "\n";
    }
    data_file.close();
    std::cout << "Saved " << N << " numbers to " << name_prefix << "_data.txt\n";

    // Create histogram
    std::vector<int> bins_count(bins, 0);
    double width = (maxVal - minVal) / bins;

    for (int i = 0; i < N; i++) {
        int index = static_cast<int>((numbers[i] - minVal) / width);
        if (index >= bins) index = bins - 1; // ensure maxVal fits last bin
        bins_count[index]++;
    }

    // Save histogram to TXT
    std::ofstream hist_file(name_prefix + "_histogram.txt");
    if (!hist_file.is_open()) {
        std::cout << "Error: cannot open " << name_prefix << "_histogram.txt\n";
        return;
    }
    for (int i = 0; i < bins; i++) {
        double lower = minVal + i * width;
        double upper = lower + width;
        hist_file << lower << "\t" << upper << "\t" << bins_count[i] << "\n";
    }
    hist_file.close();
    std::cout << "Saved histogram to " << name_prefix << "_histogram.txt\n";

    // Save histogram to CSV (count + percent)
    std::ofstream csv_file("histogram_1.csv", std::ios::app);
    if (!csv_file.is_open()) {
        std::cout << "Error: cannot open histogram_1.csv\n";
        return;
    }
    for (int i = 0; i < bins; i++) {
        double lower = minVal + i * width;
        double upper = lower + width;
        double percent = (static_cast<double>(bins_count[i]) / N) * 100.0;
        csv_file << N << "," << lower << "-" << upper << "," << bins_count[i] << "," << percent << "\n";
    }
    csv_file.close();
}

int main() {
    // Add header to CSV
    std::ofstream csv_file("histogram_1.csv");
    csv_file << "samples,bin_range,count,percent\n";
    csv_file.close();

    const int bins = 20;

    // Generate datasets
    generate_dataset(10000, bins, "data_10k_1");
    generate_dataset(100000, bins, "data_100k_1");
    generate_dataset(1000000, bins, "data_1M_1");

    std::cout << "All datasets generated. CSV contains all histograms with count and percent.\n";
    return 0;
}
