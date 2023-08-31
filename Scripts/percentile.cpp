#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm> // For std::transform
#include <cctype>    // For std::tolower

double calculateAverage(const std::vector<double>& numbers) {
    double sum = 0.0;

    // Calculate the sum of all elements in the vector
    for (const auto& num : numbers) {
        sum += num;
    }

    // Calculate the average
    double average = sum / numbers.size();

    return average;
}

bool containsSubstringInLowercase(const std::string& mainString, const std::string& subString) {
    std::string mainStringLower = mainString;
    std::string subStringLower = subString;

    // Convert both strings to lowercase
    std::transform(mainStringLower.begin(), mainStringLower.end(), mainStringLower.begin(), ::tolower);
    std::transform(subStringLower.begin(), subStringLower.end(), subStringLower.begin(), ::tolower);

    // Check if the lowercase substring exists in the lowercase main string
    return mainStringLower.find(subStringLower) != std::string::npos;
}


double percentileOfScore(std::vector<double> dataset, double score) {
    // Sort the dataset in ascending order
    std::sort(dataset.begin(), dataset.end());

    // Count how many elements are less than or equal to the score
    int count = 0;
    for (int i = 0; i < dataset.size(); i++) {
        if (dataset[i] <= score) {
            count++;
        }
    }

    // Calculate the percentile
    double percentile = (static_cast<double>(count - 1) / static_cast<double>(dataset.size() - 1)) * 100.0;

    return percentile;
}

std::unordered_map<std::string,double> getThrHash(std::string thrfile) {

    std::unordered_map<std::string, double> dataMap; // HashMap to store the data
    std::ifstream file(thrfile);
    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return dataMap;
    }


    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key, value;

        if (std::getline(ss, key, '\t') && std::getline(ss, value, '\t')) { // Assuming tab as the delimiter
            //std::cout << key << " " << value << std::endl;
            dataMap[key] = std::stod(value);
        }
    }

    file.close();

    return dataMap;
}

std::unordered_map<std::string,int> getNeutrophilHash(std::string sfile) {

    std::unordered_map<std::string, int> dataMap; // HashMap to store the data
    std::ifstream file(sfile);
    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return dataMap;
    }


    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key, value;

        if (std::getline(ss, key, '\t') && std::getline(ss, value, '\t')) {
            // Assuming tab as the delimiter
            //std::cout << key << " " << value << std::endl;
            if (containsSubstringInLowercase(value, "neutrophil")) {
              dataMap[key] = 1;
            } else  {
              dataMap[key] = 0;
            }
        }
    }

    file.close();

    return dataMap;
}

int main(int argc, char **argv) {
    char * pre = argv[1];
    std::string efile = std::string(pre) + "-expr.txt";
    std::string thrfile = std::string(pre) + "-thr.txt";
    std::unordered_map<std::string, double> thrhash = getThrHash(thrfile);
    std::string sfile = "/booleanfs2/sahoo/Data/BooleanLab/RxCovea/human-gpl570-tissue.txt";
    std::unordered_map<std::string, int> shash = getNeutrophilHash(sfile);

    std::ifstream file(efile);
    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return 1;
    }

    std::string head;
    std::getline(file, head);
    std::vector<std::string> headers;
    std::stringstream ss(head);
    std::string cell;
    while (std::getline(ss, cell, '\t')) { // Assuming tab as the delimiter
        headers.push_back(cell);
    }
    std::vector<int> selected;
    for (std::string key: headers) {
        if (shash.find(key) != shash.end()) {
            selected.push_back(shash[key]);
        } else {
            selected.push_back(0);
        }
    }

    std::string line;
    int index = 0;
    std::cout << headers[0] <<"\t" <<"Percentile" << std::endl;
    while (std::getline(file, line)) {
        std::vector<double> row_hi;
        std::vector<double> row_lo;
        std::stringstream ss(line);
        std::string key, name;
        std::getline(ss, key, '\t');
        std::getline(ss, name, '\t');
        double thr = thrhash[key];
        std::vector<double> values;

        int idx = 2;
        while (std::getline(ss, cell, '\t')) { // Assuming tab as the delimiter
            double value = std::stod(cell);
            if (selected[idx] == 1) {
                values.push_back(value);
            }
            if (value > thr) {
                row_hi.push_back(value);
            }
            else {
                row_lo.push_back(value);
            }
            idx += 1;
        }
        double m1 = calculateAverage(values);
        double percentile;
        if (m1 > thr) {
            percentile = percentileOfScore(row_hi, m1);
        } else {
            percentile = -percentileOfScore(row_lo, m1);
        }
        std::cout << key <<"\t" <<percentile << std::endl;
 
        index += 1;
        if ( (index % 1000) == 0) {
            std::cerr << index << std::endl;
        }
    }

    file.close();

    return 0;
}

