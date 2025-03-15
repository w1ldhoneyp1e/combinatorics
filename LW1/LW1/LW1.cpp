#include <algorithm>
#include <iostream>
#include <vector>
#include <chrono>
#include <string>

size_t factorialFromTo(size_t from, size_t to) {
    size_t result = 1;
    for (int i = from; i <= to; i++) {
        result *= i;
    }
    return result;
}

std::string getValueAndMeasure(size_t milliseconds) {
    size_t duration = milliseconds;
    std::string measure;
    if (milliseconds < 1000) {
        measure = "ms";
    } else if (milliseconds < 60000) {
        duration = milliseconds / 1000;
        measure = "s";
    } else if (milliseconds < 3600000) {
        duration = milliseconds / 60000;
        measure = "min";
    } else if (milliseconds < 86400000) {
        duration = milliseconds / 3600000;
        measure = "h";
    } else if (milliseconds < 31536000000) {
        duration = milliseconds / 86400000;
        measure = "d";
    } else {
        duration = milliseconds / 31536000000;
        measure = "y";
    }

    return std::to_string(duration) + " " + measure;
}

size_t findTheAssignmentQuadraticProblemResult(std::vector<size_t> v, std::vector<std::vector<size_t>> paths, std::vector<std::vector<size_t>> mass) {
    size_t result = INT_MAX;
    auto begin = std::chrono::steady_clock::now();

    do {
        size_t curr = 0;
        for (size_t i = 0; i < v.size(); i++) {
            for (size_t j = 0; j < v.size(); j++) {
                curr += paths[i][j] * mass[v[i]][v[j]];
            }
        }
        result = std::min(result, curr);
    } while (std::next_permutation(v.begin(), v.end()));

    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Result: " << result << std::endl;

    return elapsed_ms.count();
}

int main()  
{
    std::vector<size_t> v10 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    
    std::vector<std::vector<size_t>> paths10 = {
        {0,  3,  4,  5,  2,  7,  3,  4,  6,  2},
        {3,  0,  5,  2,  4,  8,  5,  3,  7,  4},
        {4,  5,  0,  6,  3,  4,  2,  5,  8,  6},
        {5,  2,  6,  0,  5,  3,  4,  7,  4,  5},
        {2,  4,  3,  5,  0,  6,  5,  2,  3,  7},
        {7,  8,  4,  3,  6,  0,  3,  4,  5,  2},
        {3,  5,  2,  4,  5,  3,  0,  6,  2,  4},
        {4,  3,  5,  7,  2,  4,  6,  0,  5,  3},
        {6,  7,  8,  4,  3,  5,  2,  5,  0,  4},
        {2,  4,  6,  5,  7,  2,  4,  3,  4,  0}
    };

    std::vector<std::vector<size_t>> mass10 = {
        {0,  4,  8,  3,  7,  9,  5,  2,  6,  1},
        {4,  0,  5,  9,  2,  8,  6,  3,  7,  4},
        {8,  5,  0,  6,  4,  3,  9,  7,  1,  8},
        {3,  9,  6,  0,  8,  5,  2,  4,  7,  3},
        {7,  2,  4,  8,  0,  1,  6,  9,  5,  7},
        {9,  8,  3,  5,  1,  0,  4,  2,  8,  6},
        {5,  6,  9,  2,  6,  4,  0,  8,  3,  5},
        {2,  3,  7,  4,  9,  2,  8,  0,  1,  9},
        {6,  7,  1,  7,  5,  8,  3,  1,  0,  2},
        {1,  4,  8,  3,  7,  6,  5,  9,  2,  0}
    };

    std::vector<size_t> v12 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    
    std::vector<std::vector<size_t>> paths12 = {
        {0,  3,  4,  5,  2,  7,  3,  4,  6,  2,  5,  3,  4,  2,  6},
        {3,  0,  5,  2,  4,  8,  5,  3,  7,  4,  2,  6,  3,  5,  4},
        {4,  5,  0,  6,  3,  4,  2,  5,  8,  6,  4,  2,  5,  7,  3},
        {5,  2,  6,  0,  5,  3,  4,  7,  4,  5,  3,  4,  6,  2,  5},
        {2,  4,  3,  5,  0,  6,  5,  2,  3,  7,  5,  3,  4,  6,  2},
        {7,  8,  4,  3,  6,  0,  3,  4,  5,  2,  6,  5,  2,  4,  7},
        {3,  5,  2,  4,  5,  3,  0,  6,  2,  4,  7,  4,  5,  3,  6},
        {4,  3,  5,  7,  2,  4,  6,  0,  5,  3,  2,  7,  4,  5,  2},
        {6,  7,  8,  4,  3,  5,  2,  5,  0,  4,  3,  2,  6,  4,  5},
        {2,  4,  6,  5,  7,  2,  4,  3,  4,  0,  5,  6,  3,  7,  4},
        {5,  2,  4,  3,  5,  6,  7,  2,  3,  5,  0,  4,  5,  2,  6},
        {3,  6,  2,  4,  3,  5,  4,  7,  2,  6,  4,  0,  3,  5,  7},
        {4,  3,  5,  6,  4,  2,  5,  4,  6,  3,  5,  3,  0,  4,  2},
        {2,  5,  7,  2,  6,  4,  3,  5,  4,  7,  2,  5,  4,  0,  3},
        {6,  4,  3,  5,  2,  7,  6,  2,  5,  4,  6,  7,  2,  3,  0}
    };

    std::vector<std::vector<size_t>> mass12 = {
        {0,  4,  8,  3,  7,  9,  5,  2,  6,  1,  8,  3,  7,  4,  2},
        {4,  0,  5,  9,  2,  8,  6,  3,  7,  4,  5,  9,  1,  6,  3},
        {8,  5,  0,  6,  4,  3,  9,  7,  1,  8,  2,  4,  8,  5,  7},
        {3,  9,  6,  0,  8,  5,  2,  4,  7,  3,  6,  1,  9,  2,  4},
        {7,  2,  4,  8,  0,  1,  6,  9,  5,  7,  3,  8,  2,  7,  5},
        {9,  8,  3,  5,  1,  0,  4,  2,  8,  6,  7,  5,  3,  9,  1},
        {5,  6,  9,  2,  6,  4,  0,  8,  3,  5,  9,  2,  6,  1,  8},
        {2,  3,  7,  4,  9,  2,  8,  0,  1,  9,  4,  7,  5,  3,  6},
        {6,  7,  1,  7,  5,  8,  3,  1,  0,  2,  8,  6,  4,  8,  2},
        {1,  4,  8,  3,  7,  6,  5,  9,  2,  0,  3,  8,  7,  5,  9},
        {8,  5,  2,  6,  3,  7,  9,  4,  8,  3,  0,  5,  2,  7,  4},
        {3,  9,  4,  1,  8,  5,  2,  7,  6,  8,  5,  0,  9,  3,  1},
        {7,  1,  8,  9,  2,  3,  6,  5,  4,  7,  2,  9,  0,  6,  8},
        {4,  6,  5,  2,  7,  9,  1,  3,  8,  5,  7,  3,  6,  0,  5},
        {2,  3,  7,  4,  5,  1,  8,  6,  2,  9,  4,  1,  8,  5,  0}
    };

    auto ms = findTheAssignmentQuadraticProblemResult(v10, paths10, mass10);
    std::cout << "The time for 10 elements: " << getValueAndMeasure(ms) << std::endl;

    ms = findTheAssignmentQuadraticProblemResult(v12, paths12, mass12);
    std::cout << "The time for 12 elements: " << getValueAndMeasure(ms) << std::endl;
    std::cout << "Prediction for 20 elements: " << getValueAndMeasure(ms * factorialFromTo(v12.size(), 20)) << std::endl;
    std::cout << "Prediction for 50 elements: " << getValueAndMeasure(ms * factorialFromTo(v12.size(), 50)) << std::endl;

    std::cout << "--------------------------------" << std::endl;

    return 0;
}