#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <climits>
#include <set>
#include <string>

std::string getValueAndMeasure(size_t milliseconds) {
    auto duration = milliseconds;
    std::string measure;
    if (milliseconds < 1000) {
        measure = "ms";
    }
    else if (milliseconds < 60000) {
        duration = milliseconds / 1000;
        measure = "s";
    }
    else if (milliseconds < 3600000) {
        duration = milliseconds / 60000;
        measure = "min";
    }
    else if (milliseconds < 86400000) {
        duration = milliseconds / 3600000;
        measure = "h";
    }
    else {
        duration = milliseconds / 86400000;
        measure = "d";
    }

    return std::to_string(duration) + " " + measure;
}

std::vector<std::vector<std::vector<size_t>>> generateCombinations(const std::vector<std::vector<size_t>>& workers, size_t k) {
    std::vector<std::vector<std::vector<size_t>>> combinations;
    size_t n = workers.size();
    std::vector<bool> bitmask(n, false);
    
    std::fill(bitmask.begin(), bitmask.begin() + k, true);
    
    do {
        std::vector<std::vector<size_t>> combination;
        for (size_t i = 0; i < n; ++i) {
            if (bitmask[i]) {
                combination.push_back(workers[i]);
            }
        }
        combinations.push_back(combination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    
    return combinations;
}

bool canGroupOfWorkersHandleAllWork(const std::vector<std::vector<size_t>>& group, const std::vector<size_t>& work) {
    std::set<size_t> tasks;
    for (const auto& worker : group) {
        for (const auto& task : worker) {
            tasks.insert(task);
        }
    }
    return tasks.size() == work.size();
}

void findMinimalGroup(const std::vector<std::vector<size_t>>& workers, const std::vector<size_t>& work) {
    size_t amountOfAllWorkers = workers.size();
    
    for (size_t amountOfWorkers = 1; amountOfWorkers <= amountOfAllWorkers; ++amountOfWorkers) {
        auto combinations = generateCombinations(workers, amountOfWorkers);
        for (const auto& combination : combinations) {
            if (canGroupOfWorkersHandleAllWork(combination, work)) {
                std::cout << "Достаточно " << amountOfWorkers << " работников" << std::endl;
                return;
            }
        }
    }
}

bool NextCombinations(size_t size, std::vector<size_t>& state) {
    size_t k = state.size() - 1;
    size_t m = k;
    while (state[m] == size - k + m) m = m - 1;
    state[m] = state[m] + 1;
    for (size_t i = m + 1; i <= k; i++) state[i] = state[i - 1] + 1;
    return m != 0;
}

void findMinimalGroup_Alternative(const std::vector<std::vector<size_t>>& workers, const std::vector<size_t>& work) {
    size_t n = workers.size();

    for (size_t k = 1; k <= n; k++) {
        std::vector<size_t> a(k + 1);

        for (size_t i = 1; i <= k; i++) a[i] = i;

         do {
            for (size_t i = 1; i <= k; i++) std::cout << a[i] << " ";
            std::vector<std::vector<size_t>> currWorkers;
            for (size_t i = 1; i <= k; i++) currWorkers.push_back(workers[a[i] - 1]);
            if (canGroupOfWorkersHandleAllWork(currWorkers, work)) {
                std::cout << "Достаточно " << a.size() - 1 << " работников" << std::endl;
                return;
            }
            std::cout << std::endl;
         } while (NextCombinations(n, a));
    }
}

int main() {
    setlocale(LC_ALL, "RU");
    const std::vector<size_t> work = {1, 2, 3, 4, 5};
    
    const std::vector<std::vector<size_t>> workers = {
        {1, 4, 5},
        {2, 1},
        {1, 2},
        {3, 4},
        {1},
        {2},
        {3},
        {4},
        {5}
    };
    auto begin = std::chrono::steady_clock::now();

    // findMinimalGroup(workers, work);
    findMinimalGroup_Alternative(workers, work);

    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "The time for " << workers.size() << " workers: " << getValueAndMeasure(elapsed_ms.count()) << std::endl;

    return 0;
}        
