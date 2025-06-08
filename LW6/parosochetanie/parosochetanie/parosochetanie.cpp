#include <iostream>
#include <vector>
#include "HopcroftKarp.h"

using namespace std;

int main() 
{
    vector<pair<int, int>> edges;
    int size = 0;
    int u, v;

    while (cin >> u >> v) 
    {
        edges.emplace_back(u, v);
        size = max(size, max(u, v));
    }
    
    HopcroftKarpAlgorithm algorithm(size, edges);
    int matching = algorithm.Solve();
    
    cout << "Max matching: " << matching << endl;
    
    auto matchingPairs = algorithm.GetMatching();
    for (auto& pair : matchingPairs) 
    {
        cout << pair.first << " - " << pair.second << endl;
    }
    
    return 0;
}