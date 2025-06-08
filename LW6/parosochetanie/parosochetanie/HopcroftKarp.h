#pragma once

#include <vector>
#include <queue>

class HopcroftKarpAlgorithm 
{
private:
    std::vector<std::vector<int>> adjacencyList;
    std::vector<int> pairU, pairV, dist;
    int size;
    
    const int nil = 0;
    const int inf = 1e9;

    bool BuildLayers();
    bool FindPath(int u);

public:
    HopcroftKarpAlgorithm(int size, const std::vector<std::pair<int, int>>& edges);
    
    int Solve();
    std::vector<std::pair<int, int>> GetMatching();
}; 