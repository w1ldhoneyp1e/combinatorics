#include "HopcroftKarp.h"

HopcroftKarpAlgorithm::HopcroftKarpAlgorithm(int size, const std::vector<std::pair<int, int>>& edges) : size(size) 
{
    adjacencyList.assign(size + 1, std::vector<int>());
    for (auto& edge : edges) 
    {
        adjacencyList[edge.first].push_back(edge.second);
    }
    pairU.assign(size + 1, nil);
    pairV.assign(size + 1, nil);
    dist.assign(size + 1, 0);
}

bool HopcroftKarpAlgorithm::BuildLayers() 
{
    std::queue<int> bfsQueue;
    for (int u = 1; u <= size; ++u) 
    {
        if (pairU[u] == nil) 
        {
            dist[u] = 0;
            bfsQueue.push(u);
        } 
        else 
        {
            dist[u] = inf;
        }
    }
    dist[nil] = inf;

    while (!bfsQueue.empty()) 
    {
        int u = bfsQueue.front();
        bfsQueue.pop();

        if (dist[u] == inf) 
        {
            continue;
        }

        for (int v : adjacencyList[u]) 
        {
            if (dist[pairV[v]] != inf) 
            {
                continue;
            }
            dist[pairV[v]] = dist[u] + 1;
            bfsQueue.push(pairV[v]);
        }
    }
    
    return dist[nil] != inf;
}

bool HopcroftKarpAlgorithm::FindPath(int u) 
{
    if (u != nil) 
    {
        for (int v : adjacencyList[u]) 
        {
            if (dist[pairV[v]] == dist[u] + 1) 
            {
                if (FindPath(pairV[v])) 
                {
                    pairU[u] = v;
                    pairV[v] = u;
                    return true;
                }
            }
        }
        dist[u] = inf;

        return false;
    }

    return true;
}

int HopcroftKarpAlgorithm::Solve() 
{
    int matching = 0;
    while (BuildLayers()) 
    {
        for (int u = 1; u <= size; ++u) 
        {
            if (pairU[u] == nil && FindPath(u))
                ++matching;
        }
    }
    
    return matching;
}

std::vector<std::pair<int, int>> HopcroftKarpAlgorithm::GetMatching() 
{
    std::vector<std::pair<int, int>> result;
    for (int u = 1; u <= size; ++u) 
    {
        if (pairU[u] != nil) 
        {
            result.emplace_back(u, pairU[u]);
        }
    }
    return result;
} 