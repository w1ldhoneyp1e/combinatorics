#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <set>
#include <algorithm>

class Graph {
private:
    struct Edge 
    {
        int from;
        int to;
        Edge(int f, int t) : from(f), to(t) {}
    };

    std::vector<std::vector<bool>> adjacencyMatrix;
    std::vector<std::pair<int, int>> times;
    std::set<int> articulationPoints;
    std::vector<Edge> edgeList;
    std::vector<bool> visited;
    std::vector<int> up;
    int time;

    void readEdgeList(const std::string& filename) 
    {
        std::ifstream file(filename);
        if (!file.is_open()) 
        {
            throw std::runtime_error("Failed to open file");
        }

        std::string line;
        if (!std::getline(file, line)) 
        {
            throw std::runtime_error("File is empty");
        }
        
        try 
        {
            int vertexCount = std::stoi(line);
            if (vertexCount < 0) {
                throw std::runtime_error("Number of vertices cannot be negative");
            }
            adjacencyMatrix.resize(vertexCount, std::vector<bool>(vertexCount, false));
        } 
        catch (const std::invalid_argument&) 
        {
            throw std::runtime_error("Invalid number of vertices in the first line");
        }

        while (std::getline(file, line)) 
        {
            std::istringstream iss(line);
            int v1, v2;
            
            if (!(iss >> v1 >> v2)) 
            {
                continue;
            }
            
            if (v1 < 0 || v1 >= adjacencyMatrix.size() || 
                v2 < 0 || v2 >= adjacencyMatrix.size()) 
            {
                throw std::runtime_error("Vertex number is out of bounds");
            }
            
            edgeList.emplace_back(v1, v2);
        }
    }

    void convertEdgeListToMatrix() 
    {
        for (const auto& edge : edgeList) 
        {
            adjacencyMatrix[edge.from][edge.to] = true;
            adjacencyMatrix[edge.to][edge.from] = true;
        }
    }

    void findArticulationPointsDFS(int v, int parent, int& timer) 
    {
        visited[v] = true;
        times[v].first = up[v] = timer++;
        int children = 0;

        for (size_t to = 0; to < adjacencyMatrix[v].size(); to++) 
        {
            if (!adjacencyMatrix[v][to]) continue;
            
            if (to == parent) continue;
            
            if (visited[to]) 
            {
                up[v] = std::min(up[v], times[to].first);
            } 
            else 
            {
                findArticulationPointsDFS(to, v, timer);
                up[v] = std::min(up[v], up[to]);
                
                if (up[to] >= times[v].first && parent != -1) 
                {
                    articulationPoints.insert(v);
                }
                children++;
            }
        }
        
        times[v].second = timer++;
        
        if (parent == -1 && children > 1) 
        {
            articulationPoints.insert(v);
        }
    }

public:
    explicit Graph(const std::string& filename) 
    {
        readEdgeList(filename);
        convertEdgeListToMatrix();
    }

    std::vector<Edge> getEdgeList() 
    {
        return edgeList;
    }

    std::vector<Edge> matrixToEdgeList() 
    {
        std::vector<Edge> result;
        for (size_t i = 0; i < adjacencyMatrix.size(); i++) 
        {
            for (size_t j = i; j < adjacencyMatrix[i].size(); j++) 
            {
                if (adjacencyMatrix[i][j]) 
                {
                    result.emplace_back(i, j);
                }
            }
        }
        return result;
    }

    std::set<int> findArticulationPoints() 
    {
        int n = adjacencyMatrix.size();
        times.assign(n, {-1, -1});
        up.assign(n, -1);
        visited.assign(n, false);
        articulationPoints.clear();
        
        int timer = 0;
        for (int i = 0; i < n; i++) 
        {
            if (!visited[i]) 
            {
                findArticulationPointsDFS(i, -1, timer);
            }
        }
        
        return articulationPoints;
    }

    void printArticulationPoints(const std::string& outputFile) 
    {
        try 
        {
            std::ofstream outFile(outputFile);
            outFile << "Searching for articulation points:";
            outFile << "\nArticulation points:";
            if (articulationPoints.empty()) 
            {
                outFile << " none\n";
            } else {
                for (int point : articulationPoints) 
                {
                    outFile << " " << point;
                }
                outFile << "\n";
            }
        }
        catch (const std::exception& e) 
        {
            throw std::runtime_error("Failed to write to output file: " + std::string(e.what()));
        }
    }
};

int main(int argc, char* argv[]) 
{    
    if (argc != 3) 
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        return 1;
    }
    
    try 
    {
        Graph g(argv[1]);
        auto articulationPoints = g.findArticulationPoints();
        g.printArticulationPoints(argv[2]);
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}