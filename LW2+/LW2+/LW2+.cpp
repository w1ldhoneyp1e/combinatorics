#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <windows.h>

class GraphException : public std::runtime_error {
public:
    explicit GraphException(const std::string& message) : std::runtime_error(message) {}
};

struct Vertex {
    int x, y;
    Vertex(int _x = 0, int _y = 0) : x(_x), y(_y) {}
};

struct Edge {
    int from, to;
    Edge(int _from = 0, int _to = 0) : from(_from), to(_to) {}
};

class GeometricGraph {
private:
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    int vertexCount;

public:
    GeometricGraph() : vertexCount(0) {}

    void readFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw GraphException("Ошибка открытия файла: " + filename);
        }

        file >> vertexCount;
        vertices.resize(vertexCount);

        char ch;
        int x1, y1, x2, y2;
        while (
            file >> 
            ch >> x1 >> ch >> y1 >> ch >>
            ch >> x2 >> ch >> y2 >> ch
        ) {
            int index1 = -1, index2 = -1;
            
            for (int i = 0; i < vertexCount; ++i) {
                if (vertices[i].x == x1 && vertices[i].y == y1) {
                    index1 = i;
                    break;
                }
                if (vertices[i].x == 0 && vertices[i].y == 0) {
                    vertices[i].x = x1;
                    vertices[i].y = y1;
                    index1 = i;
                    break;
                }
            }

            for (int i = 0; i < vertexCount; ++i) {
                if (vertices[i].x == x2 && vertices[i].y == y2) {
                    index2 = i;
                    break;
                }
                if (vertices[i].x == 0 && vertices[i].y == 0) {
                    vertices[i].x = x2;
                    vertices[i].y = y2;
                    index2 = i;
                    break;
                }
            }

            if (index1 >= 0 && index2 >= 0) {
                edges.push_back(Edge(index1, index2));
            }
        }
    }

    int calculateDistance(int v1, int v2) const {
        return std::abs(vertices[v1].x - vertices[v2].x) + 
               std::abs(vertices[v1].y - vertices[v2].y);
    }

    void printGraph() const {
        if (vertexCount == 0) {
            throw GraphException("Граф пуст");
        }

        std::cout << "Вершины:\n";
        for (int i = 0; i < vertexCount; ++i) {
            std::cout << i << ": (" << vertices[i].x << "," << vertices[i].y << ")\n";
        }
        
        std::cout << "\nРёбра:\n";
        for (const auto& edge : edges) {
            std::cout << edge.from << " -> " << edge.to << "\n";
        }
    }

    void printAdjacencyMatrix() const {
        if (vertexCount == 0) {
            throw GraphException("Граф пуст");
        }

        std::cout << "\nМатрица смежности с длинами путей:\n";
        std::cout << "   ";
        
        for (int i = 0; i < vertexCount; ++i) {
            std::cout << i << "  ";
        }
        std::cout << "\n";

        std::vector<std::vector<int>> adjacencyMatrix(vertexCount, std::vector<int>(vertexCount, 0));
        
        for (const auto& edge : edges) {
            adjacencyMatrix[edge.from][edge.to] = calculateDistance(edge.from, edge.to);
        }

        for (int i = 0; i < vertexCount; ++i) {
            std::cout << i << "  ";
            for (int j = 0; j < vertexCount; ++j) {
                std::cout << adjacencyMatrix[i][j] << "  ";
            }
            std::cout << "\n";
        }
    }
};

int main() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale(LC_ALL, "RU");
    try {
        GeometricGraph graph;
        graph.readFromFile("graph.txt");
        graph.printGraph();
        graph.printAdjacencyMatrix();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}