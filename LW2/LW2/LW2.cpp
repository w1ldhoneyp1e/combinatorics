#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>

class Graph {
private:
    struct Edge {
        int from;
        int to;
        Edge(int f, int t) : from(f), to(t) {}
    };

    std::vector<std::vector<bool>> adjacencyMatrix;
    std::vector<std::pair<int, int>> times;
    std::vector<Edge> edgeList;
    int time;

    const std::string RED = "\033[31m";
    const std::string GREEN = "\033[32m";
    const std::string RESET = "\033[0m";

    void initDFS() {
        times.resize(adjacencyMatrix.size(), {-1, -1});
        time = 0;
    }

    void printMatrixLine(std::ostream& out, size_t i, int width, bool useColor = true) const {
        const std::string divider = useColor 
            ? width == 4 
                ? "   " 
                : "  " 
            : width == 4 
                ? "  " 
                : " ";
        out << std::setw(width) << i << divider;
        for (size_t j = 0; j < adjacencyMatrix[i].size(); j++) {
            if (adjacencyMatrix[i][j]) {
                if (useColor) {
                    out << std::setw(width-1) << GREEN << "1" << RESET << divider;
                } else {
                    out << std::setw(width-1) << "1" << divider;
                }
            } else {
                if (useColor) {
                    out << std::setw(width-1) << RED << "0" << RESET << divider;
                } else {
                    out << std::setw(width-1) << "0" << divider;
                }
            }
        }
        out << "\n";
    }

    void printTimesLine(std::ostream& out, size_t i) const {
        out << std::setw(7) << i << " | ";
        out << std::setw(4) << times[i].first << " | ";
        out << std::setw(5) << times[i].second << "\n";
    }

    void printAdjacencyMatrixToStream(std::ostream& out, bool useColor) const {
        if (adjacencyMatrix.empty()) {
            out << "Граф пуст" << std::endl;
            return;
        }

        int width = 3;
        if (adjacencyMatrix.size() > 9) {
            width = 4;
        }

        out << "Матрица смежности:\n";
        out << std::setw(width) << " ";
        for (size_t i = 0; i < adjacencyMatrix.size(); i++) {
            out << std::setw(width) << i;
        }
        out << "\n";
        for (size_t i = 0; i < adjacencyMatrix.size(); i++) {
            printMatrixLine(out, i, width, useColor);
        }
    }

    void printDFSTimesToStream(std::ostream& out) const {
        if (times.empty()) {
            out << "DFS еще не был выполнен" << std::endl;
            return;
        }

        out << "Временные метки DFS:\n";
        out << "Вершина | Вход | Выход\n";
        out << "---------------------\n";
        for (size_t i = 0; i < times.size(); i++) {
            printTimesLine(out, i);
        }
    }

    void readEdgeList(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Не удалось открыть файл");
        }

        std::string line;
        if (!std::getline(file, line)) {
            throw std::runtime_error("Файл пуст");
        }
        
        try {
            int vertexCount = std::stoi(line);
            if (vertexCount < 0) {
                throw std::runtime_error("Количество вершин не может быть отрицательным");
            }
            adjacencyMatrix.resize(vertexCount, std::vector<bool>(vertexCount, false));
        } catch (const std::invalid_argument&) {
            throw std::runtime_error("Некорректное количество вершин в первой строке");
        }

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int v1, v2;
            
            if (!(iss >> v1 >> v2)) {
                continue;
            }
            
            if (v1 < 0 || v1 >= adjacencyMatrix.size() || 
                v2 < 0 || v2 >= adjacencyMatrix.size()) {
                throw std::runtime_error("Номер вершины выходит за допустимые пределы");
            }
            
            edgeList.emplace_back(v1, v2);
        }
    }

    void convertEdgeListToMatrix() {
        for (const auto& edge : edgeList) {
            adjacencyMatrix[edge.from][edge.to] = true;
            adjacencyMatrix[edge.to][edge.from] = true;
        }
    }

public:
    explicit Graph(const std::string& filename) {
        readEdgeList(filename);
        convertEdgeListToMatrix();
    }

    std::vector<Edge> getEdgeList() const {
        return edgeList;
    }

    std::vector<Edge> matrixToEdgeList() const {
        std::vector<Edge> result;
        for (size_t i = 0; i < adjacencyMatrix.size(); i++) {
            for (size_t j = i; j < adjacencyMatrix[i].size(); j++) {
                if (adjacencyMatrix[i][j]) {
                    result.emplace_back(i, j);
                }
            }
        }
        return result;
    }

    void DFS() {
        initDFS();
        std::vector<bool> visited(adjacencyMatrix.size(), false);
        
        for (size_t v = 0; v < adjacencyMatrix.size(); v++) {
            if (!visited[v]) {
                DFSRecursive(v, visited);
            }
        }
    }

    void printAdjacencyMatrix() const {
        printAdjacencyMatrixToStream(std::cout, true);

        std::ofstream outFile("output.txt");
        if (outFile.is_open()) {
            printAdjacencyMatrixToStream(outFile, false);
        }
    }

    void printDFSTimes() const {
        printDFSTimesToStream(std::cout);

        std::ofstream outFile("output.txt", std::ios::app);
        if (outFile.is_open()) {
            outFile << "\n";
            printDFSTimesToStream(outFile);
        }
    }

private:
    void DFSRecursive(int v, std::vector<bool>& visited) {
        visited[v] = true;
        times[v].first = ++time;

        for (size_t i = 0; i < adjacencyMatrix[v].size(); i++) {
            if (adjacencyMatrix[v][i] && !visited[i]) {
                DFSRecursive(i, visited);
            }
        }

        times[v].second = ++time;
    }
};

int main() {
    setlocale(LC_ALL, "RU");
    try {
        Graph g("input.txt");
        g.printAdjacencyMatrix();
        g.DFS();
        g.printDFSTimes();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}