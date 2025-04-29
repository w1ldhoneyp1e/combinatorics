#include <iostream>
#include <vector>
#include "DelaunayTriangulation.h"
#include "Draw.h"
#include <fstream>
#include <iostream>

int main() 
{
    std::ifstream input("input.txt");
    if (!input.is_open()) 
    {
        std::cout << "Cannot open input.txt" << std::endl;
        return 1;
    }

    int n;
    input >> n;

    DelaunayTriangulation triangulation;
    
    for (int i = 0; i < n; i++) 
    {
        double x, y;
        input >> x >> y;
        triangulation.AddVertex(x, y);
    }

    input.close();
    triangulation.Triangulate();
    
    Draw draw(triangulation);
    draw.Print();

    return 0;
}