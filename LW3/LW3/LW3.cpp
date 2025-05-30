﻿#include <iostream>
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
    
    // for (int i = 0; i < n; i++) 
    // {
    //     double x, y;
    //     input >> x >> y;
    //     triangulation.AddVertex(x, y);
    // }

    input.close();
    triangulation.GenerateRandomPoints(6, 10.0f, 400.0f, 300.0f, 800, 600);
    triangulation.Triangulate();
    
    Draw draw(triangulation);
    draw.Print();

    return 0;
}