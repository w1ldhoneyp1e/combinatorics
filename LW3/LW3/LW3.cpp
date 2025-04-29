#include <iostream>
#include <vector>
#include "DelaunayTriangulation.h"
#include "Draw.h"

int main() {
    DelaunayTriangulation dt;
    int n;
    double x, y;
    
    std::cin >> n;
    
    for (int i = 0; i < n; i++) {
        std::cin >> x >> y;
        dt.AddVertex(x, y);
    }
    
    dt.Triangulate();
    
    Draw draw(dt);
    draw.Print();
    
    return 0;
}