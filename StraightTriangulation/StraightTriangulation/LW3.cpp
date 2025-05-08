#include <iostream>
#include "StraightTriangulation.h"
#include "Draw.h"

int main()
{
    StraightTriangulation triangulation;
    
    triangulation.GenerateRandomPoints(20, 10.0f, 400.0f, 300.0f, 800, 600);
    
    triangulation.Triangulate();
    
    Draw draw(triangulation);
    draw.Print();
    
    return 0;
}