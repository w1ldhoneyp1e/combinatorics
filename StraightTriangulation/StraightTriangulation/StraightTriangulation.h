#pragma once
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include <vector>
#include <set>
#include <map>
#include <random>
#include <utility>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <stack>

class StraightTriangulation
{
public:

private:
    std::vector<Vertex> vertices;
    Face::Faces faces;

    std::vector<Vertex*> ConvexHull(std::vector<Vertex*> vertexPtrs);
    double Cross(const Vertex& O, const Vertex& A, const Vertex& B);
    bool InCircle(const Vertex& p, Vertex* a, Vertex* b, Vertex* c);

public:
    void Triangulate();
    void GenerateRandomPoints(int count, float scale, float offsetX, float offsetY,
        int windowWidth, int windowHeight);
    void AddVertex(double x, double y);
    Face::Faces GetFaces() const;
    std::vector<Vertex> GetVertices() const;
};
