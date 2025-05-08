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
private:
    std::vector<Vertex> vertices;
    Face::Faces faces;

    std::vector<Vertex*> ConvexHull(std::vector<Vertex*> vertexPtrs);
    bool InCircle(const Vertex& p, Vertex* a, Vertex* b, Vertex* c);
    double GetOrientation(const Vertex& O, const Vertex& A, const Vertex& B);
    bool IsTriangleUsed(const std::set<Face>& usedTriangles, Vertex* a, Vertex* b, Vertex* c);
    void AddTriangle(std::set<Face>& usedTriangles, Vertex* a, Vertex* b, Vertex* c);
    bool IsDelaunay(Vertex* a, Vertex* b, Vertex* c, const std::vector<Vertex*>& vertexPtrs);
    Vertex* FindBestPointForEdge(const Edge& edge, const std::vector<Vertex*>& vertexPtrs, const std::set<Face>& usedTriangles);

public:
    void Triangulate();
    void GenerateRandomPoints(int count, float scale, float offsetX, float offsetY,
        int windowWidth, int windowHeight);
    void AddVertex(double x, double y);
    Face::Faces GetFaces() const;
    std::vector<Vertex> GetVertices() const;
};
