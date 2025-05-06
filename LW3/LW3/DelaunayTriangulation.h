#pragma once
#include "Vertex.h"
#include "Edge.h"
#include <vector>
#include <iostream>
#include "Face.h"
#include <cmath>
#include <set>
#include <random>

class DelaunayTriangulation 
{
private:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    std::vector<Face> CreateBaseTriangulation(std::vector<Vertex>& points);
    std::vector<Face> HandleFourPoints(std::vector<Vertex>& points);
    std::vector<Face> Merge(const std::vector<Face>& left, const std::vector<Face>& right);
    bool IsLowerPoint(Vertex* a, Vertex* b) const;
    std::vector<Face> DivideAndConquer(std::vector<Vertex>& points);
    Vertex* FindDelaunayNeighbor(Edge* baseLine, const std::vector<Face>& triangulation);
    void RemoveConflictingTriangles(std::vector<Face>& triangulation, Edge* baseLine, Vertex* newVertex);
    bool IsLeftOfLine(Vertex* a, Vertex* b, Vertex* c);
    bool IsHigherTangent(Vertex* a, Vertex* b, Vertex* c);
    bool IsLowerTangent(Vertex* a, Vertex* b, Vertex* c);
    double CalculateCircumradius(Vertex* a, Vertex* b, Vertex* c);
    double CalculateBaseLineAngle(Vertex* v1, Vertex* v2, Vertex* p);
    bool EdgesIntersect(Vertex* a1, Vertex* a2, Vertex* b1, Vertex* b2) const;

public:
    const std::vector<Vertex>& GetVertices() const;
    const std::vector<Face>& GetFaces() const;
    void AddVertex(double x, double y);
    void Triangulate();
    void GenerateRandomPoints(int count, float scale, float offsetX, float offsetY, 
                             int windowWidth, int windowHeight);
};

