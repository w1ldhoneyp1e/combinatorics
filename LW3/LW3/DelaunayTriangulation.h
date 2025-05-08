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
    Face::Faces faces;

    Face::Faces CreateBaseTriangulation(std::vector<Vertex>& points);
    Face::Faces HandleFourPoints(std::vector<Vertex>& points);
    Face::Faces Merge(const Face::Faces& left, const Face::Faces& right);
    bool IsLowerPoint(Vertex* a, Vertex* b) const;
    Face::Faces DivideAndConquer(std::vector<Vertex>& points);
    Vertex* FindDelaunayNeighbor(Edge* baseLine, const std::set<Vertex*, Vertex::VertexPtrCompare>& vertices);
    void RemoveConflictingTriangles(Face::Faces& triangulation,
        Edge* baseLine,
        Vertex* newVertex,
        const std::set<Vertex*, 
        Vertex::VertexPtrCompare>& leftVertices
    );
    double CalculateBaseLineAngle(Vertex* v1, Vertex* v2, Vertex* p);
    void FindTangents(
        const Face::Faces& left, 
        const Face::Faces& right,
        Vertex*& P0, 
        Vertex*& P1, 
        Vertex*& P2, 
        Vertex*& P3
    );
    bool IsValidUpperTangent(Vertex* leftPoint, Vertex* rightPoint, const std::vector<Vertex*>& vertices);
    bool IsValidLowerTangent(Vertex* leftPoint, Vertex* rightPoint, const std::vector<Vertex*>& vertices);

public:
    const std::vector<Vertex>& GetVertices() const;
    const Face::Faces& GetFaces() const;
    void AddVertex(double x, double y);
    void Triangulate();
    void GenerateRandomPoints(
        int count, 
        float scale, 
        float offsetX, 
        float offsetY, 
        int windowWidth, 
        int windowHeight
    );
};

