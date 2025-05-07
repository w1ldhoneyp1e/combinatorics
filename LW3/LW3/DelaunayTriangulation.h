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

    struct EdgeInfo {
        Vertex* v1;
        Vertex* v2;
        
        EdgeInfo(Vertex* v1, Vertex* v2) : v1(v1), v2(v2) {
            if (Vertex::VertexPtrCompare()(v2, v1)) {
                std::swap(this->v1, this->v2);
            }
        }

        bool operator==(const EdgeInfo& other) const {
            return (v1 == other.v1 && v2 == other.v2) ||
                   (v1 == other.v2 && v2 == other.v1);
        }
    };

    struct EdgeInfoCompare {
        bool operator()(const EdgeInfo& a, const EdgeInfo& b) const {
            Vertex::VertexPtrCompare comp;
            if (a.v1 != b.v1) return comp(a.v1, b.v1);
            return comp(a.v2, b.v2);
        }
    };

    std::vector<Face> CreateBaseTriangulation(std::vector<Vertex>& points);
    std::vector<Face> HandleFourPoints(std::vector<Vertex>& points);
    std::vector<Face> Merge(const std::vector<Face>& left, const std::vector<Face>& right);
    bool IsLowerPoint(Vertex* a, Vertex* b) const;
    std::vector<Face> DivideAndConquer(std::vector<Vertex>& points);
    Vertex* FindDelaunayNeighbor(Edge* baseLine, const std::set<Vertex*, Vertex::VertexPtrCompare>& vertices);
    void RemoveConflictingTriangles(std::vector<Face>& triangulation,
        Edge* baseLine,
        Vertex* newVertex,
        const std::set<Vertex*, Vertex::VertexPtrCompare>& leftVertices);
    bool IsLeftOfLine(Vertex* a, Vertex* b, Vertex* c);
    bool IsHigherTangent(Vertex* a, Vertex* b, Vertex* c);
    bool IsLowerTangent(Vertex* a, Vertex* b, Vertex* c);
    double CalculateCircumradius(Vertex* a, Vertex* b, Vertex* c);
    double CalculateBaseLineAngle(Vertex* v1, Vertex* v2, Vertex* p);
    bool EdgesIntersect(Vertex* a1, Vertex* a2, Vertex* b1, Vertex* b2) const;
    void FindTangents(const std::vector<Face>& left, const std::vector<Face>& right,
        Vertex*& P0, Vertex*& P1, Vertex*& P2, Vertex*& P3);
    bool IsValidUpperTangent(Vertex* leftPoint, Vertex* rightPoint,
        const std::vector<Vertex*>& vertices) const;
    bool IsValidLowerTangent(Vertex* leftPoint, Vertex* rightPoint,
        const std::vector<Vertex*>& vertices) const;
        
public:
    const std::vector<Vertex>& GetVertices() const;
    const std::vector<Face>& GetFaces() const;
    void AddVertex(double x, double y);
    void Triangulate();
    void GenerateRandomPoints(int count, float scale, float offsetX, float offsetY, 
                             int windowWidth, int windowHeight);
};

