#pragma once
#include "Vertex.h"
#include "Edge.h"
#include <vector>
#include "Face.h"

class DelaunayTriangulation 
{
private:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    std::vector<Face> CreateBaseTriangulation(std::vector<Vertex>& points);
    std::vector<Face> Merge(const std::vector<Face>& left, const std::vector<Face>& right);
    Edge* FindBaseEdge(const std::vector<Face>& left, const std::vector<Face>& right);
    bool IsDelaunayEdge(Edge* edge, Vertex* v, const std::vector<Face>& triangulation);
    void FlipEdge(Edge* edge, Vertex* v, std::vector<Face>& triangulation);
    Vertex* FindTopCandidate(Edge* edge, const std::vector<Face>& triangulation);
    Vertex* FindBottomCandidate(Edge* edge, const std::vector<Face>& triangulation);
    double CalculateAngle(Vertex* v1, Vertex* v2, Vertex* v3);
    std::vector<Face> DivideAndConquer(std::vector<Vertex>& points);
    bool IsPointInTriangle(const Vertex& p, const Face& face);
    void LegalizeEdge(Vertex* p, Edge* edge, std::vector<Face>& triangulation);
    Vertex* GetOppositeVertex(const Face& face, Vertex* v1, Vertex* v2);
    Face* FindAdjacentTriangle(const Face& face, Vertex* v1, Vertex* v2, std::vector<Face>& triangulation);

public:
    void AddVertex(double x, double y);

    void Triangulate();

    const std::vector<Vertex>& GetVertices() const;
    const std::vector<Face>& GetFaces() const;
};

