#pragma once
#include "Vertex.h"
#include "Edge.h"

class Face 
{
public:
    Vertex* v1, * v2, * v3;
    Edge* e1, * e2, * e3;

    Face(Vertex* v1, Vertex* v2, Vertex* v3);

    bool ContainsVertex(const Vertex* v) const;

    bool IsDelaunay(const Vertex& p) const;

    std::pair<Face, Face> SplitAtEdge(Vertex& point, Vertex* edgeV1, Vertex* edgeV2) const;

    bool ContainsPoint(const Vertex& p) const;

    Vertex* GetOppositeVertex(Vertex* v1, Vertex* v2) const;
};

