#pragma once
#include "Vertex.h"
#include "Edge.h"
#include <cmath>
#include <utility>

class Face 
{
public:
    Vertex* v1;
    Vertex* v2;
    Vertex* v3;
    Edge* e1;
    Edge* e2;
    Edge* e3;

    using Faces = std::vector<Face>;

    Face(Vertex* v1, Vertex* v2, Vertex* v3);

    bool ContainsVertex(const Vertex* v) const;

    bool IsDelaunay(const Vertex& p) const;

    std::pair<Face, Face> SplitAtEdge(Vertex& point, Vertex* edgeV1, Vertex* edgeV2) const;

    Vertex* GetOppositeVertex(Vertex* v1, Vertex* v2) const;
    
    static bool PointInCircumcircle(Vertex* a, Vertex* b, Vertex* c, Vertex* d);
    
    static double CalculateCircumradius(Vertex* a, Vertex* b, Vertex* c);
};

