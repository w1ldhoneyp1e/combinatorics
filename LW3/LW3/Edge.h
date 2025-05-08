#pragma once
#include "Vertex.h"
#include <cmath>
#include <vector>

class Edge 
{
public:
    Vertex* v1;
    Vertex* v2;

    Edge(Vertex* v1, Vertex* v2);

    void OrderVertices();

    bool operator==(const Edge& other) const;

    bool ContainsPoint(const Vertex& p) const;

    struct EdgeCompare {
        bool operator()(const Edge& a, const Edge& b) const {
            Vertex::VertexPtrCompare comp;
            if (a.v1 != b.v1) return comp(a.v1, b.v1);
            return comp(a.v2, b.v2);
        }
    };
    
    static bool EdgesIntersect(Vertex* a1, Vertex* a2, Vertex* b1, Vertex* b2);
    
    static double OrientationCheck(Vertex* a, Vertex* b, Vertex* c);
};

