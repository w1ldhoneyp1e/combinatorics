#pragma once
#include "Vertex.h"
#include <cmath>

class Edge 
{
public:
    Vertex* v1;
    Vertex* v2;

    Edge(Vertex* v1, Vertex* v2);

    bool operator==(const Edge& other) const;

    bool ContainsPoint(const Vertex& p) const;
};

