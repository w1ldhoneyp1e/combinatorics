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

    bool operator<(const Edge& other) const {
        if (*v1 < *other.v1) return true;
        if (*other.v1 < *v1) return false;
        return *v2 < *other.v2;
    }

    bool ContainsPoint(const Vertex& p) const;
};

