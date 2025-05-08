#pragma once
#include "Vertex.h"
#include "Edge.h"
#include <vector>
#include <cmath>
#include <utility>
#include <set>

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

    bool operator==(const Face& other) const {
        std::set<Vertex*> s1 = {v1, v2, v3};
        std::set<Vertex*> s2 = {other.v1, other.v2, other.v3};
        return s1 == s2;
    }

    bool ContainsVertex(const Vertex* v) const;

    bool IsDelaunay(const Vertex& p) const;

    std::pair<Face, Face> SplitAtEdge(Vertex& point, Vertex* edgeV1, Vertex* edgeV2) const;

    Vertex* GetOppositeVertex(Vertex* v1, Vertex* v2) const;
};

