#pragma once
#include "Vertex.h"
#include "Edge.h"
#include <vector>
#include <cmath>
#include <utility>
#include <set>
#include <algorithm>

class Face
{
public:
    Vertex* v1;
    Vertex* v2;
    Vertex* v3;

    using Faces = std::vector<Face>;

    Face(Vertex* vertex1, Vertex* vertex2, Vertex* vertex3) : v1(vertex1), v2(vertex2), v3(vertex3) {};

    bool operator==(const Face& other) const {
        std::set<Vertex*> s1 = {v1, v2, v3};
        std::set<Vertex*> s2 = {other.v1, other.v2, other.v3};
        return s1 == s2;
    }

    bool operator<(const Face& other) const {
        std::vector<Vertex*> t1 = {v1, v2, v3};
        std::vector<Vertex*> t2 = {other.v1, other.v2, other.v3};
        std::sort(t1.begin(), t1.end());
        std::sort(t2.begin(), t2.end());
        return t1 < t2;
    }
};

