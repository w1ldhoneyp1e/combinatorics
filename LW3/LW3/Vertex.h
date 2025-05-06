#pragma once

class Vertex 
{
public:
    double x, y;
    size_t index;

    Vertex(double x, double y, size_t idx);

    bool operator<(const Vertex& other) const;

    struct VertexPtrCompare 
    {
        bool operator()(const Vertex* a, const Vertex* b) const {
            if (a->x != b->x) return a->x < b->x;
            return a->y < b->y;
        }
    };
};
