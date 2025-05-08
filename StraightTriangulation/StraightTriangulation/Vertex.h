#pragma once
#include <cmath>

class Vertex
{
public:
    double x, y;
    size_t index;

    Vertex(double x, double y, size_t idx) : x(x), y(y), index(idx) {}

    bool operator<(const Vertex& other) const
    {
        return x < other.x || (x == other.x && y < other.y);
    }

    bool operator==(const Vertex& other) const 
    {
        return std::abs(x - other.x) < 1e-10 &&
            std::abs(y - other.y) < 1e-10;
    }
};

