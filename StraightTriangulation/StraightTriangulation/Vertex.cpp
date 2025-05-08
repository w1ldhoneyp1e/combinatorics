#include "Vertex.h"

Vertex::Vertex(double x, double y, size_t idx) : x(x), y(y), index(idx) {}

bool Vertex::operator<(const Vertex& other) const
{
    return x < other.x || (x == other.x && y < other.y);
}