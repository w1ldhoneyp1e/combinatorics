#include "Vertex.h"

Vertex::Vertex(double x = 0, double y = 0, int idx = 0) : x(x), y(y), index(idx) {}

bool Vertex::operator<(const Vertex& other) const {
    return x < other.x || (x == other.x && y < other.y);
}