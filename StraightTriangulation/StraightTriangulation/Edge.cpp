#include "Edge.h"

Edge::Edge(Vertex* v1, Vertex* v2) : v1(v1), v2(v2) {}

bool Edge::operator==(const Edge& other) const
{
    return (v1 == other.v1 && v2 == other.v2) || (v1 == other.v2 && v2 == other.v1);
}

bool Edge::ContainsPoint(const Vertex& p) const
{
    double crossProduct = (p.y - v1->y) * (v2->x - v1->x) -
        (p.x - v1->x) * (v2->y - v1->y);

    if (std::abs(crossProduct) > 1e-10) return false;

    if (v1->x != v2->x)
    {
        return (v1->x <= p.x && p.x <= v2->x) || (v2->x <= p.x && p.x <= v1->x);
    }
    else
    {
        return (v1->y <= p.y && p.y <= v2->y) || (v2->y <= p.y && p.y <= v1->y);
    }
}