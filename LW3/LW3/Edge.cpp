#include "Edge.h"
#include <iostream>

Edge::Edge(Vertex* v1, Vertex* v2) : v1(v1), v2(v2) 
{
    OrderVertices();
}

void Edge::OrderVertices() 
{
    if (Vertex::VertexPtrCompare()(v2, v1)) {
        std::swap(v1, v2);
    }
}

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

bool Edge::EdgesIntersect(Vertex* a1, Vertex* a2, Vertex* b1, Vertex* b2)
{
    if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) {
        return false;
    }

    double d1 = ((b2->x - b1->x) * (a1->y - b1->y) - (b2->y - b1->y) * (a1->x - b1->x));
    double d2 = ((b2->x - b1->x) * (a2->y - b1->y) - (b2->y - b1->y) * (a2->x - b1->x));
    double d3 = ((a2->x - a1->x) * (b1->y - a1->y) - (a2->y - a1->y) * (b1->x - a1->x));
    double d4 = ((a2->x - a1->x) * (b2->y - a1->y) - (a2->y - a1->y) * (b2->x - a1->x));

    if (std::abs(d1) < 1e-10) return false;
    if (std::abs(d2) < 1e-10) return false;
    if (std::abs(d3) < 1e-10) return false;
    if (std::abs(d4) < 1e-10) return false;

    return (d1 * d2 < 0) && (d3 * d4 < 0);
}

double Edge::OrientationCheck(Vertex* a, Vertex* b, Vertex* c)
{
    return (b->x - a->x) * (c->y - a->y) - (b->y - a->y) * (c->x - a->x);
}
