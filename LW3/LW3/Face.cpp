#include "Face.h"
#include <iostream>

Face::Face(Vertex* v1, Vertex* v2, Vertex* v3) : v1(v1), v2(v2), v3(v3) 
{
    e1 = new Edge(v1, v2);
    e2 = new Edge(v2, v3);
    e3 = new Edge(v3, v1);
}

bool Face::ContainsVertex(const Vertex* v) const 
{
    return v == v1 || v == v2 || v == v3;
}

bool Face::IsDelaunay(const Vertex& p) const 
{
    return !PointInCircumcircle(v1, v2, v3, const_cast<Vertex*>(&p));
}

std::pair<Face, Face> Face::SplitAtEdge(Vertex& point, Vertex* edgeV1, Vertex* edgeV2) const 
{
    Vertex* oppositeVertex = GetOppositeVertex(edgeV1, edgeV2);
    
    Face newFace1(&point, edgeV1, oppositeVertex);
    Face newFace2(&point, edgeV2, oppositeVertex);
    
    return std::make_pair(newFace1, newFace2);
}

Vertex* Face::GetOppositeVertex(Vertex* v1, Vertex* v2) const 
{
    if (this->v1 != v1 && this->v1 != v2) return this->v1;
    if (this->v2 != v1 && this->v2 != v2) return this->v2;
    return this->v3;
}

bool Face::PointInCircumcircle(Vertex* p1, Vertex* p2, Vertex* p3, Vertex* pTest)
{
    if (!p1 || !p2 || !p3 || !pTest) return false;

    double orientation = Edge::OrientationCheck(p1, p2, p3);
    if (std::abs(orientation) < 1e-12)
    {
        return false;
    }

    double ax = p1->x; double ay = p1->y;
    double bx = p2->x; double by = p2->y;
    double cx = p3->x; double cy = p3->y;
    double dx = pTest->x; double dy = pTest->y;

    double adx = ax - dx; double ady = ay - dy;
    double bdx = bx - dx; double bdy = by - dy;
    double cdx = cx - dx; double cdy = cy - dy;

    double asq = adx * adx + ady * ady;
    double bsq = bdx * bdx + bdy * bdy;
    double csq = cdx * cdx + cdy * cdy;

    double det = adx * (bdy * csq - cdy * bsq) -
                 ady * (bdx * csq - cdx * bsq) +
                 asq * (bdx * cdy - cdx * bdy);

    double epsilon = 1e-9;

    if (orientation > 0)
    {
        return det > epsilon;
    }
    else
    {
        return det < -epsilon;
    }
}

double Face::CalculateCircumradius(Vertex* a, Vertex* b, Vertex* c) 
{
    double ax = a->x, ay = a->y;
    double bx = b->x, by = b->y;
    double cx = c->x, cy = c->y;

    double D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (std::abs(D) < 1e-10) return std::numeric_limits<double>::infinity();

    double x = ((ax * ax + ay * ay) * (by - cy) + 
                (bx * bx + by * by) * (cy - ay) + 
                (cx * cx + cy * cy) * (ay - by)) / D;
    double y = ((ax * ax + ay * ay) * (cx - bx) + 
                (bx * bx + by * by) * (ax - cx) + 
                (cx * cx + cy * cy) * (bx - ax)) / D;

    return std::hypot(x - ax, y - ay);
}