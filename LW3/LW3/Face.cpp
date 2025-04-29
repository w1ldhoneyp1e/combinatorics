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
    double ax = v1->x - p.x;
    double ay = v1->y - p.y;
    double bx = v2->x - p.x;
    double by = v2->y - p.y;
    double cx = v3->x - p.x;
    double cy = v3->y - p.y;

    double det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
        (bx * bx + by * by) * (ax * cy - cx * ay) +
        (cx * cx + cy * cy) * (ax * by - bx * ay);

    return det <= 0;
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