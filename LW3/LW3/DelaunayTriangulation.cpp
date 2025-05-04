#include <algorithm>
#include <cmath>
#include <iostream>
#include "DelaunayTriangulation.h"
#include "Draw.h"

const double PI = 4 * std::atan(1);

const std::vector<Vertex>& DelaunayTriangulation::GetVertices() const 
{
    return vertices;
}

const std::vector<Face>& DelaunayTriangulation::GetFaces() const 
{
    return faces;
}

std::vector<Face> DelaunayTriangulation::CreateBaseTriangulation(std::vector<Vertex>& points) 
{
    std::vector<Face> result;

    if (points.size() < 3) 
    {
        return result;
    }

    Vertex* v1 = nullptr;
    Vertex* v2 = nullptr;
    Vertex* v3 = nullptr;

    for (auto& vertex : vertices)
    {
        if (vertex.x == points[0].x && vertex.y == points[0].y) v1 = &vertex;
        if (vertex.x == points[1].x && vertex.y == points[1].y) v2 = &vertex;
        if (vertex.x == points[2].x && vertex.y == points[2].y) v3 = &vertex;
    }

    if (!v1 || !v2 || !v3) {
        std::cout << "Error: Could not find vertices in the main vector" << std::endl;
        return result;
    }

    double orientation = (v2->x - v1->x) * (v3->y - v1->y) -
        (v3->x - v1->x) * (v2->y - v1->y);

    if (orientation < 0) 
    {
        std::swap(v2, v3);
    }

    result.emplace_back(v1, v2, v3);
    return result;
}

std::vector<Face> DelaunayTriangulation::Merge(const std::vector<Face>& left, const std::vector<Face>& right) 
{
    std::vector<Face> result = left;
    result.insert(result.end(), right.begin(), right.end());

    Vertex* P0 = nullptr;
    Vertex* P1 = nullptr;
    Vertex* P2 = nullptr;
    Vertex* P3 = nullptr;

    double minY_left = std::numeric_limits<double>::infinity();
    double minY_right = std::numeric_limits<double>::infinity();

    for (const Face& face : left) {
        if (face.v1->y < minY_left) { P0 = face.v1; minY_left = face.v1->y; }
        if (face.v2->y < minY_left) { P0 = face.v2; minY_left = face.v2->y; }
        if (face.v3->y < minY_left) { P0 = face.v3; minY_left = face.v3->y; }
    }

    for (const Face& face : right) {
        if (face.v1->y < minY_right) { P1 = face.v1; minY_right = face.v1->y; }
        if (face.v2->y < minY_right) { P1 = face.v2; minY_right = face.v2->y; }
        if (face.v3->y < minY_right) { P1 = face.v3; minY_right = face.v3->y; }
    }

    double maxY_left = -std::numeric_limits<double>::infinity();
    double maxY_right = -std::numeric_limits<double>::infinity();

    for (const Face& face : left) {
        if (face.v1->y > maxY_left) { P2 = face.v1; maxY_left = face.v1->y; }
        if (face.v2->y > maxY_left) { P2 = face.v2; maxY_left = face.v2->y; }
        if (face.v3->y > maxY_left) { P2 = face.v3; maxY_left = face.v3->y; }
    }

    for (const Face& face : right) {
        if (face.v1->y > maxY_right) { P3 = face.v1; maxY_right = face.v1->y; }
        if (face.v2->y > maxY_right) { P3 = face.v2; maxY_right = face.v2->y; }
        if (face.v3->y > maxY_right) { P3 = face.v3; maxY_right = face.v3->y; }
    }

    std::cout << "P0: " << P0->x << " " << P0->y << std::endl;
    std::cout << "P1: " << P1->x << " " << P1->y << std::endl;
    std::cout << "P2: " << P2->x << " " << P2->y << std::endl;
    std::cout << "P3: " << P3->x << " " << P3->y << std::endl;
    Edge* baseLine = new Edge(P0, P1);
    
    while (true) {
        Vertex* delaunayNeighbor = FindDelaunayNeighbor(baseLine, result);
        std::cout << "delaunayNeighbor: " << delaunayNeighbor->x << " " << delaunayNeighbor->y << std::endl;
        if (!delaunayNeighbor) break;

        Face newFace(baseLine->v1, baseLine->v2, delaunayNeighbor);
        result.push_back(newFace);

        if ((baseLine->v1 == P2 && baseLine->v2 == P3) ||
            (baseLine->v1 == P3 && baseLine->v2 == P2)) {
            break;
        }

        Vertex* oldV1 = baseLine->v1;
        Vertex* oldV2 = baseLine->v2;

        if (IsLowerPoint(delaunayNeighbor, baseLine->v1) && IsLowerPoint(delaunayNeighbor, baseLine->v2)) 
        {
            if (IsLowerPoint(baseLine->v1, baseLine->v2)) {
                baseLine->v2 = delaunayNeighbor;
            } else {
                baseLine->v1 = delaunayNeighbor;
            }
        }

        if (oldV1 == baseLine->v1 && oldV2 == baseLine->v2) {
                break;
        }
    }

    delete baseLine;
    return result;
}

bool DelaunayTriangulation::IsLowerPoint(Vertex* a, Vertex* b) const {
    return a->y > b->y;
}

Vertex* DelaunayTriangulation::FindDelaunayNeighbor(Edge* baseLine, const std::vector<Face>& triangulation) {
    Vertex* bestVertex = nullptr;
    double maxAngle = -std::numeric_limits<double>::infinity();

    for (const Face& face : triangulation) {
        for (Vertex* v : {face.v1, face.v2, face.v3}) {
            if (v != baseLine->v1 && v != baseLine->v2) {
                double crossProduct = (baseLine->v2->x - baseLine->v1->x) * (v->y - baseLine->v1->y) -
                                    (baseLine->v2->y - baseLine->v1->y) * (v->x - baseLine->v1->x);
                
                if (crossProduct > 0) {
                    double angle = CalculateBaseLineAngle(baseLine->v1, baseLine->v2, v);
                    
                    if (angle > maxAngle) {
                        maxAngle = angle;
                        bestVertex = v;
                    }
                }
            }
        }
    }

    return bestVertex;
}

double DelaunayTriangulation::CalculateBaseLineAngle(Vertex* v1, Vertex* v2, Vertex* p) {
    double v1x = v1->x - p->x;
    double v1y = v1->y - p->y;
    
    double v2x = v2->x - p->x;
    double v2y = v2->y - p->y;
    
    double dot = v1x * v2x + v1y * v2y;
    double len1 = std::sqrt(v1x * v1x + v1y * v1y);
    double len2 = std::sqrt(v2x * v2x + v2y * v2y);
    
    return std::acos(dot / (len1 * len2));
}

void DelaunayTriangulation::RemoveConflictingTriangles(std::vector<Face>& triangulation, 
                                                      Edge* baseLine, Vertex* newVertex) {
    triangulation.erase(
        std::remove_if(triangulation.begin(), triangulation.end(),
            [baseLine, newVertex](const Face& face) {
                return face.ContainsVertex(baseLine->v1) && 
                       face.ContainsVertex(baseLine->v2) &&
                       face.ContainsVertex(newVertex);
            }
        ),
        triangulation.end()
    );
}

bool DelaunayTriangulation::IsLeftOfLine(Vertex* a, Vertex* b, Vertex* c) {
    return ((b->x - a->x) * (c->y - a->y) - (b->y - a->y) * (c->x - a->x)) > 0;
}

bool DelaunayTriangulation::IsHigherTangent(Vertex* a, Vertex* b, Vertex* c) {
    return IsLeftOfLine(a, b, c);
}

bool DelaunayTriangulation::IsLowerTangent(Vertex* a, Vertex* b, Vertex* c) {
    return !IsLeftOfLine(a, b, c);
}

double DelaunayTriangulation::CalculateCircumradius(Vertex* a, Vertex* b, Vertex* c) {
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

std::vector<Face> DelaunayTriangulation::HandleFourPoints(std::vector<Vertex>& points) 
{
    std::vector<Face> result;
    
    Vertex* v1 = nullptr;
    Vertex* v2 = nullptr;
    Vertex* v3 = nullptr;
    Vertex* v4 = nullptr;

    for (auto& vertex : vertices)
    {
        if (std::abs(vertex.x - points[0].x) < 1e-10 && std::abs(vertex.y - points[0].y) < 1e-10) 
            v1 = &vertex;
        if (std::abs(vertex.x - points[1].x) < 1e-10 && std::abs(vertex.y - points[1].y) < 1e-10) 
            v2 = &vertex;
        if (std::abs(vertex.x - points[2].x) < 1e-10 && std::abs(vertex.y - points[2].y) < 1e-10) 
            v3 = &vertex;
        if (std::abs(vertex.x - points[3].x) < 1e-10 && std::abs(vertex.y - points[3].y) < 1e-10) 
            v4 = &vertex;
    }

    if (!v1 || !v2 || !v3 || !v4) 
    {
        std::cout << "Error: Could not find vertices in the main vector" << std::endl;
        return result;
    }
    
    result.emplace_back(v1, v2, v3);
    
    Face& firstTriangle = result[0];
    if (firstTriangle.IsDelaunay(*v4))
    {
        result.emplace_back(v2, v3, v4);
    } 
    else 
    {
        result.emplace_back(v1, v3, v4);
        result.emplace_back(v1, v2, v4);
    }
    
    return result;
}

std::vector<Face> DelaunayTriangulation::DivideAndConquer(std::vector<Vertex>& points) 
{
    size_t N = points.size();
    
    if (N == 3) 
    {
        return CreateBaseTriangulation(points);
    }

    if (N == 4) 
    {
        return HandleFourPoints(points);
    }
    
    if (N == 8) 
    {
        std::vector<Vertex> left(points.begin(), points.begin() + 4);
        std::vector<Vertex> right(points.begin() + 4, points.end());
        
        auto leftTri = DivideAndConquer(left);
        auto rightTri = DivideAndConquer(right);
        
        return Merge(leftTri, rightTri);
    }
    
    if (N < 12) 
    {
        std::vector<Vertex> left(points.begin(), points.begin() + 3);
        std::vector<Vertex> right(points.begin() + 3, points.end());
        
        auto leftTri = DivideAndConquer(left);
        auto rightTri = DivideAndConquer(right);
        
        return Merge(leftTri, rightTri);
    }
    
    size_t mid = N / 2;
    std::vector<Vertex> left(points.begin(), points.begin() + mid);
    std::vector<Vertex> right(points.begin() + mid, points.end());

    auto leftTri = DivideAndConquer(left);
    auto rightTri = DivideAndConquer(right);

    return Merge(leftTri, rightTri);
}

void DelaunayTriangulation::AddVertex(double x, double y) 
{
    vertices.emplace_back(x, y, vertices.size());
}

void DelaunayTriangulation::Triangulate() 
{
    std::sort(vertices.begin(), vertices.end());
    faces = DivideAndConquer(vertices);
}
