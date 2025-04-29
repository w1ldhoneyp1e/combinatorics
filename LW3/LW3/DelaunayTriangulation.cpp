#include <algorithm>
#include <cmath>
#include "DelaunayTriangulation.h"
#include "Draw.h"

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

    if (points.size() < 3) {
        return result;
    }

    Vertex* v1 = &points[0];
    Vertex* v2 = &points[1];
    Vertex* v3 = &points[2];

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

    Edge* baseEdge = FindBaseEdge(left, right);
    if (!baseEdge) return result;

    while (true) {
        Vertex* topCandidate = FindTopCandidate(baseEdge, result);
        if (!topCandidate) break;

        if (!IsDelaunayEdge(baseEdge, topCandidate, result)) 
        {
            FlipEdge(baseEdge, topCandidate, result);
        }
        else 
        {
            break;
        }
    }

    while (true) {
        Vertex* bottomCandidate = FindBottomCandidate(baseEdge, result);
        if (!bottomCandidate) break;

        if (!IsDelaunayEdge(baseEdge, bottomCandidate, result)) 
        {
            FlipEdge(baseEdge, bottomCandidate, result);
        }
        else {
            break;
        }
    }

    return result;
}

Edge* DelaunayTriangulation::FindBaseEdge(const std::vector<Face>& left, const std::vector<Face>& right) 
{
    Vertex* rightmostLeft = nullptr;
    double maxX = -std::numeric_limits<double>::infinity();

    for (const Face& face : left) {
        if (face.v1->x > maxX) { rightmostLeft = face.v1; maxX = face.v1->x; }
        if (face.v2->x > maxX) { rightmostLeft = face.v2; maxX = face.v2->x; }
        if (face.v3->x > maxX) { rightmostLeft = face.v3; maxX = face.v3->x; }
    }

    Vertex* leftmostRight = nullptr;
    double minX = std::numeric_limits<double>::infinity();

    for (const Face& face : right) {
        if (face.v1->x < minX) { leftmostRight = face.v1; minX = face.v1->x; }
        if (face.v2->x < minX) { leftmostRight = face.v2; minX = face.v2->x; }
        if (face.v3->x < minX) { leftmostRight = face.v3; minX = face.v3->x; }
    }

    return new Edge(rightmostLeft, leftmostRight);
}

bool DelaunayTriangulation::IsDelaunayEdge(Edge* edge, Vertex* v, const std::vector<Face>& triangulation) 
{
    for (const Face& face : triangulation) 
    {
        if ((face.e1 == edge || face.e2 == edge || face.e3 == edge) &&
            face.ContainsVertex(v)) 
        {
            return face.IsDelaunay(*v);
        }
    }
    return true;
}

void DelaunayTriangulation::FlipEdge(Edge* edge, Vertex* v, std::vector<Face>& triangulation) 
{
    Face* t1 = nullptr;
    Face* t2 = nullptr;
    size_t t1_idx = 0, t2_idx = 0;

    for (size_t i = 0; i < triangulation.size(); ++i) 
    {
        if (triangulation[i].e1 == edge || triangulation[i].e2 == edge || triangulation[i].e3 == edge) 
        {
            if (!t1) 
            {
                t1 = &triangulation[i];
                t1_idx = i;
            }
            else 
            {
                t2 = &triangulation[i];
                t2_idx = i;
                break;
            }
        }
    }

    if (!t1 || !t2) return;

    Vertex* v1 = edge->v1;
    Vertex* v2 = edge->v2;
    Vertex* v3 = v;
    Vertex* v4 = nullptr;

    if (!t2->ContainsVertex(v)) 
    {
        if (t2->v1 != v1 && t2->v1 != v2) v4 = t2->v1;
        else if (t2->v2 != v1 && t2->v2 != v2) v4 = t2->v2;
        else v4 = t2->v3;
    }

    Face newFace1(v3, v4, v1);
    Face newFace2(v3, v2, v4);

    triangulation[t1_idx] = newFace1;
    triangulation[t2_idx] = newFace2;

    edge->v1 = v3;
    edge->v2 = v4;
}

Vertex* DelaunayTriangulation::FindTopCandidate(Edge* edge, const std::vector<Face>& triangulation) 
{
    Vertex* bestCandidate = nullptr;
    double maxAngle = -std::numeric_limits<double>::infinity();

    for (const Face& face : triangulation) 
    {
        if (face.e1 == edge || face.e2 == edge || face.e3 == edge) 
        {
            Vertex* candidate = nullptr;
            if (!face.ContainsVertex(edge->v1) && !face.ContainsVertex(edge->v2)) 
            {
                candidate = face.v1;
            }
            else if (!face.ContainsVertex(edge->v1) && !face.ContainsVertex(edge->v2)) 
            {
                candidate = face.v2;
            }
            else 
            {
                candidate = face.v3;
            }

            double angle = CalculateAngle(edge->v1, edge->v2, candidate);
            if (angle > maxAngle) 
            {
                maxAngle = angle;
                bestCandidate = candidate;
            }
        }
    }

    return bestCandidate;
}

Vertex* DelaunayTriangulation::FindBottomCandidate(Edge* edge, const std::vector<Face>& triangulation) 
{
    Vertex* bestCandidate = nullptr;
    double minAngle = std::numeric_limits<double>::infinity();

    for (const Face& face : triangulation) 
    {
        if (face.e1 == edge || face.e2 == edge || face.e3 == edge) 
        {
            Vertex* candidate = nullptr;
            if (!face.ContainsVertex(edge->v1) && !face.ContainsVertex(edge->v2)) 
            {
                candidate = face.v1;
            }
            else if (!face.ContainsVertex(edge->v1) && !face.ContainsVertex(edge->v2)) 
            {
                candidate = face.v2;
            }
            else 
            {
                candidate = face.v3;
            }

            double angle = CalculateAngle(edge->v1, edge->v2, candidate);
            if (angle < minAngle) 
            {
                minAngle = angle;
                bestCandidate = candidate;
            }
        }
    }

    return bestCandidate;
}

double DelaunayTriangulation::CalculateAngle(Vertex* v1, Vertex* v2, Vertex* v3) 
{
    double dx1 = v2->x - v1->x;
    double dy1 = v2->y - v1->y;
    double dx2 = v3->x - v1->x;
    double dy2 = v3->y - v1->y;

    return std::atan2(dx1 * dy2 - dy1 * dx2, dx1 * dx2 + dy1 * dy2);
}

std::vector<Face> HandleFourPoints(std::vector<Vertex>& points) 
{
    std::vector<Face> result;
    
    Vertex* v1 = &points[0];
    Vertex* v2 = &points[1];
    Vertex* v3 = &points[2];
    Vertex* v4 = &points[3];
    
    result.emplace_back(v1, v2, v3);
    
    Face& firstTriangle = result[0];
    if (firstTriangle.IsDelaunay(*v4))
    {
        result.emplace_back(v2, v3, v4);
    } else 
    {
        result.emplace_back(v1, v3, v4);
        result.emplace_back(v1, v2, v4);
    }
    
    return result;
}

// void DelaunayTriangulation::AddPointOutsideTriangulation(Vertex& point, std::vector<Face>& triangulation) {
//     double minDist = std::numeric_limits<double>::infinity();
//     Face* closestFace = nullptr;
//     Vertex *v1 = nullptr, *v2 = nullptr;
    
//     for (Face& face : triangulation) {
//         double d1 = DistanceToEdge(point, face.v1, face.v2);
//         double d2 = DistanceToEdge(point, face.v2, face.v3);
//         double d3 = DistanceToEdge(point, face.v3, face.v1);
        
//         if (d1 < minDist) {
//             minDist = d1;
//             closestFace = &face;
//             v1 = face.v1;
//             v2 = face.v2;
//         }
//         if (d2 < minDist) {
//             minDist = d2;
//             closestFace = &face;
//             v1 = face.v2;
//             v2 = face.v3;
//         }
//         if (d3 < minDist) {
//             minDist = d3;
//             closestFace = &face;
//             v1 = face.v3;
//             v2 = face.v1;
//         }
//     }
    
//     Face newFace(&point, v1, v2);
//     triangulation.push_back(newFace);
    
//     LegalizeEdge(&point, newFace.e1, triangulation);
// }

// void DelaunayTriangulation::LegalizeEdge(Vertex* p, Edge* edge, std::vector<Face>& triangulation) {
//     Face* currentFace = nullptr;
//     Face* adjacentFace = nullptr;
    
//     for (Face& face : triangulation) {
//         if (face.ContainsVertex(p) && 
//             (face.e1 == edge || face.e2 == edge || face.e3 == edge)) {
//             currentFace = &face;
//             break;
//         }
//     }
    
//     if (!currentFace) return;
    
//     adjacentFace = FindAdjacentTriangle(*currentFace, edge->v1, edge->v2, triangulation);
//     if (!adjacentFace) return;
    
//     Vertex* oppositeVertex = adjacentFace->GetOppositeVertex(edge->v1, edge->v2);
    
//     if (!currentFace->IsDelaunay(*oppositeVertex)) {
//         Face newFace1(p, oppositeVertex, edge->v1);
//         Face newFace2(p, oppositeVertex, edge->v2);
        
//         *currentFace = newFace1;
//         *adjacentFace = newFace2;
        
//         Edge* edge1 = newFace1.e2; 
//         Edge* edge2 = newFace1.e3; 
//         Edge* edge3 = newFace2.e3; 
        
//         LegalizeEdge(p, edge2, triangulation);
//         LegalizeEdge(p, edge3, triangulation);
//     }
// }

// Face* DelaunayTriangulation::FindAdjacentTriangle(const Face& face, Vertex* v1, Vertex* v2, std::vector<Face>& triangulation) {
//     // Ищем треугольник, который содержит оба указанных vertex и не является исходным треугольником
//     for (Face& triangle : triangulation) {
//         if (&triangle != &face && 
//             triangle.ContainsVertex(v1) && 
//             triangle.ContainsVertex(v2)) {
//             return &triangle;
//         }
//     }
//     return nullptr;
// }

// void DelaunayTriangulation::AddPointToTriangulation(Vertex& point, std::vector<Face>& triangulation) {
//     bool pointFound = false;
    
//     for (size_t i = 0; i < triangulation.size(); i++) {
//         Face& face = triangulation[i];
//         if (face.ContainsPoint(point)) {
//             Face newFace1(&point, face.v1, face.v2);
//             Face newFace2(&point, face.v2, face.v3);
//             Face newFace3(&point, face.v3, face.v1);
            
//             triangulation[i] = newFace1;
//             triangulation.push_back(newFace2);
//             triangulation.push_back(newFace3);
            
//             LegalizeEdge(&point, newFace1.e1, triangulation);
//             LegalizeEdge(&point, newFace2.e1, triangulation);
//             LegalizeEdge(&point, newFace3.e1, triangulation);
//             pointFound = true;
//             break;
//         }
//     }
    
//     if (!pointFound) {
//         for (size_t i = 0; i < triangulation.size(); i++) {
//             Face& face = triangulation[i];
//             if (face.e1->ContainsPoint(point)) {
//                 auto [newFace1, newFace2] = face.SplitAtEdge(point, face.v1, face.v2);
//                 triangulation[i] = newFace1;
//                 triangulation.push_back(newFace2);
//                 return;
//             }
//             if (face.e2->ContainsPoint(point)) {
//                 auto [newFace1, newFace2] = face.SplitAtEdge(point, face.v2, face.v3);
//                 triangulation[i] = newFace1;
//                 triangulation.push_back(newFace2);
//                 return;
//             }
//             if (face.e3->ContainsPoint(point)) {
//                 auto [newFace1, newFace2] = face.SplitAtEdge(point, face.v3, face.v1);
//                 triangulation[i] = newFace1;
//                 triangulation.push_back(newFace2);
//                 return;
//             }
//         }
        
//         AddPointOutsideTriangulation(point, triangulation);
//     }
// }

std::vector<Face> DelaunayTriangulation::DivideAndConquer(std::vector<Vertex>& points) 
{
    int N = points.size();
    
    if (N == 3) 
    {
        return CreateBaseTriangulation(points);
    }
    
    if (N == 4) 
    {
        return HandleFourPoints(points);
    }
    
    // if (N == 5) {
    //     std::vector<Vertex> first(points.begin(), points.begin() + 3);
    //     auto result = CreateBaseTriangulation(first);
        
    //     for (int i = 3; i < 5; i++) {
    //         AddPointToTriangulation(points[i], result);
    //     }
    //     return result;
    // }
    
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
    
    int mid = N / 2;
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
