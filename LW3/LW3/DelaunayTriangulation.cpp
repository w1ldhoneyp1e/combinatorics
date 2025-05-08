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
    std::cout << "\n=== Starting Merge ===\n";
    std::cout << "Left triangulation has " << left.size() << " faces\n";
    std::cout << "Right triangulation has " << right.size() << " faces\n";

    std::vector<Face> result = left;
    result.insert(result.end(), right.begin(), right.end());

    Vertex *P0 = nullptr, *P1 = nullptr, *P2 = nullptr, *P3 = nullptr;
    FindTangents(left, right, P0, P1, P2, P3);

    Edge* baseLine = new Edge(P0, P1);
    int iterationCount = 0;

    std::set<Vertex*, Vertex::VertexPtrCompare> leftVertices;
    std::set<Vertex*, Vertex::VertexPtrCompare> rightVertices;
    
    for (const Face& face : left) {
        leftVertices.insert(face.v1);
        leftVertices.insert(face.v2);
        leftVertices.insert(face.v3);
    }
    for (const Face& face : right) {
        rightVertices.insert(face.v1);
        rightVertices.insert(face.v2);
        rightVertices.insert(face.v3);
    }

    while (true) {
        iterationCount++;
        std::cout << "\n--- Iteration " << iterationCount << " ---";
        std::cout << "\nCurrent baseline: [" << baseLine->v1->index << "](" << baseLine->v1->x << "," << baseLine->v1->y 
                  << ") -> [" << baseLine->v2->index << "](" << baseLine->v2->x << "," << baseLine->v2->y << ")" << std::endl;
        
        std::set<Vertex*, Vertex::VertexPtrCompare> allVertices;
        allVertices.insert(leftVertices.begin(), leftVertices.end());
        allVertices.insert(rightVertices.begin(), rightVertices.end());
        
        Vertex* delaunayNeighbor = FindDelaunayNeighbor(baseLine, allVertices);
        if (!delaunayNeighbor) {
            std::cout << "No Delaunay neighbor found, stopping" << std::endl;
            break;
        }

        std::cout << "Found Delaunay neighbor: [" << delaunayNeighbor->index << "](" 
                  << delaunayNeighbor->x << "," << delaunayNeighbor->y << ")" << std::endl;

        RemoveConflictingTriangles(result, baseLine, delaunayNeighbor, leftVertices);

        Face newFace(baseLine->v1, baseLine->v2, delaunayNeighbor);
        result.push_back(newFace);
        std::cout << "Added new triangle: [" << newFace.v1->index << "](" 
                  << newFace.v1->x << "," << newFace.v1->y << ") - ["
                  << newFace.v2->index << "](" << newFace.v2->x << "," << newFace.v2->y << ") - ["
                  << newFace.v3->index << "](" << newFace.v3->x << "," << newFace.v3->y << ")" << std::endl;

        if ((baseLine->v1 == P2 && baseLine->v2 == P3) ||
            (baseLine->v1 == P3 && baseLine->v2 == P2)) {
            std::cout << "Reached lower base line, stopping" << std::endl;
                break;
        }

        Vertex* oldV1 = baseLine->v1;
        Vertex* oldV2 = baseLine->v2;

        bool isInLeft = leftVertices.find(delaunayNeighbor) != leftVertices.end();
        bool isInRight = rightVertices.find(delaunayNeighbor) != rightVertices.end();

        if (isInLeft) {
            std::cout << "Setting new baseline: replacing left vertex v1[" << baseLine->v1->index 
                      << "] with [" << delaunayNeighbor->index << "] (from left part)" << std::endl;
            baseLine->v1 = delaunayNeighbor;
        } else if (isInRight) {
            std::cout << "Setting new baseline: replacing right vertex v2[" << baseLine->v2->index 
                      << "] with [" << delaunayNeighbor->index << "] (from right part)" << std::endl;
            baseLine->v2 = delaunayNeighbor;
        } else {
            std::cout << "Warning: Delaunay neighbor not found in either part!" << std::endl;
        }

        if (oldV1 == baseLine->v1 && oldV2 == baseLine->v2) {
            std::cout << "Baseline didn't change, stopping" << std::endl;
            break;
        }

        std::cout << "Updated baseline to: [" << baseLine->v1->index << "](" << baseLine->v1->x << "," << baseLine->v1->y 
                  << ") -> [" << baseLine->v2->index << "](" << baseLine->v2->x << "," << baseLine->v2->y << ")" << std::endl;
    }

    std::cout << "\n=== Merge completed ===\n";
    std::cout << "Final triangulation has " << result.size() << " faces\n";

    delete baseLine;
    return result;
}

bool DelaunayTriangulation::IsLowerPoint(Vertex* a, Vertex* b) const 
{
    return a->y > b->y;
}

Vertex* DelaunayTriangulation::FindDelaunayNeighbor(Edge* baseLine, 
    const std::set<Vertex*, Vertex::VertexPtrCompare>& allVertices) 
{
    Vertex* bestVertex = nullptr;
    double maxAngle = -std::numeric_limits<double>::infinity();

    std::cout << "\nSearching for Delaunay neighbor for baseline: ["
              << baseLine->v1->index << "]" << " -> ["
              << baseLine->v2->index << "]";

    std::cout << "Checking " << allVertices.size() << " vertices\n";

    int candidateCount = 0;
    for (Vertex* v : allVertices) {
        if (v != baseLine->v1 && v != baseLine->v2) {
            double direction = (baseLine->v2->x - baseLine->v1->x) * (v->y - baseLine->v1->y) -
                             (baseLine->v2->y - baseLine->v1->y) * (v->x - baseLine->v1->x);
            
            std::cout << "Checking vertex [" << v->index << "] (" << v->x << "," << v->y << "):\n";
            std::cout << "  Direction: " << (direction > 0 ? "right" : "left") << "\n";
            
            if (direction > 0) {
                candidateCount++;
                double angle = CalculateBaseLineAngle(baseLine->v1, baseLine->v2, v);
                std::cout << "  Valid candidate! Angle: " << angle * 180.0 / PI << " degrees\n";
                
                if (angle > maxAngle) {
                    maxAngle = angle;
                    bestVertex = v;
                    std::cout << "  New best vertex found!\n";
                }
            } else {
                std::cout << "  Skipped: point is not on the correct side of baseline\n";
            }
        }
    }

    std::cout << "\nChecked " << allVertices.size() << " vertices, found " 
              << candidateCount << " valid candidates\n";
    
    if (bestVertex) {
        std::cout << "Selected best vertex: [" << bestVertex->index << "](" << bestVertex->x << "," << bestVertex->y 
                  << ") with angle " << maxAngle * 180.0 / PI << " degrees\n";
    } else {
        std::cout << "No valid vertex found!\n";
    }

    return bestVertex;
}

double DelaunayTriangulation::CalculateBaseLineAngle(Vertex* v1, Vertex* v2, Vertex* p) 
{
    double v1x = v1->x - p->x;
    double v1y = v1->y - p->y;
    
    double v2x = v2->x - p->x;
    double v2y = v2->y - p->y;
    
    double dot = v1x * v2x + v1y * v2y;
    double len1 = std::sqrt(v1x * v1x + v1y * v1y);
    double len2 = std::sqrt(v2x * v2x + v2y * v2y);
    
    return std::acos(dot / (len1 * len2));
}

void DelaunayTriangulation::RemoveConflictingTriangles(
    std::vector<Face>& triangulation, 
    Edge* baseLine, 
    Vertex* newVertex,
    const std::set<Vertex*, Vertex::VertexPtrCompare>& leftVertices)
{
    bool isInLeft = leftVertices.find(newVertex) != leftVertices.end();
    Vertex* baseVertex = isInLeft ? baseLine->v2 : baseLine->v1;
    
    std::cout << "\nChecking for conflicting edges with new edge: ["
              << baseVertex->index << "](" << baseVertex->x << "," << baseVertex->y << ") -> ["
              << newVertex->index << "](" << newVertex->x << "," << newVertex->y << ")" << std::endl;

    std::vector<Face> trianglesToRemove;
    std::set<EdgeInfo, EdgeInfoCompare> conflictingEdges;
    
    for (const Face& face : triangulation) {
        bool hasConflict = 
            EdgesIntersect(face.v1, face.v2, baseVertex, newVertex) ||
            EdgesIntersect(face.v2, face.v3, baseVertex, newVertex) ||
            EdgesIntersect(face.v3, face.v1, baseVertex, newVertex);
            
        if (hasConflict) {
            trianglesToRemove.push_back(face);
            
            if (EdgesIntersect(face.v1, face.v2, baseVertex, newVertex)) {
                std::cout << "Edge [" << face.v1->index << "]->[" << face.v2->index 
                          << "] intersects with new edge" << std::endl;
                conflictingEdges.insert(EdgeInfo(face.v1, face.v2));
            }
            if (EdgesIntersect(face.v2, face.v3, baseVertex, newVertex)) {
                std::cout << "Edge [" << face.v2->index << "]->[" << face.v3->index 
                          << "] intersects with new edge" << std::endl;
                conflictingEdges.insert(EdgeInfo(face.v2, face.v3));
            }
            if (EdgesIntersect(face.v3, face.v1, baseVertex, newVertex)) {
                std::cout << "Edge [" << face.v3->index << "]->[" << face.v1->index 
                          << "] intersects with new edge" << std::endl;
                conflictingEdges.insert(EdgeInfo(face.v3, face.v1));
            }
        }
    }
    
    for (const Face& faceToRemove : trianglesToRemove) {
        auto it = std::find_if(triangulation.begin(), triangulation.end(), 
            [&faceToRemove](const Face& f) {
                return (f.v1 == faceToRemove.v1 && f.v2 == faceToRemove.v2 && f.v3 == faceToRemove.v3) ||
                       (f.v1 == faceToRemove.v1 && f.v2 == faceToRemove.v3 && f.v3 == faceToRemove.v2) ||
                       (f.v1 == faceToRemove.v2 && f.v2 == faceToRemove.v1 && f.v3 == faceToRemove.v3) ||
                       (f.v1 == faceToRemove.v2 && f.v2 == faceToRemove.v3 && f.v3 == faceToRemove.v1) ||
                       (f.v1 == faceToRemove.v3 && f.v2 == faceToRemove.v1 && f.v3 == faceToRemove.v2) ||
                       (f.v1 == faceToRemove.v3 && f.v2 == faceToRemove.v2 && f.v3 == faceToRemove.v1);
            });
            
        if (it != triangulation.end()) {
            triangulation.erase(it);
        }
    }
    
    std::set<EdgeInfo, EdgeInfoCompare> boundaryEdges;
    for (const Face& face : trianglesToRemove) {
        EdgeInfo e1(face.v1, face.v2);
        EdgeInfo e2(face.v2, face.v3);
        EdgeInfo e3(face.v3, face.v1);
        
        if (conflictingEdges.find(e1) == conflictingEdges.end()) {
            boundaryEdges.insert(e1);
        }
        if (conflictingEdges.find(e2) == conflictingEdges.end()) {
            boundaryEdges.insert(e2);
        }
        if (conflictingEdges.find(e3) == conflictingEdges.end()) {
            boundaryEdges.insert(e3);
        }
    }
    
    for (const EdgeInfo& edge : boundaryEdges) {
        triangulation.emplace_back(edge.v1, edge.v2, newVertex);
    }
    
    std::cout << "Remaining triangles after edge removal: " << triangulation.size() << std::endl;
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

double DelaunayTriangulation::Orientation(Vertex* p1, Vertex* p2, Vertex* p3) const
{
    if (!p1 || !p2 || !p3) return 0.0;
    return (p2->x - p1->x) * (p3->y - p1->y) - (p2->y - p1->y) * (p3->x - p1->x);
}

bool DelaunayTriangulation::PointInCircumcircle(Vertex* p1, Vertex* p2, Vertex* p3, Vertex* pTest) const
{
    if (!p1 || !p2 || !p3 || !pTest) return false;

    double ax = p1->x; double ay = p1->y;
    double bx = p2->x; double by = p2->y;
    double cx = p3->x; double cy = p3->y;
    double dx = pTest->x; double dy = pTest->y;

    double orientationP1P2P3 = Orientation(p1, p2, p3);
    if (std::abs(orientationP1P2P3) < 1e-12)
    {
        return false;
    }

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

    if (orientationP1P2P3 > 0)
    {
        return det > epsilon;
    }
    else
    {
        return det < -epsilon;
    }
}

std::vector<Face> DelaunayTriangulation::HandleFourPoints(std::vector<Vertex>& pointsParam)
{
    std::vector<Face> result;
    if (pointsParam.size() != 4)
    {
        std::cerr << "HandleFourPoints был вызван с " << pointsParam.size() << " точками, когда ожидалось 4." << std::endl;
        return result;
    }

    Vertex* pOrig[4];
    for (int i = 0; i < 4; ++i)
    {
        bool found = false;
        for (auto& masterV : this->vertices)
        {
            if (std::abs(masterV.x - pointsParam[i].x) < 1e-9 && std::abs(masterV.y - pointsParam[i].y) < 1e-9)
            {
                pOrig[i] = &masterV;
                found = true;
                break;
            }
        }
        if (!found)
        {
            std::cerr << "HandleFourPoints: Не удалось найти главную вершину для точки с индексом " << pointsParam[i].index << std::endl;
            return result;
        }
    }

    Vertex* pSorted[4];
    std::copy(pOrig, pOrig + 4, pSorted);

    int minIdx = 0;
    for (int i = 1; i < 4; ++i)
    {
        if (pSorted[i]->y < pSorted[minIdx]->y || (pSorted[i]->y == pSorted[minIdx]->y && pSorted[i]->x < pSorted[minIdx]->x))
        {
            minIdx = i;
        }
    }
    std::swap(pSorted[0], pSorted[minIdx]);

    std::sort(pSorted + 1, pSorted + 4, [&](Vertex* a, Vertex* b)
    {
        double angleA = std::atan2(a->y - pSorted[0]->y, a->x - pSorted[0]->x);
        double angleB = std::atan2(b->y - pSorted[0]->y, b->x - pSorted[0]->x);
        if (std::abs(angleA - angleB) > 1e-9)
        {
            return angleA < angleB;
        }
        return ((a->x - pSorted[0]->x) * (a->x - pSorted[0]->x) + (a->y - pSorted[0]->y) * (a->y - pSorted[0]->y) <
                (b->x - pSorted[0]->x) * (b->x - pSorted[0]->x) + (b->y - pSorted[0]->y) * (b->y - pSorted[0]->y));
    });

    std::cout << "HandleFourPoints, упорядоченные вершины: "
              << "p0[" << pSorted[0]->index << "], p1[" << pSorted[1]->index
              << "], p2[" << pSorted[2]->index << "], p3[" << pSorted[3]->index << "]" << std::endl;

    auto createOrientedFace = [&](Vertex* v1, Vertex* v2, Vertex* v3)
    {
        if (Orientation(v1, v2, v3) < 0)
        {
            return Face(v1, v3, v2);
        }
        return Face(v1, v2, v3);
    };

    bool diag02IsDelaunay = false;
    Face tempFace1 = createOrientedFace(pSorted[0], pSorted[1], pSorted[2]);
    Face tempFace2 = createOrientedFace(pSorted[0], pSorted[2], pSorted[3]);

    if (!PointInCircumcircle(tempFace1.v1, tempFace1.v2, tempFace1.v3, pSorted[3]) &&
        !PointInCircumcircle(tempFace2.v1, tempFace2.v2, tempFace2.v3, pSorted[1]))
    {
        diag02IsDelaunay = true;
    }

    if (diag02IsDelaunay)
    {
        std::cout << "  HandleFourPoints: Выбрана Делоне диагональ " << pSorted[0]->index << "-" << pSorted[2]->index << std::endl;
        result.push_back(tempFace1);
        result.push_back(tempFace2);
    }
    else
    {
        std::cout << "  HandleFourPoints: Выбрана Делоне диагональ " << pSorted[1]->index << "-" << pSorted[3]->index << std::endl;
        result.push_back(createOrientedFace(pSorted[1], pSorted[0], pSorted[3]));
        result.push_back(createOrientedFace(pSorted[1], pSorted[3], pSorted[2]));
    }

    if (result.size() != 2)
    {
         std::cout << "  ПРЕДУПРЕЖДЕНИЕ HandleFourPoints: Создано " << result.size() << " граней вместо 2." << std::endl;
    }
    return result;
}

std::vector<Face> DelaunayTriangulation::DivideAndConquer(std::vector<Vertex>& points) 
{
    size_t N = points.size();
    
    std::cout << "\n=== Dividing " << N << " points ===\n";
    std::cout << "Points in set: ";
    for (const auto& p : points) {
        std::cout << "[" << p.index << "] ";
    }
    std::cout << std::endl;

    if (N == 3) 
    {
        std::cout << "Base case (3 points)\n";
        return CreateBaseTriangulation(points);
    }

    if (N == 4) 
    {
        std::cout << "Special case (4 points)\n";
        return HandleFourPoints(points);
    }
    
    if (N == 8) 
    {
        std::vector<Vertex> left(points.begin(), points.begin() + 4);
        std::vector<Vertex> right(points.begin() + 4, points.end());
        
        std::cout << "Splitting 8 points into:\nLeft: ";
        for (const auto& p : left) std::cout << "[" << p.index << "] ";
        std::cout << "\nRight: ";
        for (const auto& p : right) std::cout << "[" << p.index << "] ";
        std::cout << std::endl;

        auto leftTri = DivideAndConquer(left);
        auto rightTri = DivideAndConquer(right);
        
        return Merge(leftTri, rightTri);
    }
    
    if (N < 12) 
    {
        std::vector<Vertex> left(points.begin(), points.begin() + 3);
        std::vector<Vertex> right(points.begin() + 3, points.end());
        
        std::cout << "Splitting " << N << " points into:\nLeft: ";
        for (const auto& p : left) std::cout << "[" << p.index << "] ";
        std::cout << "\nRight: ";
        for (const auto& p : right) std::cout << "[" << p.index << "] ";
        std::cout << std::endl;

    auto leftTri = DivideAndConquer(left);
    auto rightTri = DivideAndConquer(right);

    return Merge(leftTri, rightTri);
}

    size_t mid = N / 2;
    std::vector<Vertex> left(points.begin(), points.begin() + mid);
    std::vector<Vertex> right(points.begin() + mid, points.end());

    std::cout << "Splitting " << N << " points into:\nLeft: ";
    for (const auto& p : left) std::cout << "[" << p.index << "] ";
    std::cout << "\nRight: ";
    for (const auto& p : right) std::cout << "[" << p.index << "] ";
    std::cout << std::endl;

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

void DelaunayTriangulation::GenerateRandomPoints(int count, float scale, float offsetX, float offsetY, 
                                               int windowWidth, int windowHeight) 
{
    vertices.clear();
    faces.clear();

    std::mt19937 rng(std::time(nullptr));
    std::uniform_real_distribution<double> distX(-40.0, 40.0);
    std::uniform_real_distribution<double> distY(-30.0, 30.0);

    std::set<std::pair<double, double>> uniquePoints;

    while (uniquePoints.size() < count) {
        double x = std::round(distX(rng) * 10.0) / 10.0;
        double y = std::round(distY(rng) * 10.0) / 10.0;

        float screenX = x * scale + offsetX;
        float screenY = -y * scale + offsetY;

        if (screenX >= 0 && screenX <= windowWidth &&
            screenY >= 0 && screenY <= windowHeight) {
            if (uniquePoints.insert({x, y}).second) {
                size_t idx = vertices.size();
                AddVertex(x, y);
            }
        }
    }

    std::cout << "\nGenerated " << vertices.size() << " points:\n";
    for (const auto& v : vertices) {
        std::cout << "[" << v.index << "]: " << v.x << "," << v.y << std::endl;
    }
}

bool DelaunayTriangulation::EdgesIntersect(Vertex* a1, Vertex* a2, Vertex* b1, Vertex* b2) const 
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

void DelaunayTriangulation::FindTangents(const std::vector<Face>& left, const std::vector<Face>& right,
    Vertex*& P0, Vertex*& P1, Vertex*& P2, Vertex*& P3) 
{
    std::vector<Vertex*> leftVertices;
    std::vector<Vertex*> rightVertices;
    
    for (const Face& face : left) {
        if (std::find(leftVertices.begin(), leftVertices.end(), face.v1) == leftVertices.end())
            leftVertices.push_back(face.v1);
        if (std::find(leftVertices.begin(), leftVertices.end(), face.v2) == leftVertices.end())
            leftVertices.push_back(face.v2);
        if (std::find(leftVertices.begin(), leftVertices.end(), face.v3) == leftVertices.end())
            leftVertices.push_back(face.v3);
    }
    
    for (const Face& face : right) {
        if (std::find(rightVertices.begin(), rightVertices.end(), face.v1) == rightVertices.end())
            rightVertices.push_back(face.v1);
        if (std::find(rightVertices.begin(), rightVertices.end(), face.v2) == rightVertices.end())
            rightVertices.push_back(face.v2);
        if (std::find(rightVertices.begin(), rightVertices.end(), face.v3) == rightVertices.end())
            rightVertices.push_back(face.v3);
    }

    std::sort(leftVertices.begin(), leftVertices.end(),
        [](Vertex* a, Vertex* b) { return a->y < b->y; });
    std::sort(rightVertices.begin(), rightVertices.end(),
        [](Vertex* a, Vertex* b) { return a->y < b->y; });

    std::vector<Vertex*> allVertices;
    allVertices.insert(allVertices.end(), leftVertices.begin(), leftVertices.end());
    allVertices.insert(allVertices.end(), rightVertices.begin(), rightVertices.end());

    std::cout << "\nSearching for upper tangent...\n";
    bool foundUpper = false;
    for (Vertex* leftIt : leftVertices) {
        for (Vertex* rightIt : rightVertices) {
            if (IsValidUpperTangent(leftIt, rightIt, allVertices)) {
                P0 = leftIt;
                P1 = rightIt;
                std::cout << "Found upper tangent: [" << P0->index << "] -> [" << P1->index << "]\n";
                foundUpper = true;
                break;
            }
        }
        if (foundUpper) break;
    }

    std::cout << "\nSearching for lower tangent...\n";
    bool foundLower = false;
    for (Vertex* leftV : leftVertices) {
        for (Vertex* rightV : rightVertices) {
            if (IsValidLowerTangent(leftV, rightV, allVertices)) {
                P2 = leftV;
                P3 = rightV;
                std::cout << "Found lower tangent: [" << P2->index << "] -> [" << P3->index << "]\n";
                foundLower = true;
                break;
            }
        }
        if (foundLower) break;
    }

    if (!foundUpper || !foundLower) {
        std::cout << "Error: Could not find valid tangents!\n";
        return;
    }

    std::cout << "\nFinal tangent lines:";
    std::cout << "\nUpper tangent: [" << P0->index << "] -> [" << P1->index << "]";
    std::cout << "\nLower tangent: [" << P2->index << "] -> [" << P3->index << "]\n";
}

bool DelaunayTriangulation::IsValidUpperTangent(Vertex* leftPoint, Vertex* rightPoint, const std::vector<Vertex*>& vertices) const 
{
    for (Vertex* v : vertices) {
        if (v == leftPoint || v == rightPoint) continue;
        
        double cross = (rightPoint->x - leftPoint->x) * (v->y - leftPoint->y) -
                      (rightPoint->y - leftPoint->y) * (v->x - leftPoint->x);
        
        if (cross < 0) {
            std::cout << "Point [" << v->index << "] is on the wrong side of upper tangent\n";
            return false;
        }
    }
    
    std::cout << "Valid upper tangent: [" << leftPoint->index << "] -> [" << rightPoint->index << "]\n";
    return true;
}

bool DelaunayTriangulation::IsValidLowerTangent(Vertex* leftPoint, Vertex* rightPoint, const std::vector<Vertex*>& vertices) const 
{
    for (Vertex* v : vertices) {
        if (v == leftPoint || v == rightPoint) continue;
        
        double cross = (rightPoint->x - leftPoint->x) * (v->y - leftPoint->y) -
                      (rightPoint->y - leftPoint->y) * (v->x - leftPoint->x);
        
        if (cross > 0) {
            std::cout << "Point [" << v->index << "] is on the wrong side of lower tangent\n";
            return false;
        }
    }
    
    std::cout << "Valid lower tangent: [" << leftPoint->index << "] -> [" << rightPoint->index << "]\n";
    return true;
}
