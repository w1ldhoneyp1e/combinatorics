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

    std::set<EdgeInfo, EdgeInfoCompare> edges;
    for (const Face& face : triangulation) {
        edges.insert(EdgeInfo(face.v1, face.v2));
        edges.insert(EdgeInfo(face.v2, face.v3));
        edges.insert(EdgeInfo(face.v3, face.v1));
    }

    std::vector<EdgeInfo> edgesToRemove;
    for (const EdgeInfo& edge : edges) {
        if (EdgesIntersect(edge.v1, edge.v2, baseVertex, newVertex)) {
            std::cout << "Edge [" << edge.v1->index << "]->[" << edge.v2->index 
                      << "] intersects with new edge" << std::endl;
            edgesToRemove.push_back(edge);
        }
    }

    for (const auto& edge : edgesToRemove) {
        edges.erase(edge);
    }

    triangulation.clear();
    for (const EdgeInfo& edge : edges) {
        for (const EdgeInfo& secondEdge : edges) {
            if (edge == secondEdge) continue;
            
            if (edge.v1 == secondEdge.v1 || edge.v1 == secondEdge.v2 ||
                edge.v2 == secondEdge.v1 || edge.v2 == secondEdge.v2) {
                
                Vertex* commonVertex = nullptr;
                Vertex* otherVertex1 = nullptr;
                Vertex* otherVertex2 = nullptr;

                if (edge.v1 == secondEdge.v1) {
                    commonVertex = edge.v1;
                    otherVertex1 = edge.v2;
                    otherVertex2 = secondEdge.v2;
                } else if (edge.v1 == secondEdge.v2) {
                    commonVertex = edge.v1;
                    otherVertex1 = edge.v2;
                    otherVertex2 = secondEdge.v1;
                } else if (edge.v2 == secondEdge.v1) {
                    commonVertex = edge.v2;
                    otherVertex1 = edge.v1;
                    otherVertex2 = secondEdge.v2;
                } else if (edge.v2 == secondEdge.v2) {
                    commonVertex = edge.v2;
                    otherVertex1 = edge.v1;
                    otherVertex2 = secondEdge.v1;
                }

                EdgeInfo thirdEdge(otherVertex1, otherVertex2);
                if (edges.find(thirdEdge) != edges.end()) {
                    triangulation.emplace_back(commonVertex, otherVertex1, otherVertex2);
                }
            }
        }
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

std::vector<Face> DelaunayTriangulation::HandleFourPoints(std::vector<Vertex>& points) 
{
    std::vector<Face> result;
    
    Vertex* v1 = nullptr;
    Vertex* v2 = nullptr;
    Vertex* v3 = nullptr;
    Vertex* v4 = nullptr;

    double minY = std::numeric_limits<double>::infinity();
    int topIndex = 0;
    for (int i = 0; i < 4; i++) {
        if (points[i].y < minY) {
            minY = points[i].y;
            topIndex = i;
        }
    }

    for (auto& vertex : vertices)
    {
        if (std::abs(vertex.x - points[topIndex].x) < 1e-10 && 
            std::abs(vertex.y - points[topIndex].y) < 1e-10) 
            v1 = &vertex;
        if (std::abs(vertex.x - points[(topIndex + 1) % 4].x) < 1e-10 && 
            std::abs(vertex.y - points[(topIndex + 1) % 4].y) < 1e-10) 
            v2 = &vertex;
        if (std::abs(vertex.x - points[(topIndex + 2) % 4].x) < 1e-10 && 
            std::abs(vertex.y - points[(topIndex + 2) % 4].y) < 1e-10) 
            v3 = &vertex;
        if (std::abs(vertex.x - points[(topIndex + 3) % 4].x) < 1e-10 && 
            std::abs(vertex.y - points[(topIndex + 3) % 4].y) < 1e-10) 
            v4 = &vertex;
    }

    if (!v1 || !v2 || !v3 || !v4) 
    {
        std::cout << "Error: Could not find vertices in the main vector" << std::endl;
        return result;
    }

    std::cout << "Found four vertices for special case:\n";
    std::cout << "v1[" << v1->index << "](" << v1->x << "," << v1->y << ")\n";
    std::cout << "v2[" << v2->index << "](" << v2->x << "," << v2->y << ")\n";
    std::cout << "v3[" << v3->index << "](" << v3->x << "," << v3->y << ")\n";
    std::cout << "v4[" << v4->index << "](" << v4->x << "," << v4->y << ")\n";

    std::vector<Vertex*> others = {v2, v3, v4};
    std::sort(others.begin(), others.end(),
        [v1](Vertex* a, Vertex* b) {
            double angleA = std::atan2(a->y - v1->y, a->x - v1->x);
            double angleB = std::atan2(b->y - v1->y, b->x - v1->x);
            return angleA < angleB;
        });
    
    v2 = others[0];
    v3 = others[1];
    v4 = others[2];

    double cross1 = (v3->x - v1->x) * (v4->y - v2->y) - (v3->y - v1->y) * (v4->x - v2->x);
    double cross2 = (v4->x - v2->x) * (v1->y - v3->y) - (v4->y - v2->y) * (v1->x - v3->x);

    if (cross1 * cross2 > 0) 
    {
        result.emplace_back(v1, v2, v3);
        result.emplace_back(v3, v4, v1);
    } 
    else 
    {
        result.emplace_back(v2, v3, v4);
        result.emplace_back(v4, v1, v2);
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
