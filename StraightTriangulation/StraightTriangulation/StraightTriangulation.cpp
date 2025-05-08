#include "StraightTriangulation.h"
#include <set>
#include <map>
#include <tuple>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

bool StraightTriangulation::IsTriangleUsed(const std::set<std::tuple<Vertex*, Vertex*, Vertex*>>& usedTriangles, Vertex* a, Vertex* b, Vertex* c) {
    return usedTriangles.count({a, b, c}) || usedTriangles.count({a, c, b}) ||
           usedTriangles.count({b, a, c}) || usedTriangles.count({b, c, a}) ||
           usedTriangles.count({c, a, b}) || usedTriangles.count({c, b, a});
}

void StraightTriangulation::AddTriangle(std::set<std::tuple<Vertex*, Vertex*, Vertex*>>& usedTriangles, Vertex* a, Vertex* b, Vertex* c) {
    usedTriangles.insert({a, b, c});
    usedTriangles.insert({a, c, b});
    usedTriangles.insert({b, a, c});
    usedTriangles.insert({b, c, a});
    usedTriangles.insert({c, a, b});
    usedTriangles.insert({c, b, a});
    faces.emplace_back(a, b, c);
}

bool StraightTriangulation::IsDelaunay(Vertex* a, Vertex* b, Vertex* c, const std::vector<Vertex*>& vertexPtrs) {
    for (Vertex* p : vertexPtrs) {
        if (p == a || p == b || p == c) continue;
        if (InCircle(*p, a, b, c)) {
            return false;
        }
    }
    return true;
}

Vertex* StraightTriangulation::FindBestPointForEdge(const Edge& edge, const std::vector<Vertex*>& vertexPtrs, const std::set<std::tuple<Vertex*, Vertex*, Vertex*>>& usedTriangles) {
    Vertex* a = edge.v1;
    Vertex* b = edge.v2;
    Vertex* best = nullptr;
    double minRadius = 1e100;

    for (Vertex* c : vertexPtrs) {
        if (c == a || c == b) continue;
        double orient = GetOrientation(*a, *b, *c);
        if (orient <= 0) continue;

        if (IsTriangleUsed(usedTriangles, a, b, c)) continue;

        double A = hypot(a->x - b->x, a->y - b->y);
        double B = hypot(b->x - c->x, b->y - c->y);
        double C = hypot(c->x - a->x, c->y - a->y);
        double s = (A + B + C) / 2.0;
        double area = sqrt(std::max(0.0, s * (s - A) * (s - B) * (s - C)));
        if (area < 1e-10) continue;

        double R = (A * B * C) / (4.0 * area);

        if (!IsDelaunay(a, b, c, vertexPtrs)) continue;

        if (R < minRadius) {
            minRadius = R;
            best = c;
        }
    }
    return best;
}

void StraightTriangulation::Triangulate()
{
    faces.clear();

    std::vector<Vertex*> vertexPtrs;
    for (auto& v : vertices) vertexPtrs.push_back(&v);

    std::vector<Vertex*> hull = ConvexHull(vertexPtrs);

    std::set<std::tuple<Vertex*, Vertex*, Vertex*>> usedTriangles;
    std::map<Edge, int> edgeSide;

    std::vector<Edge> activeEdges;
    for (size_t i = 0; i < hull.size(); ++i) {
        Vertex* a = hull[i];
        Vertex* b = hull[(i + 1) % hull.size()];
        Edge e(a, b);
        activeEdges.push_back(e);
        edgeSide[e] = 1;
    }

    while (!activeEdges.empty()) {
        Edge edge = activeEdges.back();
        activeEdges.pop_back();

        if (edgeSide[edge] == 2) continue;

        Vertex* best = FindBestPointForEdge(edge, vertexPtrs, usedTriangles);

        if (best) {
            AddTriangle(usedTriangles, edge.v1, edge.v2, best);

            for (Edge e2 : {Edge(edge.v1, best), Edge(best, edge.v2)}) {
                if (edgeSide[e2] < 2) {
                    activeEdges.push_back(e2);
                    edgeSide[e2]++;
                }
            }
            edgeSide[edge]++;
        }
    }
}

bool StraightTriangulation::InCircle(const Vertex& p, Vertex* a, Vertex* b, Vertex* c)
{
    double ax = a->x - p.x, ay = a->y - p.y;
    double bx = b->x - p.x, by = b->y - p.y;
    double cx = c->x - p.x, cy = c->y - p.y;

    double det = (ax * ax + ay * ay) * (bx * cy - by * cx)
               - (bx * bx + by * by) * (ax * cy - ay * cx)
               + (cx * cx + cy * cy) * (ax * by - ay * bx);
    return det > 0;
}

double StraightTriangulation::GetOrientation(const Vertex& O, const Vertex& A, const Vertex& B)
{
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

std::vector<Vertex*> StraightTriangulation::ConvexHull(std::vector<Vertex*> vertexPtrs)
{
    if (vertexPtrs.size() <= 1) return vertexPtrs;

    std::sort(vertexPtrs.begin(), vertexPtrs.end(), [](Vertex* a, Vertex* b) {
        return *a < *b;
    });

    std::vector<Vertex*> hull;
    for (int i = 0; i < 2; i++) {
        auto start = hull.size();
        for (Vertex* p : (i == 0 ? vertexPtrs : std::vector<Vertex*>(vertexPtrs.rbegin(), vertexPtrs.rend()))) {
            while (hull.size() >= start + 2 && GetOrientation(*hull[hull.size() - 2], *hull.back(), *p) <= 0)
                hull.pop_back();
            hull.push_back(p);
        }
        hull.pop_back();
    }

    return hull;
}

void StraightTriangulation::GenerateRandomPoints(int count, float scale, float offsetX, float offsetY,
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
            if (uniquePoints.insert({ x, y }).second) {
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

void StraightTriangulation::AddVertex(double x, double y)
{
    vertices.emplace_back(x, y, vertices.size());
}

Face::Faces StraightTriangulation::GetFaces() const
{
    return faces;
}
std::vector<Vertex> StraightTriangulation::GetVertices() const
{
    return vertices;
}