#include "StraightTriangulation.h"
#include <set>
#include <map>
#include <tuple>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

void StraightTriangulation::Triangulate()
{
    faces.clear();

    std::vector<Vertex*> vertexPtrs;
    for (auto& v : vertices) vertexPtrs.push_back(&v);

    std::vector<Vertex*> hull = ConvexHull(vertexPtrs);

    using EdgeKey = std::pair<Vertex*, Vertex*>;
    std::set<std::tuple<Vertex*, Vertex*, Vertex*>> usedTriangles;
    std::map<EdgeKey, int> edgeSide;

    std::vector<EdgeKey> activeEdges;
    for (size_t i = 0; i < hull.size(); ++i) {
        Vertex* a = hull[i];
        Vertex* b = hull[(i + 1) % hull.size()];
        activeEdges.emplace_back(a, b);
        edgeSide[{a, b}] = 1;
    }

    while (!activeEdges.empty()) {
        auto [a, b] = activeEdges.back();
        activeEdges.pop_back();

        if (edgeSide[{a, b}] == 2) continue;

        Vertex* best = nullptr;
        double minRadius = 1e100;

        for (Vertex* c : vertexPtrs) {
            if (c == a || c == b) continue;
            double orient = Cross(*a, *b, *c);
            if (orient <= 0) continue;

            std::set<Vertex*> tri = {a, b, c};
            if (usedTriangles.count({a, b, c}) || usedTriangles.count({a, c, b}) ||
                usedTriangles.count({b, a, c}) || usedTriangles.count({b, c, a}) ||
                usedTriangles.count({c, a, b}) || usedTriangles.count({c, b, a}))
                continue;

            double A = hypot(a->x - b->x, a->y - b->y);
            double B = hypot(b->x - c->x, b->y - c->y);
            double C = hypot(c->x - a->x, c->y - a->y);
            double s = (A + B + C) / 2.0;
            double area = sqrt(std::max(0.0, s * (s - A) * (s - B) * (s - C)));
            if (area < 1e-10) continue;

            double R = (A * B * C) / (4.0 * area);

            bool delaunay = true;
            for (Vertex* p : vertexPtrs) {
                if (p == a || p == b || p == c) continue;
                if (InCircle(*p, a, b, c)) {
                    delaunay = false;
                    break;
                }
            }
            if (!delaunay) continue;

            if (R < minRadius) {
                minRadius = R;
                best = c;
            }
        }

        if (best) {
            faces.emplace_back(a, b, best);
            usedTriangles.insert({a, b, best});
            usedTriangles.insert({a, best, b});
            usedTriangles.insert({b, a, best});
            usedTriangles.insert({b, best, a});
            usedTriangles.insert({best, a, b});
            usedTriangles.insert({best, b, a});

            for (auto e : {EdgeKey(a, best), EdgeKey(best, b)}) {
                if (edgeSide[e] < 2) {
                    activeEdges.push_back(e);
                    edgeSide[e]++;
                }
            }
            edgeSide[{a, b}]++;
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

double StraightTriangulation::Cross(const Vertex& O, const Vertex& A, const Vertex& B)
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
            while (hull.size() >= start + 2 && Cross(*hull[hull.size() - 2], *hull.back(), *p) <= 0)
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