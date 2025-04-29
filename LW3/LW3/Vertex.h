#pragma once

class Vertex {
public:
    double x, y;
    int index;

    Vertex(double x = 0, double y = 0, int idx = 0);

    bool operator<(const Vertex& other) const;
};

