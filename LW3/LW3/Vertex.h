#pragma once

class Vertex 
{
public:
    double x, y;
    size_t index;

    Vertex(double x, double y, size_t idx);

    bool operator<(const Vertex& other) const;
};

