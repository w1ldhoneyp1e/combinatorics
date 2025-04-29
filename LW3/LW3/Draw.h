#pragma once
#include <SFML/Graphics.hpp>
#include "DelaunayTriangulation.h"

class Draw 
{
public:
    Draw(const DelaunayTriangulation& triangulation);
    void Print() const;

private:
    const DelaunayTriangulation& triangulation;
    
    sf::Vector2f ScalePoint(const Vertex& point, const sf::RenderWindow& window) const;
    void FitToScreen(std::vector<sf::Vector2f>& points, const sf::RenderWindow& window) const;
}; 