#pragma once
#include <SFML/Graphics.hpp>
#include "StraightTriangulation.h"

class Draw
{
private:
    struct DrawEdge
    {
        sf::Vector2f p1;
        sf::Vector2f p2;

        DrawEdge(sf::Vector2f p1, sf::Vector2f p2) : p1(p1), p2(p2)
        {
            if (p1.x > p2.x || (p1.x == p2.x && p1.y > p2.y))
            {
                std::swap(this->p1, this->p2);
            }
        }

        bool operator==(const DrawEdge& other) const
        {
            return (p1.x == other.p1.x && p1.y == other.p1.y &&
                p2.x == other.p2.x && p2.y == other.p2.y);
        }
    };

    const StraightTriangulation& triangulation;
    sf::Vector2f ScalePoint(const Vertex& point, const sf::RenderWindow& window) const;
    void FitToScreen(std::vector<sf::Vector2f>& points, const sf::RenderWindow& window) const;

public:
    Draw(const StraightTriangulation& triangulation);
    void Print() const;
};