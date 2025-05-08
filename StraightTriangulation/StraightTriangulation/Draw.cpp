#include "Draw.h"
#include <iostream>

Draw::Draw(const StraightTriangulation& triangulation)
    : triangulation(triangulation) {
}

void Draw::Print() const
{
    sf::VideoMode videoMode;
    videoMode.width = 800;
    videoMode.height = 600;
    sf::RenderWindow window(videoMode, "Straight Delaunay Triangulation");
    window.setFramerateLimit(60);

    while (window.isOpen())
    {
        sf::Event event;
        if (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }

        window.clear(sf::Color::White);

        std::vector<DrawEdge> edges;

        for (const auto& face : triangulation.GetFaces())
        {
            sf::Vector2f p1 = ScalePoint(*face.v1, window);
            sf::Vector2f p2 = ScalePoint(*face.v2, window);
            sf::Vector2f p3 = ScalePoint(*face.v3, window);

            DrawEdge e1(p1, p2);
            DrawEdge e2(p2, p3);
            DrawEdge e3(p3, p1);

            auto addEdgeIfNotExists = [&edges](const DrawEdge& edge)
                {
                    if (std::find(edges.begin(), edges.end(), edge) == edges.end())
                    {
                        edges.push_back(edge);
                    }
                };

            addEdgeIfNotExists(e1);
            addEdgeIfNotExists(e2);
            addEdgeIfNotExists(e3);
        }

        for (const auto& edge : edges)
        {
            sf::RectangleShape line(sf::Vector2f(
                std::sqrt(std::pow(edge.p2.x - edge.p1.x, 2) +
                    std::pow(edge.p2.y - edge.p1.y, 2)), 2.0f));

            line.setPosition(edge.p1);
            line.setFillColor(sf::Color::Black);

            float angle = std::atan2(edge.p2.y - edge.p1.y, edge.p2.x - edge.p1.x);
            line.setRotation(angle * 180 / 3.14159f);

            window.draw(line);
        }

        for (const auto& vertex : triangulation.GetVertices())
        {
            sf::CircleShape point(5.0f);
            point.setFillColor(sf::Color::Red);

            sf::Vector2f pos = ScalePoint(vertex, window);
            point.setPosition(pos.x - point.getRadius(), pos.y - point.getRadius());

            window.draw(point);
        }

        window.display();
        sf::sleep(sf::milliseconds(16));
    }
}

sf::Vector2f Draw::ScalePoint(const Vertex& point, const sf::RenderWindow& window) const
{
    float minX = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::lowest();
    float minY = std::numeric_limits<float>::max();
    float maxY = std::numeric_limits<float>::lowest();

    for (const auto& v : triangulation.GetVertices())
    {
        minX = std::min(minX, static_cast<float>(v.x));
        maxX = std::max(maxX, static_cast<float>(v.x));
        minY = std::min(minY, static_cast<float>(v.y));
        maxY = std::max(maxY, static_cast<float>(v.y));
    }

    float rangeX = maxX - minX;
    float rangeY = maxY - minY;

    if (rangeX < 0.000001f) rangeX = 1.0f;
    if (rangeY < 0.000001f) rangeY = 1.0f;

    float padding = 50.0f;
    float width = static_cast<float>(window.getSize().x);
    float height = static_cast<float>(window.getSize().y);

    float scale = std::min((width - 2 * padding) / rangeX,
        (height - 2 * padding) / rangeY);

    float offsetX = (width - scale * rangeX) / 2.0f;
    float offsetY = (height - scale * rangeY) / 2.0f;

    float scaledX = offsetX + scale * (static_cast<float>(point.x) - minX);
    float scaledY = offsetY + scale * (static_cast<float>(point.y) - minY);

    return sf::Vector2f(scaledX, scaledY);
}

void Draw::FitToScreen(std::vector<sf::Vector2f>& points, const sf::RenderWindow& window) const
{
    if (points.empty()) return;

    float minX = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::lowest();
    float minY = std::numeric_limits<float>::max();
    float maxY = std::numeric_limits<float>::lowest();

    for (const auto& p : points)
    {
        minX = std::min(minX, p.x);
        maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y);
        maxY = std::max(maxY, p.y);
    }

    float padding = 50.0f;
    float scaleX = (window.getSize().x - 2 * padding) / (maxX - minX);
    float scaleY = (window.getSize().y - 2 * padding) / (maxY - minY);

    for (auto& p : points)
    {
        p.x = padding + (p.x - minX) * scaleX;
        p.y = padding + (p.y - minY) * scaleY;
    }
}