#include "Draw.h"

Draw::Draw(const DelaunayTriangulation& triangulation) 
    : triangulation(triangulation) {}

void Draw::Print() const {
    sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay Triangulation");
    
    std::vector<sf::Vector2f> points;
    for (const auto& vertex : triangulation.GetVertices()) 
    {
        points.emplace_back(vertex.x, vertex.y);
    }
    
    FitToScreen(points, window);
    
    while (window.isOpen()) 
    {
        sf::Event event;
        while (window.pollEvent(event)) 
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        window.clear(sf::Color::White);
        
        for (const auto& face : triangulation.GetFaces()) 
        {
            sf::ConvexShape triangle;
            triangle.setPointCount(3);
            
            sf::Vector2f p1 = ScalePoint(*face.v1, window);
            sf::Vector2f p2 = ScalePoint(*face.v2, window);
            sf::Vector2f p3 = ScalePoint(*face.v3, window);
            
            triangle.setPoint(0, p1);
            triangle.setPoint(1, p2);
            triangle.setPoint(2, p3);
            
            triangle.setFillColor(sf::Color::Transparent);
            triangle.setOutlineColor(sf::Color::Black);
            triangle.setOutlineThickness(1.0f);
            
            window.draw(triangle);
        }
        
        for (const auto& vertex : triangulation.GetVertices()) 
        {
            sf::CircleShape point(3.0f);
            point.setFillColor(sf::Color::Red);
            
            sf::Vector2f pos = ScalePoint(vertex, window);
            point.setPosition(pos.x - point.getRadius(), pos.y - point.getRadius());
            
            window.draw(point);
        }
        
        window.display();
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
        minX = std::min(minX, (float)v.x);
        maxX = std::max(maxX, (float)v.x);
        minY = std::min(minY, (float)v.y);
        maxY = std::max(maxY, (float)v.y);
    }
    
    float padding = 50.0f;
    float scaleX = (window.getSize().x - 2 * padding) / (maxX - minX);
    float scaleY = (window.getSize().y - 2 * padding) / (maxY - minY);
    
    return sf::Vector2f(
        padding + (point.x - minX) * scaleX,
        padding + (point.y - minY) * scaleY
    );
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