#include "Web.h"
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>

Web::Vertex* Web::GetVertexById(uint32_t id) 
{
	auto it = std::find_if(m_vertices.begin(), m_vertices.end(),
						   [id](const Vertex& v) { return v.id == id; });
	return (it != m_vertices.end()) ? &(*it) : nullptr;
}

Web::Edge* Web::GetReverseEdge(Edge* forwardEdge) 
{
	if (!forwardEdge) 
	{
		throw std::logic_error("GetReverseEdge called with a null forwardEdge pointer.");
	}
	if (!forwardEdge->reverseEdge) 
	{
		throw std::logic_error("GetReverseEdge called, but reverseEdge is null.");
	}
	return forwardEdge->reverseEdge;
}

void Web::ReadGraph(std::istream& input)
{
	if (!(input >> m_size)) 
	{
		throw std::runtime_error("Invalid number of vertices.");
	}

	m_vertices.resize(m_size);
	SINK_ID = m_size - 1;

	for (uint32_t i = 0; i < m_size; ++i) 
	{
		m_vertices[i].id = i;
		m_vertices[i].level = 0;
		m_vertices[i].excess = 0;
	}

	Vertex* sourceVertex = GetVertexById(SOURCE_ID);
	sourceVertex->level = m_size;

	std::string line;
	if (input.peek() == '\n') 
	{
		std::getline(input, line);
	}

	while (std::getline(input, line))
	{
		uint32_t id1, id2;
		int capacity;

		std::istringstream iss(line);
		if (!(iss >> id1 >> id2 >> capacity)) 
		{
			throw std::invalid_argument("Invalid arguments in graph file line.");
		}

		if (id1 >= m_size || id2 >= m_size) 
		{
			throw std::out_of_range("Vertex ID out of range.");
		}

		size_t forwardEdgeIdx = m_edges.size();
		m_edges.emplace_back(id1, id2, capacity, 0);

		size_t reverseEdgeIdx = m_edges.size();
		m_edges.emplace_back(id2, id1, 0, 0);

		Edge* fe = &m_edges[forwardEdgeIdx];
		Edge* re = &m_edges[reverseEdgeIdx];

		fe->reverseEdge = re;
		re->reverseEdge = fe;

		if (id1 == SOURCE_ID && capacity > 0) 
		{
			Vertex* sourceV = GetVertexById(SOURCE_ID);
			Vertex* neighbor = GetVertexById(id2);

			int pushAmount = capacity;
			
			fe->flow = pushAmount;
			re->flow = -pushAmount;

			neighbor->excess += pushAmount;
			sourceV->excess -= pushAmount;
		}
	}

	m_activeVertices.clear();
	for (auto& vertex : m_vertices)
	{
		if (vertex.id != SOURCE_ID && vertex.id != SINK_ID && vertex.excess > 0) 
		{
			m_activeVertices.push_back(&vertex);
		}
	}
}

void Web::Push(Vertex& u, Vertex& v, Edge& edge) 
{
	int residualCapacity = edge.capacity - edge.flow;
	int delta = std::min(u.excess, residualCapacity);

	if (delta <= 0) return;

	edge.flow += delta;
	u.excess -= delta;
	v.excess += delta;

	Edge* reverseEdge = GetReverseEdge(&edge);
	reverseEdge->flow -= delta;
}

void Web::Relabel(Vertex& u) {
    int minNeighborLevel = std::numeric_limits<int>::max();
    bool foundValidNeighbor = false;

    for (const auto& edge : m_edges) 
	{
        if (edge.vStartIndex != u.id) continue;

        int residualCapacity = edge.capacity - edge.flow;
        if (residualCapacity > 0) 
		{
            Vertex* v = GetVertexById(edge.vEndIndex);
            if (v) 
			{
                minNeighborLevel = std::min(minNeighborLevel, v->level);
                foundValidNeighbor = true;
            }
        }
    }

    if (foundValidNeighbor) 
	{
        u.level = minNeighborLevel + 1;
    } else 
	{
        u.level = m_size;
    }
}

void Web::Deactivate(Vertex& u) 
{
    while (u.excess > 0) 
	{
        bool pushed = false;
        
        for (auto& edge : m_edges) 
		{
            if (edge.vStartIndex != u.id) continue;

            Vertex* v = GetVertexById(edge.vEndIndex);
            if (!v) continue;

            int residualCapacity = edge.capacity - edge.flow;
            if (residualCapacity > 0 && u.level == v->level + 1) 
			{
                Push(u, *v, edge);
                pushed = true;

                if (v->id != SOURCE_ID && v->id != SINK_ID && v->excess > 0 && v->level < m_size) 
				{
                    auto it = std::find_if(m_activeVertices.begin(), m_activeVertices.end(),
                                           [&](const Vertex* activeV) { return activeV->id == v->id; });
                    if (it == m_activeVertices.end()) 
					{
                        m_activeVertices.push_back(v);
                    }
                }
                
                if (u.excess == 0) break;
            }
        }

        if (!pushed) 
		{
            int oldLevel = u.level;
            Relabel(u);
            
            if (u.level == oldLevel) 
			{
                break;
            }
        }
    }
}

int Web::FindMaximalCapacity() 
{
    while (!m_activeVertices.empty()) 
	{
        Vertex* u = m_activeVertices.front();
        m_activeVertices.pop_front();

		Deactivate(*u);

        if (u->excess > 0) 
		{
            if (u->level < m_size) 
			{ 
                m_activeVertices.push_back(u);
            }
        }
    }

    Vertex* sinkVertex = GetVertexById(SINK_ID);
    return (sinkVertex) ? std::max(0, sinkVertex->excess) : 0;
}