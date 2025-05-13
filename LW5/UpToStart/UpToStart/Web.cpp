#include "Web.h"
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>

Web::Vertex* Web::GetVertexById(uint32_t id) {
	auto it = std::find_if(m_vertices.begin(), m_vertices.end(),
						   [id](const Vertex& v) { return v.id == id; });
	return (it != m_vertices.end()) ? &(*it) : nullptr;
}

Web::Edge* Web::GetReverseEdge(Edge* forwardEdge) {
	if (!forwardEdge) {
		throw std::logic_error("GetReverseEdge called with a null forwardEdge pointer.");
	}
	if (!forwardEdge->reverseEdge) {
		throw std::logic_error("GetReverseEdge called with a null reverseEdge pointer.");
	}
	return forwardEdge->reverseEdge;
}

void Web::ReadGraph(std::istream& input)
{
	if (!(input >> m_size)) {
		throw std::runtime_error("Invalid number of vertices.");
	}

	m_vertices.resize(m_size);
	SINK_ID = m_size - 1;

	for (uint32_t i = 0; i < m_size; ++i) {
		m_vertices[i].id = i;
		m_vertices[i].level = 0;
		m_vertices[i].excess = 0;
	}

	Vertex* sourceVertex = GetVertexById(SOURCE_ID);
	sourceVertex->level = m_size;

	std::string line;
	std::getline(input, line); 

	while (std::getline(input, line))
	{
		uint32_t id1, id2;
		int capability;

		std::istringstream iss(line);
		if (!(iss >> id1 >> id2 >> capability)) {
			throw std::invalid_argument("Invalid arguments");
		}

		size_t forwardEdgeIdx = m_edges.size();
		m_edges.emplace_back(id1, id2, capability, 0);

		size_t reverseEdgeIdx = m_edges.size();
		m_edges.emplace_back(id2, id1, 0, 0);

		Edge* fe = &m_edges[forwardEdgeIdx];
		Edge* re = &m_edges[reverseEdgeIdx];

		fe->reverseEdge = re;
		re->reverseEdge = fe;

		if (id1 == SOURCE_ID && capability > 0) {
			Vertex* neighbor = GetVertexById(id2);
			int pushAmount = capability;
			
			fe->flow = pushAmount;
			re->flow = -pushAmount;

			neighbor->excess += pushAmount;
			sourceVertex->excess -= pushAmount;
		}
	}
}

void Web::Push(Vertex& u, Vertex& v, Edge& edge) {
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
	bool foundNeighbor = false;

	for (const auto& edge : m_edges) {
		int residualCapacity = 0;
		uint32_t neighborId = std::numeric_limits<uint32_t>::max();

		if (edge.vStartIndex == u.id) {
			residualCapacity = edge.capacity - edge.flow;
			neighborId = edge.vEndIndex; 
		} 

		if (residualCapacity > 0) { 
			Vertex* v = GetVertexById(neighborId);
			if (v) {
				minNeighborLevel = std::min(minNeighborLevel, v->level);
				foundNeighbor = true;
			}
		}
	}

	if (foundNeighbor) {
		u.level = minNeighborLevel + 1;
	} else {
		u.level = m_size * 2; 
	}
}

void Web::Deactivate(Vertex& u) {
	while (u.excess > 0) {
		bool pushed = false;
		for (auto& edge : m_edges) {
			if (u.excess == 0) break;

			if (edge.vStartIndex == u.id) {
				int residualCapacity = edge.capacity - edge.flow;
				if (residualCapacity > 0) {
					Vertex* v = GetVertexById(edge.vEndIndex);
					if (v && u.level == v->level + 1) {
						Push(u, *v, edge);
						pushed = true;
					}
				}
			}
		}

		if (!pushed) {
			Relabel(u);
			continue;
		}
	}
}

int Web::FindMaximalCapacity()
{
	bool activeVertexExists;

	do {
		activeVertexExists = false;
		for (auto& vertex : m_vertices) { 
			if (vertex.excess > 0 && vertex.id != SOURCE_ID && vertex.id != SINK_ID) {
				activeVertexExists = true;
				Deactivate(vertex); 
				break; 
			}
		}
	} while (activeVertexExists);

	Vertex* sinkVertex = GetVertexById(SINK_ID);
	return std::max(0, sinkVertex->excess); 
}
