#pragma once
#include <cstdint>
#include <vector>
#include <istream>
#include <limits>
#include <list>

class Web
{

public:
	const uint32_t SOURCE_ID = 0;
	uint32_t SINK_ID = std::numeric_limits<uint32_t>::max();

	struct Vertex
	{
		uint32_t id;
		int level;
		int excess;
	};

	struct Edge
	{
		Edge(uint32_t id1, uint32_t id2, int _capacity, int _flow = 0) : vStartIndex(id1), vEndIndex(id2), capacity(_capacity), flow(_flow) {};
		uint32_t vStartIndex;
		uint32_t vEndIndex;
		int capacity;
		int flow;
		Edge* reverseEdge = nullptr;
	};

private:
	using id = uint32_t;
	using level = uint32_t;

	std::vector<Vertex> m_vertices;
	std::vector<Edge> m_edges;
	std::list<Vertex*> m_activeVertices;
	uint32_t m_size = 0;

public:
	void ReadGraph(std::istream& input);
	int FindMaximalCapacity();

private:
	void Push(Vertex& u, Vertex& v, Edge& edge);
	void Relabel(Vertex& u);
	void Deactivate(Vertex& u);
	Web::Vertex* GetVertexById(uint32_t id);
	Web::Edge* GetReverseEdge(Edge* forwardEdge);
};

