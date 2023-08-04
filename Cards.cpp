#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::unordered_map;
using std::vector;

const int64_t kInf = 2009000999;

/*!
 * Default edge with weight for working with graphs.
 *
 * @tparam T
 * @tparam W
 */
template <typename T = int, typename W = int>
class DefaultEdge {
 public:
  DefaultEdge(const T& first, const T& second, W weight)
      : source_(first), destination_(second), weight_(weight) {}

  const T& From() const { return source_; }
  const T& To() const { return destination_; }
  const W& Weight() const { return weight_; }

 private:
  T source_ = 0;
  T destination_ = 0;
  W weight_ = 0;
};

template <typename T = int, typename W = int>
bool operator<(const DefaultEdge<T, W>& lhs, const DefaultEdge<T, W>& rhs) {
  if (lhs.Weight() == rhs.Weight()) {
    return lhs.To() < rhs.To();
  }
  return lhs.Weight() < rhs.Weight();
}

template <typename Iterator>
struct IteratorRange {
  Iterator begin() { return start; }  // NOLINT
  Iterator end() { return finish; }   // NOLINT

  Iterator start;
  Iterator finish;
};

/*!
 * @brief class for working with applied tasks needs mathematical objects
 * graphs. It can be created like adjacency matrix or list.
 * @tparam Vertex
 * @tparam Edge
 * @tparam Weight
 * @param vertex_number_ - number of vertices
 * @param edges_number_ - number of edges
 */
template <typename Vertex = int, typename Weight = int,
          typename Edge = DefaultEdge<Vertex, Weight>>
class Graph {
 protected:
  size_t vertex_number_ = 0;
  size_t edges_number_ = 0;

 public:
  using VertexType = Vertex;
  using EdgeType = Edge;
  using WeightType = Weight;

  explicit Graph(size_t vert_number, size_t edges_number = 0)
      : vertex_number_(vert_number), edges_number_(edges_number) {}

  virtual size_t GetVertexNumber() const { return vertex_number_; }

  virtual size_t GetEdgesNumber() const { return edges_number_; }

  /*!
   * @param vertex
   * @return set<pair<VertexType, WeightType>> with vertices connected with
   * given vertex by edge.
   */
  virtual typename set<pair<VertexType, WeightType>>::iterator GetNeighboursIt(
      const Vertex& vertex) = 0;

  /*!
   * @param vertex
   * @return set<pair<VertexType, WeightType>> with vertices connected with
   * given vertex by edge.
   */
  virtual set<pair<VertexType, WeightType>> GetNeighbours(
      const Vertex& vertex) = 0;
};

/*!
 * Most convenient representation of graph.
 * @tparam VertexType
 * @tparam EdgeType
 * @tparam WeightType
 */
template <typename VertexType = int, typename WeightType = int,
          typename EdgeType = DefaultEdge<VertexType, WeightType>>
class AdjList : public Graph<VertexType, WeightType, EdgeType> {
  map<VertexType, set<pair<VertexType, WeightType>>> list_;

  friend bool operator<(const pair<VertexType, WeightType>& lhs,
                        const pair<VertexType, WeightType>& rhs) {
    return lhs.second == rhs.second ? lhs.first < rhs.first
                                    : lhs.second < rhs.second;
  }

 public:
  AdjList(size_t vert_number, const vector<EdgeType>& edges)
      : Graph<VertexType, WeightType, EdgeType>(vert_number, edges.size()) {
    for (const auto& e : edges) {
      list_[e.From()].insert({e.To(), e.Weight()});
      list_[e.To()].insert({e.From(), e.Weight()});
    }
  }

  set<pair<VertexType, WeightType>> GetNeighbours(
      const VertexType& vertex) final {
    return list_[vertex];
  }

  typename set<pair<VertexType, WeightType>>::iterator GetNeighboursIt(
      const VertexType& vertex) final {
    return list_[vertex].begin();
  }
};

/*!
 * Alternative and less popular representation of graph, but may be useful in
 * some tasks.
 * @tparam VertexType
 * @tparam EdgeType
 */
template <typename VertexType, typename EdgeType = DefaultEdge<VertexType>>
class AdjMatrix : public Graph<VertexType, EdgeType> {
  vector<vector<VertexType>> data_;

 public:
  AdjMatrix(size_t vert_number, const vector<EdgeType>& edges)
      : Graph<VertexType, EdgeType>(vert_number, edges.size()),
        data_(vector<vector<VertexType>>(vert_number)) {
    for (size_t i = 0; i < vert_number; ++i) {
      data_[i] = vector<VertexType>(vert_number);
    }

    for (const auto& e : edges) {
      data_[e.From()][e.To] = 1;
      data_[e.To()][e.From] = 1;
    }
  }

  /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
  vector<VertexType> GetNeighbours(VertexType vertex) final;

  /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
  typename vector<VertexType>::iterator GetNeighboursIt(
      const VertexType& vertex) final {
    return data_[vertex].begin();
  }
};

template <typename VertexType, typename EdgeType>
vector<VertexType> AdjMatrix<VertexType, EdgeType>::GetNeighbours(
    VertexType vertex) {
  vector<VertexType> result;

  for (VertexType i = 0; i < this->vertex_number_; ++i) {
    if (data_[vertex][i]) {
      result.push_back(i);
    }
  }

  return result;
}

/*!
 * @brief Special class for saving result of Dijkstra's algorithm.
 * @tparam VertexType
 * @tparam EdgeType
 * @tparam WeightType
 * @param ancestors_ - unordered map with ancestor for every vertex.
 * @param dist_ - distances between src and other vertices
 */
template <typename VertexType = int, typename WeightType = int,
          typename EdgeType = DefaultEdge<VertexType, WeightType>>
class DijkstraVisitor {
 public:
  DijkstraVisitor(size_t size, WeightType border, VertexType src)
      : dist_(size, border) {
    dist_[src] = 0;
    ancestors_[src] = -1;
  }

  /*!
   * @brief Function prints distances between source and other vertices.
   * @tparam VertexType
   * @tparam WeightType
   * @tparam EdgeType
   */
  void PrintRo();

  /*!
   * @brief Function checks edge between vert and i and updates result.
   * @param vert - src vertex
   * @param i - dest vertex and weight
   * @param sizes
   */
  void TreeEdge(VertexType vert, pair<VertexType, WeightType> i,
                set<pair<int64_t, int64_t>>& sizes);

  ~DijkstraVisitor() = default;

 private:
  std::vector<WeightType> dist_;
  std::unordered_map<VertexType, VertexType> ancestors_;
};

/*!
 * @brief Dijkstra's algorithm
 * @tparam Graph
 * @tparam Visitor
 * @param graph - current graph
 * @param visitor - current visitor
 * @param src - start
 */
template <typename Graph, class Visitor>
void Dijkstra(Graph& graph, Visitor& visitor, int64_t src);

/*!
 * @brief Default finction for initialize vector of edges by keyboard.
 * @tparam EdgeType
 * @param edges - vector with edges
 * @param number_of_rooms - number of vertices
 * @param number_of_connections - number of edges
 */
template <typename EdgeType>
void Initialize(vector<EdgeType>& edges, int64_t number_of_rooms,
                int64_t number_of_connections);

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  std::cout.tie(nullptr);

  int64_t number_of_graphs = 0;
  std::cin >> number_of_graphs;

  for (int64_t i = 0; i < number_of_graphs; ++i) {
    size_t number_of_rooms = 0;
    int64_t number_of_connections = 0;
    int64_t start = 0;

    std::cin >> number_of_rooms >> number_of_connections;

    std::vector<DefaultEdge<int64_t, int64_t>> edges;

    Initialize(edges, number_of_rooms, number_of_connections);

    AdjList<int64_t, int64_t> graph(number_of_rooms, edges);

    std::cin >> start;

    DijkstraVisitor<int64_t, int64_t, int64_t> visitor(graph.GetVertexNumber(),
                                                       kInf, start);
    Dijkstra(graph, visitor, start);
    visitor.PrintRo();
  }

  return 0;
}

template <typename EdgeType>
void Initialize(vector<EdgeType>& edges, int64_t number_of_rooms,
                int64_t number_of_connections) {
  for (int64_t i = 0; i < number_of_rooms; i++) {
    edges[i];
  }

  for (int64_t i = 0; i < number_of_connections; i++) {
    int64_t source = 0;
    int64_t destination = 0;
    int64_t weight = 0;
    std::cin >> source >> destination >> weight;

    edges.push_back({source, destination, weight});
  }
}

template <typename Graph, class Visitor>
void Dijkstra(Graph& graph, Visitor& visitor, int64_t src) {
  set<pair<int64_t, int64_t>> sizes;
  sizes.insert({0, src});

  while (!sizes.empty()) {
    auto top = sizes.begin();
    int64_t vert = top->second;
    sizes.erase(top);

    for (auto i : graph.GetNeighbours(vert)) {
      visitor.TreeEdge(vert, i, sizes);
    }
  }
}

template <typename VertexType, typename WeightType, typename EdgeType>
void DijkstraVisitor<VertexType, WeightType, EdgeType>::PrintRo() {
  for (auto& i : dist_) {
    std::cout << i << ' ';
  }
  std::cout << '\n';
}

template <typename VertexType, typename WeightType, typename EdgeType>
void DijkstraVisitor<VertexType, WeightType, EdgeType>::TreeEdge(
    VertexType vert, pair<VertexType, WeightType> i,
    set<pair<int64_t, int64_t>>& sizes) {
  if (dist_[i.first] > dist_[vert] + i.second) {
    if (dist_[i.first] < kInf) {
      auto it = sizes.find({dist_[i.first], i.first});
      if (it != sizes.end()) {
        sizes.erase(it);
      }
    }
    dist_[i.first] = dist_[vert] + i.second;
    ancestors_[i.first] = vert;
    sizes.insert({dist_[i.first], i.first});
  }
}
