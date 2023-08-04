#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::unordered_map;
using std::vector;

/*!
 * Default edge with weight for working with graphs.
 *
 * @tparam VertexType
 * @tparam WeightType
 */
template <typename VertexType = int, typename WeightType = int>
class WeightedEdge {
 public:
  WeightedEdge() = default;
  WeightedEdge(const VertexType& first, const VertexType& second, const WeightType& weight)
      : source_(first), destination_(second), weight_(weight) {}

  [[nodiscard]] const VertexType& From() const { return source_; }
  [[nodiscard]] const VertexType& To() const { return destination_; }
  [[nodiscard]] const WeightType& Weight() const { return weight_; }

 private:
  VertexType source_ = 0;
  VertexType destination_ = 0;
  WeightType weight_ = 0;
};

template <typename VertexType = int, typename WeightType = int>
bool operator<(const WeightedEdge<VertexType, WeightType>& lhs, const WeightedEdge<VertexType, WeightType>& rhs) {
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
 * @tparam VertexType
 * @tparam EdgeType
 * @tparam WeightType
 * @param vertex_number_ - number of vertices
 * @param edges_number_ - number of edges
 */
template <typename VertexType = int, typename WeightType = int,
          typename EdgeType = WeightedEdge<VertexType>>
class WeightedGraph {
 protected:
  size_t vertex_number_ = 0;
  size_t edges_number_ = 0;

 public:
  explicit WeightedGraph(size_t vert_number, size_t edges_number = 0)
      : vertex_number_(vert_number), edges_number_(edges_number) {}

  [[nodiscard]] virtual size_t GetVertexNumber() const { return vertex_number_; }
  [[nodiscard]] virtual size_t GetEdgesNumber() const { return edges_number_; }

  /*!
   * @param vertex
   * @return set<pair<VertexType, WeightType>> with vertices connected with
   * given vertex by edge.
   */
  virtual typename set<pair<VertexType, WeightType>>::iterator GetNeighboursIt(
          const VertexType& vertex) = 0;

  /*!
   * @param vertex
   * @return set<pair<VertexType, WeightType>> with vertices connected with
   * given vertex by edge.
   */
  virtual set<pair<VertexType, WeightType>> GetNeighbours(
          const VertexType& vertex) = 0;
};

/*!
 * Most convenient representation of graph.
 * @tparam VertexType
 * @tparam EdgeType
 * @tparam WeightType
 */
template <typename VertexType = int, typename WeightType = int,
          typename EdgeType = WeightedEdge<VertexType, WeightType>>
class WeightedAdjList : public WeightedGraph<VertexType, WeightType, EdgeType> {
 public:
  WeightedAdjList(size_t vert_number, const vector<EdgeType>& edges)
      : WeightedGraph<VertexType, WeightType, EdgeType>(vert_number, edges.size()) {
    for (const auto& e : edges) {
      list_[e.From()].insert({e.To(), e.Weight()});
      list_[e.To()].insert({e.From(), e.Weight()});
    }
  }

  set<pair<VertexType, WeightType>> GetNeighbours(
          const VertexType& vertex) override {
    return list_[vertex];
  }

  typename set<pair<VertexType, WeightType>>::iterator GetNeighboursIt(
          const VertexType& vertex) override {
    return list_[vertex].begin();
  }

  vector<EdgeType> GetEdges() {
    vector<EdgeType> result;
    for (const auto& i : list_) {
      for (const auto& j : i.second) {
        result.push_back({i.first, j.first, j.second});
      }
    }
    return result;
  }

 private:
  map<VertexType, set<pair<VertexType, WeightType>>> list_;

  friend bool operator<(const pair<VertexType, WeightType>& lhs,
                        const pair<VertexType, WeightType>& rhs) {
    return lhs.second == rhs.second ? lhs.first < rhs.first
                                    : lhs.second < rhs.second;
  }
};

/*!
 * Alternative and less popular representation of graph, but may be useful in
 * some tasks.
 * @tparam VertexType
 * @tparam EdgeType
 */
template <typename VertexType, typename EdgeType = WeightedEdge<VertexType>>
class WeightedAdjMatrix : public WeightedGraph<VertexType, EdgeType> {
 public:
  WeightedAdjMatrix(size_t vert_number, const vector<EdgeType>& edges)
      : WeightedGraph<VertexType, EdgeType>(vert_number, edges.size()),
        data_(vector<vector<bool>>(vert_number)) {
    for (size_t i = 0; i < vert_number; ++i) {
      data_[i] = vector<bool>(vert_number, false);
    }

    for (const auto& e : edges) {
      data_[e.From()][e.To] = true;
      data_[e.To()][e.From] = true;
    }
  }

  /*!
   * @param vertex
   * @return vector<VertexType> with vertices connected with given vertex by edge.
   */
  vector<VertexType> GetNeighbours(const VertexType& vertex) final;

  /*!
   * @param vertex
   * @return vector<VertexType> with vertices connected with given vertex by edge.
   */
  typename vector<VertexType>::iterator GetNeighboursIt(
          const VertexType& vertex) final {
    return data_[vertex].begin();
  }
 private:
  vector<vector<bool>> data_;
};

template <typename VertexType, typename EdgeType>
vector<VertexType> WeightedAdjMatrix<VertexType, EdgeType>::GetNeighbours(
        const VertexType& vertex) {
  vector<VertexType> result;

  for (VertexType i = 0; i < this->vertex_number_; ++i) {
    if (data_[vertex][i]) {
      result.push_back(i);
    }
  }

  return result;
}

/*!
 * Default weighted tree.
 * @tparam VertexType
 * @tparam EdgeType
 * @tparam WeightType
 */
template <typename VertexType = int, typename WeightType = int,
          typename EdgeType = WeightedEdge<VertexType>>
class WeightedTree : public WeightedAdjList<VertexType, WeightType, EdgeType> {
 public:
  WeightedTree(size_t vert_number, const vector<EdgeType>& edges, WeightType weight)
      : WeightedAdjList<VertexType, WeightType, EdgeType>(vert_number, edges), weight_(weight) {}

  [[nodiscard]] WeightType Weight() const { return weight_; }

 private:
  WeightType weight_;
};

template <typename EdgeType>
void Initialize(vector<EdgeType>& edges, int64_t number_of_rooms,
                int64_t number_of_connections, std::istream& input = std::cin) {
  for (int64_t i = 0; i < number_of_connections; i++) {
    int64_t source = 0;
    int64_t destination = 0;
    int64_t weight = 0;
    input >> source >> destination >> weight;

    edges.push_back({source, destination, weight});
  }
}

template <typename VertexType = int, typename WeightType = int,
          typename EdgeType = WeightedEdge<VertexType, WeightType>>
WeightedTree<VertexType, WeightType, EdgeType> MstPrim(WeightedAdjList<VertexType, WeightType, EdgeType>& graph);

template <typename VertexType, typename WeightType,
          typename EdgeType>
WeightedTree<VertexType, WeightType, EdgeType> MstKruskal(WeightedAdjList<VertexType, WeightType, EdgeType>& graph);

template <typename ElemType>
class DSU {
 public:
  DSU(size_t size) {
    data_ = vector<ElemType>(size);
  }

  ElemType& operator[] (size_t index) {
    return data_[index];
  }

  const ElemType& at(size_t index) const {
    return data_.at(index);
  }

  ElemType& Get(const ElemType& vert) {
    return (vert == data_[vert]) ? vert : (data_[vert] = Get(data_[vert]));
  }

  void Unite(const ElemType& lhs, const ElemType& rhs) {
    lhs = Get(lhs);
    rhs = Get(rhs);

    if (lhs != rhs) data_[lhs] = rhs;
  }

 private:
  vector<ElemType> data_;
};

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  std::cout.tie(nullptr);

  size_t number_of_vertices = 0;
  int64_t number_of_edges = 0;

  std::cin >> number_of_vertices >> number_of_edges;

  std::vector<WeightedEdge<int64_t, int64_t>> edges;

  Initialize(edges, number_of_vertices, number_of_edges);

  WeightedAdjList<int64_t, int64_t> graph(number_of_vertices, edges);

  auto tree = MstPrim(graph);
  //auto tree = MstKruskal(graph);
  std::cout << tree.Weight();

  return 0;
}

/*!
 * @function MstPrim(WeightedAdjList<VertexType, WeightType, EdgeType>& graph)
 * @brief Creates MST for graph
 * @tparam VertexType
 * @tparam WeightType
 * @tparam EdgeType
 * @param graph
 * @return WeightedTree - MST for the given graph
 */
template <typename VertexType, typename WeightType,
          typename EdgeType>
WeightedTree<VertexType, EdgeType, WeightType> MstPrim(WeightedAdjList<VertexType, WeightType, EdgeType>& graph) {
  const int64_t kInf = std::numeric_limits<WeightType>::max();
  int64_t size = graph.GetVertexNumber();
  vector<int64_t> dist(size + 2, kInf);
  vector<EdgeType> prev(size + 2);
  dist[0] = 0;

  set<EdgeType> queue;
  for (size_t i = 0; i <= size; ++i) {
    queue.insert({0, i, dist[i]});
  }

  while (!queue.empty()) {
    EdgeType start = *(queue.begin());
    int64_t v = start.To();
    queue.erase(start);

    for (auto& i : graph.GetNeighbours(v)) {
      EdgeType current = {v, i.first, dist[i.first]};
      if (queue.find(current) != queue.end() && i.second < dist[i.first]) {
        queue.erase(current);
        prev[i.first] = {v, i.first, i.second};
        dist[i.first] = i.second;
        queue.insert({v, i.first, dist[i.first]});
      }
    }
  }

  auto res = std::accumulate(prev.begin(), prev.end(), 0,
                             [] (int64_t sum, const auto& item) { return sum + item.Weight(); });

  return WeightedTree<VertexType, EdgeType>(size, prev, res);
}

/*!
 * @function MstKruskal(WeightedAdjList<VertexType, WeightType, EdgeType>& graph)
 * @brief Creates MST for graph
 * @tparam VertexType
 * @tparam WeightType
 * @tparam EdgeType
 * @param graph
 * @return WeightedTree - MST for the given graph
 */
template <typename VertexType, typename WeightType,
          typename EdgeType>
WeightedTree<VertexType, WeightType, EdgeType> MstKruskal(WeightedAdjList<VertexType, WeightType, EdgeType>& graph) {
  auto edges = graph.GetEdges();
  std::sort(edges.begin(), edges.end());

  WeightType result;
  if constexpr(std::is_arithmetic_v<WeightType>) { result = 0; }

  vector<EdgeType> result_edges(graph.GetVertexNumber() + 1);
  DSU<VertexType> sets(graph.GetVertexNumber() + 1);
  for (size_t i = 0; i <= graph.GetVertexNumber(); ++i) {
    sets[i] = i;
  }

  if (graph.GetEdgesNumber() != 0) {
    for (size_t i = 0; i <= graph.GetEdgesNumber(); ++i) {
      auto from = edges[i].From();
      auto to = edges[i].To();
      auto weight = edges[i].Weight();

      if (sets.Get(from) != sets.Get(to)) {
        result += weight;
        result_edges[i] = {from, to, weight};
        sets.Unite(from, to);
      }
    }
  }

  return WeightedTree<VertexType, WeightType, EdgeType>(graph.GetVertexNumber(), result_edges, result);
}
