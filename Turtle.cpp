#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::pair;
using std::unordered_map;
using std::vector;

/*!
 * Default edge for working with graphs.
 * Inherited by std::pair
 *
 * @tparam T
 */
template <typename T = int>
struct DefaultEdge : pair<T, T> {
  DefaultEdge(const T& first, const T& second)
      : std::pair<T, T>(first, second) {}
  using BaseClass = pair<T, T>;
  [[nodiscard]] const T& From() const { return BaseClass::first; }
  [[nodiscard]] const T& To() const { return BaseClass::second; }
};

/*!
 * @brief class for working with applied tasks needs mathematical objects
 * graphs. It can be created like adjacency matrix or list.
 * @tparam Vertex
 * @tparam Edge
 * @param vertex_number_ - number of vertices
 * @param edges_number_ - number of edges
 */
template <typename Vertex = int, typename Edge = DefaultEdge<Vertex>>
class Graph {
 public:
  using VertexType = Vertex;
  using EdgeType = Edge;

  explicit Graph(size_t vert_number, size_t edges_number = 0)
      : vertex_number_(vert_number),
        edges_number_(edges_number){}

  [[nodiscard]] virtual size_t GetVertexNumber() const { return vertex_number_; }
  [[nodiscard]] virtual size_t GetEdgesNumber() const { return edges_number_; }

  /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
  [[nodiscard]] virtual typename vector<Vertex>::iterator GetNeighboursIt(
      const Vertex& vertex) = 0;
  [[nodiscard]] virtual typename vector<Vertex>::const_iterator GetNeighboursIt(
          const Vertex& vertex) const = 0;


  /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
  virtual vector<Vertex> GetNeighbours(const Vertex& vertex) = 0;
  virtual const vector<Vertex>& GetNeighbours(const Vertex& vertex) const = 0;

protected:
  size_t vertex_number_ = 0;
  size_t edges_number_ = 0;
};

/*!
 * Most convenient representation of graph.
 * @tparam VertexType
 * @tparam EdgeType
 */
template <typename VertexType = int,
          typename EdgeType = DefaultEdge<VertexType>>
class AdjList : public Graph<VertexType, EdgeType> {
 public:
  AdjList(size_t vert_number, const vector<EdgeType>& edges)
      : Graph<VertexType, EdgeType>(vert_number, edges.size()) {
    for (const auto& edge: edges) {
      list_[edge.From()].push_back(edge.To());
      list_[edge.To()].push_back(edge.From());
    }
  }

  vector<VertexType> GetNeighbours(const VertexType& vertex) final {
    return list_[vertex];
  }

  const vector<VertexType>& GetNeighbours(const VertexType& vertex) const final {
    try {
      list_.at(vertex);
    } catch (const std::out_of_range& e) {
      std::cout << "Vertex does not exists!" << std::endl;
    }
    return list_.at(vertex);
  }

  [[nodiscard]] typename vector<VertexType>::iterator GetNeighboursIt(
      const VertexType& vertex) final {
    return list_[vertex].begin();
  }

  [[nodiscard]] typename vector<VertexType>::const_iterator GetNeighboursIt(
          const VertexType& vertex) const final {
    try {
      list_.at(vertex).begin();
    } catch (const std::out_of_range& e) {
      std::cout << "Vertex does not exists!" << std::endl;
    }
    return list_.at(vertex).begin();
  }

 private:
  unordered_map<VertexType, vector<VertexType>> list_;
};

/*!
 * Alternative and less popular representation of graph, but may be useful in
 * some tasks.
 * @tparam VertexType
 * @tparam EdgeType
 */
template <typename VertexType, typename EdgeType = DefaultEdge<VertexType>>
class AdjMatrix : public Graph<VertexType, EdgeType> {
 public:
  AdjMatrix(size_t vert_number, const vector<EdgeType>& edges)
      : Graph<VertexType, EdgeType>(vert_number, edges.size()),
        data_(vector<vector<bool>>(vert_number)) {
    for (size_t i = 0; i < vert_number; ++i) {
      data_[i] = vector<bool>(vert_number);
    }

    for (const auto& edge: edges) {
      data_[edge.From()][edge.To] = true;
      data_[edge.To()][edge.From] = true;
    }
  }

  /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
  vector<VertexType> GetNeighbours(const VertexType& vertex) final;
  const vector<VertexType>& GetNeighbours(const VertexType& vertex) const final;

  /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
  [[nodiscard]] typename vector<VertexType>::iterator GetNeighboursIt(
      const VertexType& vertex) final {
    return data_[vertex].begin();
  }
  [[nodiscard]] typename vector<VertexType>::const_iterator GetNeighboursIt(
          const VertexType& vertex) const final {
    try {
      data_.at(vertex).begin();
    } catch (const std::out_of_range& e) {
      std::cout << "Vertex does not exists!" << std::endl;
    }

    return data_.at(vertex).begin();
  }

 private:
  vector<vector<bool>> data_;
};

template <typename VertexType, typename EdgeType>
vector<VertexType> AdjMatrix<VertexType, EdgeType>::GetNeighbours(
    const VertexType& vertex) {
  vector<bool> result;

  for (size_t i = 0; i < this->vertex_number_; ++i) {
    if (data_[vertex][i]) {
      result.push_back(true);
    }
  }

  return result;
}

template <typename VertexType, typename EdgeType>
const vector<VertexType>& AdjMatrix<VertexType, EdgeType>::GetNeighbours(
        const VertexType& vertex) const {
  vector<VertexType> result;

  for (size_t i = 0; i < this->vertex_number_; ++i) {
    if (data_.at(vertex)[i]) {
      result.push_back(i);
    }
  }

  return result;
}

/*!
 * @brief Breadth First Search algorithm for graph using Visitor class.
 * @tparam Vertex
 * @tparam Graph
 * @tparam Visitor
 * @param origin_vertex - start
 * @param graph - graph for round
 * @param visitor - helper class saving the result of traversal.
 */
template <class Vertex, class Graph, class Visitor>
void BreadthFirstSearch(Vertex origin_vertex, Graph& graph, Visitor& visitor) {
  std::queue<Vertex> bfs_queue;
  std::unordered_set<Vertex> visited_vertices;

  bfs_queue.push(origin_vertex);
  visited_vertices.insert(origin_vertex);

  while (!bfs_queue.empty()) {
    auto cur_vertex = bfs_queue.front();
    bfs_queue.pop();
    for (auto& neighbour : graph.GetNeighbours(cur_vertex)) {
      if (visited_vertices.find(neighbour) == visited_vertices.end()) {
        visitor.TreeEdge({cur_vertex, neighbour});
        bfs_queue.push(neighbour);
        visited_vertices.insert(neighbour);
      }
    }
  }
}

/*!
 * @brief Special class for saving result of BFS round.
 * @tparam VertexType
 * @tparam EdgeType
 * @param ancestors_ - unordored map with results.
 */
template <typename VertexType, typename EdgeType>
class BfsVisitor {
 public:
  virtual void DiscoverVertex(VertexType) /*unused*/ = 0;    // NOLINT
  virtual void ExamineEdge(const EdgeType&) /*unused*/ = 0;  // NOLINT
  virtual void ExamineVertex(VertexType) /*unused*/ = 0;     // NOLINT
  virtual void TreeEdge(const EdgeType& edge) = 0;
  virtual ~BfsVisitor() = default;
};

/*!
 * @brief Special class for saving result of BFS round.
 * @tparam VertexType
 * @tparam EdgeType
 * @param ancestors_ - unordered map with BFS results.
 */
template <typename VertexType, typename EdgeType>
class AncestorBfsVisitor : public BfsVisitor<VertexType, EdgeType> {
 public:
  virtual void DiscoverVertex(VertexType) /*unused*/ {}    // NOLINT
  virtual void ExamineEdge(const EdgeType&) /*unused*/ {}  // NOLINT
  virtual void ExamineVertex(VertexType) /*unused*/ {}     // NOLINT
  void TreeEdge(const EdgeType& edge) override {
    ancestors_[edge.To()] = edge.From();
  }

  const unordered_map<VertexType, VertexType>& GetAncMap() const { return ancestors_; }

  ~AncestorBfsVisitor() override = default;

 private:
  unordered_map<VertexType, VertexType> ancestors_;
};

template<class Container>
void Print(const Container& container, const std::string& sep = " ") {
  for (const auto& elem : container) {
    std::cout << elem << sep;
  }
}

namespace CurrentTask {
  namespace CurrentTaskConstants {
    const int kMaxPlanets = 1e5;
  }
  
  /*!
 * @brief Default function filling vector of edges by keyboard.
 * @param edges
 * @param number_of_pairs
 */
  void Initialize(vector<DefaultEdge<int>>& edges, int number_of_pairs) {
    for (int i = 0; i < number_of_pairs; ++i) {
      int from = 0;
      int to = 0;
      std::cin >> from >> to;
      
      edges.emplace_back(from, to);
    }
  }
  
  void PrintResult(vector<int>& result, int planet_counter) {
    if (planet_counter == CurrentTaskConstants::kMaxPlanets) {
      std::cout << -1;
      return;
    }
    std::cout << result.size() - 1 << '\n';
    
    Print(result);
  }
  
  vector<int> GetPath(const unordered_map<int, int>& ancestors_map, int start, int finish, int& counter) {
    vector<int> result;
    
    if (!ancestors_map.count(finish)) {
      counter = CurrentTaskConstants::kMaxPlanets;
      return result;
    }

    result.push_back(finish);
    counter = 0;
    while (finish != start) {
      try {
        finish = ancestors_map.at(finish);
      } catch (const std::out_of_range& e) {
        std::cout << "Vertex does not exists!" << std::endl;
      }

      result.push_back(finish);
      ++counter;
    }
    std::reverse(result.begin(), result.end());

    return result;
  }
}

int main() {
  int number_of_planets = 0;
  int number_of_pairs = 0;
  std::cin >> number_of_planets >> number_of_pairs;

  int start = 0;
  int finish = 0;
  std::cin >> start >> finish;

  vector<DefaultEdge<int>> edges;
  CurrentTask::Initialize(edges, number_of_pairs);
  AdjList<int, DefaultEdge<int>> graph(number_of_planets, edges);

  AncestorBfsVisitor<int, DefaultEdge<int>> visitor;
  BreadthFirstSearch(start, graph, visitor);
  auto ancestors_map = visitor.GetAncMap();

  int planet_counter = 0;
  vector<int> result = CurrentTask::GetPath(ancestors_map, start, finish, planet_counter);

  CurrentTask::PrintResult(result, planet_counter);

  return 0;
}
