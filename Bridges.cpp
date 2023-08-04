#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <stack>
#include <unordered_map>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::stack;
using std::unordered_map;
using std::vector;

namespace graph {
  enum Colors { WHITE, GREY, BLACK };

  /*!
 * Default edge for working with graphs.
 * Inherited by std::pair
 *
 * @tparam T
 */
  template <typename T = int>
  struct DefaultEdge : pair<T, T> {
    using BaseClass = pair<T, T>;
    DefaultEdge(const T& from, const T& to)
        : std::pair<T, T>(from, to) {}
    [[nodiscard]] const T& From() const { return BaseClass::first; }
    [[nodiscard]] const T& To() const { return BaseClass::second; }
  };

  /*!
 * @brief Special class for saving result of BFS round.
 * @tparam VertexType
 * @tparam EdgeType
 * @param ancestors_ - unordored map with results.
 */
  template <typename VertexType, typename EdgeType>
  class DfsVisitor {
  public:
    virtual void Update(const VertexType&, const VertexType&, bool flag = false) = 0;
    virtual void ExamineEdge(const EdgeType&) = 0;
    virtual void ExamineVertex(VertexType) = 0;
    virtual ~DfsVisitor() = default;
  };

  /*!
 * @brief Special class for saving result of DFS round.
 * @tparam VertexType
 * @tparam EdgeType
 * @param ancestors_ - unordored map with results.
 */
  template <typename VertexType = int, typename EdgeType = DefaultEdge<int>>
  class BridgeDfsVisitor final : public DfsVisitor<VertexType, EdgeType> {
  public:
    explicit BridgeDfsVisitor(int number_of_vertices)
        : bridges({}),
          t_in(number_of_vertices + 1),
          ret(number_of_vertices + 1) {}

    void ExamineEdge(const EdgeType& edge) override {
      if (ret[edge.To()] > t_in[edge.From()]) {
        auto edge_num = GetEdgeNumber(edge.From(), edge.To());
        if (edge_num != -1) {
          bridges.insert(GetEdgeNumber(edge.To(), edge.From()));
        }
      }
    }

    void ExamineVertex(VertexType current) final {
      ret[current] = timer;
      t_in[current] = timer++;
    }

    void Update(const VertexType& from, const VertexType& to, bool flag) override {
      if (!flag) {
        ret[from] = std::min(ret[from], t_in[to]);
      } else {
        ret[from] = std::min(ret[from], ret[to]);
      }
    }

    [[nodiscard]] set<int> GetBridges() const { return bridges; }

    bool IsBridge(const EdgeType& edge) {
      return bridges.count(edges[edge.first][edge.second].first);
    }

    [[nodiscard]] int GetEdgeNumber(const VertexType& from, const VertexType& to) {
      if ((*edges)[from].count(to) && (*edges)[from][to].second > 1) {
        return -1;
      }
      return (*edges)[from][to].first;
    }

    ~BridgeDfsVisitor() final {
      delete edges;
    };

    vector<unordered_map<VertexType, pair<int, int>>>* edges;

   private:
    int timer = 1;
    vector<int> t_in{};
    vector<int> ret{};
    set<int> bridges{};
  };

  /*!
 * @brief class for working with applied tasks needs mathematical objects
 * graphs. It can be created like adjacency matrix or list.
 * @tparam Vertex
 * @tparam Edge
 * @param vertex_number_ - number of vertices
 * @param edges_number_ - number of edges
 */
  template <typename VertexType = int, typename EdgeType = DefaultEdge<VertexType>>
  class Graph {
  public:
    explicit Graph(size_t vert_number, size_t edges_number = 0)
        : vertex_number_(vert_number),
          edges_number_(edges_number){}

    [[nodiscard]] virtual size_t GetVertexNumber() const {
      return vertex_number_;
    }

    /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
    [[nodiscard]] virtual typename vector<VertexType>::iterator GetNeighboursIt(
            const VertexType& vertex) = 0;
    [[nodiscard]] virtual typename vector<VertexType>::const_iterator GetNeighboursIt(
            const VertexType& vertex) const = 0;

    /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
    virtual vector<VertexType> GetNeighbours(const VertexType& vertex) = 0;
    virtual const vector<VertexType>& GetNeighbours(const VertexType& vertex) const = 0;

  private:
    size_t vertex_number_;
    size_t edges_number_;
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
    AdjList(size_t vert_number, vector<EdgeType> edges,
            BridgeDfsVisitor<VertexType, EdgeType>& bdv)
        : Graph<VertexType, EdgeType>(vert_number, edges.size()) {
      list_ = vector<vector<VertexType>>(vert_number);
      bdv.edges = new vector<unordered_map<VertexType, pair<int, int>>>(vert_number);

      auto number_of_edges = edges.size();
      for (int i = 0; i < number_of_edges; ++i) {
        auto [from, to] = edges[i];
        list_[from].push_back(to);
        list_[to].push_back(from);

        if (bdv.edges[0][from].count(to)) {
          if (bdv.edges[0][from][to].second <= 1) {
            bdv.edges[0][from][to].second++;
            bdv.edges[0][to][from].second++;
          }
          ++edge_number_;
        } else {
          bdv.edges[0][from].insert({to, {edge_number_, 1}});
          bdv.edges[0][to].insert({from, {edge_number_++, 1}});
        }
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

  protected:
    size_t vertex_number_ = 0;
    size_t edges_number_ = 0;

  private:
    vector<vector<VertexType>> list_;
    int edge_number_ = 1;
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

      for (const auto& e : edges) {
        data_[e.From()][e.To] = true;
        data_[e.To()][e.From] = true;
      }
    }

    /*!
   * @param vertex
   * @return vector<Vertex> with vertices connected with given vertex by edge.
   */
    vector<VertexType> GetNeighbours(VertexType vertex) final;
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
          VertexType vertex) {
    vector<VertexType> result;

    for (VertexType i = 0; i < this->vertex_number_; ++i) {
      if (data_[vertex][i]) {
        result.push_back(static_cast<VertexType>(i));
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
        result.push_back(static_cast<VertexType>(i));
      }
    }

    return result;
  }

  /*!
 * @brief Default finction for initialize vector of edges by keyboard.
 * @param edges
 * @param number_of_pairs
 */
  void Initialize(vector<DefaultEdge<int>>& edges, int number_of_pairs) {
    for (int i = 0; i < number_of_pairs; i++) {
      int from = 0;
      int to = 0;
      std::cin >> from >> to;

      edges.emplace_back(from, to);
    }
  }
  
  namespace traverses {
    template <typename VertexType = int,
              typename EdgeType = DefaultEdge<VertexType>>
    void Dfs(const Graph<VertexType>& graph, vector<char>& used,
             DfsVisitor<VertexType, DefaultEdge<VertexType>>& visitor,
             VertexType current, VertexType prev) {
      used[current] = GREY;
      visitor.ExamineVertex(current);

      for (auto& i : graph.GetNeighbours(current)) {
        if (i == prev) {
          continue;
        }

        if (used[i] != WHITE) {
          visitor.Update(current, i, false);
        } else {
          Dfs(graph, used, visitor, i, current);
          visitor.Update(current, i, true);

          visitor.ExamineEdge({current, i});
        }
      }
      used[current] = BLACK;
    }
  }

  template <typename VertexType = int,
            typename EdgeType = DefaultEdge<VertexType>>
  void FindBridges(const AdjList<VertexType>& graph,
                    BridgeDfsVisitor<VertexType>& visitor) {
    auto size = graph.GetVertexNumber();
    vector<char> used(size + 1, WHITE);

    for (int i = 0; i < static_cast<int>(size); ++i) {
      traverses::Dfs(graph, used, visitor, i, -1);
    }
  }

  template <typename VertexType = int,
            typename EdgeType = DefaultEdge<VertexType>>
  void PrintBridges(const BridgeDfsVisitor<VertexType>& visitor) {
    std::cout << visitor.GetBridges().size() << '\n';
    for (auto& bridge: visitor.GetBridges()) {
      std::cout << bridge << '\n';
    }
  }
}

int main() {
  using namespace graph;

  int number_of_vertices = 0;
  int number_of_edges = 0;
  std::cin >> number_of_vertices >> number_of_edges;
  vector<DefaultEdge<int>> edges;
  Initialize(edges, number_of_edges);

  BridgeDfsVisitor<int, DefaultEdge<int>> visitor(++number_of_vertices);
  AdjList<int, DefaultEdge<int>> graph(number_of_vertices, edges,
                                       visitor);

  FindBridges(graph, visitor);
  PrintBridges(visitor);

  return 0;
}
