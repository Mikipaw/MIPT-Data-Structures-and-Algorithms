#include <iostream>
#include <map>
#include <vector>

using std::map;
using std::vector;

size_t Sum(size_t lhs, size_t rhs, size_t friendship) {
  return lhs + rhs + friendship;
}

template <typename FunctionType = size_t>
class DSU_visitor {
 public:
  DSU_visitor(FunctionType (*func)(FunctionType lhs, FunctionType rhs, FunctionType), size_t size) : function_(func) {
    for (size_t i = 0; i < size; ++i) mp_[i];
  };

  void CallFunction(size_t lhs, size_t rhs, const FunctionType& item) {
    mp_[lhs] = (lhs == rhs ? function_(mp_[lhs], mp_[0], item) : function_(mp_[lhs], mp_[rhs], item));
  }

  [[nodiscard]] FunctionType Get (size_t idx) const { return mp_.at(idx); }

 private:
  FunctionType (*function_)(FunctionType lhs, FunctionType rhs, FunctionType);
  map<size_t, FunctionType> mp_;
};

template <typename VisitorType = size_t, typename Visitor = DSU_visitor<VisitorType>>
class DSU {
 public:
  explicit DSU(int64_t size, VisitorType (*func)(VisitorType lhs, VisitorType rhs, VisitorType)) :
                               parent_(vector<int64_t>(size + 1)),
                               size_(vector<int64_t>(size + 1)),
                               visitor_(func, size + 1) {
    for (size_t i = 1; i <= size; i++) {
      parent_[i] = i;
      size_[i] = 1;
    }
  };

  VisitorType Get(const VisitorType& student);
  [[nodiscard]] VisitorType GetWeight(const VisitorType& student) const { return visitor_.Get(student); }

  void Unite(size_t lhs, size_t rhs, VisitorType& friendship_lvl);

 private:
  vector<int64_t> parent_;
  vector<int64_t> size_;
  Visitor visitor_;
};

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  std::cout.tie(nullptr);

  size_t number_of_vertices = 0;
  size_t number_of_edges = 0;
  std::cin >> number_of_vertices >> number_of_edges;

  DSU<size_t> DSU(++number_of_vertices, Sum);

  for (size_t i = 0; i < number_of_edges; i++) {
    char command = 0;
    std::cin >> command;

    if (command == '1') {
      size_t first_student = 0;
      size_t second_student = 0;
      size_t friendship_lvl = 0;
      std::cin >> first_student >> second_student >> friendship_lvl;

      DSU.Unite(first_student, second_student, friendship_lvl);
    } else {
      size_t student = 0;
      std::cin >> student;
      std::cout << DSU.GetWeight(DSU.Get(student)) << '\n';
    }
  }

  return 0;
}

template <typename VisitorType, typename Visitor>
VisitorType DSU<VisitorType, Visitor>::Get(const VisitorType& student) {
  if (student == parent_[student]) {
    return student;
  }
  return parent_[student] = Get(parent_[student]);
}

template <typename VisitorType, typename Visitor>
void DSU<VisitorType, Visitor>::Unite(size_t lhs, size_t rhs, VisitorType& friendship_lvl) {
  lhs = Get(lhs);
  rhs = Get(rhs);

  if (lhs != rhs) {
    if (size_[lhs] < size_[rhs]) {
      std::swap(lhs, rhs);
    }

    size_[lhs] += size_[rhs];
    parent_[rhs] = lhs;
  }
  visitor_.CallFunction(lhs, rhs, friendship_lvl);
}