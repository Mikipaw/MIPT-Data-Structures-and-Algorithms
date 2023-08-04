#include <algorithm>
#include <iostream>
#include <list>
#include <queue>
#include <string>
#include <vector>

const int kAlphabetSize = 26;

class PrefixNode {
 private:
  PrefixNode** next = nullptr;
  PrefixNode* suff_ref = nullptr;
  PrefixNode* end_ref = nullptr;

  PrefixNode* parent = nullptr;

  friend class PrefixTree;
  std::list<int> end_line_num;
  char symbol = 0;

 public:
  PrefixNode() = default;

  PrefixNode* AddSymbol(char sym);
  PrefixNode* GetSymbol(char sym);

  [[nodiscard]] bool HasSymbol(char sym) const;
  [[nodiscard]] static int GetNextNodesSize() { return kAlphabetSize + 1; }
  [[nodiscard]] bool IsEnd() const { return HasSymbol('\0'); }
  [[nodiscard]] bool IsRoot() const { return this->parent == this; }

  ~PrefixNode() {
    if (next == nullptr) return;

    for (int i = 0; i < GetNextNodesSize(); ++i) {
      delete this->next[i];
    }

    delete[] next;
  }
};

class PrefixTree {
  PrefixNode root;
  PrefixNode* current{};

 public:
  PrefixTree() {
    root.parent = &root;
    root.suff_ref = &root;
    root.end_ref = &root;
  }

  void AddWord(const std::string& str, int line_num);
  void CountSuffixRefs();
  void CountEndRefs();

  std::list<size_t> ProcessSymbol(char sym);

  void Preparing();
  void ProcessString(const std::string& str, std::vector<std::list<size_t>>& pat,
                     std::vector<size_t>& pat_sizes);
};

template<class Container>
void Print(const Container& container, const std::string& sep = " ");
void InputPatterns(PrefixTree& tree, std::vector<size_t>& sizes, int number_of_words);

int main() {
  std::string input_text;
  std::cin >> input_text;

  int number_of_words = 0;
  std::cin >> number_of_words;

  std::vector<size_t> sizes(number_of_words);
  PrefixTree tree;
  InputPatterns(tree, sizes, number_of_words);

  std::vector<std::list<size_t>> result(number_of_words);
  tree.Preparing();
  tree.ProcessString(input_text, result, sizes);

  for (const auto& item: result) {
    std::cout << item.size() << ' ';
    Print(item);
    std::cout << '\n';
  }

  return 0;
}

PrefixNode* PrefixNode::AddSymbol(char sym) {
  if (next == nullptr) {
    next = new PrefixNode*[kAlphabetSize + 1];

    for (int i = 0; i < kAlphabetSize + 1; ++i)
      next[i] = nullptr;
  }

  auto symbol_num = sym ? sym - 'a' : kAlphabetSize;

  if (next[symbol_num] == nullptr) next[symbol_num] = new PrefixNode;
  next[symbol_num]->symbol = sym;
  next[symbol_num]->parent = this;

  return next[symbol_num];
}

PrefixNode* PrefixNode::GetSymbol(char sym) {
  if (!sym) return this->next[kAlphabetSize];
  else return this->next[sym - 'a'];
}

bool PrefixNode::HasSymbol(char sym) const {
  if (next == nullptr) return false;

  if (!sym) return next[kAlphabetSize] != nullptr;

  return next[sym - 'a'] != nullptr;
}


void PrefixTree::AddWord(const std::string& str, int line_num) {
  current = &root;

  for (char sym : str) current = current->AddSymbol(sym);

  current->AddSymbol('\0');
  current->end_line_num.push_back(line_num);
  current = &root;
}

void PrefixTree::CountSuffixRefs() {
  std::queue<PrefixNode*> bfs_queue;
  bfs_queue.push(&root);

  while (!bfs_queue.empty()) {
    auto node = bfs_queue.front();
    bfs_queue.pop();

    if (node->parent->IsRoot())
      node->suff_ref = node->parent;
    else {
      auto supref_ref = node->parent->suff_ref;

      while (!supref_ref->HasSymbol(node->symbol) && !supref_ref->IsRoot())
        supref_ref = supref_ref->suff_ref;

      if (supref_ref->HasSymbol(node->symbol))
        node->suff_ref = supref_ref->GetSymbol(node->symbol);
      else
        node->suff_ref = supref_ref;
    }

    for (int i = 0; i < node->GetNextNodesSize(); ++i)
      if (node->next != nullptr && node->next[i] != nullptr)
        bfs_queue.push(node->next[i]);
  }
}

void PrefixTree::CountEndRefs() {
  std::queue<PrefixNode*> bfs_queue;
  bfs_queue.push(&root);
  while (!bfs_queue.empty()) {
    auto node = bfs_queue.front();
    bfs_queue.pop();

    if (node->suff_ref->IsEnd())
      node->end_ref = node->suff_ref;
    else
      node->end_ref = node->suff_ref->end_ref;

    for (int i = 0; i < node->GetNextNodesSize(); ++i)
      if (node->next != nullptr && node->next[i] != nullptr)
        bfs_queue.push(node->next[i]);
  }
}

std::list<size_t> PrefixTree::ProcessSymbol(char sym) {
  while (!current->HasSymbol(sym)) {
    if (current->IsRoot())
      return {};

    current = current->suff_ref;
  }

  current = current->GetSymbol(sym);
  std::list<size_t> ended_lines;

  if (current->IsEnd()) {
    for (const auto& item : current->end_line_num) {
      ended_lines.push_back(item);
    }
  }

  auto end_ref = current->end_ref;
  while (!end_ref->IsRoot()) {
    for (const auto& item : end_ref->end_line_num) {
      ended_lines.push_back(item);
    }
    end_ref = end_ref->end_ref;
  }

  return ended_lines;
}

void PrefixTree::Preparing() {
  CountSuffixRefs();
  CountEndRefs();
}

void PrefixTree::ProcessString(const std::string& str, std::vector<std::list<size_t>>& pat,
                               std::vector<size_t>& pat_sizes) {
  this->current = &root;

  auto size = str.size();
  for (size_t i = 0; i < size; ++i) {
    auto res = ProcessSymbol(str[i]);
    for (const auto& it: res) {
      pat[it].push_back(i - pat_sizes[it] + 2);
    }
  }
}

template<class Container>
void Print(const Container& container, const std::string& sep) {
  for (const auto& elem : container) {
    std::cout << elem << sep;
  }
}

void InputPatterns(PrefixTree& tree, std::vector<size_t>& sizes, int number_of_words) {
  for (int i = 0; i < number_of_words; ++i) {
    std::string pattern;
    std::cin >> pattern;
    tree.AddWord(pattern, i);
    sizes[i] = pattern.size();
  }
}
