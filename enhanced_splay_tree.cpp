#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <unordered_map>
#include <vector>

#define normal false
#define modified true

static uint64_t prime = ((uint64_t)1 << 63) - 59;
static uint64_t base = 255;

struct node {
  node *parent;
  node *left;
  node *right;
  uint32_t subtree_size;
  unsigned char character;
  uint64_t kr;
  uint64_t base_exp;

  node(char c)
      : parent(nullptr), left(nullptr), right(nullptr), subtree_size(1),
        character(c), kr(c), base_exp(base) {}

  node(node *n) {
    parent = n->parent;
    left = n->left;
    right = n->right;
    subtree_size = n->subtree_size;
    character = n->character;
    kr = n->kr;
    base_exp = n->base_exp;
  }

  // Destructor for node
  ~node() {
    if (left)
      delete left;
    if (right)
      delete right;
  }
};

// overload << operator for node
std::ostream &operator<<(std::ostream &os, const node *n) {
  if (!n) {
    os << "NULL NODE";
    return os;
  }
  os << &n << " " << n->character << " : " << n->subtree_size << " - " << n->kr
     << " p " << n->base_exp;
  return os;
}

class SplayTree {
  node *root;

public:
  SplayTree() : root(nullptr) {}

  SplayTree(const std::string &str) {
    root = buildBalancedTree(str, 0, str.size() - 1);
  }

  SplayTree(node *n) : root(n) {}

  // Destructor
  ~SplayTree() {
    if (root)
      delete root;
  }

  uint32_t getSubtreeSize(node *n) { return n->subtree_size; }

  // Insert a new character at a given position using only split operations
  void insert(const unsigned char c, const uint32_t position) {
    if (position == 0) {
      node *n = new node(c);
      n->right = root;
      if (root)
        root->parent = n;
      root = n;
    } else if (position == root->subtree_size) {
      node *n = new node(c);
      n->left = root;
      if (root)
        root->parent = n;
      root = n;
    } else {
      SplayTree rightTree = nullptr;
      extract(position, getSubtreeSize(root) - 1, rightTree);
      node *n = new node(c);
      SplayTree newTree(n);
      introduce(getSubtreeSize(root), newTree);
      introduce(getSubtreeSize(root), rightTree);
    }
    updateNodeInfo(root);
  }

  void edit(const unsigned char c, const uint32_t position) {
    if (position > root->subtree_size) {
      std::cout << "Invalid position for edit" << std::endl;
      return;
    }
    get<normal>(position);
    root->character = c;
    updateNodeInfo(root);
  }

  // get a node at a given position
  template <auto isModified> node *get(const int position) {
    node *z = root;
    int currentPos = position;

    while (z) {
      int leftSize = z->left ? z->left->subtree_size : 0;
      if (currentPos == leftSize) {
        splay<isModified>(z);
        return root;
      } else if (currentPos < leftSize) {
        z = z->left;
      } else {
        currentPos -= leftSize + 1;
        z = z->right;
      }
    }
    return nullptr;
  }

  // Delete implementation using split and join
  void delete_node(const int position) {
    if (position == 0) {
      get<normal>(0);
      root = root->right;
      root->parent = nullptr;
      updateNodeInfo(root);
      return;
    }

    SplayTree leftTree, rightTree;
    extract(0, position - 1, leftTree);
    extract(1, getSubtreeSize(root) - 1, rightTree);
    introduce(0, leftTree);
    introduce(getSubtreeSize(root), rightTree);
    /* split(position - 1, leftTree, rightTree); */

    /* rightTree.split(0, deletedNode, rightTree); */

    /* leftTree.join(rightTree); */
    /* root = leftTree.root; */
  }

  // Introduce a substring from another tree at a given position
  void introduce(const uint32_t pos, SplayTree &other) {
    if (!root) {
      if (other.root)
        root = other.root;
      other.root = nullptr;
    } else if (pos == root->subtree_size) {
      get<normal>(pos - 1);
      root->right = other.root;
      root->right->parent = root;
      other = nullptr;
      updateNodeInfo(root);
    } else if (pos < root->subtree_size && pos > 0) {
      isolate(pos, pos - 1);

      if (root->right) {
        root->right->left = other.root;
        root->right->left->parent = root->right;
        other.root = nullptr;

        updateNodeInfo(root->right);
        updateNodeInfo(root);
      } else {
        root->right = other.root;
        root->right->parent = root;
        other.root = nullptr;
        updateNodeInfo(root);
      }
    } else if (pos == 0) {
      get<normal>(other.root->subtree_size - 1);
      other.root->right = root;
      other.root->right->parent = other.root;
      root = other.root;
      other.root = nullptr;
    }
  }

  /* // Introduce a substring from another tree at a given position */
  /* void introduce(const uint32_t pos, node *other) { */
  /*   if (!root) { */
  /*     if (other->character) */
  /*       root = other; */
  /*   } else if (pos == root->subtree_size) { */
  /*     get<normal>(pos - 1); */
  /*     root->right = other; */
  /*     root->right->parent = root; */
  /*     other = nullptr; */
  /*     updateNodeInfo(root); */
  /*   } else if (pos < root->subtree_size && pos > 0) { */
  /*     isolate(pos, pos - 1); */

  /*     if (root->right) { */
  /*       root->right->left = other; */
  /*       root->right->left->parent = root->right; */

  /*       updateNodeInfo(root->right); */
  /*       updateNodeInfo(root); */
  /*     } else { */
  /*       root->right = other; */
  /*       root->right->parent = root; */
  /*       updateNodeInfo(root); */
  /*     } */
  /*   } else if (pos == 0) { */
  /*     other->right = root; */
  /*     other->right->parent = other; */
  /*     root = other; */
  /*   } */
  /* } */

  // Extract a substring from the tree and return the root of the extracted
  // subtree
  void extract(const uint32_t i, const uint32_t j, SplayTree &other) {
    if (i < 0 || j >= root->subtree_size || i > j) {
      std::cout << "Invalid range for extract" << std::endl;
      return;
    }
    if (i == 0 && j == root->subtree_size - 1) {
      other.root = root;
      root = nullptr;
      return;
    }
    if (i == 0) {
      get<normal>(j + 1);
      other.root = root->left;
      other.root->parent = nullptr;
      root->left = nullptr;
      updateNodeInfo(root);
      return;
    }

    other.root = isolate(i, j); // newRoot;
    other.root->parent = nullptr;

    if (root->right && j != root->subtree_size - 1) {
      root->right->left = nullptr;
      updateNodeInfo(root->right);
    } else {
      root->right = nullptr;
    }

    updateNodeInfo(root);
    return;
  }

  // Retrieve a substring from the tree
  std::string retrieve(const uint32_t i, const uint32_t j) {
    if (!root) {
      return "";
    }
    node *subtree = isolate(i, j);
    std::string result;
    inorder(subtree, result);
    return result;
  }

  // Substring equality check
  bool equal(const uint32_t i, SplayTree &other, const uint32_t j,
             const uint32_t length) {
    if (i + length - 1 >= root->subtree_size ||
        j + length - 1 >= getSubtreeSize(other.root)) {
      return false;
    }

    uint64_t kr1 = 0, kr2 = 0;
    node *leftSubstring = isolate(i, i + length - 1);
    if (leftSubstring) {
      kr1 = leftSubstring->kr;
    } else {
      return false;
    }

    node *rightSubstring = other.isolate(j, j + length - 1);
    if (rightSubstring) {
      kr2 = rightSubstring->kr;
    } else {
      return false;
    }

    return kr1 == kr2;
  }

  // Compute the LCP of two substrings
  uint32_t LCP(uint32_t i, SplayTree &other, uint32_t j) {
    if (i >= root->subtree_size || j >= getSubtreeSize(other.root)) {
      std::cout << "Invalid range for LCP computation" << std::endl;
      return 0;
    }
    if (&other == this && i > j) {
      uint32_t tmp = i;
      i = j;
      j = tmp;
    }

    int n_prime =
        std::min(getSubtreeSize(root) - i, getSubtreeSize(other.root) - j);

    // Check boundary cases
    if (get<normal>(i)->character != other.get<normal>(j)->character) {
      return 0;
    }
    if (!equal(i, other, j, 2)) {
      return 1;
    }
    if (equal(i, other, j, n_prime)) {
      return n_prime;
    }

    // Check at predetermined threshold (2^{log_2^{2/3}(n)}) if substrings are
    // equal
    int threshold = std::min(
        n_prime,
        (int)pow(2, pow(log2(getSubtreeSize(root) + getSubtreeSize(other.root)),
                        2.0 / 3)));

    bool overlap = false;
    // if yes, just do exponential search, extract the substring and
    // do doubling search on the extracted substring until LCP is computed
    if (equal(i, other, j, threshold)) {
      SplayTree firstSubstr = nullptr, secondSubstr = nullptr;
      return _LCP_routine(i, other, j, n_prime, &other == this, firstSubstr,
                          secondSubstr);
    } else { // first extract then do LCP computation
      SplayTree firstSubstrThresh = nullptr, secondSubstrThresh = nullptr;
      SplayTree firstSubstr = nullptr, secondSubstr = nullptr;

      int i_extracted = 0, j_extracted = 0;
      if (i + threshold - 1 >= j && &other == this) {
        extract(i, j + threshold - 1, firstSubstrThresh);
        j_extracted = j - i;
        overlap = true;
      } else {
        extract(i, i + threshold - 1, firstSubstrThresh);
        other.extract(j - ((&other == this) * threshold),
                      j - ((&other == this) * threshold) + threshold - 1,
                      secondSubstrThresh);
      }

      uint32_t LCP_value = firstSubstrThresh._LCP_routine(
          i_extracted, (overlap) ? firstSubstrThresh : secondSubstrThresh,
          j_extracted, getSubtreeSize(firstSubstrThresh.root) - j_extracted,
          (&other == this) & overlap, firstSubstr, secondSubstr);

      if (i + threshold - 1 >= j && &other == this) {
        introduce(i, firstSubstrThresh);
      } else {
        introduce(i, firstSubstrThresh);
        other.introduce(j, secondSubstrThresh);
      }
      return LCP_value;
    }
  }
  // Visualize the tree with at most 10 levels of depth
  void visualize(node *n, int depth = 0) {
    if (depth > 5)
      return;
    if (n) {
      visualize(n->right, depth + 1);
      for (int i = 0; i < depth; i++)
        std::cout << "   ";
      std::cout << n << std::endl;
      visualize(n->left, depth + 1);
    }
  }

  uint32_t size() { return root->subtree_size; }

private:
  // Construct a fully balanced binary tree from an array of characters
  node *buildBalancedTree(const std::string &chars, const int32_t start,
                          const int32_t end) {
    if (start > end)
      return nullptr;

    int32_t mid = start + (end - start) / 2;
    node *n = new node(chars[mid]);
    n->left = buildBalancedTree(chars, start, mid - 1);
    n->right = buildBalancedTree(chars, mid + 1, end);

    if (n->left)
      n->left->parent = n;
    if (n->right)
      n->right->parent = n;

    updateNodeInfo(n);
    return n;
  }

  // Isolate a substring of the tree and return the root of the isolated tree
  node *isolate(const uint32_t i, const uint32_t j) {
    if (i < 0 || j >= root->subtree_size) {
      return nullptr;
    }
    if (i == 0 && j == root->subtree_size - 1) {
      return root;
    }
    if (i == 0) {
      get<normal>(j + 1);
      return root->left;
    }
    if (j == root->subtree_size - 1) {
      get<normal>(i - 1);
      return root->right;
    }

    get<normal>(j + 1);
    get<modified>(i - 1);

    return root->right->left;
  }

  // LCP subroutine
  uint32_t _LCP_routine(const uint32_t i, SplayTree &other, const uint32_t j,
                        const uint32_t n_prime, const bool sameString,
                        SplayTree &firstSubstr, SplayTree &secondSubstr) {

    int l_prime = exponentialSearch(i, other, j, n_prime);
    int i_extracted = 0, j_extracted = 0;
    if (i + l_prime - 1 >= j && sameString) { // same string and overlapping
      this->extract(i, j + l_prime - 1, firstSubstr);
      j_extracted = j - i;
    } else {
      this->extract(i, i + l_prime - 1, firstSubstr);
      other.extract(j, j + l_prime - 1, secondSubstr);
    }
    int range = doublingSearch(i_extracted, firstSubstr, j_extracted,
                               j_extracted ? firstSubstr : secondSubstr);
    uint32_t LCP = binarySearch(firstSubstr, i_extracted + range / 2,
                                j_extracted ? firstSubstr : secondSubstr,
                                j_extracted + range / 2, range / 2) +
                   (range / 2);

    if (i + l_prime - 1 >= j && sameString) {
      introduce(i, firstSubstr);
    } else {
      this->introduce(i, firstSubstr);
      other.introduce(j, secondSubstr);
    }
    return LCP;
  }

  // Exponential search for LCP computation
  // Returns the first length at which the exponentiation equality fails
  uint32_t exponentialSearch(const uint32_t i, SplayTree &other,
                             const uint32_t j, const uint32_t n_prime) {
    uint32_t p = 2;

    while (p < n_prime && equal(i, other, j, p)) {
      if ((uint64_t)pow((uint64_t)p, 2) >= (uint64_t)n_prime) {
        p = n_prime;
        break;
      }
      p = pow(p, 2);
    }

    return std::min(p, n_prime);
  }

  // Doubling search for LCP computation
  // Returns the first length at which the doubling equality fails
  uint32_t doublingSearch(const uint32_t i, SplayTree &firstSubstring,
                          const uint32_t j, SplayTree &secondSubstring) {
    uint32_t length = 1;
    while (firstSubstring.equal(i, secondSubstring, j, length) &&
           length < getSubtreeSize(firstSubstring.root)) {
      length *= 2;
    }
    return length;
  }

  // Binary search between two splay trees in a range
  uint32_t binarySearch(SplayTree &firstSubstring, const uint32_t i,
                        SplayTree &secondSubstring, const uint32_t j,
                        const uint32_t size_range) {
    int low = 0, high = size_range;
    int result = 0;
    while (low <= high) {
      int mid = (high + low) / 2;
      if (firstSubstring.equal(i, secondSubstring, j, mid)) {
        result = mid;
        low = mid + 1;
      } else {
        high = mid - 1;
      }
    }
    return result;
  }

  // Update the subtree size of a node
  void updateSubtreeSize(node *n) {
    if (n) {
      n->subtree_size = 1;
      if (n->left)
        n->subtree_size += n->left->subtree_size;
      if (n->right)
        n->subtree_size += n->right->subtree_size;
    }
  }

  // Update the KR has of a node */
  void updateKR(node *n) {
    if (n) {
      n->kr = n->character;
      n->base_exp = base;
      if (n->left) {
        n->kr = ((__uint128_t)n->kr + (__uint128_t)n->left->kr * base) % prime;
        n->base_exp =
            ((__uint128_t)n->left->base_exp * (__uint128_t)n->base_exp) % prime;
      }
      if (n->right) {
        n->kr = ((__uint128_t)n->kr * (__uint128_t)n->right->base_exp +
                 (__uint128_t)n->right->kr) %
                prime;
        n->base_exp =
            ((__uint128_t)n->base_exp * (__uint128_t)n->right->base_exp) %
            prime;
      }
    }
    assert(n->base_exp > 0);
  }

  void updateNodeInfo(node *n) {
    updateSubtreeSize(n);
    updateKR(n);
  }

  // Splay operation
  template <auto isModified> void splay(node *x) {
    while (x->parent) {
      if (!x->parent->parent) {
        if (x->parent->left == x) {
          rightRotate(x->parent);
        } else {
          leftRotate(x->parent);
        }
      } else if (x->parent->left == x && x->parent->parent->left == x->parent) {
        if constexpr (isModified == modified) {
          if (x->parent->parent == root) {
            rightRotate(x->parent);
            rightRotate(x->parent);
          } else {
            rightRotate(x->parent->parent);
            rightRotate(x->parent);
          }
        } else {
          rightRotate(x->parent->parent);
          rightRotate(x->parent);
        }
      } else if (x->parent->right == x &&
                 x->parent->parent->right == x->parent) {
        leftRotate(x->parent->parent);
        leftRotate(x->parent);
      } else if (x->parent->left == x &&
                 x->parent->parent->right == x->parent) {
        rightRotate(x->parent);
        leftRotate(x->parent);
      } else {
        leftRotate(x->parent);
        rightRotate(x->parent);
      }
    }
    root = x;

    updateNodeInfo(root);
  }

  // Right rotate
  void rightRotate(node *x) {
    node *y = x->left;
    x->left = y->right;
    if (y->right)
      y->right->parent = x;
    y->parent = x->parent;
    if (!x->parent)
      root = y;
    else if (x == x->parent->left)
      x->parent->left = y;
    else
      x->parent->right = y;
    y->right = x;
    x->parent = y;

    updateNodeInfo(x);
    updateNodeInfo(y);
  }

  // Left rotate
  void leftRotate(node *x) {
    node *y = x->right;
    x->right = y->left;
    if (y->left)
      y->left->parent = x;
    y->parent = x->parent;
    if (!x->parent)
      root = y;
    else if (x == x->parent->left)
      x->parent->left = y;
    else
      x->parent->right = y;
    y->left = x;
    x->parent = y;

    updateNodeInfo(x);
    updateNodeInfo(y);
  }

  // Join current tree with another tree
  node *join(SplayTree &rightTree) {
    if (!root) {
      root = rightTree.root;
      rightTree.root = nullptr;
      return root;
    } else if (rightTree.root) {
      get<normal>(root->subtree_size - 1);
      root->right = rightTree.root;
      root->right->parent = root;
      rightTree.root = nullptr;
      updateNodeInfo(root);
      return root;
    }
    return nullptr;
  }

  // Split the tree into two trees based on a position
  void split(const uint32_t position, SplayTree &leftTree,
             SplayTree &rightTree) {
    node *z = get<normal>(position);
    if (!z) {
      leftTree.root = root;
      rightTree.root = nullptr;
    } else {
      leftTree.root = z;
      rightTree.root = z->right;
      z->right = nullptr;
      if (rightTree.root)
        rightTree.root->parent = nullptr;
      updateNodeInfo(rightTree.root);
      updateNodeInfo(leftTree.root);
    }
  }

  // Inorder traversal of the tree
  void inorder(node *n, std::string &result) {
    if (n) {
      inorder(n->left, result);
      result.push_back(n->character);
      inorder(n->right, result);
    }
  }
};

class StringHandler {

public:
  uint64_t handler;

  StringHandler() : handler(0) {}

  StringHandler(const SplayTree *tree) { handler = (uint64_t)tree; }

  bool operator==(const StringHandler &other) const {
    return handler == other.handler;
  }
};
// define hash function for StringHandler
template <> struct std::hash<StringHandler> {
  std::size_t operator()(const StringHandler &str) const {
    return std::hash<uint64_t>()(str.handler);
  }
};

class FeST {
  std::unordered_map<StringHandler, SplayTree> trees;

public:
  FeST() { trees.reserve(1000); }

  FeST(const std::string &str) {
    SplayTree tree(str);
    trees[&tree] = tree;
  }

  FeST(FeST &other) { trees = other.trees; }

  void delete_string(StringHandler tree) { trees.erase(tree); }

  // Merge two strings, first+second, invalidate the key for the second string
  void merge_strings(StringHandler first, StringHandler second) {
    trees[first].introduce(trees[first].size(), trees[second]);
    trees.erase(second);
  }

  StringHandler split_string(StringHandler tree, uint32_t position) {
    SplayTree *rightTree = new SplayTree();
    trees[tree].extract(position, trees[tree].size() - 1, *rightTree);
    trees[rightTree] = *rightTree;
    return rightTree;
  }

  StringHandler add_new_string(SplayTree &tree) {
    SplayTree *newTree = new SplayTree(tree);
    trees[newTree] = *newTree;
    return newTree;
  }

  StringHandler add_new_string(const std::string &str) {
    SplayTree *newTree = new SplayTree(str);
    trees[newTree] = *newTree;
    return newTree;
  }

  SplayTree *get_string(StringHandler sh) { return &trees[sh]; }

  SplayTree *get_string(const std::string &str) {
    for (auto &tree : trees) {
      if (tree.second.retrieve(0, tree.second.size() - 1) == str) {
        return &tree.second;
      }
    }
    return nullptr;
  }

  uint64_t size() { return trees.size(); }
};

// Do some tests on LCP computation on two randomly generated strings and
// check with a scan if the results are correct
int main() {

  const std::string firstString = "banana";
  const std::string secondString = "bananana";
  const std::string thirdString = "abracadabra";
  const std::string fourthString = "bandana";
  const std::string fifthString = "mississippi";

  FeST forest;

  std::vector<StringHandler> keys;
  keys.push_back(forest.add_new_string(firstString));
  keys.push_back(forest.add_new_string(secondString));
  keys.push_back(forest.add_new_string(thirdString));
  keys.push_back(forest.add_new_string(fourthString));
  keys.push_back(forest.add_new_string(fifthString));

  // print StringHandler and the corresponding string
  for (auto &key : keys) {
    std::cout << key.handler << " : "
              << forest.get_string(key)->retrieve(
                     0, forest.get_string(key)->size() - 1)
              << std::endl;
  }

  // edit third string by inserting a character
  forest.get_string(keys[2])->insert('c', 5);
  std::cout << "Edited string: "
            << forest.get_string(keys[2])->retrieve(
                   0, forest.get_string(keys[2])->size() - 1)
            << std::endl;

  // edit third string by changing a character
  forest.get_string(keys[2])->edit('r', 5);
  std::cout << "Edited string: "
            << forest.get_string(keys[2])->retrieve(
                   0, forest.get_string(keys[2])->size() - 1)
            << std::endl;

  // delete inserted character
  forest.get_string(keys[2])->delete_node(5);
  std::cout << "Reverted string: "
            << forest.get_string(keys[2])->retrieve(
                   0, forest.get_string(keys[2])->size() - 1)
            << std::endl;

  // merge first and second string
  // CAREFUL: the key for the second string is invalidated
  forest.merge_strings(keys[0], keys[1]);
  std::cout << "Merged string: "
            << forest.get_string(keys[0])->retrieve(
                   0, forest.get_string(keys[0])->size() - 1)
            << std::endl;

  StringHandler newString = forest.split_string(keys[0], 6);
  std::cout << "Splitted string "
            << forest.get_string(newString)->retrieve(
                   0, forest.get_string(newString)->size() - 1)
            << " from "
            << forest.get_string(keys[0])->retrieve(
                   0, forest.get_string(keys[0])->size() - 1)
            << std::endl;

  /* exit(0); */

  double average_time_LCP = 0;
  double average_time_scan = 0;
  // set the seed for random number generation
  srand(time(0));
  for (int i = 0; i < 10; i++) {
    std::cout << "Round " << i << std::endl;
    uint32_t string_length = rand() % 100000000;
    std::string chars;
    chars.reserve(string_length);
    for (uint32_t j = 0; j < string_length; j++) {
      chars.push_back('a' + rand() % 26);
    }

    std::string secondChars;
    secondChars = chars;
    // pop a random character from the second string
    uint32_t position = rand() % string_length;
    std::cout << "Going to remove character at position: " << position
              << std::endl;
    secondChars.erase(secondChars.begin() + position);

    std::cout << "Building first tree" << std::endl;
    auto buildTimeStart = std::chrono::high_resolution_clock::now();
    SplayTree tree(chars);
    auto buildTimeEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for tree construction: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     buildTimeEnd - buildTimeStart)
                     .count()
              << " microseconds" << std::endl;
    std::cout << "Building second tree" << std::endl;
    buildTimeStart = std::chrono::high_resolution_clock::now();
    SplayTree secondTree(secondChars);
    buildTimeEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for tree construction: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     buildTimeEnd - buildTimeStart)
                     .count()
              << " microseconds" << std::endl;

    int start = rand() % 100;
    int secondStart = start; // rand() % 100;
    std::cout << "Start: " << start << " Second start: " << secondStart
              << std::endl;

    std::cout << "The LCP should be either " << position - start << " or "
              << string_length - start << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    uint32_t LCP_value = tree.LCP(start, secondTree, secondStart);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for LCP computation: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     end_time - start_time)
                     .count()
              << " microseconds" << std::endl;
    average_time_LCP += std::chrono::duration_cast<std::chrono::microseconds>(
                            end_time - start_time)
                            .count();

    // chek if the LCP value is correct
    start_time = std::chrono::high_resolution_clock::now();
    uint32_t scan = 0;
    while (start + scan < string_length && secondStart + scan < string_length &&
           chars[start + scan] == secondChars[secondStart + scan]) {
      scan++;
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for scan: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     end_time - start_time)
                     .count()
              << " microseconds" << std::endl;
    average_time_scan += std::chrono::duration_cast<std::chrono::microseconds>(
                             end_time - start_time)
                             .count();

    if (scan != LCP_value) {
      std::cout << "Error: LCP value is not correct" << std::endl;
      std::cout << "Computed LCP: " << LCP_value << " Scan LCP: " << scan
                << std::endl;

      return 1;
    } else {
      std::cout << "LCP computation is correct" << std::endl;
      std::cout << "Computed LCP: " << LCP_value << " Scan LCP: " << scan
                << std::endl;
    }
  }

  std::cout << "Average time for LCP computation: " << average_time_LCP / 100
            << " microseconds" << std::endl;
  std::cout << "Average time for scan computation: " << average_time_scan / 100
            << " microseconds" << std::endl;

  {
    // test LCP computation on the same string
    std::string chars = {'m', 'i', 's', 's', 'i', 's', 'a'};
    SplayTree tree(chars);
    uint32_t LCP_value = tree.LCP(1, tree, 4);
    std::cout << "LCP: " << LCP_value << std::endl;
    assert(LCP_value == 2);

    // test LCP computation on the same string and overlapping LCP
    std::string secondChars = {'m', 'i', 's', 's', 'i', 's', 's', 'i', 'a'};
    SplayTree secondTree(secondChars);
    LCP_value = secondTree.LCP(1, secondTree, 4);
    std::cout << "LCP: " << LCP_value << std::endl;
    assert(LCP_value == 4);

    LCP_value = secondTree.LCP(0, secondTree, 0);
    std::cout << "LCP: " << LCP_value << std::endl;
    assert(LCP_value == 9);

    LCP_value = secondTree.LCP(5, secondTree, 1);
    std::cout << "LCP: " << LCP_value << std::endl;
    assert(LCP_value == 0);
  }

  // random test for LCP computation on same random string
  average_time_LCP = 0;
  average_time_scan = 0;

  srand(time(0));
  for (int i = 0; i < 10; i++) {
    std::cout << "Round (same string) " << i << std::endl;
    uint32_t string_length = rand() % 100000000;
    std::string chars;
    chars.reserve(string_length);
    for (uint32_t j = 0; j < string_length; j++) {
      chars.push_back('a' + rand() % 26);
    }

    // select a region of the string and make it pseudo-periodic with period
    // 30
    int beginning = rand() % 100;
    int end = beginning + (rand() % 1000); // rand() % 100;
    for (int j = beginning; j < end; j++) {
      chars[j] = 'a';
    }

    SplayTree tree(chars);
    int start = beginning;                    // rand() % 100;
    int secondStart = start + (rand() % 100); // rand() % 100;
    std::cout << "Start: " << start << " Second start: " << secondStart
              << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    uint32_t LCP_value = tree.LCP(start, tree, secondStart);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for LCP computation: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     end_time - start_time)
                     .count()
              << " microseconds" << std::endl;
    average_time_LCP += std::chrono::duration_cast<std::chrono::microseconds>(
                            end_time - start_time)
                            .count();
    // chek if the LCP value is correct
    start_time = std::chrono::high_resolution_clock::now();
    uint32_t scan = 0;
    while (start + scan < string_length && secondStart + scan < string_length &&
           chars[start + scan] == chars[secondStart + scan]) {
      scan++;
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for scan: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     end_time - start_time)
                     .count()
              << " microseconds" << std::endl;
    average_time_scan += std::chrono::duration_cast<std::chrono::microseconds>(
                             end_time - start_time)
                             .count();
    if (scan != LCP_value) {
      std::cout << "Error: LCP value is not correct" << std::endl;
      std::cout << "Computed LCP: " << LCP_value << " Scan LCP: " << scan
                << std::endl;
      return 1;
    } else {
      std::cout << "LCP computation is correct" << std::endl;
      std::cout << "Computed LCP: " << LCP_value << " Scan LCP: " << scan
                << std::endl;
    }
  }

  return 0;
}
