#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>

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
  bool reversed;
  uint64_t kr_rev;

  node(char c)
      : parent(nullptr), left(nullptr), right(nullptr), subtree_size(1),
        character(c), kr(c), base_exp(base), reversed(false), kr_rev(c) {}

  node(node *n) {
    parent = n->parent;
    left = n->left;
    right = n->right;
    subtree_size = n->subtree_size;
    character = n->character;
    kr = n->kr;
    base_exp = n->base_exp;
    reversed = false;
    kr_rev = n->kr;
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
     << " p " << n->base_exp << " rev " << n->reversed << " kr_rev "
     << n->kr_rev;

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
      SplayTree leftTree, rightTree;
      split(position - 1, leftTree, rightTree);
      node *n = new node(c);
      root = n;
      root->left = leftTree.root;
      root->right = rightTree.root;
      if (leftTree.root)
        leftTree.root->parent = root;
      if (rightTree.root)
        rightTree.root->parent = root;
    }
    updateNodeInfo(root);
  }

  // Find a node at a given position
  template <auto isModified> node *find(const int position) {
    node *z = root;
    int currentPos = position;

    while (z) {
      updateReversed(z);
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
      find<normal>(0);
      root = root->right;
      root->parent = nullptr;
      updateNodeInfo(root);
      return;
    }

    SplayTree leftTree, rightTree, deletedNode;
    split(position - 1, leftTree, rightTree);

    rightTree.split(0, deletedNode, rightTree);

    leftTree.join(rightTree);
    root = leftTree.root;
  }

  // Introduce a substring from another tree at a given position
  void introduce(const uint32_t pos, SplayTree &other) {
    if (!root) {
      if (other.root)
        root = other.root;
      other.root = nullptr;
    } else if (pos == root->subtree_size) {
      find<normal>(pos - 1);
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
      find<normal>(other.root->subtree_size - 1);
      other.root->right = root;
      other.root->right->parent = other.root;
      root = other.root;
      other.root = nullptr;
    }
  }

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
      find<normal>(j + 1);
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
    if (find<normal>(i)->character != other.find<normal>(j)->character) {
      return 0;
    }
    if (!equal(i, other, j, 2)) {
      return 1;
    }
    if (equal(i, other, j, n_prime)) {
      return n_prime;
    }

    // Check at predetermined threshold (2^log_2(n)) if substrings are equal
    int threshold = std::min(
        n_prime,
        (int)pow(2, pow(log2(getSubtreeSize(root) + getSubtreeSize(other.root)),
                        2.0 / 3)));

    bool overlap = false;
    // if yes, just do exponential search, extract the substring and
    // doubling search the extracted substring until LCP is computed
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

  // Reverse substring from i to j
  void reverse(const uint32_t i, const uint32_t j) {
    if (i < 0 || j >= root->subtree_size || i > j) {
      std::cout << "Invalid range for reverse" << std::endl;
      return;
    }
    node *subtree = isolate(i, j);
    subtree->reversed = !subtree->reversed;
    updateNodeInfo(subtree);
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
      find<normal>(j + 1);
      return root->left;
    }
    if (j == root->subtree_size - 1) {
      find<normal>(i - 1);
      return root->right;
    }

    find<normal>(j + 1);
    find<modified>(i - 1);

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

  // Update the reverse KR hash of a node
  void updateKR_rev(node *n) {
    if (n) {
      n->kr_rev = n->character;
      if (n->right) {
        n->kr_rev =
            ((__uint128_t)n->kr_rev * n->right->base_exp + n->right->kr_rev) %
            prime;
      }
      if (n->left) {
        n->kr_rev =
            ((__uint128_t)n->left->kr_rev * n->base_exp + n->kr_rev) % prime;
      }
    }
  }

  // Update reversed
  void updateReversed(node *n) {
    if (n) {
      if (n->reversed) {
        n->reversed = false;
        if (n->left)
          n->left->reversed = !n->left->reversed;
        if (n->right)
          n->right->reversed = !n->right->reversed;
        node *tmp = n->left;
        n->left = n->right;
        n->right = tmp;
        tmp = nullptr;

        uint64_t tmp_kr = n->kr;
        n->kr = n->kr_rev;
        n->kr_rev = tmp_kr;
      }
    }
  }

  void updateNodeInfo(node *n) {
    updateSubtreeSize(n);
    updateKR(n);
    updateKR_rev(n);
    updateReversed(n);
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
      find<normal>(root->subtree_size - 1);
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
    node *z = find<normal>(position);
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

// Do some tests on LCP computation on two randomly generate strings and check
// with a scan if the results are correct
int main() {

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
    SplayTree tree(chars);
    std::cout << "Building second tree" << std::endl;
    SplayTree secondTree(secondChars);

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
              << " ms" << std::endl;
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
              << " ms" << std::endl;
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
            << " ms" << std::endl;
  std::cout << "Average time for scan computation: " << average_time_scan / 100
            << " ms" << std::endl;

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
    uint32_t string_length = rand() % 10000000;
    std::string chars;
    chars.reserve(string_length);
    for (uint32_t j = 0; j < string_length; j++) {
      chars.push_back('a' + rand() % 26);
    }

    // select a region of the string and make it pseudo-periodic with period 30
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
              << " ms" << std::endl;
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
              << " ms" << std::endl;
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

/* std::vector<char> chars = {'m', 'i', 's', 's', 'i', 's', */
/*                            's', 'i', 'p', 'p', 'i'}; */
/* SplayTree tree; */
/* tree.root = tree.buildBalancedTree(chars, 0, chars.size() - 1); */

/* std::cout << "Root character: " << tree.root->character << std::endl; */
/* tree.visualize(tree.getRoot()); */

/* std::cout << "Find node at position 5" << std::endl; */
/* tree.find<normal>(8); */
/* tree.visualize(tree.getRoot()); */

/* tree.insert('h', 3); // Insert 'h' at position 3 */
/* std::cout << "Root after insertion: " << tree.root->character << std::endl;
 */
/* tree.visualize(tree.getRoot()); */

/* tree.insert('i', 0); // Insert 'i' at position 0 */
/* std::cout << "Root after insertion: " << tree.root->character << std::endl;
 */
/* tree.visualize(tree.getRoot()); */

/* tree.deleteNode(3); // Delete node at position 3 */
/* std::cout << "Root after deletion: " << tree.root->character << std::endl;
 */
/* tree.visualize(tree.getRoot()); */

/* tree.find<normal>(3); // Find node at position 3 */
/* std::cout << "Found node: " << tree.root->character << std::endl; */
/* tree.visualize(tree.getRoot()); */

/* tree.deleteNode(0); // Delete node at position 1 */
/* std::cout << "Root after deletion: " << tree.root->character << std::endl;
 */
/* tree.visualize(tree.getRoot()); */

/* node *subtring = tree.isolate(3, 6); // Isolate substring from position 3
 * to
 * 6 */
/* std::cout << "Isolated substring: " << subtring->character << std::endl; */
/* tree.visualize(subtring); */

/* std::string firstSubstring = */
/*     tree.retrieve(3, 4); // Retrieve substring from position 3 to 4 */
/* std::string secondSubstring = */
/*     tree.retrieve(6, 7); // Retrieve substring from position 6 to 7 */
/* bool equal = tree.equal(3, tree, 6, 2); // Check if two substrings are
 * equal
 */
/* std::cout << "Are two substrings (" << firstSubstring << ", " */
/*           << secondSubstring << ") equal: " << (equal ? "TRUE" : "FALSE")
 */
/*           << std::endl; */

/* firstSubstring = tree.retrieve(0, 6); */
/* secondSubstring = tree.retrieve(3, 9); */
/* equal = tree.equal(3, tree, 0, 7); */
/* std::cout << "Are two substrings (" << firstSubstring << ", " */
/*           << secondSubstring << ") equal: " << (equal ? "TRUE" : "FALSE")
 */
/*           << std::endl; */

/* node *strangeSubstring = tree.isolate(6, 5); */
/* if (strangeSubstring == nullptr) { */
/*   std::cout << "There is no root->right->left node" << std::endl; */
/* } else { */
/*   tree.visualize(strangeSubstring); */
/* } */

/* std::vector<char> secondChars = {'h', 'e', 'l', 'l', 'o'}; */
/* SplayTree secondTree; */
/* secondTree.root = */
/*     secondTree.buildBalancedTree(secondChars, 0, secondChars.size() - 1);
 */

/* tree.introduce(3, secondTree); */
/* std::cout << "Root after introducing substring: " << tree.root->character
 */
/*           << std::endl; */
/* tree.visualize(tree.getRoot()); */

/* SplayTree extractedTree = tree.extract(2, 6); */
/* std::cout << "Extracted tree: " << extractedTree.root->character <<
 * std::endl; */
/* extractedTree.visualize(extractedTree.getRoot()); */
/* std::cout << "Root after extraction: " << tree.root->character <<
 * std::endl;
 */
/* tree.visualize(tree.getRoot()); */

/* extractedTree.insert('h', extractedTree.getRoot()->subtree_size); */
/* std::cout << "Root after insertion: " << extractedTree.root->character */
/*           << std::endl; */
/* extractedTree.visualize(extractedTree.getRoot()); */
/* std::vector<char> thirdChars = {'h', 'e', 'l', 'l', 'o'}; */
/* SplayTree thirdTree; */
/* thirdTree.root = */
/*     thirdTree.buildBalancedTree(thirdChars, 0, thirdChars.size() - 1); */

/* uint32_t LCP_value = extractedTree.LCP(2, thirdTree, 1); */
/* std::cout << "Root after LCP computation: " <<
 * extractedTree.root->character
 */
/*           << std::endl; */
/* extractedTree.visualize(extractedTree.getRoot()); */
/* std::cout << "Root after LCP computation: " << thirdTree.root->character */
/*           << std::endl; */
/* thirdTree.visualize(thirdTree.getRoot()); */
/* std::cout << "LCP: " << LCP_value << std::endl; */

/* // Insert a new character at a given position */
/*   void insert(char c, int position) { */

/*     if (position == 0) { */
/*       node *n = new node(c); */
/*       n->right = root; */
/*       if (root) */
/*         root->parent = n; */
/*       root = n; */
/*       return; */
/*     } */

/*     node *z = root; */
/*     int currentPos = position - 1; */

/*     while (z) { */
/*       int leftSize = z->left ? z->left->subtree_size : 0; */
/*       if (currentPos == leftSize + 1) { */
/*         break; */
/*       } else if (currentPos < leftSize + 1) { */
/*         z = z->left; */
/*       } else { */
/*         currentPos -= leftSize + 1; */
/*         z = z->right; */
/*       } */
/*     } */

/*     node *n = new node(c); */
/*     if (!z) { */
/*       root = n; */
/*     } else if (currentPos == 0) { */
/*       n->left = z->left; */
/*       if (z->left) */
/*         z->left->parent = n; */
/*       z->left = n; */
/*       n->parent = z; */
/*     } else { */
/*       n->left = z->right; */
/*       if (z->right) */
/*         z->right->parent = n; */
/*       z->right = n; */
/*       n->parent = z; */
/*     } */

/*     splay(n); */
