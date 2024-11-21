#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <vector>

#define normal false
#define modified true

static uint64_t prime = ((uint64_t)1 << 63) - 59;
static uint64_t base = 255;

struct node {
  node *parent;
  node *left;
  node *right;
  int subtree_size;
  char character;
  uint64_t kr;
  uint64_t base_exp;

  node(char c)
      : parent(nullptr), left(nullptr), right(nullptr), subtree_size(1),
        character(c), kr(c), base_exp(base) {}
};

class SplayTree {
public:
  node *root;

  // Getter for root
  node *getRoot() { return root; }

  SplayTree() : root(nullptr) {}

  // Construct a fully balanced binary tree from an array of characters
  node *buildBalancedTree(const std::vector<char> &chars, int start, int end) {
    if (start > end)
      return nullptr;

    int mid = start + (end - start) / 2;
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

  /* void updateKR(node *n) { */
  /*   if (n) { */
  /*     n->kr = n->character; */
  /*     n->base_exp = base; */
  /*     if (n->left) { */
  /*       n->kr = ((__uint128_t)n->kr + (__uint128_t)n->left->kr * base) %
   * prime; */
  /*       n->base_exp = */
  /*           ((__uint128_t)n->left->base_exp * (__uint128_t)n->base_exp) %
   * prime; */
  /*     } */
  /*     if (n->right) { */
  /*       n->kr = ((__uint128_t)n->kr * (__uint128_t)n->right->base_exp +
   * (__uint128_t)n->right->kr) % */
  /*               prime; */
  /*       n->base_exp = */
  /*           ((__uint128_t)n->base_exp * (__uint128_t)n->right->base_exp) % */
  /*           prime; */
  /*     } */
  /*   } */
  /*   assert(n->base_exp > 0); */
  /* } */

  // Insert a new character at a given position using only split operations
  void insert(const char c, const int position) {
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
  void deleteNode(const int position) {
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

  // Isolate a substring of the tree and return the root of the isolated tree
  node *isolate(int i, int j) {
    if (i < 0 || j >= root->subtree_size) {
      return nullptr;
    }
    if (i == 0 && j == root->subtree_size - 1) {
      return root;
    }
    if (i == 0) {
      find<normal>(j + 1);
      assert(root->left != nullptr);
      return root->left;
    }
    if (j == root->subtree_size - 1) {
      find<normal>(i - 1);
      assert(root->right != nullptr);
      return root->right;
    }

    find<normal>(j + 1);
    find<modified>(i - 1);

    return root->right->left;
  }

  // Introduce a substring from another tree at a given position
  void introduce(int pos, SplayTree &other) {
    if (!root) {
      if (other.getRoot())
        root = other.getRoot();
    } else if (pos == root->subtree_size) {
      root->right = other.getRoot();
      if (other.getRoot())
        other.getRoot()->parent = root;
      updateNodeInfo(root);
    } else if (pos <= root->subtree_size && pos > 0) {
      isolate(pos, pos - 1);

      if (root->right) {
        node *toConnect = root->right;
        toConnect->left = other.getRoot();
        if (other.getRoot())
          other.getRoot()->parent = toConnect;

        updateNodeInfo(toConnect);
        updateNodeInfo(root);
      }
    } else if (pos == 0) {
      other.join(*this);
    }
  }

  // Extract a substring from the tree and return the root of the extracted
  // subtree
  SplayTree extract(int i, int j) {
    if (i < 0 || j >= root->subtree_size || i > j) {
      std::cout << "Invalid range" << std::endl;
      return SplayTree();
    }
    if (i == 0 && j == root->subtree_size - 1) {
      SplayTree extractedTree;
      extractedTree.root = root;
      root = nullptr;
      return extractedTree;
    }

    node *newRoot = isolate(i, j);
    newRoot->parent = nullptr;
    SplayTree extractedTree;
    extractedTree.root = newRoot;

    if (root->right && j != root->subtree_size - 1) {
      root->right->left = nullptr;
      updateNodeInfo(root->right);
    } else {
      root->right = nullptr;
    }

    updateNodeInfo(root);

    return extractedTree;
  }

  // Retrieve a substring from the tree
  std::string retrieve(int i, int j) {
    node *subtree = isolate(i, j);
    std::string result;
    inorder(subtree, result);
    return result;
  }

  // Substring equality check
  bool equal(int i, SplayTree &other, int j, int length) {
    if (i + length - 1 >= root->subtree_size ||
        j + length - 1 >= other.getRoot()->subtree_size) {
      std::cout << "Invalid range" << std::endl;
      return false;
    }

    node *leftSubstring = isolate(i, i + length - 1);
    node *rightSubstring = other.isolate(j, j + length - 1);

    if (!leftSubstring || !rightSubstring)
      return false;

    std::cout << "Left substring: " << leftSubstring->character << std::endl;
    std::cout << "Right substring: " << rightSubstring->character << std::endl;
    std::cout << "KR left: " << leftSubstring->kr << std::endl;
    std::cout << "KR right: " << rightSubstring->kr << std::endl;
    std::cout << "Left tree: " << std::endl;
    visualize(root);
    std::cout << "Right tree: " << std::endl;
    visualize(other.getRoot());
    return leftSubstring->kr == rightSubstring->kr;
  }

  // Compute the LCP of two substrings
  uint32_t LCP(int i, SplayTree &other, int j) {
    if (i >= root->subtree_size || j >= other.getRoot()->subtree_size) {
      std::cout << "Invalid range" << std::endl;
      return 0;
    }

    int n_prime =
        std::min(root->subtree_size - i, other.getRoot()->subtree_size - j);
    std::cout << "N prime: " << n_prime << std::endl;
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
    int threshold =
        std::min(n_prime, (int)pow(2, pow(log2(root->subtree_size +
                                               other.getRoot()->subtree_size),
                                          2.0 / 3)));

    std::cout << "Threshold: " << threshold << std::endl;
    // if yes, just do doubling search, extract the substring and
    // exponential-search the extracted substring until LCP is computed
    if (equal(i, other, j, threshold)) {
      std::cout << "Substrings are equal at threshold" << std::endl;
      return _LCP_routine(i, other, j, n_prime);
    } else {
      std::cout << "Substrings are not equal at threshold" << std::endl;
      SplayTree firstSubstrThresh = extract(i, i + threshold - 1);
      /* std::cout << "First substring of length threshold." << std::endl; */
      /* visualize(firstSubstrThresh.getRoot()); */
      SplayTree secondSubstrThresh = other.extract(j, j + threshold - 1);
      /* std::cout << "Second substring of length threshold." << std::endl; */
      /* visualize(secondSubstrThresh.getRoot()); */
      uint32_t LCP_value =
          firstSubstrThresh._LCP_routine(0, secondSubstrThresh, 0, n_prime);
      introduce(i, firstSubstrThresh);
      other.introduce(j, secondSubstrThresh);
      return LCP_value;
    }
  }

  uint32_t _LCP_routine(int i, SplayTree &other, int j, int n_prime) {
    int l_prime = exponentialSearch(i, other, j, n_prime);
    std::cout << "Exponential search length: " << l_prime << std::endl;
    /* std::cout << "Going to extract between " << i << " and " << i + l_prime -
     * 1 */
    /*           << std::endl; */
    SplayTree firstSubstr = this->extract(i, i + l_prime - 1);
    /* std::cout << "First substring" << std::endl; */
    /* firstSubstr.visualize(firstSubstr.getRoot()); */
    /* std::cout << "Going to extract between " << j << " and " << j + l_prime -
     * 1 */
    /*           << std::endl; */
    SplayTree secondSubstr = other.extract(j, j + l_prime - 1);
    /* std::cout << "Second substring" << std::endl; */
    /* secondSubstr.visualize(secondSubstr.getRoot()); */
    int range = doublingSearch(firstSubstr, secondSubstr);
    std::cout << "Doubling search length: " << range << std::endl;
    uint32_t LCP = binarySearch(firstSubstr, range / 2, secondSubstr, range / 2,
                                range / 2) +
                   (range / 2);

    /* std::cout << "LCP: " << LCP << std::endl; */
    this->introduce(i, firstSubstr);
    /* std::cout << "First substring after introducing" << std::endl; */
    /* visualize(this->getRoot()); */
    other.introduce(j, secondSubstr);
    /* std::cout << "Second substring after introducing" << std::endl; */
    /* visualize(other.getRoot()); */
    return LCP;
  }

  // Visualize the tree with at most 10 levels of depth
  void visualize(node *n, int depth = 0) {
    if (depth > 5)
      return;
    if (n) {
      visualize(n->right, depth + 1);
      for (int i = 0; i < depth; i++)
        std::cout << "   ";
      std::cout << n->character;
      std::cout << " : " << n->subtree_size << " - " << n->kr << " p "
                << n->base_exp << std::endl;
      visualize(n->left, depth + 1);
    }
  }

private:
  // Exponential search for LCP computation
  // Returns the first length at which the exponentiation equality fails
  int exponentialSearch(int x, SplayTree &other, int i, int n_prime) {
    int p = 2;

    while (p < n_prime && equal(x, other, i, p)) {
      p = pow(p, 2);
    }
    return p;
  }

  // Doubling search for LCP computation
  // Returns the first length at which the doubling equality fails
  int doublingSearch(SplayTree &firstSubstring, SplayTree &secondSubstring) {
    int length = 1;
    while (firstSubstring.equal(0, secondSubstring, 0, length) &&
           length < firstSubstring.getRoot()->subtree_size) {
      length *= 2;
    }
    return length;
  }

  // Binary search between two splay trees in a range
  int binarySearch(SplayTree &firstSubstring, int i, SplayTree &secondSubstring,
                   int j, int size_range) {
    assert(firstSubstring.getRoot()->subtree_size ==
           secondSubstring.getRoot()->subtree_size);
    int low = 0, high = size_range;
    /* std::cout << "Starting index: " << i << " " << j << std::endl; */
    /* std::cout << "Range: " << low << " " << high << std::endl; */
    int result = 0;
    while (low <= high) {
      int mid = (high + low) / 2;
      /* std::cout << "Mid: " << mid << std::endl; */
      if (firstSubstring.equal(i, secondSubstring, j, mid)) {
        result = mid;
        low = mid + 1;
      } else {
        high = mid - 1;
      }
    }
    return result;
  }

  // Log any base
  double logb(double x, double base) { return log(x) / log(base); }

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

  // Update the KR has of a node
  void updateKR(node *n) {
    if (n) {
      n->kr = n->character;
      n->base_exp = base;
      if (n->left) {
        n->kr = (n->kr + n->left->kr * base) % prime;
        n->base_exp = (n->left->base_exp * n->base_exp) % prime;
      }
      if (n->right) {
        n->kr = (n->kr * n->right->base_exp + n->right->kr) % prime;
        n->base_exp = (n->base_exp * n->right->base_exp) % prime;
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
          std::cout << "Modified" << std::endl;
          rightRotate(x->parent);
          rightRotate(x->parent);
        } else {
          rightRotate(x->parent->parent);
          rightRotate(x->parent);
        }
      } else if (x->parent->right == x &&
                 x->parent->parent->right == x->parent) {
        if constexpr (isModified == modified) {
          std::cout << "Modified WARNING" << std::endl;
          leftRotate(x->parent);
          leftRotate(x->parent);
        } else {
          leftRotate(x->parent->parent);
          leftRotate(x->parent);
        }
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
  void join(SplayTree &rightTree) {
    if (!root)
      root = rightTree.getRoot();
    else if (rightTree.getRoot()) {
      node *maxNode = root;
      while (maxNode->right)
        maxNode = maxNode->right;
      splay<normal>(maxNode);
      maxNode->right = rightTree.getRoot();
      rightTree.getRoot()->parent = maxNode;
      updateNodeInfo(maxNode);
    }
  }

  // Split the tree into two trees based on a position
  void split(int position, SplayTree &leftTree, SplayTree &rightTree) {
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

  for (int i = 0; i < 10; i++) {
    uint32_t string_length = rand() % 10000;
    std::vector<char> chars;
    chars.reserve(string_length);
    for (int j = 0; j < string_length; j++) {
      chars.push_back('a' + rand() % 26);
    }

    std::vector<char> secondChars;
    secondChars = chars;
    // pop a random character from the second string
    uint32_t position = rand() % string_length;
    std::cout << "Going to remove character at position: " << position
              << std::endl;
    secondChars.erase(secondChars.begin() + position);

    SplayTree tree;
    tree.root = tree.buildBalancedTree(chars, 0, chars.size() - 1);
    SplayTree secondTree;
    secondTree.root =
        secondTree.buildBalancedTree(secondChars, 0, secondChars.size() - 1);

    int start = rand() % 100;
    int secondStart = start; // rand() % 100;
    std::cout << "Start: " << start << " Second start: " << secondStart
              << std::endl;

    std::cout << "The LCP should be either " << position - start << " or "
              << string_length - start << std::endl;

    uint32_t LCP_value = tree.LCP(start, secondTree, secondStart);

    // chek if the LCP value is correct
    uint32_t scan = 0;
    while (start + scan < string_length && secondStart + scan < string_length &&
           chars[start + scan] == secondChars[secondStart + scan]) {
      scan++;
    }

    if (scan != LCP_value) {
      std::cout << "Error: LCP value is not correct" << std::endl;
      std::cout << "Computed LCP: " << LCP_value << " Scan LCP: " << scan
                << std::endl;

      // print the two strings for verification
      std::string firstString =
          std::string(chars.begin(), chars.end()).substr(start, scan + 1);
      std::string secondString =
          std::string(secondChars.begin(), secondChars.end())
              .substr(secondStart, scan + 1);
      std::cout << "First string: \t" << firstString << std::endl;
      std::cout << "Second string: \t" << secondString << std::endl;

      firstString =
          std::string(chars.begin(), chars.end()).substr(start, LCP_value);
      secondString = std::string(secondChars.begin(), secondChars.end())
                         .substr(secondStart, LCP_value);
      std::cout << "First string: \t" << firstString << std::endl;
      std::cout << "Second string: \t" << secondString << std::endl;
      return 1;
    } else {
      std::cout << "LCP computation is correct" << std::endl;
      std::cout << "Computed LCP: " << LCP_value << " Scan LCP: " << scan
                << std::endl;
      // print the two strings for verification
      std::string firstString = tree.retrieve(start, start + scan);
      std::string secondString =
          secondTree.retrieve(secondStart, secondStart + scan);
      std::cout << "First string: " << firstString << std::endl;
      std::cout << "Second string: " << secondString << std::endl;
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
/*   } */
