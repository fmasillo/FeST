#include <cassert>
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

  // Construct a fully balanced binary tree from a sorted array of characters
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

    updateSubtreeSize(n);
    updateKR(n);
    return n;
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

  // Splay operation
  template <auto Modified> void splay(node *x) {
    while (x->parent) {
      if (!x->parent->parent) {
        if (x->parent->left == x) {
          rightRotate(x->parent);
        } else {
          leftRotate(x->parent);
        }
      } else if (x->parent->left == x && x->parent->parent->left == x->parent) {
        if constexpr (Modified == modified) {
          rightRotate(x->parent);
          rightRotate(x->parent);
        } else {
          rightRotate(x->parent->parent);
          rightRotate(x->parent);
        }
      } else if (x->parent->right == x &&
                 x->parent->parent->right == x->parent) {
        if constexpr (Modified == modified) {
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

    updateSubtreeSize(root);
    updateKR(root);
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
    else if (x == x->parent->right)
      x->parent->right = y;
    else
      x->parent->left = y;
    y->right = x;
    x->parent = y;

    updateSubtreeSize(x);
    updateSubtreeSize(y);
    updateKR(x);
    updateKR(y);
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

    updateSubtreeSize(x);
    updateSubtreeSize(y);
    updateKR(x);
    updateKR(y);
  }

  // Insert a new character at a given position using only split operations
  void insert(const char c, const int position) {
    if (position == 0) {
      node *n = new node(c);
      n->right = root;
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
    updateSubtreeSize(root);
    updateKR(root);
  }

  // Find a node at a given position
  template <auto Modified> node *find(const int position) {
    node *z = root;
    int currentPos = position;

    while (z) {
      int leftSize = z->left ? z->left->subtree_size : 0;
      if (currentPos == leftSize) {
        splay<Modified>(z);
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
      updateSubtreeSize(maxNode);
      updateKR(maxNode);
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
      updateSubtreeSize(rightTree.root);
      updateSubtreeSize(leftTree.root);
    }
  }

  // Delete implementation using split and join
  void deleteNode(const int position) {
    if (position == 0) {
      find<normal>(0);
      root = root->right;
      root->parent = nullptr;
      updateSubtreeSize(root);
      updateKR(root);
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
      find<modified>(i - 1);
      assert(root->right != nullptr);
      return root->right;
    }

    find<normal>(j + 1);
    find<modified>(i - 1);

    return root->right->left;
  }

  void introduce(int pos, SplayTree &other) {
    if (!root) {
      if (other.getRoot())
        root = other.getRoot();
    } else if (pos <= root->subtree_size && pos > 0) {
      isolate(pos, pos - 1);

      if (root->right) {
        node *toConnect = root->right;
        toConnect->left = other.getRoot();
        if (other.getRoot())
          other.getRoot()->parent = toConnect;

        updateSubtreeSize(toConnect);
        updateKR(toConnect);
        updateSubtreeSize(root);
        updateKR(root);
      }
    } else if (pos == 0) {
      other.join(*this);
    }
  }

  SplayTree extract(int i, int j) {
    if (i < 0 || j >= root->subtree_size || i > j) {
      std::cout << "Invalid range" << std::endl;
      return SplayTree();
    }

    node *newRoot = isolate(i, j);
    newRoot->parent = nullptr;
    SplayTree extractedTree;
    extractedTree.root = newRoot;

    if (root->right) {
      root->right->left = nullptr;
      updateSubtreeSize(root->right);
      updateKR(root->right);
    }
    updateSubtreeSize(root);
    updateKR(root);

    return extractedTree;
  }

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

    return leftSubstring->kr == rightSubstring->kr;
  }

  std::string retrieve(int i, int j) {
    node *subtree = isolate(i, j);
    std::string result;
    inorder(subtree, result);
    return result;
  }

  void inorder(node *n, std::string &result) {
    if (n) {
      inorder(n->left, result);
      result.push_back(n->character);
      inorder(n->right, result);
    }
  }

  // Visualize the tree
  void visualize(node *n, int depth = 0) {
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
};

int main() {
  std::vector<char> chars = {'m', 'i', 's', 's', 'i', 's',
                             's', 'i', 'p', 'p', 'i'};
  SplayTree tree;
  tree.root = tree.buildBalancedTree(chars, 0, chars.size() - 1);

  std::cout << "Root character: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  std::cout << "Find node at position 5" << std::endl;
  tree.find<normal>(8);
  tree.visualize(tree.getRoot());

  tree.insert('h', 3); // Insert 'h' at position 3
  std::cout << "Root after insertion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.insert('i', 0); // Insert 'i' at position 0
  std::cout << "Root after insertion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.deleteNode(3); // Delete node at position 3
  std::cout << "Root after deletion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.find<normal>(3); // Find node at position 3
  std::cout << "Found node: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.deleteNode(0); // Delete node at position 1
  std::cout << "Root after deletion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  node *subtring = tree.isolate(3, 6); // Isolate substring from position 3 to 6
  std::cout << "Isolated substring: " << subtring->character << std::endl;
  tree.visualize(subtring);

  std::string firstSubstring =
      tree.retrieve(3, 4); // Retrieve substring from position 3 to 4
  std::string secondSubstring =
      tree.retrieve(6, 7); // Retrieve substring from position 6 to 7
  bool equal = tree.equal(3, tree, 6, 2); // Check if two substrings are equal
  std::cout << "Are two substrings (" << firstSubstring << ", "
            << secondSubstring << ") equal: " << (equal ? "TRUE" : "FALSE")
            << std::endl;

  firstSubstring = tree.retrieve(0, 6);
  secondSubstring = tree.retrieve(3, 9);
  equal = tree.equal(3, tree, 0, 7);
  std::cout << "Are two substrings (" << firstSubstring << ", "
            << secondSubstring << ") equal: " << (equal ? "TRUE" : "FALSE")
            << std::endl;

  node *strangeSubstring = tree.isolate(6, 5);
  if (strangeSubstring == nullptr) {
    std::cout << "There is no root->right->left node" << std::endl;
  } else {
    tree.visualize(strangeSubstring);
  }

  std::vector<char> secondChars = {'h', 'e', 'l', 'l', 'o'};
  SplayTree secondTree;
  secondTree.root =
      secondTree.buildBalancedTree(secondChars, 0, secondChars.size() - 1);

  tree.introduce(3, secondTree);
  std::cout << "Root after introducing substring: " << tree.root->character
            << std::endl;
  tree.visualize(tree.getRoot());

  SplayTree extractedTree = tree.extract(2, 6);
  std::cout << "Extracted tree: " << extractedTree.root->character << std::endl;
  extractedTree.visualize(extractedTree.getRoot());
  tree.visualize(tree.getRoot());

  return 0;
}

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
