#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

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
  void splay(node *x) {
    while (x->parent) {
      if (!x->parent->parent) {
        if (x->parent->left == x) {
          rightRotate(x->parent);
        } else {
          leftRotate(x->parent);
        }
      } else if (x->parent->left == x && x->parent->parent->left == x->parent) {
        rightRotate(x->parent->parent);
        rightRotate(x->parent);
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
      updateSubtreeSize(root);
      updateKR(root);
      return;
    }

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

    updateSubtreeSize(root);
    updateKR(root);
  }

  node *find(const int position) {
    node *z = root;
    int currentPos = position;

    while (z) {
      int leftSize = z->left ? z->left->subtree_size : 0;
      if (currentPos == leftSize) {
        splay(z);
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

  // Find the minimum node in a subtree
  node *minimum(node *x) {
    while (x->left)
      x = x->left;
    return x;
  }

  // Join current tree with another tree
  void join(SplayTree &rightTree) {
    if (!root)
      root = rightTree.getRoot();
    else if (rightTree.getRoot()) {
      node *maxNode = root;
      while (maxNode->right)
        maxNode = maxNode->right;
      splay(maxNode);
      maxNode->right = rightTree.getRoot();
      rightTree.getRoot()->parent = maxNode;
      updateSubtreeSize(maxNode);
      updateKR(maxNode);
    }
  }

  // Split the tree into two trees based on a poisition
  // The left tree will contain all the nodes with positions less than or equal
  // to the position, and the right tree will contain all the nodes with
  // positions greater than the position
  void split(int position, SplayTree &leftTree, SplayTree &rightTree) {
    node *z = find(position);
    if (!z) {
      std::cout << "Z is null" << std::endl;

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
      find(0);
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
  tree.find(8);
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

  tree.find(3); // Find node at position 3
  std::cout << "Found node: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.deleteNode(0); // Delete node at position 1
  std::cout << "Root after deletion: " << tree.root->character << std::endl;
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
