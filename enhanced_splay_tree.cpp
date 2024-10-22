#include <iostream>
#include <vector>

struct node {
  node *parent;
  node *left;
  node *right;
  int subtree_size;
  char character;

  node(char c)
      : parent(nullptr), left(nullptr), right(nullptr), subtree_size(1),
        character(c) {}
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
  }

  // Insert a new character at a given position
  void insert(char c, int position) {
    node *z = root;
    node *p = nullptr;
    int currentPos = 0;

    while (z) {
      p = z;
      int leftSize = z->left ? z->left->subtree_size : 0;
      if (position <= currentPos + leftSize) {
        z = z->left;
      } else {
        currentPos += leftSize + 1;
        z = z->right;
      }
    }

    z = new node(c);
    z->parent = p;
    if (!p)
      root = z;
    else if (position <= currentPos)
      p->left = z;
    else
      p->right = z;

    splay(z);
  }

  // Find a node at a given position
  node *find(int position) {
    node *z = root;
    int currentPos = 0;

    while (z) {
      int leftSize = z->left ? z->left->subtree_size : 0;
      if (position < currentPos + leftSize) {
        z = z->left;
      } else if (position > currentPos + leftSize) {
        currentPos += leftSize + 1;
        z = z->right;
      } else {
        splay(z);
        return z;
      }
    }
    return nullptr;
  }

  // Delete a node at a given position
  void deleteNode(int position) {
    node *z = find(position);
    if (!z)
      return;

    splay(z);

    if (!z->left) {
      replace(z, z->right);
    } else if (!z->right) {
      replace(z, z->left);
    } else {
      node *y = minimum(z->right);
      if (y->parent != z) {
        replace(y, y->right);
        y->right = z->right;
        y->right->parent = y;
      }
      replace(z, y);
      y->left = z->left;
      y->left->parent = y;
    }

    updateSubtreeSize(root);
    delete z;
  }

  // Replace one subtree as a child of its parent with another subtree
  void replace(node *u, node *v) {
    if (!u->parent)
      root = v;
    else if (u == u->parent->left)
      u->parent->left = v;
    else
      u->parent->right = v;
    if (v)
      v->parent = u->parent;
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
    }
  }

  // Split the tree into two trees based on a poisition
  // The left tree will contain all the nodes with positions less than or equal
  // to the position, and the right tree will contain all the nodes with
  // positions greater than the position
  void split(int position, SplayTree &leftTree, SplayTree &rightTree) {
    node *z = find(position);
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
  void delete2(int position) {
    SplayTree leftTree, rightTree;
    split(position, leftTree, rightTree);
    rightTree.split(0, rightTree, rightTree);
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
      std::cout << " : " << n->subtree_size << std::endl;
      visualize(n->left, depth + 1);
    }
  }
};

int main() {
  std::vector<char> chars = {'a', 'b', 'c', 'd', 'e', 'f', 'g'};
  SplayTree tree;
  tree.root = tree.buildBalancedTree(chars, 0, chars.size() - 1);

  std::cout << "Root character: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.insert('h', 3); // Insert 'h' at position 3
  std::cout << "Root after insertion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.delete2(3); // Delete node at position 3
  std::cout << "Root after deletion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.find(3); // Find node at position 3
  std::cout << "Found node: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  tree.delete2(0); // Delete node at position 1
  std::cout << "Root after deletion: " << tree.root->character << std::endl;
  tree.visualize(tree.getRoot());

  return 0;
}
