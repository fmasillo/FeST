#include <iostream>
using namespace std;

struct Node {
  int key;
  int subtree_size; // Dimensione del sottoalbero
  Node *parent;
  Node *left;
  Node *right;
  unsigned char character;

  Node(int value, unsigned char c = 0)
      : key(value), subtree_size(1), parent(nullptr), left(nullptr),
        right(nullptr), character(c) {}
};

class SplayTree {
private:
  Node *root;

  // Aggiorna la dimensione del sottoalbero di un nodo
  void updateSize(Node *node) {
    if (node) {
      node->subtree_size = 1;
      if (node->left) {
        node->subtree_size += node->left->subtree_size;
      }
      if (node->right) {
        node->subtree_size += node->right->subtree_size;
      }
    }
  }

  // Rotazione a destra
  void rotateRight(Node *x) {
    Node *y = x->left;
    if (y) {
      x->left = y->right;
      if (y->right) {
        y->right->parent = x;
      }
      y->parent = x->parent;
    }

    if (!x->parent) {
      root = y;
    } else if (x == x->parent->right) {
      x->parent->right = y;
    } else {
      x->parent->left = y;
    }

    if (y) {
      y->right = x;
    }
    x->parent = y;

    // Aggiorna la dimensione dei sottoalberi
    updateSize(x);
    updateSize(y);
  }

  // Rotazione a sinistra
  void rotateLeft(Node *x) {
    Node *y = x->right;
    if (y) {
      x->right = y->left;
      if (y->left) {
        y->left->parent = x;
      }
      y->parent = x->parent;
    }

    if (!x->parent) {
      root = y;
    } else if (x == x->parent->left) {
      x->parent->left = y;
    } else {
      x->parent->right = y;
    }

    if (y) {
      y->left = x;
    }
    x->parent = y;

    // Aggiorna la dimensione dei sottoalberi
    updateSize(x);
    updateSize(y);
  }

  // Operazione di splaying: porta il nodo `x` alla radice
  void splay(Node *x) {
    while (x->parent) {
      if (!x->parent->parent) { // Caso Zig
        if (x == x->parent->left) {
          rotateRight(x->parent);
        } else {
          rotateLeft(x->parent);
        }
      } else if (x == x->parent->left &&
                 x->parent == x->parent->parent->left) { // Caso Zig-Zig
        rotateRight(x->parent->parent);
        rotateRight(x->parent);
      } else if (x == x->parent->right &&
                 x->parent == x->parent->parent->right) { // Caso Zig-Zag
        rotateLeft(x->parent->parent);
        rotateLeft(x->parent);
      } else if (x == x->parent->left &&
                 x->parent == x->parent->parent->right) { // Caso Zig-Zag
        rotateRight(x->parent);
        rotateLeft(x->parent);
      } else {
        rotateLeft(x->parent);
        rotateRight(x->parent);
      }
    }
    updateSize(x); // Assicurati di aggiornare la dimensione della radice
  }

  Node *subtreeMin(Node *node) {
    while (node->left) {
      node = node->left;
    }
    return node;
  }

  void replace(Node *u, Node *v) {
    if (!u->parent) {
      root = v;
    } else if (u == u->parent->left) {
      u->parent->left = v;
    } else {
      u->parent->right = v;
    }
    if (v) {
      v->parent = u->parent;
    }
  }

public:
  SplayTree() : root(nullptr) {}

  SplayTree(Node *node) : root(node) {}

  // Getter per root
  Node *getRoot() { return root; }

  // Funzione di ricerca
  Node *search(const int key) {
    Node *x = root;
    while (x) {
      if (key < x->key) {
        x = x->left;
      } else if (key > x->key) {
        x = x->right;
      } else {
        splay(x); // Porta il nodo trovato alla radice
        return x;
      }
    }
    return nullptr; // Chiave non trovata
  }

  // Inserimento di un nodo
  void insert(const int key, const unsigned char character = 0) {
    Node *z = root;
    Node *p = nullptr;

    while (z) {
      p = z;
      if (key < z->key) {
        z = z->left;
      } else {
        z = z->right;
      }
    }

    z = new Node(key, character);
    z->parent = p;

    if (!p) {
      root = z;
    } else if (key < p->key) {
      p->left = z;
    } else {
      p->right = z;
    }

    splay(z); // Porta il nodo inserito alla radice
  }

  // Cancellazione di un nodo
  void erase(const int key) {
    Node *x = search(key);
    if (!x)
      return;

    splay(x); // Porta il nodo da cancellare alla radice

    if (!x->left) {
      replace(x, x->right);
    } else if (!x->right) {
      replace(x, x->left);
    } else {
      Node *y = subtreeMin(x->right);
      if (y->parent != x) {
        replace(y, y->right);
        y->right = x->right;
        y->right->parent = y;
      }
      replace(x, y);
      y->left = x->left;
      y->left->parent = y;
    }

    updateSize(root);
    delete x;
  }

  // Join two splay trees
  Node *join(Node *leftTree, Node *rightTree) {
    if (!leftTree)
      return rightTree;
    if (!rightTree)
      return leftTree;

    // Find the maximum Node in the left tree
    Node *maxNode = leftTree;
    while (maxNode->right) {
      maxNode = maxNode->right;
    }

    // Splay the maximum node to the root of the left tree
    splay(maxNode);

    // Attach the right tree to the right of the left tree's root
    maxNode->right = rightTree;
    if (rightTree)
      rightTree->parent = maxNode;

    // Update the size of the left tree
    updateSize(maxNode);
    return maxNode;
  }

  // Split the splay tree into two trees based on a key
  void split(int key, SplayTree &leftTree, SplayTree &rightTree) {
    splay(search(key));
    if (root->key <= key) {
      leftTree.root = root;
      rightTree.root = root->right;
      if (rightTree.root)
        rightTree.root->parent = nullptr;
      leftTree.root->right = nullptr;
    } else {
      rightTree.root = root;
      leftTree.root = root->left;
      if (leftTree.root)
        leftTree.root->parent = nullptr;
      rightTree.root->left = nullptr;
    }

    // Update the size of the trees
    leftTree.updateSize(leftTree.root);
    rightTree.updateSize(rightTree.root);
  }

  // Stampa in ordine
  void inorder(const Node *node) {
    if (!node)
      return;
    inorder(node->left);
    cout << "Chiave: " << node->key << ", |sottoalbero|: " << node->subtree_size
         << endl;
    inorder(node->right);
  }

  void inorderTraversal() {
    inorder(root);
    cout << endl;
  }

  // Stampa albero visivamente
  void show(Node *node, int spaces) {
    if (!node)
      return;
    spaces += 5;
    show(node->right, spaces);
    cout << endl;
    for (int i = 5; i < spaces; i++)
      cout << " ";
    cout << node->key << ": " << node->character << " (" << node->subtree_size
         << ")" << endl;
    show(node->left, spaces);
  }
};

int main() {
  SplayTree tree;

  tree.insert(10, 'a');
  tree.insert(20, 'b');
  tree.insert(30, 'c');
  tree.insert(40, 'd');
  tree.insert(50, 'e');

  tree.insert(60, 'f');
  tree.insert(70, 'g');
  tree.insert(80, 'h');
  tree.insert(90, 'i');
  tree.insert(100, 'j');

  cout << "Albero in ordine dopo inserimenti: " << endl;
  tree.inorderTraversal();

  tree.show(tree.getRoot(), 0);

  cout << "Cercando il nodo con chiave 30: " << endl;
  tree.search(30);
  tree.inorderTraversal();
  tree.show(tree.getRoot(), 0);

  cout << "Cancellando il nodo con chiave 20: " << endl;
  tree.erase(20);
  tree.inorderTraversal();
  tree.show(tree.getRoot(), 0);

  cout << "Splitting the splay tree at key 60: " << endl;
  SplayTree leftTree, rightTree;
  tree.split(60, leftTree, rightTree);
  cout << "Left tree: " << endl;
  tree.show(leftTree.getRoot(), 0);

  cout << "Right tree: " << endl;
  tree.show(rightTree.getRoot(), 0);

  cout << "Joining the two trees: " << endl;

  SplayTree newTree(leftTree.join(leftTree.getRoot(), rightTree.getRoot()));

  newTree.show(newTree.getRoot(), 0);

  return 0;
}
