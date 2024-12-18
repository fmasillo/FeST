#include <iostream>
using namespace std;

class Node {
  Node *left, *right;

public:
  Node(Node *l, Node *r);
  ~Node();
};

class Tree {
  Node *root;

public:
  Tree(Node *rt);
  ~Tree();
};

Tree::Tree(Node *rt) : root(rt) {
  cout << "new Tree with root node at " << rt << endl;
}
Tree::~Tree() {
  cout << "Destructor of Tree" << endl;
  if (root)
    delete root;
}
Node::Node(Node *l, Node *r) : left(l), right(r) {
  cout << "Node@" << this << "(left:" << l << ", right:" << r << ")" << endl;
}
Node::~Node() {
  cout << "~Node@" << this << endl;
  if (left)
    delete left;
  if (right)
    delete right;
}

int main() {
  Tree t(new Node(new Node(new Node(new Node(0, 0), 0), 0),
                  new Node(0, new Node(0, 0))));
}
