#include "enhanced_splay_tree.hpp"
#include <chrono>
#include <unordered_map>
#include <vector>

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
