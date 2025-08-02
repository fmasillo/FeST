# FeST

Implementation of the dynamic string solution based on a Forest of enhanced Splay Trees (FeST).

## Installation

```sh
git clone https://github.com/fmasillo/FeST.git
cd FeST
make
```

## Example usage

In the file ``FeST.cpp``, you can find an example of how to use the FeST data structure. The code in the main function includes some simple operations that can be performed on a FeST, such as inserting a new string or concatenating two existing strings into a single one. Then, some longest common prefix (LCP) queries are done on random strings.

The core component of this repository is found in file ``enhanced_splay_tree.hpp``. This file contains the functions which operate on the individual trees. The set of operations supported is the following:

```c++

// Insert a new character at a given position
void insert(const unsigned char c, const uint32_t position);

// Edit a character at a given position
void edit(const unsigned char c, const uint32_t position);

// Delete node at a given position
void delete_node(const int position);

// Introduce a substring from another tree at a given position
void introduce(const uint32_t pos, SplayTree &other);

// Extract a substring from the tree and return the root of the extracted subtree
void extract(const uint32_t i, const uint32_t j, SplayTree &other);

// Retrieve a substring from the tree
std::string retrieve(const uint32_t i, const uint32_t j);

// Substring equality check
bool equal(const uint32_t i, SplayTree &other, const uint32_t j, const uint32_t length);

// Compute the LCP of two substrings
uint32_t LCP(uint32_t i, SplayTree &other, uint32_t j);
```

## Citation

Conference paper:

Zsuzsanna Lipták, Francesco Masillo, and Gonzalo Navarro. A Textbook Solution for Dynamic Strings. In Proc. of the 32nd Annual European Symposium on Algorithms (ESA 2024), volume 308 of LIPIcs, pages 86:1-86:16. Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2024: <https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2024.86>

```bibtex
@inproceedings{LiptakMN24,
  author       = {{\relax Zs}uzsanna Lipt{\'{a}}k and
                  Francesco Masillo and
                  Gonzalo Navarro},
  editor       = {Timothy M. Chan and
                  Johannes Fischer and
                  John Iacono and
                  Grzegorz Herman},
  title        = {A Textbook Solution for Dynamic Strings},
  booktitle    = {Proc. of the 32nd Annual European Symposium on Algorithms, {ESA} 2024, September
                  2-4, 2024, Royal Holloway, London, United Kingdom},
  series       = {LIPIcs},
  volume       = {308},
  pages        = {86:1--86:16},
  publisher    = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year         = {2024},
  url          = {https://doi.org/10.4230/LIPIcs.ESA.2024.86},
  doi          = {10.4230/LIPICS.ESA.2024.86}
}
```
