#ifndef DATASTRUCT_H
#define DATASTRUCT_H
#include <inttypes.h>
#include "vector.h"



typedef struct tree_node{
  uint64_t key;
  struct tree_node *left, *right;
} tree_node;

tree_node * tree_insert(tree_node ** root, uint64_t key);
tree_node * tree_find(uint64_t key, tree_node * root);
void tree_delete(tree_node * root, uint64_t key);
void tree_print(tree_node * node);
size_t tree_count(tree_node * node);
#endif
