#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tetraeder.h"
#include "triangle.h"
#include "datastructures.h"

typedef size_t tri_int;

tree_node * tree_find(tri_int key, tree_node * root) {
  tree_node * current = root;
  while (current) {
    if (current->key == key)
      return current;
      
    if (key < current->key)
      current = current->left;
    else
      current = current->right;
  }
  return NULL;
}


tree_node * tree_insert(tree_node ** root, tri_int key) {
  tree_node ** current = root;
  while (*current) { //As long as we have not found an empty leaf, go on
    if ((*current)->key == key) 
      return NULL; //Already in the tree motherfuckers!
    if (key < (*current)->key)
      current = &((*current)->left); //Store address so we can edit this node
    else
      current = &((*current)->right);
  }
  tree_node * new = calloc(1,sizeof(tree_node));
  new->key = key;
  if (*root == NULL) //We just created the root, let root point to this
    *root = new;
  else //Edit parent to point to this new node
    *current = new; 
  return new;
}
/*
 * Returns the node with the minimum value for this (sub)tree.
 * Stores the parent of this node (if the root is not the minimum one)
 */
tree_node * tree_min(tree_node * root, tree_node **parent) { //Minimum node in this tree
  tree_node * current = root;
  *parent = NULL; //No parent found
  while (current->left) {
    *parent = current; //One level deeper, set parent
    current = current->left; //Left is smaller!
  }
  return current;
}

void tree_delete(tree_node * root, tri_int key) {
  tree_node ** current = &root;
  tree_node * parent = NULL;
  while (*current && ((*current)->key != key)) {
    parent = *current;
    if (key < (*current)->key)
      current = &(*current)->left;
    else
      current = &(*current)->right;
  }
  if (*current == NULL) //Key does not exist..
    return;
  if (parent == NULL) //We want to delete the root..
    printf("Deleting root. fuck off\n");
  
  tree_node * del = *current; //Holds the to delete node
  if (del->left && del->right) { //Two children. Find the minimum leaf on the right side, place this as new root
    tree_node * new_root = tree_min(del->right, &parent);

    del->key = new_root -> key; //Set key to the new_root key
    free(new_root); //Delete "new_root"
    if (parent) 
      parent->left = NULL; //Remove reference in parent of new_root
  } else if (del->left) {
    *current = del->left; //Parent left/right now points to the only child of our deleted node
    free(del);
  } else if (del->right) { //Right child only
    *current = del->right; //Parent left/right now points to the only child of our deleted node
    free(del);
  } else {
    free(del);
    *current = NULL; //Current is address of parent->left or parent->right, accordingly. One of parents childs is now empty
  }
}

void tree_print(tree_node * node) {
  if (node == NULL)
    return;
  tree_print(node->left);
  printf("%d, ", node->key);
  tree_print(node->right);
}

size_t tree_count(tree_node * node) {
  if (node == NULL)
    return 0;
  return tree_count(node->left) + tree_count(node->right) + 1;
}
