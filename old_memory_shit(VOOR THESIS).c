
tri_tree_list list_to_tree_list(triangle_list * list) {
  size_t dim_size = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  tri_tree_list result = {NULL,dim_size,{list->dim[0],list->dim[1],list->dim[2]}};
  tree_node *** t_arr = calloc(dim_size, sizeof(tree_node **));
  for (int i = 0; i < dim_size; i++) 
    t_arr[i] = calloc(dim_size, sizeof(tree_node *));
  //Now lets set all the nodes, exciting!
  for (size_t i = 0; i < list->len; i++) {
    for (int j=0; j < 48; j++) {
      triangle cur_triang = triangle_symmetry(&list->t_arr[i], j, list->dim[0]);
      tri_unique_int(&cur_triang);  
      int index1 = index_vertex(cur_triang.vertices[0], list->dim);
      int index2 = index_vertex(cur_triang.vertices[1], list->dim);
      tri_int index3 = index_vertex(cur_triang.vertices[2], list->dim);
      tree_insert(&t_arr[index1][index2],index3);
    }
  }
  result.t_arr = t_arr;
  return result;  
}


tri_int_list list_to_int_list(triangle_list list) {
  tri_int_list result = {NULL,NULL, list.len, {list.dim[0],list.dim[1],list.dim[2]}};
  tri_int * t_arr = malloc(list.len * sizeof(tri_int));
  for (size_t i = 0; i < list.len; i++) {
    tri_unique_int(&list.t_arr[i]);
    t_arr[i] = triangle_to_int(&list.t_arr[i]);
  }
  result.t_arr = t_arr;
  sort_int_list(&result);
  markers_int_list(&result);
  return result;
}

void dup_int_list(tri_int_list *list) {
  if (list->len == 0)
    return;
  int c = 0;
  for (int i = 1; i < list->len; i++) {
    if (list->t_arr[i] != list->t_arr[c]) {
      c++;
      list->t_arr[c] = list->t_arr[i];
    }
  }
  list->len = c + 1;
  list->t_arr = realloc(list->t_arr, (c + 1) * sizeof(tri_int));
}

#define index_3d(x,y,z,dim) (x * (dim[1] + 1) * (dim[2] + 1) +  y * (dim[2] + 1) + z)

int triangle_in_int_list(ptriangle tri, tri_int_list * list){
  int marker_start = index_3d(tri->vertices[0][0], tri->vertices[0][1], tri->vertices[0][2], list->dim);
  tri_int *ptr = list->markers[marker_start];
  tri_int *end = list->markers[marker_start + 1];
  tri_int search = triangle_to_int(tri);
  
  while (ptr != end) {
    if ((*ptr) == search)
      return ptr-list->t_arr + 1;
    ptr++;   
  }
}


int compare (const void * a, const void * b)
{ 
  int64_t tmp = *(tri_int*)a - *(tri_int*)b;
   return (tmp > 0) - (tmp < 0);
}

void sort_int_list(tri_int_list *list){
  qsort (list->t_arr, list->len, sizeof(tri_int), compare);
}

void markers_int_list(tri_int_list *list) {
  int markers_len = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  tri_int **markers = calloc(markers_len + 1, sizeof(tri_int *));
  int x,y,z;
  x = y = z = -1;
  for (unsigned int i = 0; i < list->len; i++) {
    triangle cur = int_to_triangle(list->t_arr[i]);
    if (cur.vertices[0][0] != x || cur.vertices[0][1] != y || cur.vertices[0][2] != z) {//New marker!
      x = cur.vertices[0][0];
      y = cur.vertices[0][1];
      z = cur.vertices[0][2];
      markers[index_3d(x,y,z,list->dim)] = &list->t_arr[i];
    }            
  }
  /*
   * Some markers may not be defined, as there are no triangles with such a startpoint.
   * Se the start points of these to the next one in line. Work backwards to ensure that all are
   * set. Set additional marker as the endpoint, for computational reasons.
   */
  markers[markers_len] = &list->t_arr[list->len - 1];
  for (int i = markers_len - 1; i >= 0; i--) {
    if (markers[i] == NULL)
      markers[i] = markers[i+1];
  }
  list->markers = markers; 
}



triangle_list triangle_list_symmetries(triangle_list list) {
  triangle_list result = {NULL, list.len * 48, {list.dim[0],list.dim[1],list.dim[2]}};
  ptriangle t_arr = malloc(result.len * sizeof(triangle));
  for (size_t i=0; i < list.len; i++) {
    for (int j=0; j < 48; j++) {
      t_arr[i*48 + j] = triangle_symmetry(&list.t_arr[i], j, list.dim[0]);
    }
  }
  result.t_arr = t_arr;
  return result;
}

tri_bit_list triangle_bit_symmetries(triangle_list * list) {
  int Rc;
  int size = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  Pvoid_t * t_arr = calloc(size*size, sizeof(Pvoid_t));
  //Pvoid_t t_arr = NULL;
  PWord_t  PValue;  
  Word_t Index;
  for (size_t i = 0; i < list->len; i++) {
    for (int j=0; j < 48; j++) {
      triangle cur_triang = triangle_symmetry(&list->t_arr[i], j, list->dim[0]);
      tri_unique_int(&cur_triang);
      
      int k = index_vertex(cur_triang.vertices[0],list->dim);
      int j = index_vertex(cur_triang.vertices[1],list->dim);
      int l = index_vertex(cur_triang.vertices[2],list->dim);
      
      J1S(Rc, t_arr[k*size+j],l);
      /*
      Index = triangle_to_int(&cur_triang);
      JLI(PValue, t_arr, Index);
      if (PValue == PJERR) {
        printf("ERROR BITCHES");
        exit(0);
      }
      *PValue =1;
      * */
    }
  }
  
  int totalc = 0;
  Word_t  Rc_word;
  //JLC(Rc_word, t_arr, 0, -1);
  //totalc = Rc_word;
  
  for (int i = 0; i < size*size; i++){
    if (t_arr[i]) {
      J1C(Rc_word, t_arr[i], 0 , -1);
      totalc += Rc_word;
    }
  }
  printf("TOTAAL AANTAL... : %d", totalc);
  
  return (tri_bit_list) {NULL, totalc};
}

size_t tree_list_count(tri_tree_list * list) {
  size_t result = 0;
  for (size_t i = 0; i < list->dim_size; i++)
    for (size_t j = 0; j < list->dim_size; j++) {
      result += tree_count(list->t_arr[i][j]);      
    }
  return result;
}


/*
 * Bitwise
 */
tri_int triangle_to_int(ptriangle tri){
  tri_int res =  tri->vertices[2][2] | (tri->vertices[2][1] << 5) | (tri->vertices[2][0] << 10) | \
                (tri->vertices[1][2] << 15) | (tri->vertices[1][1] << 20) | (tri->vertices[1][0] << 25) | \
                ((tri_int) tri->vertices[0][2] << 30) | ((tri_int) tri->vertices[0][1] << 35) | ((tri_int)tri->vertices[0][0] << 40);
  return res;
}

triangle int_to_triangle(tri_int tri) {
  triangle res;
  res.vertices[0][0] = (tri >> 40) & BIT_MASK;
  res.vertices[0][1] = (tri >> 35) & BIT_MASK;
  res.vertices[0][2] = (tri >> 30) & BIT_MASK;
  res.vertices[1][0] = (tri >> 25) & BIT_MASK;
  res.vertices[1][1] = (tri >> 20) & BIT_MASK;
  res.vertices[1][2] = (tri >> 15) & BIT_MASK;
  res.vertices[2][0] = (tri >> 10) & BIT_MASK;
  res.vertices[2][1] = (tri >>  5) & BIT_MASK;
  res.vertices[2][2] = (tri      ) & BIT_MASK;
  return res;
}

triangle bint_to_triangle(arr3 indices, arr3 dim) {
  triangle res;
  int XY = (dim[0] + 1) * (dim[1] + 1);
  int X  = (dim[0] + 1);
  int t;
  res.vertices[0][2] = indices[0] / XY;
  t = indices[0] % XY;
  res.vertices[0][1] = t / X;
  res.vertices[0][0] = t % X;
  res.vertices[1][2] = indices[1] / XY;
  t = indices[1] % XY;
  res.vertices[1][1] = t / X;
  res.vertices[1][0] = t % X;
  res.vertices[2][2] = indices[2] / XY;
  t = indices[2] % XY;
  res.vertices[2][1] = t / X;
  res.vertices[2][0] = t % X;
  return res;

}

#define index_vertex(vertex,dim) (vertex[2] * (dim[0] + 1) * (dim[1] + 1) + vertex[1] * (dim[0] + 1) + vertex[0])


/*
 * Three dimensional list. Third dimension is a bit-array
 */
tri_bit_list list_to_bit_list(triangle_list * list){
  
  int size = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  Pvoid_t * t_arr = calloc(size*size, sizeof(Pvoid_t));
  for (size_t i = 0; i < list->len; i++) {
    int k = index_vertex(list->t_arr[i].vertices[0],list->dim);
    int j = index_vertex(list->t_arr[i].vertices[1],list->dim);
    int l = index_vertex(list->t_arr[i].vertices[2],list->dim);
    int Rc;
    J1S(Rc, t_arr[k*size+j],l);   
  }
  return (tri_bit_list) {t_arr, size*size};
}


int compare_pt(arr3 A, arr3 B) {
  if (A[0] != B[0])
    return A[0] - B[0];
  else if (A[1] != B[1])
    return A[1] - B[1];
  else
    return A[2] - B[2];
}

tri_int tri_unique_int(ptriangle triang) {
  int swapped = 1;
  int t;
  while (swapped) {
    swapped = 0;
    for (int i = 1; i < 3; i++) {
      int cmp = compare_pt(triang->vertices[i-1], triang->vertices[i]);
      if (cmp > 0) {
        swap(triang->vertices[i-1][0],triang->vertices[i][0]);
        swap(triang->vertices[i-1][1],triang->vertices[i][1]);
        swap(triang->vertices[i-1][2],triang->vertices[i][2]);
        swapped = 1;
      }
    }
  }
}

tri_int_list list_to_int_list(triangle_list list);
void sort_int_list(tri_int_list *list);
void markers_int_list(tri_int_list *list);
int triangle_in_int_list(ptriangle tri, tri_int_list * list);
void dup_int_list(tri_int_list *list);

tri_tree_list list_to_tree_list(triangle_list * list);
size_t tree_list_count(tri_tree_list * list);
tri_bit_list list_to_bit_list(triangle_list * list);
triangle_list triangle_list_symmetries(triangle_list list);
tri_bit_list triangle_bit_symmetries(triangle_list * list);


typedef uint64_t tri_int;


#define BIT_MASK 31

typedef struct tri_int_list
{
  tri_int * t_arr;
  tri_int **markers;
  size_t len;
  arr3 dim;
} tri_int_list;

typedef struct tri_tree_list
{
  tree_node *** t_arr; //t_arr[vertex[0]][vertex[1]], then it's a BST
  size_t dim_size;
  arr3 dim;
} tri_tree_list;

typedef struct tri_bit_list
{
  Pvoid_t * t_arr;
  size_t len;
} tri_bit_list;



void test_trees(void) {
  tri_int test_list[] = {0,5,1,7,12,3,6,7,8,1,23,505,555,125346,43634,3253,2342,34234};
  tree_node * root = NULL, * find;
  for (int i= 0; i < sizeof(test_list)/sizeof(test_list[0]); i++) {
    if (!tree_insert(&root, test_list[i]))
      printf("Value was already inserted: %d\n", test_list[i]);
    else
      printf("Succesfully inserted: %d\n", test_list[i]);
  }
  tree_print(root);
  
  printf("\nTotal length of our tree:%d\n",tree_count(root));
  find = tree_find(8, root);
  printf("Find: %d,%d,%d\n", find->key, find->left, find->right);
  tree_delete(root, 8);
  tree_print(root);
  printf("\n");  
}

void test_lists() {
    arr3 dim =  {4,4,4};
  triangle_list tri_list;
  tri_int_list int_list;
  tri_list = acute_triangle_recur(dim);
  int_list = list_to_int_list(tri_list);
  int idx =  triangle_in_int_list(&tri_list.t_arr[1000],&int_list);
  printf("%d\n",idx);
  triangle cur = int_to_triangle(int_list.t_arr[idx-1]);
  print_triangle(&cur);
  print_triangle(&tri_list.t_arr[1000]);  
}

tri_int triangle_to_int(ptriangle tri);
triangle int_to_triangle(tri_int tri);

/*
 * Convert a memory list from the triangle list. 
 * 
 * NOT USED ANYMORE?
 */
tri_mem_list mem_list_from_triangle_list(triangle_list * list) {
  printf("Not implemented, mem_list_from_triangle_list");
  /*
  arr3 dim_size;
  dim_size[0] = (list->dim[0] + 1) * (list->dim[1] + 1) * (list->dim[2] + 1);
  dim_size[1] = dim_size[0];
  dim_size[2] = dim_size[1]/8 + 1;  
  tri_mem_list result = mem_list_init(list->dim,dim_size);
  //Now lets set all the nodes, exciting!
  for (int i = 0; i < list->len; i++) 
    mem_list_set_sym(&result, &list->t_arr[i]);
  return result;  
  */
}
