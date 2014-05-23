
/*
 * Returns whether this tetrahedron only has facets in the mem_list
 */
 
int tetra_facets_list(ptriangle triang, arr3 new_vertex, tri_mem_list * acute_list) {
  //triang itself is included in the list.. Otherwise dumb call!
  tri_index indices[3];
  
  if (acute_list->fund) //We need to symmetry the triangle to a fundamental-area representation first
    return (mem_list_get_fund(acute_list, triang->vertices[0], triang->vertices[1], new_vertex) &&
            mem_list_get_fund(acute_list, triang->vertices[0], triang->vertices[2], new_vertex) &&
            mem_list_get_fund(acute_list, triang->vertices[1], triang->vertices[2], new_vertex)); 
  
  else {
    vertices_to_index(triang->vertices[0], triang->vertices[1], new_vertex, acute_list->dim_mult, indices[0]);
    vertices_to_index(triang->vertices[0], triang->vertices[2], new_vertex, acute_list->dim_mult, indices[1]);
    vertices_to_index(triang->vertices[1], triang->vertices[2], new_vertex, acute_list->dim_mult, indices[2]);
    for (int i = 0; i < 3; i++) {
      if (!EMI(acute_list->t_arr, indices[i]) ||  //Point does not exist or is not acute
          !GMI(acute_list->t_arr, indices[i]))
        return 0;
     } 
    return 1;
  }
}



/*
 * Returns whether this triangle has an acute tetrahedron above and below. More information
 * is stored in result. Convention: normal points to above (inwards for triangles on the boundary).
 * 
 * If acute_list is set, only tetrahedrons with facets in this list may be included.
 */
int triangle_tetra_acute(ptriangle triang, cube_points * cube, ptriang_tetra_result res  , tri_mem_list *acute_list ) {
  //Check if triangle is acute..
  if (!mat3_triangle_acute(triang->vertices)) //Can replace with code below
    return 0;
    
 /* 
  * Initalize variables
  */
  res->boundary_above = 0;
  res->boundary_below = 0;
  res->acute_above    = 0;
  res->acute_below    = 0;
  arr3 * cube_pts = cube->points;
    
  mat3 P;
  //Store sides in matrix P
  //P[0] = V1 - V0
  //P[1] = V2 - V0
  //P[2] = V2 - V1
  triangle_sides(triang->vertices[0], triang->vertices[1], triang->vertices[2],P);
  
  /*
   * Calculate the normals for this triangle
   * - The normal on the triangle plane, to determine on which side a point is
   * - The three normals on the "side"-planes, used to check if a point lies in the toblerone
   *
   * Now calculate the plane equations normal . p = d
   */
  //Convention P[1] cross P[0]
  arr3 tri_normal;
  crossArr3(P[1], P[0], tri_normal);
  int  tri_d = dotArr3(tri_normal, triang->vertices[0]); //For the plane equation of the triangle
  //triangle_normal(triang, normal);
  
  
  mat3 side_normals;
  arr3 side_d; //The constant expression for the plane equations of the sides
  //Convention, need to explain why it works?? Third point must be on the other side!!
  
  crossArr3(tri_normal, P[0], side_normals[0]);
  crossArr3(P[1], tri_normal, side_normals[1]);
  crossArr3(tri_normal, P[2], side_normals[2]);
  
  side_d[0] = dotArr3(side_normals[0], triang->vertices[0]); 
  side_d[1] = dotArr3(side_normals[1], triang->vertices[0]);
  side_d[2] = dotArr3(side_normals[2], triang->vertices[2]);
  
  int dotprod;
  //Determine whether triangle lies entirely on boundary plane.
  if (triangle_boundary_cube(triang,cube->dim)) {
    arr3 testpt = {1,1,1}; //Point inside the cube
    /*
     * Calculate angle between normal and vector from vertex to point inside cube (testpt)
     * If both point the same way, the angle < 90 degrees, and the dotproduct must be positive.
     */
    dotprod = dotArr3(testpt, tri_normal);
    if (dotprod > tri_d) {
    /*
     * Normal points inward. So the boundary lies "below" the triangle. And therefore we don't need to check below
     */
      res->boundary_below = 1;
      res->acute_below    = 1;
    } else {
    /*
     * Normal points outside the cube. Boundary lies "above". We only need to check "below" the triangle"
     */
      res->boundary_above  = 1;
      res->acute_above     = 1;
    } 
  }
  
  for (size_t i = 0; i < cube->len; i++){
    dotprod = dotArr3(cube_pts[i], tri_normal); //Calculate angle between normal and vector from plane to cube_pt
    
    
    //Check if current point lies on side of the triangle which isn't sharp yet and check if point lies in the
    //prism of possible sharp tetrahedron
    
    if (((dotprod > tri_d && !res->acute_above) ||  //Dotprod > tri_d implies that the point lies "above" the triangle (same side as normal) 
         (dotprod < tri_d && !res->acute_below))     //Dotprod < tri_d implies lies below
          && //Check if point lies in the prism
          dotArr3(cube_pts[i],side_normals[0]) < side_d[0] &&
          dotArr3(cube_pts[i],side_normals[1]) < side_d[1] &&
          dotArr3(cube_pts[i],side_normals[2]) < side_d[2]) 
      {
       if (acute_list && !tetra_facets_list(triang, cube_pts[i], acute_list))
            continue;
      tetra test_tetra;
      memcpy(test_tetra.vertices, &cube_pts[i], sizeof(arr3));
      memcpy(test_tetra.vertices + 1, triang->vertices, 3 * sizeof(arr3)); 
      //test_tetra now holds the four vertices.
      if (tetra_acute(&test_tetra)) {//Yes this one is acute!
        if (dotprod > tri_d) {
          res->acute_above = 1;
        } else {
          res->acute_below  = 1;
        }
        if (res->acute_above && res->acute_below) {
          return 1;
        }
      }
    }
  }
  return 0;  
}



/*
 * Dindex * comb_list = combinations_list(cube.len, 2,&len); //All possible pairs
    cube_points fund = gen_fund_points(dim); 
    index_list = malloc(sizeof(arr3) * len * fund.len);
    c = 0;
    for (size_t i = 0; i < len; i++) {
      triangle cur_triangle;
      Dindex v1,v2;
      v1 = comb_list[i * 2];
      v2 = comb_list[i * 2 + 1];
      memcpy(cur_triangle.vertices, cube.points[v1], sizeof(arr3));
      memcpy(cur_triangle.vertices + 1, cube.points[v2], sizeof(arr3));
      //First two points of cur_triangle now hold the combination index i
      for (int j = 0; j < fund.len; j++)
      {
        memcpy(cur_triangle.vertices + 2, fund.points[j], sizeof(arr3));
        arr3 index;
        triangle_to_index(&cur_triangle, acute_list->dim, index);
        if (EMI(acute_list->t_arr,index) && GMI(acute_list->t_arr,index)) {//Triangle is acute
          memcpy(index_list + c, index, sizeof(arr3));
          c++;
        }
      }
    }
    free(comb_list);

*/


/*
 * Filter list on face2face property. 
 * If sym is true, it only checks for triangles in the fundamental domain (and then filters the 48 symmetries).
 * If sym is false, it checks all triangles that are set in the acute_list
 * BROKEN AS FUCK
 */
void mem_list_face2face_old(tri_mem_list * acute_list, int sym){
  cube_points cube = gen_cube_points(acute_list->dim); 
  
  int *  fund_index = NULL;
  size_t fund_len   = 0;
  //Calculate the indices of the fundamental domain
  if (sym) {
    cube_points fund = gen_fund_points(acute_list->dim[0]);
    fund_len = fund.len;
    fund_index = malloc(sizeof(int) * fund_len);
    for (size_t i = 0; i < fund_len; i++) 
      fund_index[i] = vertex_to_index_cube(fund.points[i], acute_list->dim_mult);
    free(fund.points);
  }
  
  int changed = 1;
  triang_tetra_result tri_result;
  while (changed) {
    changed = 0;
    printf("Start filter loop, amount of acute facets: %zu\n", mem_list_count(acute_list));
    //Loop over all indices
    int vertex_fund = 0; //Indicates whether we have a vertex in the fundamental domain
    for (int i = 0; i < acute_list->dim_size[0]; i++) {
      if (!acute_list->t_arr[i])
        continue;
      if (sym)
        vertex_fund = index_in_array(i, fund_index, fund_len);
      for (int j = 0; j < acute_list->dim_size[1]; j++) { 
        if (acute_list->t_arr[i][j]) { //Third dimension exists!
          if (sym && !vertex_fund)
            vertex_fund = index_in_array(j, fund_index, fund_len);
          if (sym && !vertex_fund) // Only loop over fundamental domain!!
            for (size_t k = 0; k < fund_len; k++) {
              tri_index index = {i,j,fund_index[k]};
              if (!GMI(acute_list->t_arr, index))
                continue; //Triangle no longer in the acute list, no need to check
              triangle cur_tri = triangle_from_index(index, cube.dim);
              if (!triangle_tetra_acute(&cur_tri,&cube, &tri_result,acute_list)) { //This triangle is not acute anymore
                mem_list_clear_sym(acute_list, &cur_tri);
                changed = 1;
              }
            }
          else //Loop over entire domain
            for (int k = 0; k < acute_list->dim_size[2]*8; k++) {
              tri_index index = {i,j,k};
              if (!GMI(acute_list->t_arr, index))
                continue; //Triangle no longer in the acute list, no need to check
              triangle cur_tri = triangle_from_index(index, cube.dim);
              if (!triangle_tetra_acute(&cur_tri,&cube, &tri_result,acute_list)) { //This triangle is not acute anymore
                mem_list_clear_sym(acute_list, &cur_tri);
                changed = 1;
              }
            }
        }
      }
    }
  }
  if (sym)
    free(fund_index);
  free(cube.points);
}

void mem_list_face2face(tri_mem_list * acute_list){
  size_t len = 0;
  arr3 * cube_pts;
  
  
  cube_points fund_points = gen_fund_points(acute_list->dim[0]);
  arr3 *fund_pts = fund_points.points;
  cube_points cube = gen_cube_points(acute_list->dim);
  len = cube.len;
  cube_pts = cube.points;
  Dindex	*arr;
  int changed = 1;
  while (changed) {
    changed = 0;
    printf("Starting face2face loop with %zu acute triangles\n", mem_list_count(acute_list));
    //Check every point in the fundamental domain against all other points 
    arr = revdoor_init(len, 2);  
    if (arr == 0) {
      puts("Error involving revdoor_init");
      exit(1);
    }
    triang_tetra_result tri_result;
    //Loop through all combinations
    do {
      size_t j = arr[0], k = arr[1];
      tri_index indices;
      for (size_t i = 0; i < fund_points.len; i++){
        //Check if the current triangle is in the memory_list "acute"
        if (acute_list->fund) {
          //if (!mem_list_get_fund(acute_list,fund_pts[i],cube_pts[j],cube_pts[k])) //Checks all symmetries
            //continue; 
          vertices_to_index_fund(fund_pts[i],cube_pts[j],cube_pts[k], acute_list, indices); 
        } else
          vertices_to_index(cube_pts[j],cube_pts[k], fund_pts[i], acute_list->dim_mult, indices);
        if (!EMI(acute_list->t_arr, indices) || !GMI(acute_list->t_arr,indices))
          continue;
          
        triangle cur_tri = (triangle) {{{fund_pts[i][0],fund_pts[i][1],fund_pts[i][2]},
                                  {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                                  {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]}}};
        if (!triangle_tetra_acute(&cur_tri,&cube, &tri_result,acute_list)) { //remove from list
          changed = 1;
          if (acute_list->fund)
            mem_list_clear_fund(acute_list, &cur_tri);
          else
            mem_list_clear_sym(acute_list, &cur_tri);

        }
      }
    } while (revdoor_next(arr));
    revdoor_free(arr);
  } 
  free(cube_pts);
  free(fund_pts);
}


tri_mem_list acute_triangles_tetra(arr3 dim) {
  if (dim[0] != dim[1] || dim[1] != dim[2]) {
    printf("Cannot execute acute_triangles_tetra with different dimensions\n");
    exit(0);
  }
  /*
  tri_mem_list result = {NULL, 0, {dim[0],dim[1],dim[2]}}; 
  result.index_list = malloc(1000 * sizeof(tri_index));
  result.len = 1000;
  size_t count = 0;
  */  
  tri_mem_list result = mem_list_init_fund(dim[0],MEM_LIST_FALSE);
  
  cube_points fund = gen_fund_points(dim[0]);
  cube_points cube = gen_cube_points(dim);
  
  Dindex	*arr;
  //Check every point in the fundamental domain against all other points 

  arr = revdoor_init(cube.len, 2);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
  }
  
  triang_tetra_result tri_result;
  do{  
    //Loop through all combinations
    size_t j = arr[0], k = arr[1];  
    for (size_t i = 0; i < fund.len; i++){

      triangle cur_tri = (triangle) {{{fund.points[i][0],fund.points[i][1],fund.points[i][2]},
                                {cube.points[j][0],cube.points[j][1],cube.points[j][2]},
                                {cube.points[k][0],cube.points[k][1],cube.points[k][2]}}};
      if (triangle_tetra_acute(&cur_tri,&cube, &tri_result,NULL)) {
        /*
        if (count == result.len) {
          result.index_list = realloc(result.index_list, (count + count/100) * sizeof(tri_index));
          result.len = (count + count / 100);
        }          
        if (result.index_list == NULL) {
          puts("Error allocating memory!!");
          exit(1);
        }
        triangle_to_index(&cur_tri, dim_mult, result.index_list[count]);
        count++;
        */
        mem_list_set_fund(&result, &cur_tri);
        //mem_list_set_sym(&result,&cur_tri);
      }
    } 
  } while (revdoor_next(arr));
  revdoor_free(arr);
  free(cube.points);
  free(fund.points);
  //result.len = count;
  return result;
}


triangle_list acute_triangles_tetra_old(arr3 dim) {
  triangle_list result = {NULL, 0, {dim[0],dim[1],dim[2]}};
  size_t count = 0;
  ptriangle t_arr = NULL;
  size_t len = 0;
  arr3 *cube_pts;
  
  cube_points cube = gen_cube_points(dim);
  len = cube.len;
  cube_pts = cube.points;
  
  
  Dindex	*arr;
   
  arr = revdoor_init(len, 3);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
    return result;
  }
  triang_tetra_result tri_result;
  
  //Loop through all combinations
  do {
    size_t j = arr[0], k = arr[1], i = arr[2];  
    triangle cur_tri = (triangle) {{{cube_pts[i][0],cube_pts[i][1],cube_pts[i][2]},
                              {cube_pts[j][0],cube_pts[j][1],cube_pts[j][2]},
                              {cube_pts[k][0],cube_pts[k][1],cube_pts[k][2]}}};
    if (triangle_tetra_acute(&cur_tri,&cube, &tri_result, NULL)) {
      
      count++;
      t_arr =  realloc(t_arr, count * sizeof(triangle));
      if (t_arr == NULL) {
        puts("Error allocating memory!!");
        exit(1);
      }
      t_arr[count - 1] = cur_tri;
    }
  } while (revdoor_next(arr));
  revdoor_free(arr);
  result.t_arr = t_arr;
  result.len = count;
  free(cube_pts);
  return result;
}


//For history sake!
int tetra_acute_old(ptetra tet) {
  //First check if it's facets are acute.
  //Old code
  /*
  if (!mat3_triangle_acute(tet->vertices) || //facet 0,1,2
      !mat3_triangle_acute(&tet->vertices[1]) || //facet 1,2,3
      !arr3_triangle_acute(tet->vertices[0],tet->vertices[1],tet->vertices[3]) || //facet 0,1,3
      !arr3_triangle_acute(tet->vertices[0],tet->vertices[2],tet->vertices[3])) //facet 0,2,3
    return 0;
  mat3 P;
  subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  */
  //Optimized "hacky-code". Does the same as above, but then entirely inlined. As this function is called the most, this 
  //means some performance gain!
  
  arr3 P[6]; //All the edges
  
  #define triangle_side(A,B,P) subArr3(tet->vertices[A],tet->vertices[B],P)
  
  
   // Calculate the sides of our tetrahedron check if they are acute
   
   
  //Facet 0,1,2
  triangle_side(1,0,P[0]); //P[0] = v1 - v0
  triangle_side(2,0,P[1]); //P[1] = v2 - v0
  triangle_side(2,1,P[3]); //P[3] = v2 - v1 NOTE P[3], P[2] is saved for something else
  //Convention third argument must be the side vector from v1 to v2
  if (!triangle_sides_acute(P[0],P[1],P[3])) 
    return 0;
    
  //Facet 0,1,3  
  triangle_side(3,0,P[2]); //P[2] = v3 - v0
  triangle_side(3,1,P[4]); //P[4] = v3 - v1
  if (!triangle_sides_acute(P[0], P[2], P[4]))
    return 0;
  
  //Facet 0,2,3
  triangle_side(3,2,P[5]); //P[5] = v3 - v2
  if (!triangle_sides_acute(P[1], P[2], P[5]))
    return 0;
  //Facet 1,2,3
  if (!triangle_sides_acute(P[3], P[4], P[5]))
    return 0;   
  
  
  /*
   * Get first three sides. P already holds V1-V0 and V2-V0. Outward pointing normal
   */
  //subArr3(tet->vertices[1], tet->vertices[0], P[0]);
  //subArr3(tet->vertices[2], tet->vertices[0], P[1]);
  //subArr3(tet->vertices[3], tet->vertices[0], P[2]);
  mat3 Q_cross;
  crossArr3(P[1],P[2], Q_cross[0]);
  crossArr3(P[2],P[0], Q_cross[1]);
  crossArr3(P[0],P[1], Q_cross[2]);
  
  //Normals on three facets, point outward
  //q1 loodrecht op p1, p2
  //als q1 ook loodrecht op p0, dan liggen drie vectors in hetzelfde vlak, degenerate..
  if (dotArr3(P[0], Q_cross[0]) == 0)
    return 0;
  //Calculate inner products <q_i, q_j>
  mat3 Q_cross_trans;
  transMat3(Q_cross, Q_cross_trans);
  //Q_cross has normals in rows, so transpose the right one to get inner product matrix
  mat3 Q_prod;
  dotMat3(Q_cross,Q_cross_trans, Q_prod);
  
  //Check of eerste drie facetten niet-scherpe hoek maken
    if (Q_prod[0][1] >= 0 || Q_prod[0][2] >= 0 || Q_prod[1][2] >= 0)
      return 0;
    
  //Controleer of de andere drie facetten een niet-scherpe hoek maken
  vec3 Q_last = mat3Vector(Q_prod, (vec3) {{-1,-1,-1}});
  return (Q_last.pt[0] < 0 && Q_last.pt[1] < 0 && Q_last.pt[2] < 0);
}


int index_in_array(int index, int * indices, size_t len) {
  for (size_t i = 0; i < len; i++)
    if (index == indices[i])
      return 1;
  return 0;
}
