/*
 * gencombo.c  --  Functions to generate combinations.
 * Copyright (c) 2006 by James Dow Allen
 *
 *  This program can be distributed or modified subject to
 *  the terms of the GNU General Public License (version 2
 *  or, at your option, any later version) as published by
 *  the Free Software Foundation.
 *
 *
 * Two different beautiful algorithms from Knuth are provided.
 * (See Fund. Comp. Prog. 7.2.1.3 preprint by Donald Knuth.)
 *	*  A revolving-door algorithm
 *	*  Random-access to combinations
 * Coding by James Dow Allen (jamesdowallen at gmail)
 *
 * In either case, the routine creates m-sized combinations
 *	of {0, 1, 2, ..., n-1}.  (Often you will use these
 *	as indices into an array whose elements are those
 *	the application wishes to combinate.)
 * Note: `m' must be in range:   0 < m <= n
 *
 * Random-access assigns a canonical number to each combination
 *   and allows non-successive retrieval; it uses Pascal's
 *   Triangle.
 *
 * Revolving-door is used to retrieve combinations in sequence,
 *   with successive results differing in only two positions.
 *   The order in which Revolving-door presents combinations
 *   is NOT the canonical order used by the Random-access method.
 *
 * To demonstrate the routines, compile via
 * 	gcc -Wall -DTESTING -O -o combo_test gencombo.c
 * and get usage instructions via
 *	combo_test
 * The 'B' option uses 'R' (revolving-door) algorithm, but
 * with alternate formatting support.
 *
 * Modify this module to support your needs (e.g. remove Pascal
 *   Triangle generation if you don't need it), but do not
 *   remove license notice or attribution at top.
 */

#include	<stdio.h>
#include	<stdlib.h>
#include  <string.h>
#include "combinations.h"

/*
 * Revolving-Door (Payne's algorithm,
 *	cited as Algorithm R in Knuth 7.2.1.3)
 * Three entry_points:
 *  revdoor_init(n, m) -- malloc's array for combo and sets the
 *	first combo.  For the algorithm's own convenience,
 *	arr[-1] is set to m, and arr[m] set to n.
 *	Returns pointer to the array, or NULL on failure.
 *  revdoor_next(arr, in, out) -- modifies the array passed to it
 *	and returns 1 on success, or 0 if all combinations have
 *	been presented. It also returns, by writing argument pointers
 *	`in' and `out', the values changed (this is suppressed
 *	if `in' is NULL).
 *  revdoor_free() -- free the malloc'ed array
 */

Dindex *revdoor_init(Dindex n, Dindex m)
{
	Dindex	*arr, j;

	if (m < 1 || m > n)
		return 0;
	arr = malloc((m + 2) * sizeof(Dindex));
	if (arr == 0)
		return 0;
	*arr++ = m;
	for (j = 0; j < m; j++)
		*arr++ = j;
	*arr = n;
	return arr - m;
}

void revdoor_free(Dindex *arr)
{
	free((void *)(arr - 1));
}

int revdoor_next(Dindex *arr)
{
	Dindex	j, k, a0, a1;

	k = arr[-1];
	a0 = arr[0];
	a1 = arr[1];
	if (k & 1) {
		if (a0 + 1 < a1) {
			arr[0] = a0 + 1;
			return 1;
		}
		j = 1;
	} else if (a0 > 0) {
		arr[0] = a0 - 1;
		return 1;
	} else if (a1 + 1 < arr[2]) {
		arr[0] = a1;
		arr[1] = a1 + 1;
		return 1;
	} else {
		j = 2;
		a0 = a1;
		a1 = arr[2];
	}
	while (j < k) {
		if (a1 > j) {
			arr[j] = a0;
			arr[j - 1] = j - 1;
			return 1;
		}
		a0 = a1;
		a1 = arr[++j];
		if (a1 + 1 < arr[j + 1]) {
			arr[j - 1] = a1;
			arr[j] = a1 + 1;
			return 1;
		}
		a0 = a1;
		a1 = arr[++j];
	}
	return 0;
}

long combination(long n, long r) {
    if(r==0)
        return 1;
    else {
        long num = n * combination(n - 1, r - 1);
        return num/r;
    }
}

//Returns a list of all the combinations. 
//Indexing happens with [combination * m]
//Note that the list may get VERY big!!
Dindex * combinations_list(Dindex n, Dindex m, size_t * len) {
  *len = combination(n,m);
  Dindex	*arr;
  Dindex * result = malloc(sizeof(Dindex) * m * (*len));
  
  arr = revdoor_init(n, m);  
  if (arr == 0) {
    puts("Error involving revdoor_init");
    exit(1);
  }   
  memcpy(result, arr, m * sizeof(Dindex)); //Store this number at result[i*m]
  
  
  for (size_t i = 1; i < *len; i++) {
    if (!revdoor_next(arr)) {
      printf("Weird error revdoor_next");
      exit(1);
    }
    memcpy(result + i*m, arr, m * sizeof(Dindex)); //Store this number at result[i*m]
  
  }
  revdoor_free(arr);    
  return result;
}
