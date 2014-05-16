/*
 * gencombo.c  --  Functions to generate combinations.
 * Copyright (c) 2006 by James Dow Allen
 *
 *  This program can be distributed or modified subject to
 *  the terms of the GNU General Public License (version 2
 *  or, at your option, any later version) as published by
 *  the Free Software Foundation.
 */

#ifndef COMBINATIONS_H
#define COMBINATIONS_H
typedef unsigned short Dindex;
typedef long Combnum;




Dindex *revdoor_init(Dindex n, Dindex m);
void revdoor_free(Dindex *arr);
int revdoor_next(Dindex *arr);
long combination(long n, long r);
Dindex * combinations_list(Dindex n, Dindex m, size_t * len);
#endif

