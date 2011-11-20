/* vor.c - compute the edges of a Voronoi Diagram from the Delaunay Triangles */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "math.h"
#include "stdlib.h"
#include "cti_globals.h"

#include <iostream>

using namespace std;

void init_tri_hash_table();
void free_tri_hash_table_space();
extern int outsidecontour(int,int,char);

#define HASH_TABLE_SIZE 401

typedef struct tri_hash_ele
{
    int tr;
    char edge;
    char edgeout;
    int p1,p2;
    struct tri_hash_ele *next;
} TRI_HASH_ELE;

static TRI_HASH_ELE *hash_table[HASH_TABLE_SIZE];


/* hash triangles into triangle hash table. When no match is found, we check
 * whether the edge is outside or inside here. We mark all corresponding npts
 * of outside edges as -1, and set outside field for all triangles containing 
 * outside edge(s). The outside triangles and outside edges will not be used
 * to create Tetras, and we don't need to check those created tetras are inside
 * or not.
 */
void hash_tri(CTI_TRI tri[], int tr, int edge, char cw)
{
    int p1, p2;
    int minp,maxp;
    int index, find_match;
    TRI_HASH_ELE *ptr, *ele;

    p1=tri[tr].p[(edge+1)%3];
    p2=tri[tr].p[(edge+2)%3];

    if (p1<p2)
    {
	minp=p1; maxp=p2;
    }
    else
    {
	minp=p2; maxp=p1;
    }
   
    /* Hash function is here! don't know if it is well behaved. */
    index = (13*minp + 17*maxp) % HASH_TABLE_SIZE;

    ptr = hash_table[index];
    find_match = 0;
    while (ptr != NULL) 	/* there is a collision */
    {
	if (ptr->p1==minp && ptr->p2==maxp)  /* when there is a match */
	{
	    /* we don't put the current triangle into hash table if there is
	     * a match, because only two triangles can share the edge (p1,p2)
	     */
	    tri[tr].t[edge] = ptr->tr;
	    tri[ptr->tr].t[ptr->edge] = tr;
	    /* when the edge is outside, mark triangle outside and npts = 9 */
	    if (ptr->edgeout)
	    {
		tri[tr].outside = 1;
		tri[tr].npts[edge] = 9;
	    }
	    /* else we orient the v-edge */
	    else
	    {
		tri[tr].npts[edge] = 1;
		tri[ptr->tr].npts[ptr->edge] = 3;
	    }
		
	    find_match = 1;
	    break;
	}
	else ptr = ptr->next;
    }
	    
    /* if no match, we put the triangle into hash table, and check the edge
     * to see if it is inside or not.
     */
    if (!find_match)
    {
	ele = (TRI_HASH_ELE *)malloc(sizeof(TRI_HASH_ELE));
	ele->next = hash_table[index];
	hash_table[index] = ele;
	ele->tr = tr;
	ele->edge = edge;
	ele->p1 = minp;
	ele->p2 = maxp;
	ele->edgeout = 0;
	if (outsidecontour(p1,p2,cw))
	{
	    ele->edgeout = 1;
	    tri[tr].outside = 1;
	    tri[tr].npts[edge] = 9;
	}
    }
}

void vor(CTI_TRI tri[], int ntri, char cw)
{
    int tr, edge;

    init_tri_hash_table();

    /* hash triangles using two nodes' number of each edge as key.  if
     * there is a collision in hash table, then check whether their
     * keys are really equal or not. if equal, the neighbor triangle
     * is found. This algorithm runs much faster than O(n^2) if the
     * hash function is well selected.
     */
    for (tr=0; tr < ntri; tr++){
	    for (edge=0; edge<3; edge++){
	        hash_tri(tri, tr, edge, cw);
        }
    }
        

    free_tri_hash_table_space();
}

void init_tri_hash_table()
{
    int i;
    for (i=0; i<HASH_TABLE_SIZE; i++){
	    hash_table[i] = NULL;
    }
}

void free_tri_hash_table_space()
{
    int i;
    TRI_HASH_ELE *ptr, *freeptr;

    for (i=0; i<HASH_TABLE_SIZE; i++) {
	    ptr = hash_table[i];
	    while (ptr != NULL) {
	        freeptr = ptr;
	        ptr = ptr->next;
	        free(freeptr);
	    }
    }
}



