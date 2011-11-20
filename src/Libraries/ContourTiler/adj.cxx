/* calculate adjacency of tetras: which tetras share faces? */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include <stdlib.h>
#include "math.h"
#include "cti_globals.h"

#define HASH_TABLE_SIZE 401


typedef struct tetra_hash_ele {
    int tet;
    int face;
    int p1,p2,p3;
    struct tetra_hash_ele *next;
} TETRA_HASH_ELE;

static TETRA_HASH_ELE *hash_table[HASH_TABLE_SIZE];

void hash_tetra(CTI_TETRA tetra[], int tet, int face);
void init_tetra_hash_table();
void free_tetra_hash_table_space();

void adj(CTI_TETRA tetra[], int ntetra)
{
    int tet, face;

    for (tet=0; tet<ntetra; tet++) {
	tetra[tet].t[0] =
	tetra[tet].t[1] =
	tetra[tet].t[2] =
	tetra[tet].t[3] =
	    -1;
    }

    init_tetra_hash_table();

    /* hash tetras using three nodes' numbers as key. if there is a collision,
     * we check whether the keys are equal or not. if equal, then the neighbour
     * tetra is found. this algorithm runs much faster than O(n^2) if the
     * hash function is well selected.
     */
    for (tet=0; tet<ntetra; tet++) {
	for (face=0; face<4; face++) {
	    hash_tetra(tetra,tet,face);
	}
    }

    free_tetra_hash_table_space();
}

void hash_tetra(CTI_TETRA tetra[], int tet, int face)
{
    int p1, p2, p3;
    int minp,medp,maxp;
    int index, find_match;
    TETRA_HASH_ELE *ptr, *ele;

    p1=tetra[tet].p[(face+1)%4];
    p2=tetra[tet].p[(face+2)%4];
    p3=tetra[tet].p[(face+3)%4];

    /* we don't consider horizontal faces */
    if (savept[p1].z == savept[p2].z && savept[p1].z == savept[p3].z)
	return;

    /* sort nodes' numbers here */
    if (p1<p2) {
	minp=p1; maxp=p2;
    }
    else {
	minp=p2; maxp=p1;
    }

    if (p3>maxp) {
	medp=maxp; maxp=p3;
    }
    else if (p3<minp) {
	medp=minp; minp=p3;
    }
    else 
	medp=p3;

    /* hash function is here! don't know if it's well behaved. */
    index = (minp + 3*medp + 7*maxp) % HASH_TABLE_SIZE;

    ptr = hash_table[index];
    find_match = 0;
    while (ptr != NULL) {
	/* we don't need to put the tetra into hash table if a match is
	 * found, because only two tetras share the face, thus we don't
	 * expect another match later.
	 */
	if (ptr->p1==minp && ptr->p2==medp && ptr->p3==maxp) {
	    tetra[tet].t[face] = ptr->tet;
	    tetra[ptr->tet].t[ptr->face] = tet;
	    find_match = 1;
	    break;
	}
	else ptr = ptr->next;
    }
	    
    if (!find_match) {
	ele = (TETRA_HASH_ELE *)malloc(sizeof(TETRA_HASH_ELE));
	ele->next = hash_table[index];
	hash_table[index] = ele;
	ele->tet = tet;
	ele->face = face;
	ele->p1 = minp;
	ele->p2 = medp;
	ele->p3 = maxp;
    }
}

void init_tetra_hash_table()
{
    int i;
    for (i=0; i<HASH_TABLE_SIZE; i++)
	hash_table[i] = NULL;
}

void free_tetra_hash_table_space()
{
    int i;
    TETRA_HASH_ELE *ptr, *freeptr;

    for (i=0; i<HASH_TABLE_SIZE; i++) {
	ptr = hash_table[i];
	while (ptr != NULL) {
	    freeptr = ptr;
	    ptr = ptr->next;
	    free(freeptr);
	}
    }
}

void sortpts(int p[3])
{
    int t0,t1,t2, p0=p[0],p1=p[1],p2=p[2];

    t0 = MIN(p0,p1);
    t0 = MIN(t0,p2);		/* t0 is the smallest */
    t2 = MAX(p0,p1);
    t2 = MAX(t2,p2);		/* t2 is the largest */
    if (t0 == p0)		/* t1 is what's left */
	t1 = MIN(p1,p2);
    else if (t0 == p1)
	t1 = MIN(p0,p2);
    else
	t1 =  MIN(p0,p1);

    p[0] = t0;
    p[1] = t1;
    p[2] = t2;
}
