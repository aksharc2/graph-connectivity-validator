/*
This script validates Graph Connectivity of an input graph, G(U,V,E) and the transformed/compressed/restructured graph G*(U*, V*, E*). 
The validity check is for Undirected and Unweighted graph.

What is checked?
    For all origianal vertices in G:
        u and v are the same connected component in G
    iff u and v are also the same connected component in G*.

Thus it checks it the connectivity in G* is preserved.

Supported input formats:
    1. Matrix Matket (.mtx /.mm)
    2. CSV edge list

Indexing:
    It handels 0 or 1 based indices

Build:
    cc -02 -std=c11 -Wall -Wextra -o gcval src/gcval.c

Example Usage:
    ./gcval original_graph.mtx transformed_graph.mtx

Exit Codes:
    0 = PASS (connectivity is preserved among origianl vertices)
    2 = FAIL (connectivity differes among the original vertices)
    1 = ERROR (prase / IO / memory)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <errno.h>

// Disjoint Set Union (Union-Find)

typedef struct {
    int* parent;
    int* rank;
    int n;
} DSU;

static void dsu_init(DSU *d, int n){
    d->n = n;
    d->parent = (int*)malloc((size_t)n * sizeof(int));
    d->rank = (int*)calloc((size_t)n, sizeof(int));
    if (!d->parent || !d->rank){
        fprintf(stderr, "ERROR: Out of memory allocating DSU of size %d\n", n);
        exit(1);
    }
    for (int i = 0; i <n; i++){
        d->parent[i] = i;
    }
}