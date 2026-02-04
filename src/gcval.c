/*
Connectivity validator between original graph G and transformed graph G*.

Checks (undirected/unweighted connectivity):
For original vertices (those present in G),
u connected to v in G  <=>  u connected to v in G*

Supported formats:
  - Matrix Market (.mtx/.mm) coordinate format
  - CSV edge list (.csv) with "u v" or "u,v"

Indexing:
  No special handling required: all labels are mapped to compact IDs [0..n-1].

Modes:
  - General (default): one global label space.
  - Bipartite (--bipartite/-b): each edge line is interpreted as (left_label, right_label),
    and left/right may reuse the same numeric labels (Lx != Rx).

Build:
  cc -O2 -std=c11 -Wall -Wextra -o gcval gcval.c

Usage:
  ./gcval [--general|-g | --bipartite|-b] original_graph.mtx transformed_graph.mtx

Exit codes:
  0 PASS
  2 FAIL
  1 ERROR
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <errno.h>

#define EDGEVEC_INIT_CAP 1024

/* ==================== DSU (Union-Find) =================== */

typedef struct {
    int *parent;
    int *rank;
    int n;
} DSU;

static void dsu_init(DSU *d, int n) {
    d->n = n;
    d->parent = (int*)malloc((size_t)n * sizeof(int));
    d->rank   = (int*)calloc((size_t)n, sizeof(int));
    if (!d->parent || !d->rank) {
        fprintf(stderr, "ERROR: Out of memory allocating DSU of size %d\n", n);
        exit(1);
    }
    for (int i = 0; i < n; i++) d->parent[i] = i;
}

static int dsu_find(DSU *d, int x) {
    int root = x;
    while (d->parent[root] != root) root = d->parent[root];
    while (d->parent[x] != x) {
        int nx = d->parent[x];
        d->parent[x] = root;
        x = nx;
    }
    return root;
}

static void dsu_union(DSU *d, int a, int b) {
    int ra = dsu_find(d, a);
    int rb = dsu_find(d, b);
    if (ra == rb) return;

    if (d->rank[ra] < d->rank[rb]) {
        d->parent[ra] = rb;
    } else if (d->rank[ra] > d->rank[rb]) {
        d->parent[rb] = ra;
    } else {
        d->parent[rb] = ra;
        d->rank[ra]++;
    }
}

static void dsu_free(DSU *d) {
    free(d->parent);
    free(d->rank);
    d->parent = NULL;
    d->rank = NULL;
    d->n = 0;
}

/* ==================== Edge storage =================== */

typedef struct { int64_t u, v; } Edge64;

typedef struct {
    Edge64 *a;
    size_t n;
    size_t cap;
} EdgeVec;

static void ev_init(EdgeVec *ev) {
    ev->a = NULL;
    ev->n = 0;
    ev->cap = 0;
}

static void ev_push(EdgeVec *ev, int64_t u, int64_t v) {
    if (ev->n == ev->cap) {
        size_t ncap = ev->cap ? ev->cap * 2 : EDGEVEC_INIT_CAP;
        Edge64 *na = (Edge64*)realloc(ev->a, ncap * sizeof(Edge64));
        if (!na) {
            fprintf(stderr, "ERROR: Out of memory storing edges\n");
            exit(1);
        }
        ev->a = na;
        ev->cap = ncap;
    }
    ev->a[ev->n++] = (Edge64){u, v};
}

static void ev_free(EdgeVec *ev) {
    free(ev->a);
    ev->a = NULL;
    ev->n = 0;
    ev->cap = 0;
}

/* ==================== Helpers / parsing =================== */

static int ends_with(const char *s, const char *suffix) {
    size_t ls = strlen(s), lf = strlen(suffix);
    if (lf > ls) return 0;
    return strcmp(s + (ls - lf), suffix) == 0;
}

static int is_blank_or_comment(const char *line) {
    const char *p = line;
    while (*p && isspace((unsigned char)*p)) p++;
    return (*p == '\0' || *p == '\n' || *p == '\r' || *p == '#' || *p == '%');
}

/* parse "u v" or "u,v" (ignores trailing columns) */
static int parse_two_int64(const char *line, int64_t *a, int64_t *b) {
    const char *p = line;
    while (*p && isspace((unsigned char)*p)) p++;
    if (!*p || *p == '#' || *p == '%') return 0;

    errno = 0;
    char *end1 = NULL;
    long long x = strtoll(p, &end1, 10);
    if (end1 == p || errno) return 0;

    p = end1;
    while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;

    errno = 0;
    char *end2 = NULL;
    long long y = strtoll(p, &end2, 10);
    if (end2 == p || errno) return 0;

    *a = (int64_t)x;
    *b = (int64_t)y;
    return 1;
}

/* ==================== Readers =================== */

static void read_csv_edges(const char *path, EdgeVec *edges) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open %s\n", path);
        exit(1);
    }

    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        if (is_blank_or_comment(line)) continue;
        int64_t u, v;
        if (parse_two_int64(line, &u, &v)) {
            if (u != v) ev_push(edges, u, v); /* drop self-loops */
        }
    }
    fclose(fp);
}

static void read_mtx_edges(const char *path, EdgeVec *edges) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open %s\n", path);
        exit(1);
    }

    char line[4096];

    /* Find size line: nrows ncols nnz (or at least nrows ncols) */
    int got_size = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '%') continue;
        if (is_blank_or_comment(line)) continue;

        long long r=0, c=0, nz=0;
        if (sscanf(line, "%lld %lld %lld", &r, &c, &nz) >= 2) {
            got_size = 1;
            break;
        }
    }
    if (!got_size) {
        fprintf(stderr, "ERROR: MTX parse error (missing size) in %s\n", path);
        fclose(fp);
        exit(1);
    }

    /* Read coordinate entries: "i j [value]" */
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '%') continue;
        if (is_blank_or_comment(line)) continue;

        long long i=0, j=0;
        if (sscanf(line, "%lld %lld", &i, &j) >= 2) {
            if (i != j) ev_push(edges, (int64_t)i, (int64_t)j);
        }
    }

    fclose(fp);
}

static void load_edges_auto(const char *path, EdgeVec *edges) {
    if (ends_with(path, ".mtx") || ends_with(path, ".mm")) {
        read_mtx_edges(path, edges);
    } else if (ends_with(path, ".csv")) {
        read_csv_edges(path, edges);
    } else {
        fprintf(stderr, "ERROR: Unknown file format for %s (expected .mtx/.mm/.csv)\n", path);
        exit(1);
    }
}

/* ==================== int64 -> int mapping =================== */

/* Store hash keys as uint64_t so shifting/tagging is well-defined. */
typedef struct {
    uint64_t *keys;
    int      *vals;
    unsigned char *used;
    size_t cap;  /* power of two */
    size_t n;
} Map64to32;

/* MurmurHash3 finalizer-style mixer */
static uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

static void map_init(Map64to32 *m, size_t cap_pow2) {
    m->cap = cap_pow2;
    m->n = 0;
    m->keys = (uint64_t*)malloc(m->cap * sizeof(uint64_t));
    m->vals = (int*)malloc(m->cap * sizeof(int));
    m->used = (unsigned char*)calloc(m->cap, 1);
    if (!m->keys || !m->vals || !m->used) {
        fprintf(stderr, "ERROR: Out of memory allocating label map\n");
        exit(1);
    }
}

static void map_free(Map64to32 *m) {
    free(m->keys);
    free(m->vals);
    free(m->used);
    memset(m, 0, sizeof(*m));
}

static void map_rehash(Map64to32 *m) {
    size_t oldcap = m->cap;
    uint64_t *oldk = m->keys;
    int *oldv = m->vals;
    unsigned char *oldu = m->used;

    map_init(m, oldcap * 2);

    for (size_t i = 0; i < oldcap; i++) {
        if (!oldu[i]) continue;

        size_t mask = m->cap - 1;
        size_t idx = (size_t)mix64(oldk[i]) & mask;
        while (m->used[idx]) idx = (idx + 1) & mask;

        m->used[idx] = 1;
        m->keys[idx] = oldk[i];
        m->vals[idx] = oldv[i];
        m->n++;
    }

    free(oldk);
    free(oldv);
    free(oldu);
}

static int map_get_or_put(Map64to32 *m, uint64_t key, int *next_id, int allow_new) {
    if (allow_new && (m->n + 1) * 10 >= m->cap * 7) map_rehash(m); /* load factor > 0.7 */

    size_t mask = m->cap - 1;
    size_t idx = (size_t)mix64(key) & mask;

    while (m->used[idx]) {
        if (m->keys[idx] == key) return m->vals[idx];
        idx = (idx + 1) & mask;
    }

    if (!allow_new) return -1;

    m->used[idx] = 1;
    m->keys[idx] = key;
    m->vals[idx] = (*next_id)++;
    m->n++;
    return m->vals[idx];
}

static int map_get(const Map64to32 *m, uint64_t key) {
    size_t mask = m->cap - 1;
    size_t idx = (size_t)mix64(key) & mask;

    while (m->used[idx]) {
        if (m->keys[idx] == key) return m->vals[idx];
        idx = (idx + 1) & mask;
    }
    return -1;
}

/* Key encoding:
   - general:   key = bits of label (uint64_t cast)
   - bipartite: key = (label << 1) | sidebit   (sidebit: 0 left, 1 right)
*/
typedef enum { MODE_GENERAL = 0, MODE_BIPARTITE = 1 } GraphMode;

static inline uint64_t key_general(int64_t label) {
    return (uint64_t)label;
}

static inline uint64_t key_bipartite(int64_t label, uint64_t sidebit) {
    return (((uint64_t)label) << 1) | (sidebit & 1ULL);
}

/* Build/extend one map (Design B). */
static int build_map(const EdgeVec *raw, Map64to32 *map, int next_id, int do_init, GraphMode mode) {
    if (do_init) {
        size_t need = (raw->n ? raw->n * 4 : 1024);
        size_t cap = 1;
        while (cap < need) cap <<= 1;
        map_init(map, cap);
    }

    for (size_t i = 0; i < raw->n; i++) {
        uint64_t ku, kv;
        if (mode == MODE_BIPARTITE) {
            /* interpret each edge as (left, right) */
            ku = key_bipartite(raw->a[i].u, 0);
            kv = key_bipartite(raw->a[i].v, 1);
        } else {
            ku = key_general(raw->a[i].u);
            kv = key_general(raw->a[i].v);
        }

        (void)map_get_or_put(map, ku, &next_id, 1);
        (void)map_get_or_put(map, kv, &next_id, 1);
    }
    return next_id;
}

static void union_all_edges(DSU *d, const EdgeVec *raw, const Map64to32 *map, GraphMode mode) {
    for (size_t i = 0; i < raw->n; i++) {
        uint64_t ku, kv;
        if (mode == MODE_BIPARTITE) {
            ku = key_bipartite(raw->a[i].u, 0);
            kv = key_bipartite(raw->a[i].v, 1);
        } else {
            ku = key_general(raw->a[i].u);
            kv = key_general(raw->a[i].v);
        }

        int u = map_get(map, ku);
        int v = map_get(map, kv);
        if (u < 0 || v < 0) {
            fprintf(stderr, "ERROR: Mapping failure for edge (%lld, %lld)\n",
                    (long long)raw->a[i].u, (long long)raw->a[i].v);
            exit(1);
        }
        dsu_union(d, u, v);
    }
}

/* ==================== Partition comparison (robust) =================== */

typedef struct { int root; int v; } RootV;

static int cmp_rootv(const void *a, const void *b) {
    const RootV *x = (const RootV*)a;
    const RootV *y = (const RootV*)b;
    if (x->root != y->root) return (x->root < y->root) ? -1 : 1;
    return (x->v < y->v) ? -1 : (x->v > y->v);
}

/* Exact comparison of connected-component partition on vertices [0..nA-1]. */
static int partitions_equal_exact(DSU *dA, DSU *dB, int nA) {
    RootV *a = (RootV*)malloc((size_t)nA * sizeof(RootV));
    RootV *b = (RootV*)malloc((size_t)nA * sizeof(RootV));
    if (!a || !b) { fprintf(stderr, "ERROR: OOM comparing partitions\n"); exit(1); }

    for (int i = 0; i < nA; i++) {
        a[i].root = dsu_find(dA, i);
        a[i].v = i;
        b[i].root = dsu_find(dB, i);
        b[i].v = i;
    }

    qsort(a, (size_t)nA, sizeof(RootV), cmp_rootv);
    qsort(b, (size_t)nA, sizeof(RootV), cmp_rootv);

    int ia = 0, ib = 0;
    while (ia < nA && ib < nA) {
        int ra = a[ia].root;
        int rb = b[ib].root;

        int sa = ia;
        while (ia < nA && a[ia].root == ra) ia++;
        int ea = ia;

        int sb = ib;
        while (ib < nA && b[ib].root == rb) ib++;
        int eb = ib;

        int lena = ea - sa;
        int lenb = eb - sb;
        if (lena != lenb) { free(a); free(b); return 0; }

        for (int k = 0; k < lena; k++) {
            if (a[sa + k].v != b[sb + k].v) { free(a); free(b); return 0; }
        }
    }

    free(a); free(b);
    return (ia == nA && ib == nA);
}

/* ==================== Main =================== */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [--general|-g | --bipartite|-b] <original.(mtx|mm|csv)> <transformed.(mtx|mm|csv)>\n"
        "Checks connectivity equivalence on the original vertex set.\n"
        "Modes:\n"
        "  (default) general\n"
        "  --general|-g   treat as general graph\n"
        "  --bipartite|-b treat as bipartite (each edge is left,right; labels may overlap)\n",
        prog);
}

int main(int argc, char **argv) {
    GraphMode mode = MODE_GENERAL;
    int argi = 1;

    if (argc >= 2) {
        if (strcmp(argv[1], "--bipartite") == 0 || strcmp(argv[1], "-b") == 0) {
            mode = MODE_BIPARTITE;
            argi++;
        } else if (strcmp(argv[1], "--general") == 0 || strcmp(argv[1], "-g") == 0) {
            mode = MODE_GENERAL;
            argi++;
        }
    }

    if (argc - argi != 2) { usage(argv[0]); return 1; }

    const char *pathA = argv[argi];
    const char *pathB = argv[argi + 1];

    EdgeVec Araw, Braw;
    ev_init(&Araw);
    ev_init(&Braw);

    load_edges_auto(pathA, &Araw);
    load_edges_auto(pathB, &Braw);

    if (Araw.n == 0) {
        fprintf(stderr, "ERROR: No edges parsed from original graph (%s). Format mismatch?\n", pathA);
        ev_free(&Araw);
        ev_free(&Braw);
        return 1;
    }

    Map64to32 map;
    memset(&map, 0, sizeof(map));

    int next_id = 0;
    next_id = build_map(&Araw, &map, next_id, 1, mode);
    int nA = next_id;

    next_id = build_map(&Braw, &map, next_id, 0, mode);
    int nB = next_id;

    DSU dA, dB;
    dsu_init(&dA, nA);
    dsu_init(&dB, nB);

    union_all_edges(&dA, &Araw, &map, mode);
    union_all_edges(&dB, &Braw, &map, mode);

    int ok = partitions_equal_exact(&dA, &dB, nA);

    if (ok) {
        printf("PASS: Connectivity preserved among %d original vertices.\n", nA);
        printf("      (Transformed graph has %d total mapped vertices.)\n", nB);
    } else {
        printf("FAIL: Connectivity differs among original vertices.\n");
    }

    dsu_free(&dA);
    dsu_free(&dB);
    map_free(&map);
    ev_free(&Araw);
    ev_free(&Braw);

    return ok ? 0 : 2;
}
