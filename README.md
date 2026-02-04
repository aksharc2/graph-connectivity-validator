# gcval — Graph Connectivity Validator

`gcval` is a lightweight command-line tool for validating **connectivity (reachability) preservation**
between an original graph **G** and a transformed graph **G\***.

It is intended for **researchers and practitioners working on graph transformations**, including
compression, sparsification, bipartite expansions, reductions, and preprocessing pipelines, where
preserving connectivity among original vertices is a correctness requirement.

---

## What Does gcval Verify?

`gcval` checks whether **undirected connectivity is preserved** among the vertices of the original
graph.

Formally, let `V(G)` be the set of vertices appearing in the original graph. The tool verifies:

```
∀ u, v ∈ V(G):
    u is connected to v in G  ⇔  u is connected to v in G*
```

The transformed graph **may introduce additional vertices**. These are allowed, provided they do
*not* change connectivity relationships among the original vertices.

---

## Key Properties

- Connectivity is treated as **undirected reachability**
- Edge weights (if present) are ignored
- Self-loops are ignored
- Only connectivity among **original vertices** is validated
- Suitable for use in scripts, CI pipelines, and artifact evaluation

---

## Supported Input Formats

### Matrix Market
- `.mtx`, `.mm`
- Coordinate format
- Lines of the form:
  ```
  i j [value]
  ```

### CSV Edge Lists
- `.csv`
- Supported formats:
  ```
  u v
  u,v
  ```
- Extra columns are ignored
- Comment lines starting with `#` or `%` are ignored

---

## Graph Interpretation Modes

### 1. General Graph Mode (default)

- Single global vertex label space
- Suitable for standard undirected graphs

```bash
./gcval original.mtx transformed.mtx
```

or explicitly:

```bash
./gcval --general original.mtx transformed.mtx
```

---

### 2. Bipartite Mode

- Each edge is interpreted as `(left_vertex, right_vertex)`
- Left and right label spaces are treated as **distinct**
- Numeric labels may overlap safely (e.g., `L3` ≠ `R3`)

This mode is useful for:
- Bipartite graph transformations
- Factor graphs
- Left–right expansions and reductions

```bash
./gcval --bipartite original.mtx transformed.mtx
```

---

## Build Instructions

A C11-compatible compiler is required.

```bash
cc -O2 -std=c11 -Wall -Wextra -o gcval gcval.c
```

The tool has no external dependencies.

---

## Usage

```bash
./gcval [--general|-g | --bipartite|-b] original_graph transformed_graph
```

---

## Exit Codes

| Code | Meaning |
|----:|--------|
| 0 | PASS — connectivity preserved |
| 2 | FAIL — connectivity differs |
| 1 | ERROR — parsing or runtime failure |

These exit codes make `gcval` suitable for automated testing and CI workflows.

---

## Example

### Input

Original graph:
```
1 2
2 3
```

Transformed graph:
```
1 4
4 2
2 5
5 3
```

### Output

```
PASS: Connectivity preserved among 3 original vertices.
      (Transformed graph has 5 total mapped vertices.)
```

---

## How It Works (High-Level Overview)

1. Parse edge lists from both input graphs
2. Map arbitrary vertex labels to compact integer IDs
3. Build disjoint-set (Union–Find) structures
4. Compute connected components in both graphs
5. Compare component partitions restricted to original vertices

If and only if the partitions match, connectivity is preserved.

---

## Intended Use Cases

- Validating correctness of graph compression algorithms
- Checking connectivity preservation in graph reductions
- Regression testing graph preprocessing pipelines
- Artifact evaluation and reproducibility checks
- Sanity-checking large-scale graph transformations

---

## Limitations

- Directed reachability is **not supported**
- Strongly connected components are **not supported**
- Only undirected connectivity is checked
- Graphs must fit in memory

Future extensions may include directed and weighted variants.

---

## License

MIT License (or replace with your preferred license)

---

## Citation

If you use `gcval` in academic work, please cite the repository or include a reference in your artifact
evaluation or appendix.