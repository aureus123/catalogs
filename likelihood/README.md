# Cross-Identification by Maximum Likelihood

This folder contains a cross-identification tool that, given two star catalogs, finds the matching that *maximizes the joint likelihood* over the bipartite graph of candidate pairs. It is an alternative to the greedy nearest-neighbor pipeline used by the C++ tools (`cross_gc`, `cross_south`, `compare_ppm`, etc.) whose results are stored in [results/cross](../results/cross).

The motivation is to enforce a true one-to-one matching: greedy assigns the closest counterpart star-by-star, which can map several stars from catalog *A* onto the same star in catalog *B*. The likelihood approach, instead, decomposes the problem into connected components and solves each as a minimum-cost assignment, producing an injective matching by construction.

## 1. Model

Let *A* and *B* be the two catalogs. For each candidate pair (*a*, *b*) we model the joint observation by independent Gaussians on the angular separation $\theta(a,b)$ and the magnitude difference $\Delta m = m_a - m_b$. Following my 2018 work [1] and dropping the doubles/multiples component (we only consider single-star matches here), the *log-likelihood weight* of an edge is

$$
w(a,b) = \frac{\theta(a,b)^2}{2 \sigma_{pos}^{2}} + \frac{\Delta m^{2}}{2 \sigma_{m}^{2}}
$$

i.e. minus the log of the product of two Gaussians, with constants dropped (they shift every edge equally and do not affect the optimum). The optimal matching is the one that *minimizes* the sum of $w(a,b)$ over selected edges.

The treatment differs from [1] in two ways:

- The angular distance is treated as a single isotropic parameter, not as separate $(\alpha^*, \delta)$ Gaussians. This is appropriate when no per-coordinate uncertainty is available and the catalog precision is roughly direction-independent. As a side benefit, the formulation avoids the $\cos\delta$ projection inherent in $\alpha^* = \alpha\cos\delta$, which is numerically delicate near the celestial poles: tiny tabulation errors in right ascension can blow up after dividing by a small $\cos\delta$, whereas the angular distance $\theta$ is well-defined and stable everywhere on the sphere.
- Doubles/multiples are not modelled. The bipartite graph allows only one-to-one assignments. This restriction has the additional benefit that the resulting minimum-cost assignment is solvable in polynomial time (Hungarian algorithm, $O(n^3)$), whereas the more general problem with doubles/multiples is NP-Hard, as shown in [1].

**Hard cutoffs.** No edge is created beyond an angular threshold $\theta_{max}$, nor between stars whose magnitudes are both known and differ by more than $\Delta m_{max}$. When at least one of the two stars has no reported magnitude, the magnitude term is replaced by the constant $\Delta m_{miss}^{2}/(2\sigma_{m}^{2})$, so the edge falls back to a position-only weight with a fixed photometric penalty.

**Settings used.** The current implementation uses the following hardcoded constants (see [cross_likelihood.py](cross_likelihood.py)):

| Constant | Value | Meaning |
|---|---|---|
| $\sigma_{pos}$ | 30″ | Standard deviation for angular distance |
| $\sigma_{m}$ | 0.5 mag | Standard deviation for magnitude difference |
| $\theta_{max}$ | 5′ | Hard cutoff for angular distance |
| $\Delta m_{max}$ | 3.0 mag | Hard cutoff for magnitude difference (only when both magnitudes are known) |
| $\Delta m_{miss}$ | 2.0 mag | Constant Δm assumed when at least one magnitude is missing |

*Table 1. Configuration used by the cross-identification tool.*

## 2. Algorithm

The procedure is split in three phases.

### 2.1 Edge construction

Both catalogs are read from CSV files in [cat1875/](cat1875/), where each row stores the star identifier together with its rectangular unit-vector coordinates $(x, y, z)$ at epoch B1875.0 and its magnitude.

Candidate pairs are generated with a 3-D *k*-d tree on the unit-vector coordinates of *B*, querying within the chord length $2\sin(\theta_{max}/2)$ for every star of *A*. The angular distance for each candidate is then computed via the numerically stable $\theta = 2\arcsin(\lVert v_a - v_b\rVert/2)$ rather than $\arccos(v_a \cdot v_b)$, which loses precision for nearly identical vectors. After re-normalizing the input vectors (whose 8-decimal CSV representation is not exactly unit-norm), this formula gives sub-arcsecond accuracy down to chord lengths below $10^{-4}$.

For each candidate that survives the angular and magnitude cutoffs, the edge weight $w(a,b)$ is recorded.

### 2.2 Connected components

The surviving edges induce a bipartite graph on $A \cup B$. Connected components are extracted with `scipy.sparse.csgraph.connected_components` so each component can be solved independently — a substantial complexity reduction, since the number of components is comparable to the number of stars and components are generally tiny.

### 2.3 Per-component assignment

For each non-trivial component, a *minimum-cost assignment* is computed:

- **Brute force** when $\min(|A_c|, |B_c|) \le 5$: all matchings are enumerated (at most $5! = 120$). Trivially correct, dominates the workload because most components are very small.
- **Hungarian algorithm** (Kuhn-Munkres, $O(n^3)$) for larger components, via `scipy.optimize.linear_sum_assignment`. The cost matrix is padded with a sentinel cost (above the maximum admissible weight) for absent edges so unmatched stars are explicit.

The result is an injective partial matching $A \to B$.

### 2.4 Output

The output CSV follows the schema used in [results/cross](../results/cross):

`index1, index2, mag, dist`

where `mag` is the magnitude of `index1` (omitted when missing) and `dist` is the angular separation in arcseconds. Cross-identification files generated by this tool live in [cross/](cross).

## 3. Results: GC × PPM

As a sanity check the tool was run between the *Argentine General Catalog* (GC, 32379 stars) and the *Positions and Proper Motions* catalog (PPM southern stars, 292998 stars), both at epoch B1875.0. The construction phase produced 43224 surviving edges, distributed across 29060 non-trivial connected components. Of these, 29046 were solved by brute force and only 14 required the Hungarian algorithm; the largest component had 12 stars from GC and 23 from PPM. The matching contains **31426 pairs**.

Compared against the greedy-matching file [results/cross/cross_gc_ppm.csv](../results/cross/cross_gc_ppm.csv) (after restricting it to PPM-only matches; the greedy tool falls back to GSC for unmatched GC stars), the two matchings agree on 31074 entries, i.e. **98.0%** of the greedy-PPM matches.

The remaining differences are dominated by *non-injective* assignments in the greedy file: 512 PPM stars were assigned by the greedy tool to more than one GC star. The likelihood matching cannot do this and instead routes each disputed PPM to the closest GC competitor. A typical case is

| | GC | PPM | dist (″) | mag |
|---|---|---|---|---|
| greedy | GC 10006 | PPM 369604 | 0.95 | 6.7 |
| greedy | GC 10007 | PPM 369604 | 1.94 | 7.0 |
| likelihood | GC 10006 | PPM 369604 | 0.95 | 6.7 |
| likelihood | GC 10007 | (unmatched) | — | — |

*Table 2. Example of a non-injective greedy assignment resolved by the likelihood tool: GC 10006 wins PPM 369604 because it is closer.*

This pattern accounts for the bulk of the 280-row gap between the two files. Among the 181 cases where both methods produce a match but to a different counterpart, the likelihood choice is consistently the closer one when ties are broken by mutual consistency in the surrounding component.

## 4. Limitations and notes

The position uncertainty $\sigma_{pos}$ used here should be calibrated before performing cross-identifications. Same for magnitude to a lesser extent.

The treatment of missing magnitudes uses a fixed-penalty surrogate. Stars where this surrogate "wins" against a real magnitude difference of comparable size will be matched by position alone.

Doubles and multiples are not modelled. If a star in *A* corresponds to two unresolved components in *B*, only one of them will be picked; the other will be reported as unmatched (or matched to a different star of *A* if the geometry allows it).

## References

[1] Severin, D. E. (2018). Cross-identification of stellar catalogs with multiple stars: Complexity and Resolution. *Electronic Notes in Discrete Mathematics*, Vol. 69, pp. 29–36. https://doi.org/10.1016/j.endm.2018.07.005
