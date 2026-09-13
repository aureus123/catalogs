# Cross-Identification by Maximum Likelihood

This folder contains tools that, given two star catalogs, find the matching that *maximizes the joint likelihood* over the bipartite graph of candidate pairs. They are an alternative to the greedy nearest-neighbour pipeline used by the C++ tools (`cross_gc`, `cross_south`, `compare_ppm`, etc.) whose results are stored in [results/cross](../results/cross).

Two approaches are described below. [Section 2](#2-baseline-independent-gaussians-superseded) is the original baseline, which modelled position and magnitude as independent Gaussians; it is retained for context and is no longer implemented. [Section 3](#3-current-approach-evidence-against-a-field-model) is the current method, which follows the Bayesian cross-identification framework of Budavári & Szalay (2008) and Budavári & Basu (2016), extended with the double-star formulation of Severin (2018). [Section 4](#4-tools) documents the two scripts.

For the completed cross-identification of Córdoba Durchmusterung volume I against PPM and GSC — its pipeline, fitted constants, results and machine-readable format — see **[CD_CROSS.md](CD_CROSS.md)**.

## 1. The matching problem

Let *A* and *B* be the two catalogs. The goal is an *injective* partial map $A \to B$: each historical entry gets at most one modern counterpart, and no modern object is claimed twice.

This is the central motivation. Greedy nearest-neighbour assigns the closest counterpart star by star, which can map several stars of *A* onto the same star of *B*. A likelihood formulation instead decomposes the candidate graph into connected components and solves each as a minimum-cost assignment, so exclusivity holds by construction rather than by post-hoc cleanup.

Both approaches share this skeleton:

1. **Edge construction.** A 3-D *k*-d tree on the unit vectors of *B* is queried within chord length $2\sin(R/2)$ for every star of *A*. The angular distance is computed as $\theta = 2\arcsin(\lVert v_a - v_b\rVert/2)$ rather than $\arccos(v_a \cdot v_b)$, which loses precision for nearly identical vectors. Input vectors are re-normalized first, since an 8-decimal CSV representation is not exactly unit-norm — without this, $\arccos(v \cdot v)$ yields some 17″ of phantom self-distance.
2. **Connected components.** Surviving edges induce a bipartite graph on $A \cup B$; components are extracted with `scipy.sparse.csgraph.connected_components` and solved independently. The number of components is comparable to the number of stars and most are tiny, so this is a substantial complexity reduction.
3. **Per-component assignment.** Each component is solved as a minimum-cost assignment.

Working in angular distance rather than separate $(\alpha^*, \delta)$ components is appropriate when no per-coordinate uncertainty is available and the catalog precision is roughly direction-independent. It also avoids the $\cos\delta$ projection inherent in $\alpha^* = \alpha\cos\delta$, which is numerically delicate near the poles: small tabulation errors in right ascension blow up after dividing by a small $\cos\delta$, whereas $\theta$ is well defined and stable everywhere on the sphere.

## 2. Baseline: independent Gaussians (superseded)

> This section describes the first implementation. It is **no longer the code in this folder** — `cross_likelihood.py` was rewritten to the model of section 3. The description and its results are kept because they explain what the current approach was built to fix.

### 2.1 Model

For each candidate pair, the joint observation was modelled by independent Gaussians on the angular separation $\theta(a,b)$ and the magnitude difference $\Delta m = m_a - m_b$. Following Severin (2018) [3] and dropping the doubles/multiples component, the *log-likelihood weight* of an edge was

$$
w(a,b) = \frac{\theta(a,b)^2}{2 \sigma_{pos}^{2}} + \frac{\Delta m^{2}}{2 \sigma_{m}^{2}}
$$

i.e. minus the log of the product of two Gaussians with constants dropped, since they shift every edge equally and do not affect the optimum. The optimal matching *minimizes* the sum of $w(a,b)$ over selected edges.

**Hard cutoffs.** No edge was created beyond an angular threshold $\theta_{max}$, nor between stars whose magnitudes were both known and differed by more than $\Delta m_{max}$. When at least one star had no reported magnitude, the magnitude term was replaced by the constant $\Delta m_{miss}^{2}/(2\sigma_{m}^{2})$, so the edge fell back to a position-only weight with a fixed photometric penalty.

| Constant | Value | Meaning |
|---|---|---|
| $\sigma_{pos}$ | 30″ | Standard deviation for angular distance |
| $\sigma_{m}$ | 0.5 mag | Standard deviation for magnitude difference |
| $\theta_{max}$ | 5′ | Hard cutoff for angular distance |
| $\Delta m_{max}$ | 3.0 mag | Hard cutoff for magnitude difference, when both are known |
| $\Delta m_{miss}$ | 2.0 mag | Constant Δm assumed when a magnitude is missing |

*Table 1. Hardcoded configuration of the superseded baseline.*

Components were solved by brute force when $\min(|A_c|, |B_c|) \le 5$ — at most $5! = 120$ matchings, trivially correct and dominating the workload since most components are very small — and by the Hungarian algorithm (Kuhn–Munkres, $O(n^3)$, via `scipy.optimize.linear_sum_assignment`) otherwise, with a sentinel cost above the maximum admissible weight standing in for absent edges.

Restricting to one-to-one assignments keeps the problem polynomial; the general problem with doubles and multiples is NP-hard, as shown in [3].

### 2.2 Baseline result: GC × PPM

As a sanity check the baseline was run between the *Argentine General Catalog* (GC, 32,379 stars) and *Positions and Proper Motions* (PPM southern stars, 292,998), both at epoch B1875.0. Edge construction produced 43,224 surviving edges across 29,060 non-trivial components; 29,046 were solved by brute force and only 14 needed the Hungarian algorithm, the largest component holding 12 GC and 23 PPM stars. The matching contained **31,426 pairs**.

Against the greedy file [results/cross/cross_gc_ppm.csv](../results/cross/cross_gc_ppm.csv), restricted to PPM-only matches (the greedy tool falls back to GSC for unmatched GC stars), the two agreed on 31,074 entries, or **98.0%** of the greedy-PPM matches.

The differences were dominated by *non-injective* greedy assignments: 512 PPM stars had been assigned to more than one GC star. The likelihood matching cannot do this and routes each disputed PPM to the closest competitor:

| | GC | PPM | dist (″) | mag |
|---|---|---|---|---|
| greedy | GC 10006 | PPM 369604 | 0.95 | 6.7 |
| greedy | GC 10007 | PPM 369604 | 1.94 | 7.0 |
| likelihood | GC 10006 | PPM 369604 | 0.95 | 6.7 |
| likelihood | GC 10007 | (unmatched) | — | — |

*Table 2. A non-injective greedy assignment resolved by the likelihood tool: GC 10006 wins PPM 369604 because it is closer.*

Among the 181 cases where both methods matched but to different counterparts, the likelihood choice was consistently the closer one once ties were broken by mutual consistency within the component.

### 2.3 Why it was replaced

Three limitations motivated the current approach:

- **Uncalibrated scales.** $\sigma_{pos}$ and $\sigma_{m}$ were hardcoded guesses rather than quantities fitted from the catalogs at hand.
- **Gaussian tails.** A Gaussian underestimates the frequency of large errors in nineteenth-century positions, so genuine but badly placed entries fell outside the cutoff.
- **No alternative to matching.** Every edge competed only against other edges, never against the possibility that the true counterpart is absent. A lone weak candidate won simply by being alone, and there was no notion of how crowded its neighbourhood was.

## 3. Current approach: evidence against a field model

The current model scores each candidate association by its **likelihood ratio against a field model**, and admits an explicit *unmatched* hypothesis. All angular quantities are arcseconds and densities are per square arcsecond.

### 3.1 Bayesian evidence for a common source

[Budavári & Szalay (2008), *Probabilistic Cross-Identification of Astronomical Sources*](https://arxiv.org/abs/0707.1611), DOI [10.1086/587156](https://doi.org/10.1086/587156), compare the evidence for a common source against independent sources, marginalizing the unknown true position. Their spherical Fisher model has a circular small-angle approximation, and their framework accommodates selection effects and physical information. This motivates the explicit field-density comparison and the unmatched hypothesis used here.

The departure is deliberate: rather than treating modern formal astrometric errors as a complete description of nineteenth-century errors, the implementation fits **empirical heavy-tailed residuals**. The positional kernel and its search-disk mass are

$$
f(r;s)=\frac{1}{2\pi s^2}\left(1+\frac{r^2}{4s^2}\right)^{-3},
\qquad
F(R;s)=1-\left(1+\frac{R^2}{4s^2}\right)^{-2},
$$

mixed 97% core with 3% broad kernel at broad scale $\max(6s, 60)$ so that large historical errors remain representable. The implemented Student kernel is an empirical extension of the circular approach, not the exact spherical Fisher model of the 2008 paper.

A subtlety worth stating plainly: **a single angular scale does not uniquely specify a probability density.** A half-normal is a density per unit radial distance; a circular positional kernel is a density per unit area, with a different radial distribution. The field comparison here uses the areal measure consistently. Substituting a radial density without its corresponding field measure would introduce a spurious geometric factor.

The local field density is estimated as $(k-1)/(\pi r_k^2)$ with $k = 64$ neighbours. For a singleton the likelihood ratio is

$$
LR_{ab} = \frac{f(r_{ab}) \times \text{photometric factor}}{\rho_a},
$$

where the photometric factor is a normalized Student *t* density with three degrees of freedom divided by a smoothed marginal magnitude density, clipped to $[1/20, 20]$. Missing photometry contributes a factor of one, making it neutral evidence rather than a penalty. This empirical marginal ratio is not a complete luminosity-function or selection model.

Each hypothesis carries incremental cost $-\log(\text{weight}/\text{unmatched weight})$, with the empty hypothesis at cost zero; positive-cost hypotheses are dominated by staying empty. Components are nevertheless built from **all geometric candidate edges before this pruning**.

### 3.2 Global exclusivity and orphan detections

[Budavári & Basu (2016), *Probabilistic Cross-Identification in Crowded Fields as an Assignment Problem*](https://arxiv.org/abs/1609.03065), formulate globally consistent partitions and assignment, including orphan detections. Their two-catalog single-source associations contain at most one detection from each catalog, which corresponds exactly to the single-star exclusivity problem solved here by the Hungarian algorithm with private empty alternatives.

It does *not* directly represent one historical entry associated with two modern stars — that case is section 3.3. Their discussion also motivates retaining ambiguity diagnostics: one optimal assignment does not summarize all nearly optimal alternatives, which is why the tools emit candidate and review sidecars alongside the chosen matching.

### 3.3 Visually double historical entries

[Severin (2018), *Cross-identification of stellar catalogs with multiple stars: Complexity and Resolution*](https://doi.org/10.1016/j.endm.2018.07.005) supplies the motivation for singleton and pair hypotheses and their combinatorial optimization. An unordered pair $\{b,c\}$ has evidence

$$
LR_{a,\{b,c\}}=(LR_{ab}+LR_{ac})\,
\frac{g(d_{bc})}{2\pi d_{bc}\,\rho_a},
$$

where $g$ is a truncated normal density on the separation $d_{bc}$. Both members must lie inside the search disk. The sum allows either member to explain the historical position and photometry; **fluxes are not added**, and the model does not represent blended flux or a photocentre.

A double-marked entry splits its prior between the two readings: singleton weight $Q\,p_{single}\,LR_{ab}$ against pair weight $Q(1-p_{single})\,LR_{pair}$, where $Q$ is the fitted probability that a counterpart is present at all and the unmatched weight is approximated by $1 - QF$.

Two differences from the 2018 work are worth noting. That workflow used PPMX and APASS, solving the 2-matching as an integer program, and its results for zones −22 to −24 are published as [4]; the present one uses PPM and GSC and covers zones −22 to −31, while keeping the record layout of [4]. More importantly, it normalized within each candidate class, whereas the current score contrasts every hypothesis against a field model and an explicit empty alternative — so a sole weak candidate no longer receives high probability merely because it is alone.

Components containing a double-marked entry cannot be solved by the Hungarian algorithm, since the assignment is no longer one-to-one. They are solved as a binary program (SCIP): one hypothesis per historical entry, at most one use of each modern object. This is the NP-hard case identified in [3]; solves are required to reach `optimal` with zero gap, and the extracted assignment is independently checked for integrality, objective consistency and exclusivity.

## 4. Tools

Three Python files:

| File | Role |
|---|---|
| [`cross_likelihood.py`](cross_likelihood.py) | Generic single-star matching |
| [`cross_double_likelihood.py`](cross_double_likelihood.py) | CD single/double matching, export and validation |
| [`likelihood_common.py`](likelihood_common.py) | Shared geometry, readers, calibration, likelihoods, assignment, output |

Dependencies are NumPy, pandas, SciPy and Astropy, plus PySCIPOpt for the double-star optimization.

### 4.1 Single stars

```sh
python cross_likelihood.py A.csv B.csv matches.csv
```

Catalog CSVs require `name,x,y,z,mag`, one row per star, with $(x,y,z)$ a rectangular unit vector and `mag` the magnitude. Coordinates must already refer to compatible frames and epochs — **the CLI does not infer or transform them**. Identifiers must be unique; unit vectors are normalized internally. Missing magnitudes are neutral evidence, and zero is treated as missing unless `--zero-is-valid` is given. Catalog exports in this format are in [cat1875/](cat1875/), at epoch B1875.0.

Defaults: search radius 300″, calibration anchor radius 120″, broad fraction 0.03 with a 60″ scale floor, 64 density neighbours, and $Q$ fitted from the data. Relevant options:

| Option | Effect |
|---|---|
| `--radius`, `--anchor-radius` | Override the search and anchor radii |
| `--q` | Fix the prior probability that a counterpart exists |
| `--no-systematic` | Skip the validation-gated tangent-plane offset |
| `--no-photometry` | Position-only evidence |
| `--model-in` | Reuse a frozen fitted model for reproducible comparisons |
| `--zero-is-valid` | Treat magnitude 0 as a real measurement |

A smooth tangent-plane offset is fitted on training sectors and **accepted only if it improves held-out median residuals**; otherwise the run proceeds with no positional shift.

Output retains the baseline schema `index1,index2,mag,dist`, where `dist` is the angular separation in arcseconds, alongside diagnostic, model, candidate and review sidecars. The `.secure.csv` subset applies operational margin and residual thresholds — its name is not a guarantee of correctness. Local margins are scores, not globally marginalized posterior probabilities.

### 4.2 Double stars

`cross_double_likelihood.py` adds the pair hypotheses of section 3.3 and the SCIP solver. It is currently specialized to the CD workflow — preparation, matching, export and validation — and is documented in [CD_CROSS.md](CD_CROSS.md).

### 4.3 Variable stars

Either tool can widen the photometric evidence for stars known to vary. A variable table is a CSV keyed by the identifiers of the corresponding input catalog:

```csv
name,variability_sigma
PPM 123456,1.2
PPM 234567,
```

*(A schema example, not a list of actual variable stars.)*

```sh
python cross_likelihood.py old.csv ppm.csv matches.csv --variables-b variables_ppm.csv
```

`variability_sigma` is the **per-epoch standard deviation in that catalog's magnitude units**, not peak-to-peak amplitude. Listing a name marks it variable; leaving the scatter blank disables photometric evidence for that star rather than guessing. `--variable-sigma` supplies a common assumed scatter, and `--variables-a` applies the same mechanism to the other catalog. Catalog CSVs may instead carry `variable` and `variability_sigma` columns directly.

Two cautions. Known variables are excluded from calibration anchors, so tagging them also protects the fitted scales. And a table keyed by GCVS names **cannot** be passed directly as a PPM-name table — its entries must first be identified with the input catalog, which is itself a cross-identification problem.

Contributions from the two catalogs are added in quadrature using the fitted magnitude-conversion slope.

## 5. The CD × PPM/GSC catalog

[CD_CROSS.md](CD_CROSS.md) documents the completed cross-identification of the 179,804 entries of Córdoba Durchmusterung volume I (zones −22 to −31) against PPM and GSC: the pipeline and its commands, catalog preparation, the constants fitted for that run, variable-star handling, results, validation, sensitivity controls, limitations, and the machine-readable format of the catalog [cd_ppm_gsc.txt](cd_ppm_gsc.txt).

## References

[1] Budavári, T. & Szalay, A. S. (2008). Probabilistic Cross-Identification of Astronomical Sources. *The Astrophysical Journal*, 679, 301–309. https://doi.org/10.1086/587156 — [arXiv:0707.1611](https://arxiv.org/abs/0707.1611)

[2] Budavári, T. & Basu, A. (2016). Probabilistic Cross-Identification in Crowded Fields as an Assignment Problem. *The Astronomical Journal*, 152, 86. [arXiv:1609.03065](https://arxiv.org/abs/1609.03065)

[3] Severin, D. E. (2018). Cross-identification of stellar catalogs with multiple stars: Complexity and Resolution. *Electronic Notes in Discrete Mathematics*, 69, 29–36. https://doi.org/10.1016/j.endm.2018.07.005

[4] Severin, D. E. (2018). *Cross-identification between Cordoba Durchmusterung catalog (declinations −22, −23 and −24) and PPMX catalog*. Mendeley Data, V1. https://doi.org/10.17632/5wwwtv7c8c.1 — CC BY 4.0. Companion dataset to [3]; its `cat/new_format.txt` defines the fixed-width record layout reused by [CD_CROSS.md](CD_CROSS.md).
