# Source review and provenance

Access date: September 30, 2026. The requested research branch was used throughout. Public manuscript/source text was accessed via GitHub/raw GitHub. This is a targeted source audit; no full-repository clone, complete-suite execution, or immutable upstream commit pin is claimed.

## First paper

Repository: <https://github.com/pleblond/weighted-graph-paper>

Reviewed `manuscript.md`, generated `manuscript.tex`, `header.tex`, `build.sh`, and the editorial revision audit. The new article borrows the typography and proof/counterexample/audit structure. The first paper's distinction among metric contraction, volume inflation, and fitted slope is carried through the SI and MDS analysis.

## Research branch

Repository branch: <https://github.com/pleblond/alltoall-black-holes/tree/cursor/d1-dynamics-u-2031>

| Path | Role in review |
| --- | --- |
| `README.md` | Research scope, repository organization, and attribution |
| `docs/model.md` | Microscopic state and geometric interpretation |
| `docs/DEFERRED.md` | D10/D10b/D11 status, weighted-state direction, preregistration intentions, and queued corrections |
| `src/bh_graph/weighted.py` | Arrival-event comparisons, edge-census logic, and static coarsening |
| `tests/test_weighted.py` | Price-tolerance and coarsening fixtures |
| `src/bh_graph/scrambling.py` | Persistent SI convention, finite cap, independent row estimator |
| `tests/test_mds.py` | Healthy 2D null, SI calibration, Delaunay radial mode |
| `src/bh_graph/graphs.py` | Grid and finite-box Delaunay graph construction |
| `src/bh_graph/update_rule.py` | Degree-preserving shortcut/damage rule and detour-span concept |
| `tests/test_emergent_dim.py` | D11 static tail normalizations and D10b disk-interface geometry |

### Statements requiring special care

The D11 helper `_d11_tail` returns `(E, E*r, E*r*r/100)`; omitting `/100` would introduce a factor-of-100 reporting error. The `ceff_cost_fn` in those tests is a static cost function, not the SI probability map used here.

The D10b disk is the 29 lattice sites in a radius-three disk, densified inside a larger grid. It is not a free-standing reconstructed disk. Its interface negativity is a diagnostic outcome, not a boundary defect claimed to be removed.

The source SI matrix estimates rows independently and symmetrizes afterward. The supplied manuscript's common-field estimator changes that finite-sample coupling. Numerical differences cannot be ascribed solely to a mathematical improvement in the population metric.

The static RG prototype uses coordinate division by two and drops internal block edges. A minimum crossing cost is a static path-choice rule; it is not the exponential parallel-rate sum appropriate to simultaneous SI channels.

## Literature checks

Primary publication/author sources used for the bibliography include Cambridge University Press (Richardson), JMLR (Arias-Castro/Javanmard/Pelletier), the Science publication record (Tenenbaum/de Silva/Langford), JSTOR (Schoenberg), Project Euclid (van den Berg/Kesten), and Ollivier's author manuscript on arXiv. The paper does not claim priority for these established theories.

## Scope of bundled material

The supplied scripts are self-contained independent implementations. Only the new manuscript and outputs are bundled, not either full upstream repository or any upstream font files. The PDF's normal embedded fonts are part of its rendering, not separately distributed font software.
