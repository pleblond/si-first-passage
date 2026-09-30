# Collapse of Binary Graph Embeddings and Metric Reconstruction via SI First-Passage Dynamics

**Author:** Philippe Leblond  
**Draft date:** September 30, 2026  
**Title:** *Collapse of Binary Graph Embeddings and Metric Reconstruction via SI First-Passage Dynamics*  
**Subtitle:** *Weighted audits, conditional recovery, and reproducible MDS calibration*

## Repository layout

This is the complete, standalone **Paper 2** repository package. Extract
`si-first-passage-paper/` and use its contents as the root of a new GitHub
repository. The previously delivered scientific manuscript, figures, tables,
code, and results are unchanged. No Git history, local virtual environment,
or TeX build cache is included.

```text
si-first-passage-paper/
├── README.md
├── manuscript.tex           # Canonical editable TeX
├── manuscript.pdf           # Compiled research draft
├── manuscript.txt           # Plain-text reading copy
├── references.bib           # Portable bibliography
├── build.sh                 # Compile, or recompute and compile
├── requirements.txt         # Pinned Python dependencies
├── code/                    # Benchmarks, figures, verification, designs
├── figures/                 # Vector PDF and PNG figures
├── tables/                  # Generated LaTeX tables
├── results/                 # JSON, CSV, and NumPy outputs
├── editorial/               # Claim audit and source-review ledger
├── NOTICE.md                # Attribution and draft status
├── SHA256SUMS.txt            # File-integrity checksums
├── .gitattributes
└── .gitignore
```

The TeX source is authoritative; the plain-text copy is a reading aid, not a
build input. Unlike Paper 1, no Markdown-to-TeX conversion is required.

## Read and edit

`manuscript.pdf` is the compiled research draft. `manuscript.tex` is the editable canonical source. It uses the first paper's article-scale layout: 11-point Latin Modern, letter paper, one-inch margins, compact headings, numeric references, a running header, and a two-pass pdfLaTeX build. Unlike the first paper's Markdown/Pandoc authoring pipeline, this package uses native TeX directly so that the proofs, cross-references, and mathematical environments remain editable without conversion.

The source archive contains everything needed to compile the document. Downloading the TeX file alone does **not** include its relative `figures/` and `tables/` dependencies. Upload the entire archive to an Overleaf project, or unpack it locally and run:

```sh
bash build.sh
```

Set `manuscript.tex` as the main document when using Overleaf. The build uses the inline numeric bibliography and does not require BibTeX. `references.bib` is included as a portable matching bibliography for later journal conversion; edits to it alone do not update the inline bibliography.

## What is established and what is not

The manuscript proves exact shell-counting identities, the high-rank spectrum of equidistant binary metrics, a quantitative near-equidistance bound, finite-sample metric preservation for shared-field SI estimation, a weak-shortcut recovery theorem on intact finite backbones, and an explicit MDS Gram-matrix perturbation bound. It includes counterexamples to universal rank-two collapse, universal impossibility of unweighted continuum limits, and applying deterministic mean-edge-cost certificates to mean first-passage distances.

The supplied computations show finite-dimensional improvements, not universal smooth-manifold reconstruction. The Delaunay transfer test deliberately retains the failure to remove the radial artifact. The Euclidean-edge control is an **oracle**, not a coordinate-free correction. D11 numerical values in Section 14 are **source-reported historical values**, not newly rerun simulations. Disk curvature is an **inclusion-interface** test, not an isolated flat-disk reconstruction. Static coarsening is **not** a demonstrated dynamical RG fixed point.

Read `editorial/claim_audit.md` before using the original outline's stronger claims. No priority claim is made for first-passage percolation itself or for classical MDS perturbation theory. This is a research draft for author review, not a peer-reviewed validation certificate.

## Reproduce the new computations

Tested environment: Python 3.13.5, with exact library versions in `requirements.txt`. The mathematical computations are self-contained; the external research repository is not imported.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python code/run_benchmarks.py --part all
python code/run_extension.py
python code/make_assets.py
bash build.sh
```

To check the supplied numerical outputs without rerunning the full experiments:

```sh
python code/verify_results.py
```

This checks result-group completeness, triangle inequalities of the saved SI
metrics, coupled shortcut monotonicity, the Gram perturbation bound, and table
availability. It rewrites `results/final_checks.json`; it is a check of the
supplied results, not a full rerun of every experiment.

The combined command is `bash build.sh --recompute`. Resource use is dominated by the 625-node Delaunay SI batch and full eigendecompositions. Runtime depends on the machine. Sparse all-pairs calculations use roughly quadratic storage in the graph size; this implementation is not intended for very large graphs.

Individual primary batches can be run with `--part checks`, `static`, `si`, `delaunay`, `rg`, or `disk`. For a secondary Delaunay-only single-seed rerun:

```sh
ONLY_DELAUNAY=1 GRAPH_SEED=0 python code/run_extension.py
```

This replaces the selected seed's record in an existing `results/delaunay_si.json`, or creates it when absent. Run every configured seed before regenerating the five-seed summary. A normal call without these environment variables regenerates the complete secondary batch.

## Design chronology and estimator conventions

`code/design.json` was created before the primary calculations. `code/extension_design.json` was created after the initial results but before the secondary dimension controls and Delaunay SI transfer calculations. Neither was externally registered. Their hashes lock content, not public timing; the secondary batch is explicitly diagnostic rather than confirmatory.

Each SI trial samples independent **undirected** positive geometric edge passage times, computes an all-pairs shortest-path matrix in that shared field, and averages the matrices. It is not an independently sampled source-row estimator. It has the same single-source SI marginal law under the symmetric persistent-transmission assumptions in the manuscript, but a different cross-source Monte Carlo coupling. No arrival cap or post-hoc triangle repair is used.

Reported square-grid sample standard deviations describe five estimators on the same graph. Delaunay means and error bars describe five graph seeds, with one SI estimator per graph. Pairwise distances are not treated as independent replicates.

## Result map

| File | Contents |
| --- | --- |
| `results/static.json` | New 20 x 20 damage/price diagnostics; full all-source volume audits |
| `results/shells.json` | Shell counts and spectral diagnostics |
| `results/si_baseline.json` | Five 256-field square-grid SI estimators and hop control |
| `results/weak_shortcuts.json` | Coupled 512-field intact-backbone recovery experiment |
| `results/rg.json` | Minimum/mean static coarsening, all three initial long prices |
| `results/dimension_controls.json` | Secondary 1D, 2D, and 3D controls |
| `results/delaunay.json` | Hop and coordinate-length-oracle controls on five graphs |
| `results/delaunay_si.json` | Uniform SI transfer to all five Delaunay graphs |
| `results/disk_curvature.json` | Sampled edges and individual transport curvatures |
| `results/verification.json` | Exact-case guards and the triangle Monte Carlo cross-check |
| `results/aggregates.json` | Values used in summary tables |
| `results/*_matrices.npz`, `weak_shortcuts.npz` | Saved distances, coordinates, and spectral quantities |
| `results/*edges.csv` | Realized version-sensitive damage/coarsening edges |

The JSON and NumPy files are the source of the generated tables and figures. The manuscript's rounded values should be reviewed after changing any configuration or implementation. This build does not automatically rewrite narrative sentences when numerical outputs change.

## Upstream provenance

The first manuscript and its TeX/build materials were reviewed at:

- <https://github.com/pleblond/weighted-graph-paper>

The research source was reviewed on **this specific branch**, not the default branch:

- <https://github.com/pleblond/alltoall-black-holes/tree/cursor/d1-dynamics-u-2031>

Access date: September 30, 2026. The review was targeted; the complete upstream test suite was not executed. No immutable upstream commit SHA was independently resolved. Do not treat a moving branch URL as a version pin. `editorial/source_review.md` identifies the files and their roles.

The new code independently implements the described algorithms. Its graph realizations and random-number coupling differ in the ways stated in the manuscript. Upstream code, manuscripts, and documentation retain their original attribution. No complete copy of either upstream repository is included here.

## Package integrity

`SHA256SUMS.txt` gives content hashes for the delivered package, excluding the manifest itself. Check it from the project directory with:

```sh
sha256sum -c SHA256SUMS.txt
```

Recomputation or editing will intentionally change hashes. The manifest is an integrity aid, not evidence of external preregistration or peer review.

## Licensing

No new license is selected by this packaging step. The author should choose
the applicable manuscript and code licenses before granting reuse permissions.
Existing attributions and the research-draft limitations in `NOTICE.md` and
`editorial/` remain in force.
