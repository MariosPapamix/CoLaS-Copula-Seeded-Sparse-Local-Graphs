# Decomposing Degree Assortativity in Sparse Spatial Networks — replication code

This repository contains the complete code, data, and frozen results
for the experiments in *Decomposing Degree Assortativity in Sparse
Spatial Networks*: the exact-rational identification and information
certificates, the simulation study, the eleven-city two-platform
analysis, the closure-augmented model class, the competing-model
study, and the collaboration-network two-view analysis.

## Layout

    code/     all analysis scripts, exact coefficient tables, and the
              frozen expected outputs of the three certificates
    data/     raw inputs: Brightkite and Gowalla friendship lists,
              the derived per-user home summaries, and the OpenAlex
              co-authorship dataset (see Data below)
    runs/     frozen outputs of every experiment; every number quoted
              in the paper derives from these files via
              code/make_numbers.py
    run_experiments.sh   reproduction driver (see below)

## Requirements

Python 3.11 with the pinned packages in `requirements.txt`
(`pip install -r requirements.txt`).  The three certificates use only
the standard library (`fractions.Fraction`); NumPy/SciPy/scikit-learn
are needed for the simulations, the city pipeline, and the
competing-model study.

## Quick start (minutes)

    bash run_experiments.sh

This verifies the three exact-rational certificates against their
frozen outputs byte for byte:

1. `colas_certificate_v2.py` — global injectivity and uniform Jacobian
   nondegeneracy of the three-moment map on
   B = [1,5] x [5,15] x [0,0.15], including the boundary psi = 0;
2. `colas_certificate_d2.py` — the planar hard-disc companion;
3. `ieff_face_certificate.py` — a positive uniform lower bound for the
   efficient boundary information on the null face
   [1,5] x [5,15] x {0};

and then regenerates the derived-numbers file from the frozen results
in `runs/`.

## Full reproduction (hours)

    bash run_experiments.sh --full

reruns every experiment from the raw data:

| script | output | content |
|---|---|---|
| `run_experiments.py`, `sim_extended.py`, `sim_spec_test.py` | `runs/sim_*.npz` | coverage, bias, boundary size and power, misspecification, specification-test size and power |
| `brightkite_cities.py`, `gowalla_cities.py` | `runs/*_cities.json` | eleven-city battery: observed moments, configuration nulls, direct and graph estimates, geometry baselines |
| `cities_gof.py` | `runs/cities_r8_*.json` | scaled bootstrap bands, goodness-of-fit residuals, rail flags, edge-length shares, signed counterfactuals |
| `cities_distnull.py`, `sf_extras.py` | `runs/*.json` | distance-matched null, fitted latent-distance competitor, boundary-calibrated two-view test |
| `design_functionals.py` | `runs/design_functionals.json` | first-order signed-amplitude design functionals per city |
| `ieff_profile.py` | `runs/ieff_profile.json` | numerical profile of the efficient information along the null face |
| `rich_fit_city.py` | `runs/rich_*.json` | closure-augmented two-scale class: fits, goodness of fit, counterfactual attribution |
| `competitors_sf.py` | `runs/competitors_sf.json` | degree-corrected blockmodel, spatial logistic with node effects, augmented class, geometry baseline; held-out edge prediction and posterior-predictive checks |
| `oa_analysis.py` (plus `OA_BASE=1`) | `runs/oa_descriptives.json`, `runs/oa_base_NLsub.json` | collaboration network: regional battery, direct estimates with atom-clustered errors, configuration nulls, hard-disc conditional fit |
| `oa_rich.py free\|fixdir\|free2` | `runs/oa_rich_NLsub_*.json` | collaboration network: augmented-class fits (b free, b pinned to the direct estimate, dispersed start), goodness of fit, counterfactuals, signed-b profile |

`code/make_numbers.py` maps the frozen results to every number quoted
in the paper; `code/make_figures.py` regenerates the figures.

## Data

`data/` contains the public Brightkite and Gowalla friendship lists
from the Stanford SNAP repository (Cho, Myers, and Leskovec,
*Friendship and mobility: user movement in location-based social
networks*, KDD 2011; https://snap.stanford.edu/data/loc-brightkite.html
and https://snap.stanford.edu/data/loc-gowalla.html), redistributed
here for reproducibility, together with the derived per-user home
summaries (median check-in coordinate and check-in count for users
with at least five valid check-ins).  The summaries were produced from
the full check-in files, which are not included because of their size,
by `code/summarize_checkins.py`; rerunning that step requires
downloading the check-in files from SNAP.

`data/` also contains the co-authorship dataset built from the
OpenAlex catalog (https://openalex.org) by `code/openalex_fetch.py`:
physics articles published 2018--2022 with at least one author
affiliated to a Netherlands institution, works with more than 25
authors dropped, and each author with at least two corpus works
placed at the coordinates of their most frequently listed geolocated
institution.  `oa_homes.csv.gz` holds one row per placed author
(id, institution latitude and longitude, corpus work count) and
`oa_edges.txt.gz` the co-authorship edges.  The frozen copies are
included because the live catalog updates continuously: rerunning
the fetch reproduces the construction, not necessarily the identical
corpus.

`code/openalex_fetch.py` builds an additional co-authorship dataset
(institution coordinates as positions, productivity as the observable
mark) from the OpenAlex API in the same input format; it is provided
for the extension discussed in the paper and requires direct internet
access to api.openalex.org.

## Notes

Scripts resolve all paths relative to the repository root (override
with the environment variable `COLAS_BASE`), so the repository can be
cloned anywhere.  Randomness is controlled by fixed seeds recorded in
the scripts; the certificate verifications are exact and
deterministic.

## License

The code is released under the MIT License (see `LICENSE`).  The
Brightkite and Gowalla data are due to their original authors; see the
Data section.

## Citation

If you use this code or data, please cite the accompanying paper,
*Decomposing Degree Assortativity in Sparse Spatial Networks*, by
Marios Papamichalis and Regina Ruane.
