<h1 align="center">Cost–information trade-offs partially rationalize<br>the feed-forward loop census</h1>

<p align="center">
  <b>Rafael Velásquez</b> &nbsp;·&nbsp; <b>Juan Manuel Pedraza</b><br>
  <i>Department of Physics, Universidad de los Andes, Bogotá 111711, Colombia</i>
</p>

<p align="center">
  <a href="https://github.com/Rafther0112/FFL_Research/blob/main/paper/FFL_paper.pdf">
    <img src="https://img.shields.io/badge/📄%20Paper-PDF-b31b1b?style=for-the-badge" alt="Read the paper (PDF)">
  </a>
  <a href="https://github.com/Rafther0112/FFL_Research/raw/main/paper/FFL_paper.pdf">
    <img src="https://img.shields.io/badge/⬇%20Download-PDF-4c566a?style=for-the-badge" alt="Download the PDF">
  </a>
  <a href="#citation">
    <img src="https://img.shields.io/badge/📑%20Cite-BibTeX-2e7d32?style=for-the-badge" alt="BibTeX">
  </a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.10%2B-3776ab?logo=python&logoColor=white" alt="Python 3.10+">
  <img src="https://img.shields.io/badge/NumPy-SciPy-013243?logo=numpy&logoColor=white" alt="NumPy / SciPy">
  <img src="https://img.shields.io/badge/Numba-accelerated-00a3e0" alt="Numba">
  <img src="https://img.shields.io/badge/License-MIT-yellow" alt="MIT License">
  <img src="https://img.shields.io/badge/Reproducible-fixed%20seeds-success" alt="Reproducible">
</p>

---

## Overview

Among the thirteen connected three-node directed subgraphs, only the **feed-forward loop (FFL)** is
consistently over-represented in transcriptional networks. Yet its eight structural variants occur
very unequally: in *E. coli* and *S. cerevisiae* the **C1** and **I1** variants dominate and the
other six are rare.

This repository contains the full code for the question we ask in the paper: **can the trade-off
between metabolic cost and transmitted information rationalize that asymmetry?**

Each variant is treated as a noisy channel from the environment $E = \log_{10} K_X$ to the output
expression level $Z$:

- the **channel capacity** is obtained analytically from the Lyapunov equation of the linearized
  chemical Langevin dynamics,
- the **mutual information** $\hat I_{\mathrm{SSA}}(E;Z)$ is estimated independently from Gillespie
  simulations, at steady state and after step inputs,
- the **cost** is the protein burden $B = \langle Y\rangle + \langle Z\rangle$,
- the eight variants are then compared on the **Pareto frontier** of the cost–information plane,
  under four information measures and ten promoter logics.

### Main findings

| | Result |
|---|---|
| **1** | Only **C1, C3, I1, I3** are ever Pareto sole winners within their family. **C2, C4, I2, I4** have exactly zero sole-winner mass under every measure, consistent with their rarity in the census. |
| **2** | **C1 robustly dominates the coherent family** (stationary Pareto shares 30.1%–51.7%, bootstrap intervals above every other coherent variant), in agreement with the census. |
| **3** | **No measure separates I1 from I3** by more than a few percentage points, and the measures disagree on which is ahead — while the census favours I1 over I3 by 21:1 in *S. cerevisiae*. The separating axis may be a measure of **response timing**, which our level-information measures do not probe. |

---

## The eight variants

A feed-forward loop is the directed subgraph $X \to Y$, $X \to Z$, $Y \to Z$. Each edge is an
activation or a repression, giving $2^3 = 8$ variants, split into four **coherent** (C1–C4, direct
and indirect paths reinforce) and four **incoherent** (I1–I4, the paths compete), following the
nomenclature of Mangan & Alon.

| Family | Variants | Dose-response | $I_\mathrm{loc}(E)$ profile |
|---|---|---|---|
| Coherent | C1, C2, C3, C4 | monotonic | single-peaked |
| Incoherent | I1, I2, I3, I4 | band-pass (non-monotonic) | two-peaked, valley where the gain vanishes |

At the $Z$ promoter we implement **all ten non-trivial two-input logics** (AND, OR, and eight more
drawn from a thermodynamic occupancy model of the regulated promoter), not only the two canonical
Boolean limits.

---

## Repository structure

```text
FFL_Research/
├── paper/
│   └── FFL_paper.pdf            # the manuscript (author's version)
├── src/ffl/
│   ├── gates.py                 # Hill functions + the ten logics, Eqs. (2)–(14)
│   ├── model.py                 # the six-equation CLE system, Eq. (1); the 8 variants
│   ├── lyapunov.py              # Jacobian, Lyapunov solver, capacity, Eqs. (15)–(22)
│   ├── gillespie.py             # Numba SSA kernel: steady state and step response
│   ├── information.py           # plug-in + Miller–Madow, fine-grid plug-in, KSG
│   ├── cost.py                  # burden B, Eq. (23), and alternative cost specifications
│   ├── pareto.py                # frontier, 1/k credit rule, shares, gate simplex
│   ├── sampling.py              # log-uniform prior, reweighting, Kish ESS
│   └── stats.py                 # bootstrap CIs, Wilcoxon + Holm, rank-biserial
├── scripts/
│   ├── run_stationary_scan.py   # N = 1000 samples, stationary MI
│   ├── run_capacity_scan.py     # N = 1000 samples, analytical capacity
│   ├── run_transient_scan.py    # N = 200 samples, peak + time-integrated MI
│   ├── run_robustness.py        # cost specs, priors, information thresholds
│   └── make_figures.py          # regenerates every figure and table
├── config/                      # YAML configs for each scan
├── data/                        # generated outputs (git-ignored)
├── figures/                     # generated figures
├── notebooks/                   # exploratory analyses
├── requirements.txt
├── environment.yml
└── LICENSE
```

---

## Installation

```bash
git clone https://github.com/Rafther0112/FFL_Research.git
cd FFL_Research

# with conda
conda env create -f environment.yml
conda activate ffl

# or with pip
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Requirements: Python ≥ 3.10, `numpy`, `scipy`, `numba`, `matplotlib`, `pandas`, `pyyaml`,
`scikit-learn` (KSG estimator), `tqdm`.

---

## Quick start

```python
from ffl.model import FFL
from ffl.lyapunov import channel_capacity
from ffl.cost import protein_burden

loop = FFL(variant="C1", gate="OR", n=2,
           beta_y=20, beta_z=20, K_XY=500, K_XZ=500, K_YZ=1000)

C = channel_capacity(loop)      # bits, Eq. (22)
B = protein_burden(loop)        # ⟨Y⟩ + ⟨Z⟩, Eq. (23)
print(f"C = {C:.2f} bits   B = {B:.0f} proteins")
```

Full scans:

```bash
python scripts/run_stationary_scan.py --config config/stationary.yaml   # N = 1000
python scripts/run_capacity_scan.py   --config config/capacity.yaml     # N = 1000
python scripts/run_transient_scan.py  --config config/transient.yaml    # N = 200, R = 20
python scripts/make_figures.py --all
```

---

## Reproducing the paper

Every figure and table regenerates from the configs in `config/`:

| Output | Content | Command |
|---|---|---|
| Fig. 2 | Analytical dose-response and $I_\mathrm{loc}(E)$ for the eight variants | `make_figures.py --fig 2` |
| Fig. 3 | Stratified Pareto shares under four information measures | `make_figures.py --fig 3` |
| Fig. 4 | Sole-winner mass over the ten-gate simplex | `make_figures.py --fig 4` |
| Fig. 5 | Unstratified eight-variant Pareto shares | `make_figures.py --fig 5` |
| Fig. 6 | Cross-validation: analytical capacity vs. SSA mutual information | `make_figures.py --fig 6` |
| Fig. 7 | Per-gate Pareto shares across all ten logics | `make_figures.py --fig 7` |
| Fig. 8 | Common-support filter: $D_k$ vs. information threshold | `make_figures.py --fig 8` |
| Tables I, IV | Cost-specification robustness; per-subgroup linear fits | `run_robustness.py --tables` |
| Tables II, III, V | Rank correlations; bootstrap CIs and paired tests; census comparison | `make_figures.py --tables` |

### Reproducibility

All numerical experiments use `numpy.random.default_rng` with fixed seeds: **seed 0** for the
parameter scan, seeds $1,\dots,R$ for the Gillespie replicates within each scan point, and
**seed 0** for the bootstrap analyses. The Gillespie kernel is Numba-accelerated; the Lyapunov
equation is solved with `scipy.linalg.solve_continuous_lyapunov`.

### Default parameters

| Parameter | Value / range |
|---|---|
| $\beta_y, \beta_z$ | log-uniform in $[5, 80]$ |
| $K_{XY}, K_{XZ}$ | log-uniform in $[150, 2000]$ |
| $K_{YZ}$ | log-uniform in $[300, 4000]$ |
| Hill coefficient $n$ | $\{2, 4\}$ |
| mRNA degradation $\gamma_i$ | $1$ (time unit) |
| Protein degradation $\mu_i$ | $0.05$ |
| Translation rates $K_{p i}$ | $5$ |
| Capacity integration range | $K_X \in [0.5, 10]$ |
| SSA prior | 8 log-spaced values of $K_X$ in $[1, 8]$ |

Expected runtime: the analytical capacity takes a fraction of a second per parameter sample; the
stationary Gillespie scan ($N = 1000$ samples × 8 variants × 10 gates × 2 Hill coefficients) is the
expensive step and is parallelized over samples.

---

## Validation

The two information estimators are computed at the **same** parameter samples and cross-validated:
Pearson $r = 0.82$ over all 160 000 paired values and $r \in [0.80, 0.95]$ within subgroups. The
within-family ranking by mean information is identical under both estimators in five of the eight
subgroups and differs by a single adjacent transposition in the other three. The within-family
rankings also survive seven alternative cost specifications ($\bar\rho = 0.87$–$1.00$) and two
alternative kinetic-parameter priors ($\rho > 0.85$, Kish ESS ≈ 67%–69%).

---

## Citation

<a id="citation"></a>

If you use this code or build on this work, please cite:

```bibtex
@article{Velasquez2026FFL,
  title   = {Cost-information trade-offs partially rationalize the feed-forward loop census},
  author  = {Vel{\'a}squez, Rafael and Pedraza, Juan Manuel},
  journal = {Physical Review E},
  year    = {2026},
  note    = {Submitted}
}
```

And the software:

```bibtex
@software{Velasquez2026FFLcode,
  title     = {FFL\_Research: cost-information analysis of feed-forward loop variants},
  author    = {Vel{\'a}squez, Rafael and Pedraza, Juan Manuel},
  year      = {2026},
  publisher = {GitHub},
  url       = {https://github.com/Rafther0112/FFL_Research}
}
```

---

## License

Code released under the [MIT License](LICENSE). The manuscript in `paper/` is the authors' own
version and is shared for academic use.

## Contact

Rafael Velásquez — <r.velasquez@uniandes.edu.co><br>
Juan Manuel Pedraza — <jmpedraza@uniandes.edu.co><br>
Systems Biology Group, Universidad de los Andes, Bogotá.

## Acknowledgments

We thank the Systems Biology Group at Universidad de los Andes for useful discussions.
