# An essential convex hull criterion for false confidence in location models

Author: Anonymous. Text CC0-1.0; code MIT (see `LICENSE.md`). Version 0.1.0-candidate, an unrefereed research candidate. Its review was internal and AI-assisted, not human peer review, and no proof assistant was used.

- Version DOI: [10.5281/zenodo.23948952](https://doi.org/10.5281/zenodo.23948952); concept DOI (all versions): [10.5281/zenodo.23948951](https://doi.org/10.5281/zenodo.23948951)
- Release page: <https://evidencepress.org/releases/essential-hull-false-confidence/>

## The question in one paragraph

A satellite operator computes the probability that two objects will collide. With large uncertainty that probability is small even if a collision is certain. Balch, Martin and Ferson (2019) named this *false confidence* and proved that affected hypotheses exist. Martin (2024) asked which ones. This package answers that question exactly for flat-prior posteriors (equivalently, here, fiducial distributions and, in dimension one, confidence distributions) in location models with radially decreasing errors.

## Main results

- **Hull criterion (Theorem 1).** Take the smallest closed convex shape enclosing a hypothesis $H$, ignoring pieces of zero volume. With light-tailed errors, $H$ is calibrated at all thresholds at a parameter exactly when that parameter is not strictly inside the shape.
- **Exact tail condition (Theorem 1, Lemma 0).** The criterion holds if and only if $\liminf_{R\to\infty} S(R)/g(R-b)=0$ for all $b>0$. Exponential-power errors qualify exactly for $p>1$. Tails like $\exp(-r\log r)$ qualify in dimension one and fail in every higher dimension.
- **Two-sided classification (Theorem 4, Corollary 5).** Only halfspaces are calibrated from both sides. This corrects, in these models, Martin's conjecture that all non-linear hypotheses are affected: a ball is not.
- **Exact marginalisation (Theorem 6).** Only strictly monotone functions of one linear form have exact marginal posteriors.
- **High confidence at a fixed level (Theorem 8).** With log-concave errors, if everything outside $H$ fits in a closed convex set with posterior mass at most $s<1/2$ when the data land on the truth, then $H$ gets posterior probability at least $1-s$ at least half the time. Example: with $\sigma/R=10$, a certain collision gives $P_c\leq0.01$ with probability one (Table 2).
- **Rates and heavy tails (Theorem 7, Proposition 9).** The exact logarithmic tail rate at bounded hypotheses, and what happens when $p\leq1$.

## What is prior work, and what is not claimed

The sufficient direction is elementary, and halfspace validity is in Martin (2023). The Gaussian converse combines classical Laplace-transform growth with a tail comparison. All-threshold failure may occur only at low thresholds, so the hull criterion alone says nothing about 95% confidence: $(-2,-1)\cup(1,2)$ in the standard Gaussian model never exceeds posterior mass 0.39. The model assumes known shape and exact matching of sampling and posterior error laws. The prior-art record in `audit/` reports dated searches; it is not proof of priority.

## Files

| File | Contents |
| --- | --- |
| `paper.pdf`, `paper.md` | The paper with all proofs. `code/build.sh` builds the PDF with Pandoc and Tectonic, using `header.tex`. |
| `refs.bib` | Bibliography records of the reference list. |
| `code/verify_examples.py`, `code/verify_exponential_power.py`, `code/verify_extensions.py` | Deterministic floating-point sanity checks. They write `code/*.json`. |
| `code/make_figure.py` | Draws Figure 1, writing to `code/figures/`; the copy used by the paper is in `figures/`. |
| `figures/calibration_cdf.{pdf,png}` | Figure 1. |
| `requirements.txt` | Pinned packages (numpy 2.3.5, scipy 1.17.0, matplotlib 3.10.8). |
| `audit/prior_art_record.md` | Dated prior-art record. |
| `CLAIMS.json` | Each result with its hypotheses and evidence type. |
| `AI_INDEX.md`, `RESEARCH_GATES.json`, `MANIFEST.sha256` | Package index, gate record, checksums. |
| `CITATION.cff`, `.zenodo.json`, `LICENSE.md` | Citation, archive metadata, licences. |
| `.github/workflows/replay.yml` | CI replay (with and without `-O`). |

## Reproduce

Work in a copy, because the scripts write beside themselves.

```bash
python -m pip install -r requirements.txt
cd code
python3 -I -B verify_examples.py
python3 -I -B verify_exponential_power.py
python3 -I -B verify_extensions.py
python3 -I -B make_figure.py
```

Repeat with `-I -B -O`. Each script takes seconds, and the JSON output should match the shipped files up to last-digit floating-point differences across platforms. `make_figure.py` reproduces the PNG; the figure PDF differs only in its embedded creation date. **The theorems rest on the written proofs, not on these finite checks.**

## How to cite

Anonymous (2026). *An essential convex hull criterion for false confidence in location models*. Version 0.1.0-candidate. Evidence Press, unrefereed candidate. https://doi.org/10.5281/zenodo.23948952
