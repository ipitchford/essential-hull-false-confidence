# AI index: essential convex hull criterion for false confidence (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the paper [paper.pdf](paper.pdf) (source [paper.md](paper.md)). Structured hypotheses are in [CLAIMS.json](CLAIMS.json).
Version DOI 10.5281/zenodo.23948952; concept DOI 10.5281/zenodo.23948951.

## Notation

Model (M): X = theta + AZ, Q_x = law of x + AZ', A invertible, Z with density f(z) = g(|z|), g positive, bounded, nonincreasing. q_H(x) = Q_x(H); F(t) = P_theta(q_H(X) <= t). Valid at theta (V): F(t) >= t for all t in (0,1). C_H = closed convex hull of the essential support of H. S(R) = P(|Z| >= R). (T): liminf S(R)/g(R-b) = 0 for every b > 0; (T+): the limit is 0. f_(p,d) proportional to exp(-|z|^p/p).

## Claims and where they are proved

| ID | Claim | Location in [paper.md](paper.md) | Evidence type |
|---|---|---|---|
| L0 | Lemma 0: tail asymptotics; f_(p,d) satisfies (T+) iff p > 1 and fails (T) otherwise; exp(-r log r) tails satisfy (T+) in d=1 and fail (T) for d>=2 | Section 2 | proof |
| T1 | Theorem 1: theta outside int C_H gives validity; under (T) theta inside gives failure; (T) is necessary for the equivalence | Section 3 | proof |
| C2, C3 | Corollaries 2, 3: validity at all false parameters iff int C_H is inside H; for open or regular closed H, iff convex | Section 4 | proof |
| T4 | Theorem 4: four-way classification; Exact iff almost a halfspace through theta | Section 4 | proof |
| C5 | Corollary 5: both sides valid only for almost-halfspaces (corrected form of Martin's conjecture) | Section 4 | proof |
| T6 | Theorem 6: exact marginal posteriors iff phi = psi(a^T u) a.e., psi strictly increasing | Section 4 | proof |
| T7 | Theorem 7: logarithmic tail rate with coefficient the signed hull distance, p > 1, bounded essential support | Section 5 | proof; finite tables |
| T8 | Theorem 8: Q_theta(C) <= s < 1/2 and H^c inside closed convex C give P(q_H >= 1-s) >= 1/2, under log-concavity | Section 6 | proof; tested on Table 2 |
| P9 | Proposition 9: tail trichotomy for bounded hypotheses | Section 7 | proof; constant checked |

## Verification artefacts

- [code/verify_examples.py](code/verify_examples.py) -> [code/gaussian_false_confidence_checks.json](code/gaussian_false_confidence_checks.json): 25 Gaussian cases for (-2,-1) union (1,2).
- [code/verify_exponential_power.py](code/verify_exponential_power.py) -> [code/exponential_power_checks.json](code/exponential_power_checks.json): 24 cases at p = 1.5 and p = 3.
- [code/verify_extensions.py](code/verify_extensions.py) -> [code/extension_checks.json](code/extension_checks.json): crossover 0.2865, K of Section 7, tail-ratio direction, Table 2 with Theorem 8 tested.
- [code/make_figure.py](code/make_figure.py) -> [figures/calibration_cdf.png](figures/calibration_cdf.png).
- Pinned packages: [requirements.txt](requirements.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml). PDF build: [code/build.sh](code/build.sh).
- Scope: floating-point sanity checks of formulas, not proofs. Failures exit non-zero rather than rely on assert, so they stay active under python -O. verify_examples.py records values for comparison with the shipped JSON and has no internal pass/fail test.

## Sources and novelty record

[audit/prior_art_record.md](audit/prior_art_record.md) and [RESEARCH_GATES.json](RESEARCH_GATES.json) record dated searches. They report what was found and are not proof of priority. Bibliography: [refs.bib](refs.bib). Martin (2023) is PMLR with a URL and no DOI.

## Limits

Known shape; unrestricted Euclidean parameter; exact matching of sampling and posterior error laws; all-threshold statements may concern only low thresholds, so use Theorem 8 for high confidence; unknown scale, informative priors and approximate posteriors are not covered. Status: unrefereed candidate; internal AI-assisted review only; no proof-assistant verification or independent reproduction.
