---
title: "An essential convex hull criterion for false confidence in location models"
author: "Anonymous"
date: "4 November 2026"
abstract: |
  Suppose a statistician reports the posterior probability of a fixed hypothesis about a location parameter. When can that probability be trusted not to favour the hypothesis more often than a calibrated probability would, when the hypothesis is false? We answer this exactly for flat-prior posteriors in location models with radially decreasing errors. For exponential-power errors with density proportional to $\exp(-\|z\|^p/p)$, $p>1$, including the Gaussian case $p=2$, the posterior probability is stochastically no larger than a uniform variable at every threshold exactly when the true parameter lies outside the interior of the closed convex hull of the essential support of the hypothesis. More generally, this hull criterion holds for a radial density $g$ with tail probability $S$ if and only if $\liminf_{R\to\infty} S(R)/g(R-b)=0$ for every $b>0$. This covers $p>1$, excludes every $0<p\leq1$, and reveals a dimension effect: tails like $\exp(-r\log r)$ satisfy the criterion in dimension one and fail it in every higher dimension. Applying the theorem to a hypothesis and its complement gives a four-way calibration classification, and shows that both can be free of false confidence only when the hypothesis is almost everywhere a halfspace. This corrects a conjecture of Martin (2024) in these models and yields a sharp marginalisation statement: the posterior of a real function is an exact confidence distribution at every parameter if and only if its sublevel sets are almost everywhere parallel halfspaces. Because the classification concerns all thresholds, we add a finite-level theorem for high confidence. With log-concave errors, if a closed convex set contains the complement of the hypothesis and has posterior probability at most $s<1/2$ when the observation equals the truth, then the hypothesis receives posterior probability at least $1-s$ with sampling probability at least one half. A worked satellite conjunction example shows what each result does and does not say.
---

*Status.* Version 0.1.0-candidate (4 November 2026), archived at <https://doi.org/10.5281/zenodo.23948952> (all versions: <https://doi.org/10.5281/zenodo.23948951>). AI-assisted research paper, not externally refereed; internal AI-assisted review only. Finite checks are supplied as scripts; no proof-assistant verification has been done.

# 1. Introduction

**The problem in plain terms.** A satellite operator computes the probability that two objects will collide, given a noisy estimate of their closest approach. When the noise is large the computed collision probability is small, even if a collision is certain. Balch, Martin and Ferson (2019) called the general phenomenon *false confidence*: for any additive probability produced from data, some false hypotheses tend to receive high probability. "No collision" is the standard example. Their theorem says that such hypotheses exist. It does not say which hypotheses they are.

**The question.** Martin (2024) asked which hypotheses are affected. He proved a sufficient condition in the Gaussian location model, showed that complements of convex sets are roughly the affected class, and conjectured that all non-linear hypotheses have at least a mild case of false confidence. The review by Martin (2026) still treats a characterisation as open (its §2.4 and its Open Question in §7).

**What this paper shows.** We fix a hypothesis $H$ and a parameter $\theta$, and ask whether the posterior probability $q_H(X)$ is stochastically no larger than a uniform variable. In words: at each threshold, $q_H(X)$ exceeds that threshold no more often than a perfectly calibrated probability would. This is the *all-threshold validity* condition of Martin (2024). The paper has four main results.

1. **The hull criterion (Theorem 1, §3).** Take the smallest closed convex shape that encloses $H$, ignoring pieces of zero volume. Validity at $\theta$ fails exactly when $\theta$ lies strictly inside that shape, provided the errors have sufficiently light tails. The tail condition is itself necessary and sufficient. This shape is the closed convex hull of the essential support, defined in §2.
2. **Two-sided classification (Theorem 4, §4).** Combining $H$ with its complement gives four calibration regimes. A hypothesis and its complement are both globally valid only if the hypothesis is almost everywhere a halfspace (Corollary 5). This replaces "non-linear" in Martin's conjecture by the correct geometric condition. It also yields an exact marginalisation theorem (Theorem 6).
3. **A high-confidence theorem (Theorem 8, §6).** The hull criterion can be triggered only at low thresholds, so on its own it says little about 95% confidence. Theorem 8 gives a finite, checkable condition under which a false hypothesis receives high posterior probability at least half the time.
4. **Sharp boundaries (§§5 and 7).** Theorem 7 gives the exact logarithmic tail rate at bounded hypotheses. Proposition 9 shows what happens when $p\leq1$. At $p=1$ in dimension one, hull position no longer decides anything.

**A warning about language.** In the field, "false confidence" usually suggests *high* posterior probability for a false hypothesis. The formal definition quantifies over all thresholds, and the two ideas differ. For example, $H=(-2,-1)\cup(1,2)$ in the standard Gaussian model has false confidence in the formal sense at $\theta=0$, yet its posterior probability never exceeds $0.39$. We therefore use *all-threshold validity* for the theorems, and say *high-confidence failure* when we mean the practical concern.

**Roadmap.** §2 sets up the model. §3 proves the hull criterion and the tail characterisation. §4 treats complements and marginalisation. §5 gives the tail rate. §6 proves the finite-level theorem and applies it to conjunction analysis. §7 covers heavy tails. §8 places the work in the literature and states its limits. §9 records the reproducible checks, and §10 concludes.

# 2. Setting and definitions

Let $d\geq1$ and let $A$ be an invertible real $d\times d$ matrix. Let $f(z)=g(\|z\|)$ be a probability density on $\mathbb R^d$, where $g:[0,\infty)\to(0,\infty)$ is bounded and nonincreasing. We assume
$$
X=\theta+AZ,\qquad Q_x=\mathcal L(x+AZ'),
\tag{M}
$$
where $Z,Z'$ have density $f$. Then $Q_x$ is the flat-prior posterior for the location likelihood: symmetry turns the likelihood in $x-u$ into a density in $u-x$. In these models it is also the usual fiducial distribution, and in dimension one a confidence distribution. So the results apply to those objects as well.

The *exponential-power* family is
$$
f_{p,d}(z)=c_{p,d}\exp(-\|z\|^p/p),\qquad c_{p,d}^{-1}=|\mathbb S^{d-1}|\,p^{d/p-1}\Gamma(d/p),\qquad p>0.
$$
When $p=2$, (M) is $X\sim N_d(\theta,\Sigma)$, $Q_x=N_d(x,\Sigma)$, with $\Sigma=AA^{T}$. For other $p$, $AA^T$ is a shape matrix and need not be the covariance.

For a Lebesgue-measurable $H\subseteq\mathbb R^d$ write
$$
q_H(x)=Q_x(H),\qquad F_{H,\theta}(t)=P_\theta\{q_H(X)\leq t\}.
$$
The hypothesis is *false at $\theta$* when $\theta\notin H$. We say $H$ is *valid at $\theta$* if
$$
F_{H,\theta}(t)\geq t\quad\text{for every }t\in(0,1).
\tag{V}
$$
Equivalently, $P_\theta\{q_H(X)>1-\alpha\}\leq\alpha$ for every $\alpha\in(0,1)$. So (V) says $q_H(X)$ is stochastically no larger than $\mathrm{Unif}(0,1)$. Failure of (V) at a false parameter is *false confidence* in the sense of Martin (2024). His definition uses $\geq1-\alpha$ rather than $>1-\alpha$. The choice does not matter here: the sufficient bounds below are uniform pointwise bounds, and the failures in Theorems 1 and 4 are proved for the smaller strict event. Theorem 8 is stated with $\geq$, and it implies failure of the strict form at every level $1-s'$ with $s<s'<1/2$.

Define the essential support and its closed convex hull:
$$
E_H=\{y:\lambda_d(H\cap B(y,r))>0\text{ for every }r>0\},\qquad C_H=\overline{\mathrm{conv}}(E_H).
$$
The set $E_H$ is closed, and $H\setminus E_H$ is null. Changing $H$ on a null set changes neither $q_H$ nor $E_H$. If $\lambda_d(H)>0$ then $C_H$ has nonempty interior, because a positive-measure set cannot lie almost everywhere in a hyperplane.

Let $S(R)=P(\|Z\|\geq R)$. The tail conditions are
$$
\liminf_{R\to\infty}\frac{S(R)}{g(R-b)}=0\ \text{ for every } b>0,
\tag{T}
$$
and the stronger form (T$^+$), in which the $\liminf$ is a limit. Decreasing $b$ decreases $g(R-b)$, so only small $b$ matter.

**Lemma 0 (tails).** (a) For $f_{p,d}$, as $R\to\infty$,
$$
S(R)\sim c_{p,d}|\mathbb S^{d-1}|R^{d-p}e^{-R^p/p}.
\tag{1}
$$
(b) $f_{p,d}$ satisfies (T$^+$) for every $p>1$ and fails (T) for every $0<p\leq1$.
(c) Let $g(r)\propto\exp\{-(1+r)\log(1+r)\}$. Then (T$^+$) holds in dimension $d=1$, and (T) fails in every dimension $d\geq2$.

*Proof.* (a) After $s=r^p/p$ the radial integral is a constant times $\int_{R^p/p}^\infty s^{d/p-1}e^{-s}\,ds$, and l'Hôpital's rule gives (1).

(b) By (1), $S(R)/g(R-b)\sim|\mathbb S^{d-1}|R^{d-p}\exp[-\{R^p-(R-b)^p\}/p]$. For $p>1$ the exponent is $-bR^{p-1}(1+o(1))$, and this beats any power of $R$. For $p=1$ the ratio is asymptotic to $|\mathbb S^{d-1}|R^{d-1}e^{-b}$, which stays positive. For $p<1$ the exponent tends to $0$ and the ratio tends to infinity.

(c) Write $\psi(r)=(1+r)\log(1+r)$, so $\psi'(r)=1+\log(1+r)$. L'Hôpital gives $S(R)\sim\kappa R^{d-1}e^{-\psi(R)}/\psi'(R)$, where $g=\kappa_0e^{-\psi}$ and $\kappa=\kappa_0|\mathbb S^{d-1}|$. The mean value theorem gives $\psi(R)-\psi(R-b)=b\log R+b+o(1)$. Hence $S(R)/g(R-b)$ is asymptotic to a constant times $R^{d-1-b}/\log R$. In dimension one this tends to zero for all $b>0$. For $d\geq2$ it diverges for every $b<d-1$. $\square$

# 3. The hull criterion

This section answers the basic question for a single hypothesis. Separation gives validity outside the hull. Inside the hull, finitely many pieces of $H$ lie ahead of every large observation in some direction, and that keeps $q_H$ too large in the tails.

**Theorem 1.** Assume (M) and $\lambda_d(H)>0$.

(a) If $\theta\notin\mathrm{int}\,C_H$ then $H$ is valid at $\theta$.

(b) If (T) holds and $\theta\in\mathrm{int}\,C_H$, then $\liminf_{t\downarrow0}F_{H,\theta}(t)/t=0$. So (V) fails, and if $\theta\notin H$ the hypothesis has false confidence at $\theta$. Under (T$^+$), $F_{H,\theta}(t)=o(t)$.

(c) If (T) fails, then there is a bounded nonconvex open $H$ that is valid at every parameter of $C_H$ and at every false parameter, although $\mathrm{int}\,C_H$ contains false parameters. So (T) is necessary as well as sufficient for the equivalence
$$
H\text{ valid at }\theta\iff\theta\notin\mathrm{int}\,C_H
$$
to hold for all $H$ and $\theta$.

In particular, by Lemma 0, the equivalence holds for exponential-power errors exactly when $p>1$. The classification is the same for every invertible $A$. The theorem does not require $H$ to be bounded, convex, open or smooth. If $\lambda_d(H)=0$ then $q_H\equiv0$, and (V) holds trivially.

*Proof of (a).* Here $\theta\notin\mathrm{int}\,C_H$, so a supporting or separating hyperplane gives a nonzero vector $a$ with $a^T(u-\theta)\leq0$ on $C_H$. Since $H\subseteq C_H$ almost everywhere,
$$
q_H(x)\leq Q_x\{u:a^T(u-\theta)\leq0\}=G_a\bigl(a^T(\theta-x)\bigr),
\tag{2}
$$
where $G_a$ is the distribution function of $a^TAZ$. This variable has a positive density and is symmetric about zero, so $G_a$ is continuous and strictly increasing. Under $P_\theta$, $a^T(\theta-X)=-a^TAZ$ also has distribution function $G_a$. The right side of (2) is therefore uniform on $(0,1)$. $\square$

*Proof of (b).* Applying $u\mapsto A^{-1}(u-\theta)$ preserves essential supports, hulls and interiors, so we may take $\theta=0$ and $A=I$. Choose $c>0$ with $B(0,4c)\subseteq C_H$. For each unit vector $v$, the support function of $C_H$ at $v$ is at least $4c$, so some $u_v\in E_H$ has $v^Tu_v>3c$. By continuity and compactness there are finitely many $u_1,\dots,u_m\in E_H$ with
$$
\max_{i}v^Tu_i>2c\quad\text{for every unit }v.
\tag{3}
$$
Each $A_i=H\cap B(u_i,c)$ has positive finite measure. Put $M=\max_i(\|u_i\|+c)$ and $K=\min_i\lambda_d(A_i)>0$.

Take $x=rv$ and pick $i$ from (3). Every $u\in A_i$ has $v^Tu>c$, so
$$
\|rv-u\|^2=r^2-2rv^Tu+\|u\|^2\leq(r-c)^2+M^2.
$$
Hence $\|rv-u\|\leq r-c/2$ for all $r\geq r_0$, where $r_0$ depends only on $c$ and $M$. Since $g$ is nonincreasing, with $b=c/2$,
$$
q_H(x)\geq K\,g(\|x\|-b)\qquad(\|x\|\geq r_0).
\tag{4}
$$
This lower bound uses only finitely many positive-measure pieces of $H$, so unbounded hypotheses cause no difficulty.

The function $q_H$ is continuous, by translation continuity in $L^1$, and strictly positive, so it has a positive minimum on each compact ball. Note $g(R)\to0$ because $f$ is integrable and $g$ is monotone. Given a large $R$, put $t=\tfrac12Kg(R-b)$. Once $R$ is large, $t$ lies below the minimum of $q_H$ on $B(0,r_0)$. If $r_0\leq\|x\|<R$, then (4) gives $q_H(x)\geq Kg(R-b)>t$. So $q_H(x)\leq t$ forces $\|x\|\geq R$, and
$$
\frac{F_{H,0}(t)}{t}\leq\frac{2}{K}\,\frac{S(R)}{g(R-b)}.
\tag{5}
$$
By (T) the right side tends to zero along a sequence $R_k\to\infty$, which gives thresholds $t_k\downarrow0$ with $F_{H,0}(t_k)/t_k\to0$.

Now assume (T$^+$), and take any small $t$. Let $R_t=\inf\{R:Kg(R-b)\leq t\}$, so $R_t\to\infty$ as $t\downarrow0$. If $r_0\leq\|x\|<R_t$ then (4) gives $q_H(x)\geq Kg(\|x\|-b)>t$. So, as before, $q_H(x)\leq t$ forces $\|x\|\geq R_t$. Since $R_t+b/2>R_t$ and $g$ is nonincreasing, $t\geq Kg(R_t-b/2)$. Therefore
$$
\frac{F_{H,0}(t)}{t}\leq\frac1K\,\frac{S(R_t)}{g(R_t-b/2)}\longrightarrow0,
$$
by (T$^+$) with shift $b/2$. $\square$

*Proof of (c).* Suppose (T) fails: there are $b,\kappa>0$ and $R_0\geq b$ with $S(R)\geq\kappa g(R-b)$ for $R\geq R_0$. Take $A=I$ first. Let $e_1$ be the first coordinate vector and, for $0<\varepsilon\leq b/8$,
$$
H=B(\tfrac b4e_1,\varepsilon)\cup B(-\tfrac b4e_1,\varepsilon).
$$
If $\theta\in C_H$ and $u\in H$ then $\|u-\theta\|\leq b$, so $\|x-u\|\geq r-b$ when $r=\|x-\theta\|\geq b$. For $r\geq R_0$,
$$
q_H(x)\leq\lambda_d(H)\,g(r-b)\leq\frac{\lambda_d(H)}{\kappa}S(r).
$$
For $r<R_0$, $q_H(x)\leq\lambda_d(H)g(0)\leq\lambda_d(H)g(0)S(r)/S(R_0)$. Choose $\varepsilon$ so small that both constants are at most one. Then $q_H(x)\leq S(\|x-\theta\|)$ for every $x$. Under $P_\theta$, $S(\|X-\theta\|)$ is uniform by the probability integral transform, so (V) holds at each $\theta\in C_H$. Outside $C_H$ it holds by (a). The set $H$ is open and nonconvex, and $0\notin H$ lies in the interior of its hull. Invertible linear images handle general $A$. $\square$

# 4. Complements, marginalisation and Martin's conjecture

Theorem 1 controls low values of $q_H$. Applying it to $H^c$ controls high values. Together they classify the whole distribution of $q_H(X)$. Assumption (T) is in force, so all of this applies to exponential-power errors with $p>1$.

**Corollary 2.** Assume (T). A positive-measure hypothesis $H$ satisfies (V) at every $\theta\notin H$ if and only if
$$
\mathrm{int}\,C_H\subseteq H.
\tag{6}
$$

**Corollary 3.** Assume (T). A nonempty open hypothesis is valid at every false parameter if and only if it is convex. The same holds for a nonempty regular closed set, meaning $H=\overline{\mathrm{int}\,H}$.

For open $H$, $E_H=\overline H$ and $H\subseteq\mathrm{int}\,C_H$, so (6) says $H=\mathrm{int}\,C_H$, which is convex. Conversely an open convex set is the interior of its closure. For regular closed $H$, $E_H=H$, and (6) with closure gives $C_H\subseteq H$.

The regularity restriction matters. A ball together with an isolated point is nonconvex but has the same essential hull and posterior probabilities as the ball. Conversely, removing one interior point from a ball leaves posterior probabilities unchanged, yet the hypothesis becomes false at that point, and (b) of Theorem 1 applies there. The theorem uses the essential hull but retains actual set membership as the definition of truth.

**Theorem 4 (four-way classification).** Assume (T), and that $H$ and $H^c$ both have positive measure. At each $\theta$ exactly one row applies, with $F=F_{H,\theta}$.

| Name | $\theta$ and $\mathrm{int}\,C_H$ | $\theta$ and $\mathrm{int}\,C_{H^c}$ | Behaviour of $q_H(X)$ |
|:---|:---|:---|:---|
| Conservative | outside | inside | $F(t)>t$ for all $0<t<1$ |
| Anti-conservative | inside | outside | $F(t)<t$ for all $0<t<1$ |
| Crossing | inside | inside | $F(t)<t$ for arbitrarily small $t$, and $F(1-s)>1-s$ for arbitrarily small $s$ |
| Exact | outside | outside | $q_H(X)$ is uniform |

The Exact row holds precisely when $H$ agrees almost everywhere with a halfspace whose boundary passes through $\theta$. Under (T$^+$), in the Crossing row, $F(t)=o(t)$ and $P_\theta\{q_H(X)\geq1-s\}=o(s)$.

*Proof.* Since $q_{H^c}=1-q_H$, Theorem 1 for the complement shows that $q_H(X)$ stochastically dominates a uniform variable if and only if $\theta\notin\mathrm{int}\,C_{H^c}$.

Crossing. Theorem 1(b) applied to $H^c$ gives arbitrarily small $s$ with $P_\theta\{q_H(X)\geq1-s\}=P_\theta\{q_{H^c}(X)\leq s\}<s$, so $F(1-s)>1-s$ for those $s$, with no atomlessness assumption.

Exact. The two comparisons force uniformity. Conversely, if $q_H(X)$ is uniform, Theorem 1 gives a halfspace $L$ through $\theta$ with $H\subseteq L$ almost everywhere, and $q_H\leq q_L$ with $q_L(X)$ uniform. Equal expectations force equality almost surely. If $L\setminus H$ had positive measure, the positive posterior density would give $q_H<q_L$ everywhere. So $H=L$ almost everywhere.

Strictness in the first row. Take $L$ from (2). It cannot agree with $H$ almost everywhere, since that would be the Exact row. So $q_H<q_L$ everywhere, and both are continuous. For $t\in(0,1)$ choose $x_0$ with $q_L(x_0)=t$. Continuity gives an open set, on the side $q_L>t$, where $q_H<t$. It has positive sampling probability, and $\{q_L\leq t\}\subseteq\{q_H\leq t\}$, so $F(t)>t$. For the second row, take a closed halfspace $L$ through $\theta$ containing $H^c$ almost everywhere, and put $w=q_{L^c}$. Then $w(X)$ is uniform, and $L$ does not agree with $H^c$ almost everywhere, so $q_H>w$ everywhere. Thus $\{w\geq t\}\subseteq\{q_H>t\}$. Continuity gives an open set on the side $w<t$ where $q_H>t$, and it has positive probability. Hence $P_\theta\{q_H(X)>t\}>1-t$, that is, $F(t)<t$. $\square$

The Exact row is a global form of a classical fact: in symmetric location models, a flat-prior posterior probability of a one-sided hypothesis equals a one-sided p-value. (Casella and Berger 1987 obtain the p-value as the infimum of posterior probabilities over classes of symmetric priors, under symmetry and monotone likelihood ratio, reached in the diffuse limit.) Theorem 4 says that this exact agreement occurs for no other sets.

**Corollary 5 (corrected form of Martin's conjecture).** Assume (T), and that $H,H^c$ have positive measure. If $H$ and $H^c$ are each valid at every parameter where they are false, then $H$ agrees almost everywhere with a halfspace.

*Proof.* By (6), the open convex sets $\mathrm{int}\,C_H\subseteq H$ and $\mathrm{int}\,C_{H^c}\subseteq H^c$ are disjoint and nonempty. A hyperplane separates them, and it places $H$ and $H^c$ on opposite closed sides almost everywhere. $\square$

Martin (2024, §3) conjectured that "all non-linear hypotheses" have at least a mild case of false confidence. Taken literally in the sense of (V), this is false: the hypothesis $\{u:\|u\|<r\}$ is non-linear, and by Theorem 1(a) it is valid at every false parameter. In the other direction, failure does not need a curved boundary: $(-2,-1)\cup(1,2)$ has flat boundary points only, and it fails at $0$. Martin himself suggested that non-linearity alone may produce only a milder effect. Corollary 5 states the repaired version: what cannot be avoided is failure *on one of the two sides* of a non-halfspace partition.

**Theorem 6 (exact marginalisation).** Assume (T). Let $\phi:\mathbb R^d\to\mathbb R$ be measurable. Say the posterior of $\phi$ is *exact at $\theta$* if $Q_X\{u:\phi(u)\leq\phi(\theta)\}$ is uniform under $P_\theta$. Then the posterior of $\phi$ is exact at every $\theta$ if and only if there is a unit vector $a$ such that, for every $\theta$,
$$
\{u:\phi(u)\leq\phi(\theta)\}=\{u:a^Tu\leq a^T\theta\}\quad\text{almost everywhere.}
$$
In that case $\phi(u)=\psi(a^Tu)$ for almost every $u$, with $\psi$ strictly increasing.

*Proof.* The "if" direction is (2) with equality. For "only if", uniformity rules out null or co-null sublevel sets, so the Exact row applies to each $H_\theta=\{\phi\leq\phi(\theta)\}$, and $H_\theta$ is almost everywhere a closed halfspace $L_\theta$ with $\theta$ on its boundary. Sublevel sets are nested. Two halfspaces nested up to null sets have the same outward normal, so all the $L_\theta$ share one normal $a$. For the final statement, choose $\theta_k$ with $a^T\theta_k$ dense in $\mathbb R$. If $a^T\theta_j<a^T\theta_k$ then $H_{\theta_j}$ is strictly smaller than $H_{\theta_k}$ by positive measure, so $\phi(\theta_j)<\phi(\theta_k)$. Put $\psi(s)=\inf\{\phi(\theta_k):a^T\theta_k\geq s\}$, which is strictly increasing. Off a countable union of null sets, $\phi(u)\leq\phi(\theta_k)$ if and only if $a^Tu\leq a^T\theta_k$ for all $k$, so $\phi(u)$ lies between the one-sided limits of $\psi$ at $a^Tu$. Since $\psi$ has countably many jumps, and the corresponding hyperplanes are null, $\phi=\psi(a^Tu)$ almost everywhere. $\square$

Theorem 6 identifies the obstruction to exact marginal calibration precisely. Strictly monotone reparametrisations of a linear form, such as $(a^Tu)^3$, stay exact. Everything else fails at some parameter. Earlier accounts of this effect are discussed in §8.

**Example (exact crossing).** A two-dimensional exact example shows that crossing can occur with the full range of posterior probabilities. With $H=\{u:u_1u_2>0\}$, $\theta=0$ and standard Gaussian errors, $q_H(X)=\{1+(2U_1-1)(2U_2-1)\}/2$ for independent uniforms $U_i$. Then $F(t)=t+\tfrac{1-2t}{2}\log(1-2t)<t$ for $0<t<1/2$, and $F(1-t)=1-F(t)$. The crossover is exactly $1/2$. This should not be taken as typical: for $\theta$ inside a small convex $H^c$, §6 shows that the crossover can be extremely small.

# 5. The logarithmic tail rate at bounded hypotheses

Theorem 1(b) can be sharpened for exponential-power errors, $p>1$, and bounded hypotheses.

**Theorem 7.** In model (M) with $f=f_{p,d}$, $p>1$, suppose $H$ has positive measure and bounded essential support. Define the signed distance in the shape metric,
$$
\rho=\begin{cases}\inf_{u\in\partial C_H}\|A^{-1}(u-\theta)\|,&\theta\in C_H,\\ -\inf_{u\in C_H}\|A^{-1}(u-\theta)\|,&\theta\notin C_H.\end{cases}
$$
Then
$$
\lim_{t\downarrow0}\frac{\log\{F_{H,\theta}(t)/t\}}{(p\log(1/t))^{(p-1)/p}}=-\rho.
\tag{7}
$$

At interior points $F(t)/t$ tends to zero at this exponential rate, with coefficient the distance to the hull boundary. At exterior points the ratio diverges. At boundary points the coefficient is zero and the theorem does not determine the finer behaviour. For $p=2$ the denominator is $\sqrt{2\log(1/t)}$ and $\rho$ is Mahalanobis distance. The rate is asymptotic and can be slow. It should not be read as a finite-threshold approximation when $\rho$ is small: a small $\rho$ means the limit is reached only at extremely small $t$.

*Proof.* Take $\theta=0$, $A=I$, and let $h(v)=\max_{u\in C}v^Tu$ for unit $v$, where $C$ is the transformed hull. Then $\min_vh(v)=\rho$. For interior points, $B(0,r)\subseteq C$ if and only if $r\leq h(v)$ for all $v$. On the boundary a supporting hyperplane gives minimum zero. At an exterior point let $y$ be the nearest point of $C$ to zero. The projection inequality gives $h(-y/\|y\|)=-\|y\|$, and evaluating the support at $y$ gives the reverse bound.

Uniformly over unit $v$ and $\|u\|\leq M$, for $r\geq2M$,
$$
\frac{\|rv-u\|^p}{p}=\frac{r^p}{p}-r^{p-1}v^Tu+O_{p,M}(r^{p-2}),
\tag{8}
$$
since the Hessian of $\|z\|^p/p$ has operator norm $O_p(\|z\|^{p-2})$ and the segment from $rv$ to $rv-u$ stays at distance at least $r/2$ from zero. Fix $\delta>0$. Then uniformly in $v$, for all large $r$,
$$
e^{-r^p/p+r^{p-1}(h(v)-\delta)}\leq q_H(rv)\leq e^{-r^p/p+r^{p-1}(h(v)+\delta)}.
\tag{9}
$$
For the upper bound, use (8) with $v^Tu\leq h(v)$ on the bounded essential support. For the lower bound, choose $u_v\in E_H$ with $v^Tu_v>h(v)-\delta/6$. A small ball around $u_v$ has positive $H$-measure, and continuity gives $w^Tu>h(w)-\delta/2$ on it for $w$ near $v$. A finite subcover supplies a common mass factor. Fixed factors and the remainder are $o(r^{p-1})$.

Let $L=\log(1/t)$, $r_L=(pL)^{1/p}$, and for real $\beta$ let $R_\beta(L)$ be the large root of $R^p/p-\beta R^{p-1}=L$. Expanding at $r_L+s$ with $s$ bounded gives $R_\beta(L)=r_L+\beta+o(1)$. Compact balls are excluded as before. The lower bound in (9) and $h\geq\rho$ show that $q_H(x)\leq t$ implies $\|x\|\geq R_{\rho-\delta}(L)$. With (1), and since $\log S(R)=-R^p/p+O(\log R)$, this gives
$$
\limsup_{t\downarrow0}\frac{\log\{F_{H,0}(t)/t\}}{r_L^{p-1}}\leq-(\rho-\delta).
$$
For the reverse inequality, take a set $V$ of directions of positive surface measure on which $h<\rho+\delta$. (In dimension one, a minimising direction has positive angular probability.) On that cone the upper bound in (9) gives $q_H(rv)\leq t$ for $r\geq R_{\rho+2\delta}(L)$. Direction and radius are independent, so the cone tail is a fixed fraction of $S$, which gives the matching lower bound $-(\rho+2\delta)$. Let $\delta\downarrow0$. $\square$

# 6. High confidence at a fixed level

Theorem 1 does not establish failure at conventional confidence levels. In the example $H=(-2,-1)\cup(1,2)$, $\theta=0$, Gaussian errors, the posterior mass never exceeds $0.3893500291$, and numerically (V) is violated only for $t<0.2865$. The complement is the right object for high-confidence statements. This section gives a non-asymptotic result.

![Calibration of the fixed Gaussian hypothesis $H=(-2,-1)\cup(1,2)$ at $\theta=0$. The distribution function of the posterior mass crosses the uniform line; the logarithmic panel shows the violation at low thresholds. Computed from normal distribution functions and deterministic roots.](figures/calibration_cdf.png){width=88%}

**Theorem 8.** Assume (M) with $f$ log-concave. This covers $f_{p,d}$ for all $p\geq1$. Let $C$ be a closed convex set with $H^c\subseteq C$ almost everywhere, and let $0<s<1/2$. If $\theta$ satisfies
$$
Q_\theta(C)\leq s,
\tag{10}
$$
then
$$
P_\theta\{q_H(X)\geq1-s\}\geq\tfrac12>s.
\tag{11}
$$
If in addition $C$ is symmetric about $\theta$, then $q_H(X)\geq1-s$ with probability one.

So if $\theta\notin H$, then $H$ fails at level $1-s$, by a wide margin. Here (10) has a direct meaning: it is the posterior mass of $C$ computed as though the observation had landed exactly on the truth.

*Proof.* Let $w(x)=Q_x(C)$. The function $(x,u)\mapsto\mathbf 1_C(u)f(A^{-1}(x-u))$ is log-concave on $\mathbb R^{2d}$, so $w$ is log-concave by Prékopa's theorem (Saumard and Wellner 2014, Theorem 3.3). Therefore $D=\{x:w(x)>s\}$ is convex. By (10), $\theta\notin D$, so a hyperplane through $\theta$ has $D$ on one closed side. By central symmetry, $P_\theta(X\in D)\leq1/2$. Off $D$, $q_H(x)\geq1-Q_x(C)\geq1-s$, because $H^c\subseteq C$ almost everywhere. If $C$ is symmetric about $\theta$, then $w(\theta+y)$ is even and log-concave, so $w(\theta)\geq\sqrt{w(\theta+y)w(\theta-y)}=w(\theta+y)$. Hence $w\leq s$ everywhere. $\square$

**How this fits with the classification.** Let $C$ be bounded and convex, $H=C^c$, and $\theta\in\mathrm{int}\,C$, with $f=f_{p,d}$, $p>1$. Then Theorem 4 places $H$ in the Crossing row, because $C_H=\mathbb R^d$. Theorem 1 applied to $C$ says that $P_\theta\{q_H(X)\geq1-s\}=o(s)$ as $s\to0$. Theorem 8 says that, when $Q_\theta(C)<1/2$, this probability exceeds $s$ for every $s$ in $[Q_\theta(C),1/2)$. So below level $1/2$, the behaviour near $t=1$ in the Crossing row can only appear at $s<Q_\theta(C)$. The order of quantifiers matters: the crossover depends on $H$, and it can be negligibly small. This matches Balch, Martin and Ferson (2019), whose existence theorem chooses the hypothesis after fixing the level.

**Worked example: conjunction analysis.** Take $d=2$ and isotropic Gaussian error with standard deviation $\sigma$ in the encounter plane. A collision occurs if the true miss vector $\theta$ lies in the disc $C$ of radius $R$ about the origin. So $H=C^c$ ("no collision") is false, and the computed collision probability is $P_c(x)=Q_x(C)$. Let $m=\|\theta\|<R$. Table 2 compares bound (11) with exact values from non-central chi-square distributions. Since $Q_\theta(C)\leq Q_0(C)=1-e^{-R^2/(2\sigma^2)}$, condition (10) depends mainly on $\sigma/R$. A collision is certain in every row.

Table: Probability that the computed collision probability is at most $s$, when a collision is certain.

| $\sigma/R$ | $m/R$ | $Q_\theta(C)$ | $P(P_c\leq0.1)$ | $P(P_c\leq0.01)$ | $P(P_c\leq0.001)$ |
|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 0.393 | 0.174 | 0.0101 | 0.0006 |
| 2 | 0 | 0.118 | 0.842 | 0.0728 | 0.0063 |
| 2 | 0.5 | 0.114 | 0.847 | 0.0788 | 0.0074 |
| 10 | 0 | 0.00499 | 1 | 1 | 0.200 |
| 10 | 0.9 | 0.00497 | 1 | 1 | 0.201 |
| 20 | 0.5 | 0.00125 | 1 | 1 | 0.800 |

In each cell with $s\geq Q_\theta(C)$, the theorem guarantees at least one half. The table also shows the size of the effect. When $\sigma/R=10$, a certain collision produces $P_c\leq0.01$ with probability one. This is the *dilution* effect, and Theorem 8 extends it from a disc to any hypothesis whose complement fits inside a convex set of small posterior mass. Anisotropic $\Sigma$ is covered by the shape metric. The model ignores covariance error, non-Gaussian tails and curved encounters. Elkantassi and Davison (2022) study this problem with likelihood inference on the miss distance, and note the downward tendency of the estimated collision probability.

# 7. Heavy tails

For $0<p\leq1$, Theorem 1(c) gives globally valid nonconvex sets, so the equivalence fails. The picture is sharper for bounded sets.

**Proposition 9 (tail trichotomy for bounded hypotheses).** Let $H$ have positive measure and bounded essential support, in the model with $f_{p,d}$ and $A=I$, and fix $\theta$.

(i) If $p>1$, then $F_{H,\theta}(t)/t\to0$, $\to\infty$ or has logarithmic scale coefficient zero, according as $\theta$ is interior to, exterior to, or on the boundary of $C_H$ (Theorem 7).

(ii) If $0<p<1$, or $p=1$ and $d\geq2$, then $F_{H,\theta}(t)/t\to\infty$ at every $\theta$, regardless of hull position.

(iii) If $p=1$ and $d=1$, then with $\theta=0$, $K_\pm=\tfrac12\int_He^{\pm u}\,du$, for all small $t$,
$$
F_{H,0}(t)=\frac t2\Bigl(\frac1{K_+}+\frac1{K_-}\Bigr).
\tag{12}
$$
The constant depends on exponential moments of the set, and not only on hull position.

*Proof of (ii).* Take $\theta=0$ and let $H\subseteq B(0,M)$ up to a null set. For $p\leq1$, subadditivity gives $\|x-u\|^p\geq\|x\|^p-M^p$, so $q_H(x)\leq Ce^{-\|x\|^p/p}$ with $C=c_{p,d}\lambda_d(H)e^{M^p/p}$. Define $R_t$ by $Ce^{-R_t^p/p}=t$. Then $F_{H,0}(t)\geq S(R_t)$, and (1) gives $S(R_t)/t\sim\mathrm{const}\cdot R_t^{d-p}\to\infty$, since $d>p$.

*Proof of (iii).* Outside a large interval, $q_H(x)=K_+e^{-x}$ for $x>0$ and $K_-e^{x}$ for $x<0$. Inside, $q_H$ has a positive minimum. So for small $t$ the sublevel set consists exactly of the two tails $X\geq\log(K_+/t)$ and $X\leq-\log(K_-/t)$, which have probabilities $t/(2K_+)$ and $t/(2K_-)$. $\square$

**Hull position does not decide validity at $p=1$, $d=1$.** Let $H=(-B,-\delta)\cup(\delta,B)$ with $0<\delta<B$, and $\theta=0$. Then $\theta$ is false and interior to the hull, and $K_+=K_-=\sinh B-\sinh\delta=:K$, so $F(t)=t/K$ for small $t$. If $K>1$ then $F(t)<t$ and (V) fails. For example $B=1.5$, $\delta=0.1$ gives $K\approx2.03$. If instead $a=\log2$ and $b=\log(5/2)$ with $H=(-b,-a)\cup(a,b)$, then $K=3/10$. In that case $|x-u|\geq|x|-u\,\mathrm{sign}(x)$ gives $q_H(x)\leq\tfrac3{10}e^{-|x|}$, and $e^{-|X|}$ is uniform at $\theta=0$, so (V) holds. Two sets with the same qualitative shape therefore behave oppositely.

# 8. Interpretation, prior work and limits

**Which results are classical.** The sufficient direction is elementary. Halfspace validity is a special case of Martin (2023, Theorem 5), which treats hypotheses defined through ordered-group homomorphisms in invariant models, and the extension by separation is immediate. Martin (2024, Proposition 1) gives the Gaussian dominance mechanism behind the Anti-conservative row, for a hypothesis whose complement lies in a supporting halfspace, and complementation gives the corresponding Conservative comparison. In the Gaussian case, the converse in Theorem 1(b) combines the classical growth of a Laplace transform along rays, which is governed by the support function of the convex hull of the measure's support, with a comparison to the sampling tail. The contribution claimed here is the exact calibration classification and the items built on it: the tail characterisation (T), including its dimension effect, Corollary 5, Theorem 6, Theorem 7, Theorem 8 and Proposition 9. Searches on 1 and 4 November 2026 of the false-confidence, inferential-model, confidence-distribution and problem-of-regions literatures found no earlier statement of these. That is evidence about what we found, not proof of priority.

**The problem of regions.** At $p=2$, $q_H(x)$ is the first-order confidence value of Efron and Tibshirani (1998): the bootstrap probability of a region, which they identify with a flat-prior posterior probability. They show by asymptotic expansion that for smooth boundaries it exceeds or falls below the frequentist level according to whether the boundary curves away from or towards the observation. Shimodaira (2004, 2008) develops multiscale bootstrap corrections for smooth and non-smooth boundaries, and Fraser (2011) shows the same curvature effect for the length of a normal mean. The problem-of-regions corrections use local asymptotic expansions and aim to remove the bias, while Fraser also gives exact calculations for particular non-linear functionals. Our results give an exact classification for arbitrary measurable hypotheses in the stated location models, at the cost of saying less about moderate thresholds. Theorem 6 states exactly which functions escape the shared effect.

**Confidence distributions.** Cunen, Hjort and Schweder (2020) argue that confidence distributions are free of false confidence for interval statements. In dimension one, with these models, intervals are convex, so Theorem 1 agrees, and it shows that unions of separated intervals fall outside that protection. Carmichael and Williams (2018) give an exposition of the false confidence theorem, with worked examples, and Martin (2019) develops the view that validity, which requires non-additive beliefs, protects against it.

**Fragility and limits.**

- The hull depends on all of the essential support. Adding a tiny positive-measure speck far from $H$ can move $\theta$ into the hull interior, and so changes the all-threshold verdict, while changing $q_H$ everywhere by at most $\lambda_d(\text{speck})\,g(0)/|\det A|$. This is why Theorem 8 is a necessary companion.
- All-threshold failure can coexist with fixed-threshold consistency. If $X$ is a Gaussian mean of $n$ observations and $\theta$ is at positive distance from $H$, then $E_\theta q_H(X)\to0$, and the violating thresholds move towards zero.
- The model assumes known shape, an unrestricted Euclidean parameter and exact matching of sampling and posterior error laws. With sampling law $N(\theta,\Sigma)$ and posterior $N(x,\Psi)$, the posterior probability of a halfspace through $\theta$ is $\Phi(\tau Z)$, with $\tau^2=(a^T\Sigma a)/(a^T\Psi a)$, and it is uniform only when $\tau=1$. Unknown scale, informative priors and approximate posteriors are not covered.
- Theorem 8 needs log-concavity. It does not locate the crossover exactly, and the bound one half can be far from the exact probability, as Table 2 shows.

**Open directions.** Natural next steps are the unknown-scale Gaussian case, group-invariant models with right-Haar posteriors, and an exact crossover formula for convex sets.

# 9. Reproducibility

The package contains three deterministic check scripts. They use floating-point arithmetic and are sanity checks of the formulas, not proofs.

- `verify_examples.py` covers 25 Gaussian cases for $H=(-2,-1)\cup(1,2)$, with thresholds between $10^{-2}$ and $10^{-300}$. At $t=10^{-300}$ the scaled log ratio in (7) is about $-2.03$ at $\theta=0$ (limit $-2$), and about $0.98$ at $\theta=3$ (limit $1$).
- `verify_exponential_power.py` covers 24 cases at $p=1.5$ and $p=3$. At $\theta=0$ the ratios move towards $-2$ from opposite sides.
- `verify_extensions.py` confirms the crossover $0.2865$, the constant $K$ of §7, the direction of the tail ratios in Lemma 0, and Table 2, including a test of Theorem 8 at each eligible entry.

`make_figure.py` draws Figure 1.

# 10. Conclusion

The question "which hypotheses are affected?" has an exact geometric answer in location models with light tails. A hypothesis is valid at a parameter when that parameter lies outside the interior of its essential convex hull, and only then. The light-tail condition is exact. For a practitioner, the lesson comes from the complement: if everything outside a hypothesis fits inside a convex set that carries little posterior mass even when the data land on the truth, then that hypothesis will look well supported at least half the time, whether or not it is true. A low computed collision probability under large uncertainty is the standard example. Only halfspaces are calibrated from both sides, and only essentially linear functions have exactly calibrated marginal posteriors.

# References

Balch, M. S., Martin, R., & Ferson, S. (2019). Satellite conjunction analysis and the false confidence theorem. *Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 475*(2227), 20180565. https://doi.org/10.1098/rspa.2018.0565

Carmichael, I., & Williams, J. P. (2018). An exposition of the false confidence theorem. *Stat, 7*(1), e201. https://doi.org/10.1002/sta4.201

Casella, G., & Berger, R. L. (1987). Reconciling Bayesian and frequentist evidence in the one-sided testing problem. *Journal of the American Statistical Association, 82*(397), 106–111. https://doi.org/10.1080/01621459.1987.10478396

Cunen, C., Hjort, N. L., & Schweder, T. (2020). Confidence in confidence distributions! *Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 476*(2237), 20190781. https://doi.org/10.1098/rspa.2019.0781

Efron, B., & Tibshirani, R. (1998). The problem of regions. *The Annals of Statistics, 26*(5), 1687–1718. https://doi.org/10.1214/aos/1024691353

Elkantassi, S., & Davison, A. C. (2022). Space oddity? A statistical formulation of conjunction assessment. *Journal of Guidance, Control, and Dynamics, 45*(12), 2258–2274. https://doi.org/10.2514/1.G006282

Fraser, D. A. S. (2011). Is Bayes posterior just quick and dirty confidence? *Statistical Science, 26*(3), 299–316. https://doi.org/10.1214/11-STS352

Martin, R. (2019). False confidence, non-additive beliefs, and valid statistical inference. *International Journal of Approximate Reasoning, 113*, 39–73. https://doi.org/10.1016/j.ijar.2019.06.005

Martin, R. (2023). Fiducial inference viewed through a possibility-theoretic inferential model lens. In E. Miranda, I. Montes, E. Quaeghebeur, & B. Vantaggi (Eds.), *Proceedings of the Thirteenth International Symposium on Imprecise Probability: Theories and Applications* (PMLR Vol. 215, pp. 299–310). https://proceedings.mlr.press/v215/martin23a.html

Martin, R. (2024). Which statistical hypotheses are afflicted with false confidence? In *Lecture Notes in Computer Science* (pp. 140–149). Springer. https://doi.org/10.1007/978-3-031-67977-3_15

Martin, R. (2026). Possibilistic inferential models: A review. *Journal of the American Statistical Association, 121*(553), 807–826. https://doi.org/10.1080/01621459.2025.2606127

Saumard, A., & Wellner, J. A. (2014). Log-concavity and strong log-concavity: A review. *Statistics Surveys, 8*, 45–114. https://doi.org/10.1214/14-SS107

Shimodaira, H. (2004). Approximately unbiased tests of regions using multistep-multiscale bootstrap resampling. *The Annals of Statistics, 32*(6), 2616–2641. https://doi.org/10.1214/009053604000000823

Shimodaira, H. (2008). Testing regions with nonsmooth boundaries via multiscale bootstrap. *Journal of Statistical Planning and Inference, 138*(5), 1227–1241. https://doi.org/10.1016/j.jspi.2007.04.001
