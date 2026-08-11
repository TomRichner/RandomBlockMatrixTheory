# RandomBlockMatrixTheory

A MATLAB class for sparse random connectivity matrices with **block-structured means
and block-structured standard deviations**, indexed by both the postsynaptic and the
presynaptic cell type.

## Why

Harris et al. (2023) build `W = S ∘ (A·D + M)` with `D` diagonal and `M = uvᵀ` rank-1.
Both are **presynaptic-only**: every entry in column *j* shares one mean and one standard
deviation. There is no way to say "E→E is stronger than I→E."

`RMTBlocks` generalizes this to

```
W = S .* ( A .* (g_sigma * Sigma) + g_mu * M ) + shift * eye(N)

Sigma_ij = sigma_tilde( type(i), type(j) )     % N x N, block-constant
M_ij     = mu_tilde(    type(i), type(j) )     % N x N, block-constant
A_ij     ~ N(0,1) i.i.d.,   S_ij ~ Bernoulli(alpha)
```

where `mu_tilde` and `sigma_tilde` are `D × D`. Right-multiplication by a diagonal is
pure column scaling and is *structurally incapable* of row dependence, which is why the
random term must become elementwise.

## Index convention

`(a ← b)` = row `a` **postsynaptic**, column `b` **presynaptic**, matching `A` in
`dx/dt = A·x`. For `D = 2` with type 1 = E and type 2 = I:

| entry | meaning |
|---|---|
| `mu_tilde(1,1)` = `mu_EE` | E receives from E |
| `mu_tilde(1,2)` = `mu_EI` | E receives from I |
| `mu_tilde(2,1)` = `mu_IE` | I receives from E |
| `mu_tilde(2,2)` = `mu_II` | I receives from I |

Dale's law is a **column** constraint: signs must be constant down each column of
`mu_tilde`. Magnitudes may differ freely between rows.

## Usage

```matlab
setup_paths();                       % once per session, from the repo root

rmt = RMTBlocks(300);
rmt.alpha = 1/3;
rmt.f = 0.5;                         % scalar -> [0.5 0.5]
rmt.mu_tilde    = [ 0.23 -0.31 ;     % (a <- b)
                    0.23 -0.31 ];
rmt.sigma_tilde = [ 0.08  0.08 ;
                    0.08  0.08 ];

W = rmt.W;
rmt.display_parameters();            % set / sparse-effective / measured, per block
rmt.plot_spectrum();                 % eigenvalues with the predicted radius overlaid
```

`g_mu` and `g_sigma` scale the mean and the disorder separately; setting them equal
reproduces a single global "level of chaos" multiplier on all of `W`.

`f` may be a vector of `D` fractions summing to 1. Empty populations (`f = [1 0]`) are
legal and the predictions collapse correctly to the lower-`D` answer. Use
`set_types(f, mu_tilde, sigma_tilde)` to change `D` atomically.

## Predictions

```
mu_s(a,b)       = alpha * g_mu * mu_tilde(a,b)
sigma_s_sq(a,b) = alpha(1-alpha)*(g_mu*mu_tilde(a,b))^2 + alpha*(g_sigma*sigma_tilde(a,b))^2

lambda_O = eigenvalues of  K(a,b) = N * f_b * mu_s(a,b)        (sorted by |lambda|)
R        = sqrt( Perron root of  V(a,b) = N * f_b * sigma_s_sq(a,b) )
```

With **column-uniform** blocks `K` and `V` have identical rows, hence rank 1, and these
reduce *exactly* to Harris Eqs. (17) and (18). That reduction is asserted in the tests.

Note that `sigma_s_sq` contains an `alpha(1-alpha)·mu²` term: **sparsity converts mean
into variance**, so `mu` and `sigma` are not independent contributors to the bulk radius.
For a typical E/I configuration at `alpha = 1/3` this term supplies ~86–91% of the block
variance, meaning the bulk radius is dominated by the means rather than by `sigma_tilde`.

### Status of the predictions

No published result covers sparse **and** block mean **and** block variance
simultaneously:

| | mean structure | variance structure | sparse |
|---|---|---|---|
| Harris et al. (2023) | rank-1, column-only | column-only | yes |
| Aljadeff, Renfrew & Stern (2015) | zero mean | full `D × D` blocks | no |
| Ahmadian, Fumarola & Miller (2015) | arbitrary, any rank | **separable only** (`σ_ij = l_i r_j`) | no |

Ahmadian's `LJR` ensemble cannot express four independent `sigma` blocks, since
separability forces `σ_EE·σ_II = σ_EI·σ_IE`.

The formulas above **compose** these results — Harris's blockwise sparse moment matching
feeding Aljadeff's Perron-root radius, with the finite-rank block mean contributing only
outliers (Ahmadian §V C) — and are therefore a **conjecture verified numerically**, not a
citable theorem. This is the same standard Harris meets for his own Eq. (18), which is
validated against numerical spectra in his Fig. 2(b) rather than proven.

## Tests

```matlab
test_RMTBlocks_equivalence     % 38 assertions
test_RMTBlocks_predictions     %  8 assertions
```

**`test_RMTBlocks_equivalence`** proves the generalization is strict. It builds seven
configurations (all four `zrs_mode`s, `alpha ∈ {1, 0.5}`, nonzero `shift`, and the
`f = 1` empty-population case) through the *same* alias API on `RMT`, `RMTMatrix` and
`RMTBlocks`, and asserts the resulting `W` is **bit-for-bit identical**. It also checks
the Harris Eq. (17)/(18) reduction, alias round-trips and their error identifiers, empty
populations, and genuine block behavior. Sections needing `RMT` or `RMTMatrix` skip with
a printed note if those classes are not on the path, so this repo is standalone-clonable.

**`test_RMTBlocks_predictions`** validates the conjecture empirically across five
configurations × two `alpha` × two `N`. Measured results:

- bulk edge within **2.6%** of `R` at `N = 1000` in every configuration
- error shrinks with `N` in every configuration (median 4.0% → 0.9%)
- convergence sweep: 3.7% → 0.9% → 1.2% → 0.4% for `N` = 200, 400, 800, 1600
- predicted outliers matched within **4.9%** of their modulus at `N = 1000`

The discriminating test uses `f = [0.9 0.1]` with a hyper-variable minority, where the
block Perron root and the naive population-averaged variance disagree by 2×
(Aljadeff's Fig. 1(a,b) point): the block formula lands within **4.3%** of the measured
edge while the naive average is off by **90%**.

The bulk edge is estimated from a **percentile** of `|lambda|`, not the max: Harris
(§III D 2) documents a small number of "local eigenvalue outliers" escaping the circular
support — which is exactly why his Eq. (18) is introduced as the *approximate* radius —
so `max(|lambda|)` would reject a correct formula.

## ZRS modes

All four of Harris's zero-row-sum variants are supported. With non-uniform blocks a
one-time suppressible warning (`RMTBlocks:UntestedZRSWithBlocks`) is issued, because:

- **`SZRS`** zeroes *all* row sums, which no longer removes all mean-induced outliers
  once `M` has rank > 1 — a second outlier can survive.
- **`Partial_SZRS`** corrects only the random component and so preserves block means by
  construction. It is the mode that composes correctly with block structure.

`R` remains valid under every mode: the ZRS conditions exist precisely to pull stray
eigenvalues back *inside* that radius.

## References

- Harris, Meffin, Burkitt & Peterson (2023) *Effect of sparsity on network stability in
  random neural networks obeying Dale's law.* Phys. Rev. Res. **5**, 043132.
- Aljadeff, Stern & Sharpee (2015) *Transition to chaos in random networks with
  cell-type-specific connectivity.* Phys. Rev. Lett. **114**, 088101.
- Aljadeff, Renfrew & Stern (2015) *Eigenvalues of block structured asymmetric random
  matrices.* J. Math. Phys. **56**, 103502.
- Ahmadian, Fumarola & Miller (2015) *Properties of networks with partially structured
  and partially random connectivity.* Phys. Rev. E **91**, 012820.
- Sompolinsky, Crisanti & Sommers (1988) *Chaos in random neural networks.*
  Phys. Rev. Lett. **61**, 259.
