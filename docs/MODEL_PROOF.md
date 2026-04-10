# Model Proof (Sketches)

This document provides proof sketches for core matrices used by GPBLUP.
It is intended as a practical reference, not a full mathematical treatise.

## Notation

- M: genotype matrix with entries in {0,1,2}
- p: vector of allele frequencies per SNP
- Z = M - 2p: centered genotype matrix
- c = 2 * sum_j p_j (1 - p_j): scaling constant
- G = (Z * Z^T) / c: genomic relationship matrix

## Proposition 1: G is positive semidefinite

Claim:
- For any vector x, x^T G x >= 0.

Proof sketch:
- G = (Z * Z^T) / c with c > 0.
- x^T G x = (1 / c) * x^T Z Z^T x = (1 / c) * || Z^T x ||^2.
- Norm squares are nonnegative, therefore G is PSD.

## Proposition 2: Blended G is positive definite

Define:
- G_star = (1 - w) G + w I, with 0 < w <= 1.

Claim:
- G_star is positive definite and invertible.

Proof sketch:
- For any nonzero x, x^T G_star x = (1 - w) x^T G x + w x^T x.
- x^T G x >= 0 (PSD) and x^T x > 0 for x != 0.
- Therefore x^T G_star x >= w ||x||^2 > 0.

## Proposition 3: A is positive semidefinite

Let A be the pedigree-based additive relationship matrix.
Under standard assumptions (Mendelian inheritance and nonnegative variances),
A is PSD.

Proof sketch:
- A is a covariance matrix of additive genetic values.
- Covariance matrices are PSD by construction.

## Proposition 4: A inverse via Henderson's rules

Claim:
- A inverse can be assembled from pedigree records using local rules on
  (animal, sire, dam) triplets.

Proof sketch:
- Henderson's method constructs A inverse from the inverse of the
  L D L^T decomposition of A, where L is unit lower-triangular and D is
  diagonal. Each pedigree record contributes a small update to A inverse.
- This yields the same result as inverting A directly but is far more efficient.

## Assumptions and Practical Notes

- Allele frequencies p are estimated from the available data.
- Missing genotypes are handled by program-specific rules before matrix builds.
- Scaling and weights should be reported with results for reproducibility.

For detailed implementation and parameter names, see relgeno/RELGENO_PARAMS.md
and relped/README.md.
