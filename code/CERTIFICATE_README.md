# Reusable exact-certification toolkit

`colas_certificate_v2.py` is written so that its core is independent of the
specific model and can certify sign conditions, injectivity discriminants, and
Jacobian nondegeneracy for OTHER moment maps with polynomial/analytic
parameter structure. The reusable pieces:

- `IV`: outward-rounded fixed-point dyadic interval arithmetic (PREC bits;
  exact sums, directed-rounded products/quotients; inclusion isotone).
- `PE`: error-carrying polynomials in several variables with coefficient-level
  arithmetic. This is the key device: covariance-type combinations
  (e.g. PE - S^2), derivative numerators, discriminants, and determinants are
  formed at the coefficient level, where cancellation is exact, so interval
  dependency does not blow up on large boxes. Truncation/model error rides
  along as a scalar budget with sup-norm propagation.
- `shift_at_one`: rescale-to-the-cell then Taylor-shift at exactly 1 by pure
  additions -- no multiplicative width amplification (the naive shift at
  t = lambda_0 multiplies widths by lambda_0^degree).
- adaptive bisection drivers with per-gate leaves and achieved-bound reports.

Recipe for a new problem: (1) express your moment map's ingredients as
polynomials with exact rational coefficients plus a rigorous truncation
budget; (2) build every quantity whose SIGN or BOUND you need as a single PE
via coefficient arithmetic; (3) write gates; (4) run the driver. The
monotone-triangular lemma (SI S3 v2) is likewise generic: for a triangular
moment map (one coordinate free of one parameter), certified sign conditions
on two partial derivatives plus one 2x2 discriminant upgrade pointwise rank
to global injectivity on boxes.

`ieff_face_certificate.py` applies the same devices (exact coefficient
tables, midpoint Taylor shifts, adaptive bisection) to a third gate:
a uniform positive lower bound for the efficient boundary information
on the null face.
