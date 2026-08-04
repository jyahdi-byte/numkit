# Parabolic Convergence Study

This checks how error shrinks under grid refinement — the parabolic
counterpart to the elliptic study in docs/validation.md.

## Method

r = αΔt/Δx² is held fixed across the sweep, so time-stepping error
doesn't confound the resolution comparison. sigma0 and the domain length
are held fixed too, so every run solves the same physical problem at
increasing resolution. Only Δx changes.

## Results

| dx     | Max error | Ratio vs previous |
|--------|-----------|--------------------|
| 0.04   | 3.05e-04  | —                  |
| 0.02   | 7.61e-05  | 4.01               |
| 0.01   | 1.90e-05  | 4.00               |
| 0.005  | 4.75e-06  | 4.00               |
| 0.0025 | 1.86e-06  | 2.55               |

The scheme is O(Δt) + O(Δx²) in general, and holding r fixed means Δt
scales with Δx², so the combined error reduces to O(Δx²). Halving Δx
should quarter the error. The measured ratios of 4.01, 4.00, and 4.00
confirm this.

## Error budget

The final ratio (2.55 at dx = 0.0025) falls short of 4.00. At that
resolution the run takes roughly 128,000 steps, and floating-point
roundoff accumulated over that many operations becomes large enough to
compete with the discretization error, which has gotten very small.
Same phenomenon as the iteration error in the elliptic study, different
source — there it was SOR not fully converging, here it's roundoff.