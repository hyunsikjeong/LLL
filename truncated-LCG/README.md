# Breaking Truncated LCG

A truncated LCG is a pseudo random generator that outputs the `n` leading bits `y_i` of `x_i`,
where `(x_i)` is such that `x_i + 1 = α ⋅ x_i + β mod M`, and `α` is an integer, coprime with `M`.
The first state `x_0` is unknown.

## How to use

The constructor of TLCGBreaker is `TLCGBreaker(lcg_a, lcg_b, lcg_modulus_bit, num_truncated)`, where
- `lcg_a`: `α` of LCG
- `lcg_b`: `β` of LCG
- `lcg_modulus_bit`: `m` when `M = 1 << m`
- `num_truncated`: `l` when `y_i = x_i >> l`

Several methods are given:
- `set_outputs(self, outputs)` - Set outputs `y_0`, `y_1`, ... , `y_k`
- `recover_state(self)` - After setting outputs with `set_outputs()`, use to recover the original state `x_0`.
  It will return `True` when successful to recover the state, and `False` when failed to.
- `get_output(self, idx)` - Get output `y_idx` after the state is recoverd. `idx` can be a negative number.

## Related Challenges
- DownUnderCTF 2020 - LSB||MSB Calculation Game

## References
- https://www.math.cmu.edu/~af1p/Texfiles/RECONTRUNC.pdf
- https://crypto.stackexchange.com/questions/37836/problem-with-lll-reduction-on-truncated-lcg-schemes
- https://www.josephsurin.me/posts/2020-09-20-downunderctf-2020-writeups#lsb-msb-calculation-game
