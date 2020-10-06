# Low-Density Attack

For a given set of positive integers `A = {a_1, . . . , a_n} (a_i != a_j)` and a given positive integer `s`,
determining whether there exists a subset of `A` with its sum being `s`, or finding a vector
`e = (e_1, . . . , e_n) ∈ {0, 1}^n` satisfying `Σ_{i=1}^n a_i·e_i = s`, is called the subset sum problem
(or the knapsack problem), and is known as an NP-hard problem in general.

Low-density attack is a method which works effectively against subset sum problems with low density.
The density of the subset sum problem `d` is defined by `d = n/(log2 max_i{a_i})`.

There are two well-known algorithms:
- Lagarias and Odlyzko (LO) algorithm (works on `d < 0.6463`)
- Coster, Joux, LaMacchia, Odlyzko, Schnorr, and Stern (CJLOSS) algorithm (works on `d < 0.9408`)

## How to use

If the integers (`A`) is `arr` and the target sum is `s`, use as:
```python
attack = LOAttack(arr, s)
print(attack.solve()) # Print result of LO Algorithm
attack = CJLOSSAttack(arr, s)
print(attack.solve()) # Print result of CJLOSS Algorithm
```

## Related Challenges

## Reference

- Low-Density Attack Revisited https://eprint.iacr.org/2007/066.pdf
