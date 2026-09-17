# Minimal Product Sets in Finite Metacyclic Groups

This directory implements the recursive arithmetic region $\mathcal T_G$, reproduces the numerical examples in Section 7, and performs the independent computational checks accompanying the manuscript. The calculations are consistency checks and are not used in the mathematical proofs.

## Requirements

Julia 1.12 and Python 3.10 or later. Only standard libraries are used; no package installation, solver, or external data is required. Run the following commands from this directory.

## Section 7 examples

```sh
julia --startup-file=no scripts/examples.jl
```

Evaluates the recursion and checks explicit witnesses for the three intermediate-value examples and for $K(9,57,19,4)$ at $(180,180)$.

## Small-group checks

```sh
julia --startup-file=no scripts/check_small.jl
```

Compares the recursion with exhaustive product-set enumeration for every valid Hölder presentation of order at most 12. Also checks the group law, Yang subgroups, the Berchenko–Kogan tables at orders 21, 39 and 55, and selected faithful-complement and split/nonsplit dihedral-extension tables.

## Independent dihedral lifts

```sh
python3 scripts/check_dihedral.py
```

Checks 19,202 concrete boundary lifts, including nonsplit extensions and 248 lifts with a nonabelian kernel, independently of the Julia recursion.

The Julia implementation uses Yang's convention $\eta\tau\eta^{-1}=\tau^h$. The Python checker uses its separately specified integer coordinates. All scripts print their results and stop on a failed assertion; they create no output files.
