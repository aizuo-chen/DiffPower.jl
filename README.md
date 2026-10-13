# DiffPower.jl

This repository contains the PowerEMT Julia module and examples for electromagnetic transient (EMT) simulation of power systems.

## Contents

- `src/`: network composition, power flow, device models, fault models, and numerical Jacobian utilities.
- `src/ieee_cases/`: IEEE 9-, 14-, 39-, and 118-bus case data.
- `examples/`: single-device and IEEE network simulation scripts.

## Running an example

Run from a Julia environment with the dependencies used by the selected script installed:

```julia
include("examples/single_device_test.jl")
```

The core module imports ModelingToolkit, SciMLBase, JLD2, ForwardDiff, SymbolicIndexingInterface, DifferentiationInterface, ADTypes, and SparseMatrixColorings, plus Julia standard libraries. Examples additionally use solver and plotting packages such as OrdinaryDiffEq, OrdinaryDiffEqSDIRK, CairoMakie, LinearSolve, and DelayDiffEq.

A standalone `Project.toml` and compatibility bounds have not yet been supplied. The examples have not been validated in a fresh standalone environment as part of this repository preparation.
