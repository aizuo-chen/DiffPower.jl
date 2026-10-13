# DiffPower.jl

This repository contains the DiffPower Julia module for electromagnetic transient (EMT) simulation using averaged inverter models. The tested case is a single inverter connected to an infinite bus, with simulation results consistent with MATLAB/Simulink using the same averaged model. Voltage-command delay is modeled explicitly using delay differential equations (DDEs). We plan to register DiffPower as a Julia package in the near future.
## Contents

- `src/`: network composition, power flow, device models, fault models, and numerical Jacobian utilities.
- `src/ieee_cases/`: IEEE 9-, 14-, 39-, and 118-bus case data.
- `examples/`: single-device and IEEE network simulation scripts.

## Running an example

From the repository directory, install dependencies and run the example:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("examples/single_device_test.jl")
```

The core module imports ModelingToolkit, SciMLBase, JLD2, ForwardDiff, SymbolicIndexingInterface, DifferentiationInterface, ADTypes, and SparseMatrixColorings, plus Julia standard libraries. Examples additionally use solver and plotting packages such as OrdinaryDiffEq, OrdinaryDiffEqSDIRK, CairoMakie, LinearSolve, and DelayDiffEq.

The Julia package is named `DiffPower`; its entry point is `src/DiffPower.jl`.
