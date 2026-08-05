#!/usr/bin/env julia

# Activate the project environment
using Pkg
Pkg.activate(".")

# Run the test script
include("src/julia/test_kam_torus.jl")
