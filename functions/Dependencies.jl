# init all package dependencies

using Unitful
using DataFrames
using CSV
using XLSX
import PhysicalConstants.CODATA2018: N_A, k_B # Avogadro's number and Boltzmann's constant with Unitful units
using NearestNeighbors
using LinearAlgebra
using Molly
if isaac
    using CairoMakie
else
    using GLMakie
end
# using GLMakie
# using CairoMakie
using Dates
using Distances
# using FLoops
