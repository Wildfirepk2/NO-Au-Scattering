import Pkg;

dependencies = [
                "CSV",
                "CUDA",
                "CairoMakie",
                "DataFrames",
                "Distances",
                "GLMakie",
                "LinearAlgebra",
                "Molly",
                "NearestNeighbors",
                "PhysicalConstants",
                "StatsBase",
                "Unitful",
                "XLSX",
                ]

Pkg.add(dependencies)
Pkg.update()
