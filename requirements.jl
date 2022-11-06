import Pkg;

dependencies = [
                "Molly",
                "NearestNeighbors", 
                "LinearAlgebra", 
                "DataFrames", 
                "CSV", 
                "Unitful", 
                "CUDA",
                "BenchmarkTools",
            ]

Pkg.add(dependencies)
Pkg.precompile()
