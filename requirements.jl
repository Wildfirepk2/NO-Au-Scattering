import Pkg;

dependencies = [# for Au slab
                "NearestNeighbors", 
                "LinearAlgebra", 

                # read/write files
                "DataFrames", 
                "CSV", 
                "XLSX",

                # unit integration
                "Unitful",
		"PhysicalConstants",

                # GPU integration
                "CUDA",

                # MD
                "Molly",
		"GLMakie",
		"CairoMakie",
		"Distances",
                ]

Pkg.add(dependencies)
