import Pkg;

dependencies = [# for Au slab
                "NearestNeighbors", 
                "LinearAlgebra", 

                # read/write files
                "DataFrames", 
                "CSV", 

                # unit integration
                "Unitful",
		"PhysicalConstants",

                # GPU integration
                "CUDA",

                # MD
                "Molly",
		"GLFW",
		"GLMakie",
                ]

Pkg.add(dependencies)
