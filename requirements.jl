import Pkg;

dependencies = [# for Au slab
                "NearestNeighbors", 
                "LinearAlgebra", 

                # read/write files
                "DataFrames", 
                "CSV", 

                # unit integration
                "Unitful", 

                # GPU integration
                "CUDA",

                # MD
                "Molly",
                "GLMakie",
                ]

Pkg.add(dependencies)
