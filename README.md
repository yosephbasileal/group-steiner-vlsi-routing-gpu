# GPU-Accelerated VLSI Routing Using Group Steiner Trees

A parallel approximation algorithm for the GSP based off an existing heuristic on a distributed architecture. Our implementation uses a CUDA-aware MPI-based approach to compute the approximate minimum-cost Group Steiner tree for several industry-standard VLSI graphs.


# Usage
`make`

`mpirun -n 64 --hostfile $HOME/hosts ./main -t < testdata/wrp2-11.stp`

# Dependencies
- Open MPI 1.8.5
- CUDA 6.5
