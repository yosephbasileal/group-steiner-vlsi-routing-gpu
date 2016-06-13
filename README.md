# GPU-Accelerated VLSI Routing Using Group Steiner Trees

A parallel approximation algorithm for the GSP based off an existing heuristic on a distributed architecture. Our implementation uses a CUDA-aware MPI-based approach to compute the approximate minimum-cost Group Steiner tree for several industry-standard VLSI graphs.

#Contributors
- Basileal Imana
- Venkata Suhas Maringanti
- Peter Yoon

#Publications
- [Poster and Abstract presented at SC15](http://sc15.supercomputing.org/sites/all/themes/SC15images/tech_poster/tech_poster_pages/post119.html)

# Usage
`make`

`mpirun -n 64 --hostfile $HOME/hosts ./main -t < testdata/wrp2-11.stp`

# Dependencies
- Open MPI 1.8.5
- CUDA 6.5
