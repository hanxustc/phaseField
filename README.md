## Requirements
- An MPI C/C++ compiler that supports the C++11 standard or newer  
  (e.g. `mpic++`)
- MATLAB R2023a or newer
- Python 3.9 or newer
- A unix environment (Linux recommended)

## File Structure
```text
phaseField/
├── ics/            Matlab files for generating Voronoi diagrams (generates poly_matrix.in)
├── test_poly/      Source and header files (requires poly_matrix.in and generates theta.out)
├── post/           Python and Matlab scripts for visualization and data processing (requires theta.out)

## Note
The parameter 'nsmooth' in 'parameter.c' is the interval at which a smoothing function is executed to remove wrinkles in the grain morphology for polycrystalline simulations, typically every 5000–10000 steps.
It should be turned off (set to a sufficiently large value) when computing energies using the quasi-1D grain slabs.
├── README.md
