## Requirements
- An MPI C/C++ compiler that supports the C++11 standard or newer  
  (e.g. `mpic++`)
- MATLAB R2023a or newer
- Python 3.9 or newer
- A unix-like environment (Linux or macOS recommended)

## File Structure
```text
phaseField/
├── ics/            Matlab files for generating Voronoi diagrams (generates poly_matrix.in)
├── test_poly/      Source and header files (requires poly_matrix.in and generates theta.out)
├── post/           Python and Matlab scripts for visualization and data processing (requires theta.out)
├── README.md
