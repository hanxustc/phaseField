## Requirements
- An MPI C/C++ compiler that supports the C++11 standard or newer  
  (e.g. `mpic++`)
- MATLAB R2023a or newer
- Python 3.9 or newer
- A unix-like environment (Linux or macOS recommended)

## File Structure
```text
phaseField/
├── test_poly/      Source and header files
├── ics/            Matlab files for generating Voronoi diagrams
├── post/           Python and Matlab scripts for visualization and data processing (requires theta.out file computed in test_poly)
├── README.md
