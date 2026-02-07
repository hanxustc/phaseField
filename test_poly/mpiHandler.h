#ifndef _mpi_h_
#define _mpi_h_

#include <mpi.h>

class MPIHandler {
    int rank;
    int size;
public:
    MPIHandler(int& argc, char**& argv);
    ~MPIHandler();
    const int getRank();
    const int getSize();
};

#endif