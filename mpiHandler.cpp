#include "mpiHandler.h"

MPIHandler::MPIHandler(int& argc, char**& argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

MPIHandler::~MPIHandler() {
    MPI_Finalize();
}

const int MPIHandler::getRank() {
    return rank;
}
const int MPIHandler::getSize() {
    return size;
}