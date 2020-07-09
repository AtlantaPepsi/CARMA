#include <mpi.h>
#include <iostream>
#include ".\CARMA.h"

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double *A;
    double *B;
    double *C;
    int *param;

    if (rank == 0) {
        double temp[4*4] = {10., -1., 2., 0.,
            -1., 11., -1., 3.,
            2., -1., 10., -1.,
            0.0, 3., -1., 8.};
        double temp2[4*4] = {27., 4., -8., 3.,
            -1., -6., -12., 3.,
            -9., 13., 37., 0.,
            6., -3., -1., 20.};
        A = temp;
        B = temp2;


        int paramm[3] = {4,4,4};
        param = paramm;
    }


    CARMA(&A, &B, &C, param, MPI_COMM_WORLD);

    return 0;


}