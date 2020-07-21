#include <mpi.h>
#include "mkl.h"
#include <iostream>
#include "CARMA.h"


int main(int argc, char* argv[]) {
    int result = 0;


    MPI_Init(&argc, &argv);

    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double *A;
    double *B;
    double *C;
    int *param;
    double expected_C[4*4];

    if (rank == 0) {
        double temp[4*4] = {10., -1., 2., 0.,
            -1., 11., -1., 3.,
            2., -1., 10., -1.,
            0.0, 3., -1., 8.};
        double temp2[4*4] = {27., 4., -8., 3.,
            -1., -6., -12., 3.,
            -9., 13., 37., 0.,
            6., -3., -1., 20.};
        A = (double*)malloc(sizeof(double)*16);
        B = (double*)malloc(sizeof(double)*16);
        C = (double*)malloc(sizeof(double)*16);
        std::copy(temp, temp+16, A);
        std::copy(temp2, temp2+16, B);


        int paramm[3] = {4,4,4};
        param = paramm;
    }


    CARMA(&A, &B, &C, param, MPI_COMM_WORLD);
    if(rank == 0) {
        free(A);
        free(B);
        free(C);
    }



    MPI_Finalize();


    return 0;



}