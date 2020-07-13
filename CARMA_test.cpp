#include <mpi.h>
#include <gtest/gtest.h>
#include "mkl.h"
#include <iostream>
#include "CARMA.h"

TEST(MPI_Test, CARMA) {


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
        A = temp;
        B = temp2;
        C = (double*)malloc(sizeof(double)*16);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    4, 4, 4, 1, A, 4, B, 4, 0, expected_C, 4);

        int paramm[3] = {4,4,4};
        param = paramm;
    }


    CARMA(&A, &B, &C, param, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < 16; i++) {
            EXPECT_NEAR(expected_C[i], C[i], 1e-3) << " element y[" << i <<
                "] is wrong:" << expected_C[i] << " " << C[i];
        }
        free(A);
        free(B);
        free(C);
    }
    


}

int main(int argc, char* argv[]) {
    int result = 0;

    testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);


    result = RUN_ALL_TESTS();


    MPI_Finalize();


    return 0;



}
