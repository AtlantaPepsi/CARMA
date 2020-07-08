#include <mpi.h>
#include <gtest/gtest.h>
#include "mkl.h"
#include <iostream>
#include <CARMA.h>

TEST(MPI_Test, CARMA) {


    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double A[4*4];
    double B[4*4];
    double C[4*4];
    int param[3];
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

        double expected_C[4*4];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    4, 4, 4, 1, A, 4, B, 4, 0, expected_C, 4);

        param = {4,4,4};
    }


    CARMA(&A, &B, &C, param, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            EXPECT_NEAR(expected_C[i], C[i], 1e-3) << " element y[" << i <<
                "] is wrong:" << expected_C[i] << " " << C[i];
        }
    }

}

int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);


    result = RUN_ALL_TESTS();


    MPI_Finalize();


    return 0;



}
