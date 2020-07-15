#include <mpi.h>
#include <math.h>
#include <random>
#include <vector>
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
    double expected_C[20*50];

    if (rank == 0) {
        /*
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
        */
        std::vector<double> a(20*30);
        std::vector<double> b(30*50);
        A = (double*)malloc(sizeof(double)*600);
        B = (double*)malloc(sizeof(double)*1500);
        C = (double*)malloc(sizeof(double)*1000);
        std::normal_distribution<double> distribution(200.0, 20.0);

        std::default_random_engine generator;
        for (int i = 0;i<600;i++) {
            a[i] = distribution(generator);
        }
        for (int i = 0;i<1500;i++) {
            b[i] = distribution(generator);
        }
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    20, 50, 30, 1, A, 30, B, 50, 0, expected_C, 50);

        printf("B:%p\n",&B);
        
        int paramm[3] = {20,30,50};
        param = paramm;
    }


    CARMA(&A, &B, &C, param, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < 1000; i++) {
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
