#include <mpi.h>
#include <math.h>
#include <random>
#include <vector>
#include <gtest/gtest.h>
#include "mkl.h"
#include <iostream>
#include "CARMA.h"

void test(int m, int k, int n) {


    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double *A;
    double *B;
    double *C;
    int *param=(int*)malloc(sizeof(int)*3);
    double expected_C[m*n];

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
        std::vector<double> a(m*k);
        std::vector<double> b(k*n);
        A = (double*)malloc(sizeof(double)*m*k);
        B = (double*)malloc(sizeof(double)*k*n);
        C = (double*)malloc(sizeof(double)*m*n);
	std::copy(a.begin(), a.end(), A);
        std::copy(b.begin(), b.end(), B);
        
        std::normal_distribution<double> distribution(200.0, 20.0);

        std::default_random_engine generator;
        for (int i = 0;i < a.size();i++) {
            a[i] = distribution(generator);
        }
        for (int i = 0;i < b.size();i++) {
            b[i] = distribution(generator);
        }
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, 1, A, k, B, n, 0, expected_C, n);

        printf("A:%p B:%p C:%p\n",&A,&B,&C);
       //printf("%f\n",C[0]); 
        param[0] = m;
        param[1] = k;
        param[2] = n;
        
    }

//	printf("%p\n",&A);
    CARMA(&A, &B, &C, param, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < m*n; i++) {
            if(abs(expected_C[i]-C[i]) > 1e-3) {
		printf(" element y[%d] is wrong: %f, %f\n", i, expected_C[i], C[i]);
	    }
        }
        free(A);
        free(B);
        free(C);
    }
    
    free(param);

}

int main(int argc, char* argv[]) {
    int result = 0;

    testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);


    test(60,50,10);
    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
   printf("%d ok\n",rank); 
//printf("wtf %d\n",rank);

    MPI_Finalize();

	//printf("congrats!\n");
    return 0;



}
