#include <stdlib.h>
#include <algorithm>
#include <mpi.h>
#include "CARMA.h"

using namespace std;




void CARMA(double** A, double** B, double** C, int* param, MPI_Comm comm)  //pass in &(double*)
{
    int rank, size, m, k, n, *colors;
    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &size );

    //MPI_Group group, new_group;
    //MPI_Comm comm01;
    //MPI_Comm_group( comm, &group );
        
    if (rank == 0) printf("B: %p\n", B);  ///this is different from the value passed in by caller (CARMA_test.cpp 56) 
    
    if (rank != 0) {
        MPI_Recv(param, 3, MPI_INT, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
    }
    m = param[0];
    k = param[1];
    n = param[2];
    if (rank != 0) {
        if (m==0||k==0||n==0){
            printf("k=%d\n",k );
            printf("m=%d\n",m );
            printf("n=%d\n",n );
            printf("rank=%d\n",rank );
            return;
        }
        *A = (double*) malloc(sizeof(double)*(m*k));
        *B = (double*) malloc(sizeof(double)*(k*n));
        *C = (double*) malloc(sizeof(double)*(m*n));
        MPI_Recv(*A, m*k, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
        MPI_Recv(*B, k*n, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
    }

    //calculate address of next receiver
    int temp = rank;
    int log = 1;
    while (temp >>= 1)
        log++;
    if (rank == 0)
        log--;

    //calculate cell(log(rank))
    temp = size;
    int level = 0;
    while (temp >>= 1)
        level++;
    if ((size- (1<<level)) != 0)
        level++;
    printf("rank %d: log:%d level:%d\n", rank, log, level);
    colors = (int*) malloc(sizeof(int)*level);

    //recursively split matrix
    for (int i = log; i < level; i++) {
        temp = rank + (1<<i);
        if (temp < size) {

            int maxx = (m>n) ? max(m,k) : max(n,k);

            if (maxx == n) {
                n /= 2;
                colors[i] = 3;
                int new_param[3] = {m, k, n};
                printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
                MPI_Request req; //dummy
                MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req);

                //split B vertically
                double* B_left = (double*) malloc(sizeof(double)*(k*n));
                double* B_right = (double*) malloc(sizeof(double)*(k*n));
                printf("copy begins: %d\n", rank);
                for(int j = 0; j < k; j++) {
                    copy(*B + j*2*n, *B + j*2*n + n, B_left + j*n);
                    copy(*B + j*2*n + n, *B + (j+1)*2*n, B_right + j*n);
                }
                printf("copy ends: %d\n", rank);
                MPI_Isend(*A, m*k, MPI_DOUBLE, temp, 0, comm, &req);
                MPI_Isend(B_right, k*n, MPI_DOUBLE, temp, 0, comm, &req);
                MPI_Request_free(&req);

                free(B_right);
                free(*B);
                *B = B_left;
                continue;
            }

            if (maxx == m) {
                m /= 2;
                colors[i] = 1;
                int new_param[3] = {m, k, n};
                printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
                MPI_Request req; //dummy
                MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req);

                //split A horizantally
                double* A_bot = *A + (m * k);

                MPI_Isend(A_bot, m*k, MPI_DOUBLE, temp, 0, comm, &req);
                MPI_Isend(*B, k*n, MPI_DOUBLE, temp, 0, comm, &req);
                MPI_Request_free(&req);
                continue;
            }


            if (maxx == k) {
                k /= 2;
                colors[i] = 2;
                int new_param[3] = {m, k, n};
                printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
                MPI_Request req; //dummy
                MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req);

                //split A vertically
                double* A_left = (double*) malloc(sizeof(double)*(m*k));
                double* A_right = (double*) malloc(sizeof(double)*(m*k));
                printf("copy begins: %d\n", rank);
                for(int j = 0; j < m; j++) {
                    copy(*A + j*2*k, *A + j*2*k + k, A_left + j*k);
                    copy(*A + j*2*k + k, *A + (j+1)*2*k, A_right + j*k);
                }
                printf("copy ends: %d\n", rank);
                //split B horizantally
                double* B_bot = *B + (k * n);

                MPI_Isend(A_right, m*k, MPI_DOUBLE, temp, 0, comm, &req);
                MPI_Isend(B_bot, k*n, MPI_DOUBLE, temp, 0, comm, &req);
                MPI_Request_free(&req);

                free(A_right);
                free(*A);
                *A = A_left;
                continue;
            }
        }
    }
    printf("rank %d begins\n", rank);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int x = 0; x < k; x++) {
                (*C)[i*n + j] += ((*A)[i*k + x]) * ((*B)[x*n + j]);
            }
        }
    }
    printf("rank %d finishes\n", rank);
    for (int i = level - 1; i >= log; i--) {
        printf("rank %d receiving from %d\n",rank, temp);
        temp = rank + (1<<i);

        if (temp < size) {

            if (colors[i] == 3) {
                double* C_right = (double*) malloc(sizeof(double)*(m*n));
                double* new_C = (double*) malloc(sizeof(double)*(m*2*n));
                MPI_Recv(C_right, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
                for(int j = 0; j < m; j++) {
                    copy(*C + n*j, *C + (j+1)*n, new_C + j*(2*n));
                    copy(C_right + n*j, C_right + (j+1)*n, new_C + j*(2*n)+n);
                }
                free(C_right);
                free(*C);
                *C = new_C;
                n *= 2;
                continue;
            }

            if (colors[i] == 1) {
                double* C_bot = (double*) malloc(sizeof(double)*(m*n));
                double* new_C = (double*) malloc(sizeof(double)*(m*2*n));
                MPI_Recv(C_bot, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
                copy(*C, *C + m*n, new_C);
                copy(C_bot, C_bot + m*n, new_C + m*n);
                free(C_bot);
                free(*C);
                *C = new_C;
                m *= 2;
                continue;
            }

            if (colors[i] == 2) {

                //int subcomm[2] = {rank,temp};
                //MPI_Group_incl( group, 2, subcomm, &new_group );
                //MPI_Comm_create( comm, new_group, &comm01 );

                double* new_C = (double*) malloc(sizeof(double)*(m*n));

                //MPI_Reduce(*C, new_C, m*n, MPI_DOUBLE, MPI_SUM, rank, comm01);
                MPI_Recv(new_C, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
                for (int j = 0; j < m*n; j++) {
                    (*C)[j] += new_C[j];
                }
                free(new_C);

                k *= 2;
                //MPI_Group_free( &new_group );
                //MPI_Comm_free( &comm01 );
                continue;

            }

            printf("somthing wrong\n");
            return;
        }
    }

    if (rank != 0) {
        temp = rank - (1 << (log-1));
        printf("rank %d sending back to %d\n",rank, temp);
        MPI_Request req; //dummy
        MPI_Isend(*C, m*n, MPI_DOUBLE, temp, 0, comm, &req);
        MPI_Request_free(&req);
        free(*C);
        free(*A);
        free(*B);
    }

}
