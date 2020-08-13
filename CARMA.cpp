#include <stdlib.h>
#include <algorithm>
#include "CARMA.h"

using namespace std;




void CARMA(double** A, double** B, double** C, int* param, MPI_Comm comm)  //pass in &(double*)
{
	int rank, size, m, k, n;//, *colors, *parity;
	MPI_Comm_rank( comm, &rank );
	MPI_Comm_size( comm, &size );


	if (rank != 0) {
		MPI_Recv(param, 3, MPI_INT, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
	}
	m = param[0];
	k = param[1];
	n = param[2];
	if (rank != 0) {
		if (m==0||k==0||n==0){
			printf("%d: redundant proceesor\n", rank);
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

	int colors[level];//colors = (int*) malloc(sizeof(int)*level); not freed yet
	int parity[level];//parity = (int*) malloc(sizeof(int)*level);

	//recursively split matrix
	for (int i = log; i < level; i++) {
		temp = rank + (1<<i);
		if (temp < size) {

			int maxx = (m>n) ? max(m,k) : max(n,k);

			if (maxx == 1) {       //redundant processor
				colors[i] = 4;
				int new_param[3] = {0, 0, 0};
				MPI_Request req; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req);
				MPI_Wait(&req,MPI_STATUS_IGNORE);
				MPI_Request_free(&req);
				continue;
			}

			if (maxx == m) {
				parity[i] = m%2 == 0 ? 0:1;
				int m1 = m/2;
				m = m/2 + m%2;
				
				colors[i] = 1;
				int new_param[3] = {m1, k, n};
				//printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
				MPI_Request req1, req2, req3; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req1);

				//split A horizantally
				double* A_bot = *A + (m * k);

				MPI_Isend(A_bot, m1*k, MPI_DOUBLE, temp, 0, comm, &req2);
				MPI_Isend(*B, k*n, MPI_DOUBLE, temp, 0, comm, &req3);
				MPI_Wait(&req1,MPI_STATUS_IGNORE);
				MPI_Wait(&req2,MPI_STATUS_IGNORE);
				MPI_Wait(&req3,MPI_STATUS_IGNORE);
				MPI_Request_free(&req1);
				MPI_Request_free(&req2);
				MPI_Request_free(&req3);
				continue;
			}

			if (maxx == n) {
				int N = n;
				int n1 = n/2;
				n = N - n1;
				colors[i] = 3;
				parity[i] = N%2 == 0 ? 0:1;
				int new_param[3] = {m, k, n1};
				//printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
				MPI_Request req1, req2, req3; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req1);

				//split B vertically
				double* B_left = (double*) malloc(sizeof(double)*(k*n));
				double* B_right = (double*) malloc(sizeof(double)*(k*n1));
				//printf("copy begins: %d\n", rank);
				for(int j = 0; j < k; j++) {
					copy(*B + j*N, *B + j*N + n, B_left + j*n);
					copy(*B + j*N + n, *B + (j+1)*N, B_right + j*n1);
				}
				//printf("copy ends: %d\n", rank);
				MPI_Isend(*A, m*k, MPI_DOUBLE, temp, 0, comm, &req2);
				MPI_Isend(B_right, k*n1, MPI_DOUBLE, temp, 0, comm, &req3);
				MPI_Wait(&req1,MPI_STATUS_IGNORE);
				MPI_Wait(&req2,MPI_STATUS_IGNORE);
				MPI_Wait(&req3,MPI_STATUS_IGNORE);
				MPI_Request_free(&req1);
				MPI_Request_free(&req2);
				MPI_Request_free(&req3);

				free(B_right);
				free(*B);
				*B = B_left;
				continue;
			}

			if (maxx == k) {
				int K = k;
				int k1 = k/2;
				k = K - k1;
				colors[i] = 2;
				parity[i] = K%2 == 0 ? 0:1;
				int new_param[3] = {m, k1, n};
				//printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
				MPI_Request req1, req2, req3; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req1);

				//split A vertically
				double* A_left = (double*) malloc(sizeof(double)*(m*k));
				double* A_right = (double*) malloc(sizeof(double)*(m*k1));

				for(int j = 0; j < m; j++) {
					copy(*A + j*K, *A + j*K + k, A_left + j*k);
					copy(*A + j*K + k, *A + (j+1)*K, A_right + j*k1);
				}

				//split B horizantally
				double* B_bot = *B + (k * n);

				MPI_Isend(A_right, m*k1, MPI_DOUBLE, temp, 0, comm, &req2);
				MPI_Isend(B_bot, k1*n, MPI_DOUBLE, temp, 0, comm, &req3);
				MPI_Wait(&req1,MPI_STATUS_IGNORE);
				MPI_Wait(&req2,MPI_STATUS_IGNORE);
				MPI_Wait(&req3,MPI_STATUS_IGNORE);
				MPI_Request_free(&req1);
				MPI_Request_free(&req2);
				MPI_Request_free(&req3);

				free(A_right);
				free(*A);
				*A = A_left;
				continue;
			}
		}
	}
	//printf("rank %d begins\n", rank);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					m, n, k, 1, *A, k, *B, n, 0, *C, n);
/*
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			for (int x = 0; x < k; x++) {
				(*C)[i*n + j] += ((*A)[i*k + x]) * ((*B)[x*n + j]);
			}
		}
	}
*/
	//printf("rank %d finishes\n", rank);
	for (int i = level - 1; i >= log; i--) {
		temp = rank + (1<<i);

		if (temp < size) {

		//printf("rank %d receiving from %d\n",rank, temp);
			if (colors[i] == 4) {
				continue;
			}

			if (colors[i] == 3) {
				int n2 = parity[i] == 0 ? n:n-1;
				int N = n + n2;
				double* C_right = (double*) malloc(sizeof(double)*(m*n2));
				double* new_C = (double*) malloc(sizeof(double)*(m*N));
				//if (temp == 10) printf("!!!!%d: %d,%d\n",rank,m,n);
				MPI_Recv(C_right, m*n2, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
				for(int j = 0; j < m; j++) {
					copy(*C + n*j, *C + (j+1)*n, new_C + j*N);
					copy(C_right + n2*j, C_right + (j+1)*n2, new_C + j*N + n);
				}
				free(C_right);
				free(*C);
				*C = new_C;
				n = N;
		//printf("rank %d roger %d: %d\n",rank, temp, colors[i]);
				continue;
			}

			if (colors[i] == 1) {
				int m2 = parity[i]==0 ? m:m-1;
				int M = m + m2; 
				double* C_bot = (double*) malloc(sizeof(double)*(m2*n));
				double* new_C = (double*) malloc(sizeof(double)*(M*n));
				//if (temp == 10) printf("!!!!%d: %d,%d\n",rank,m,n);
		
				MPI_Recv(C_bot, m2*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
				copy(*C, *C + m*n, new_C);
				copy(C_bot, C_bot + m2*n, new_C + m*n);
				free(C_bot);
				free(*C);
				*C = new_C;
				m = M;
		//printf("rank %d roger %d: %d\n",rank, temp, colors[i]);
				continue;
			}

			if (colors[i] == 2) {

				double* new_C = (double*) malloc(sizeof(double)*(m*n));
				//if (temp == 10) printf("!!!!%d: %d,%d\n",rank,m,n);

				//MPI_Reduce(*C, new_C, m*n, MPI_DOUBLE, MPI_SUM, rank, comm01);
				MPI_Recv(new_C, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
/*
				for (int j = 0; j < m*n; j++) {
					(*C)[j] += new_C[j];
				}
*/
				vdAdd(m*n, *C, new_C, *C);
				free(new_C);

				k = parity[i]==1 ? k*2-1:2*k;
		//printf("rank %d roger %d: %d\n",rank, temp, colors[i]);
				continue;

			}

			printf("somthing wrong\n");
			return;
		}
	}

	if (rank != 0) {
		temp = rank - (1 << (log-1));
		//printf("rank %d sending back to %d\n",rank, temp);
		MPI_Send(*C, m*n, MPI_DOUBLE, temp, 0, comm);
		free(*C);
		free(*A);
		free(*B);
	}

}
