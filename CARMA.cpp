#include <stdlib.h>
#include <algorithm>
#include "CARMA.h"

using namespace std;




void CARMA(double** A, double** B, double** C, int* param, MPI_Comm comm)  //pass in &(double*)
{
	int rank, size, m, k, n, *colors;
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

	colors = (int*) malloc(sizeof(int)*level);

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
				m /= 2;
				colors[i] = 1;
				int new_param[3] = {m, k, n};
				//printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
				MPI_Request req1, req2, req3; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req1);

				//split A horizantally
				double* A_bot = *A + (m * k);

				MPI_Isend(A_bot, m*k, MPI_DOUBLE, temp, 0, comm, &req2);
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
				n /= 2;
				colors[i] = 3;
				int new_param[3] = {m, k, n};
				//printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
				MPI_Request req1, req2, req3; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req1);

				//split B vertically
				double* B_left = (double*) malloc(sizeof(double)*(k*n));
				double* B_right = (double*) malloc(sizeof(double)*(k*n));
				//printf("copy begins: %d\n", rank);
				for(int j = 0; j < k; j++) {
					copy(*B + j*2*n, *B + j*2*n + n, B_left + j*n);
					copy(*B + j*2*n + n, *B + (j+1)*2*n, B_right + j*n);
				}
				//printf("copy ends: %d\n", rank);
				MPI_Isend(*A, m*k, MPI_DOUBLE, temp, 0, comm, &req2);
				MPI_Isend(B_right, k*n, MPI_DOUBLE, temp, 0, comm, &req3);
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
				k /= 2;
				colors[i] = 2;
				int new_param[3] = {m, k, n};
				//printf("spliting: source:%d target:%d m:%d n:%d k:%d\n", rank, temp, m, n, k);
				MPI_Request req1, req2, req3; //dummy
				MPI_Isend(new_param, 3, MPI_INT, temp, 0, comm, &req1);

				//split A vertically
				double* A_left = (double*) malloc(sizeof(double)*(m*k));
				double* A_right = (double*) malloc(sizeof(double)*(m*k));

				for(int j = 0; j < m; j++) {
					copy(*A + j*2*k, *A + j*2*k + k, A_left + j*k);
					copy(*A + j*2*k + k, *A + (j+1)*2*k, A_right + j*k);
				}

				//split B horizantally
				double* B_bot = *B + (k * n);

				MPI_Isend(A_right, m*k, MPI_DOUBLE, temp, 0, comm, &req2);
				MPI_Isend(B_bot, k*n, MPI_DOUBLE, temp, 0, comm, &req3);
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

		printf("rank %d receiving from %d\n",rank, temp);
			if (colors[i] == 4) {
				continue;
			}

			if (colors[i] == 3) {
				double* C_right = (double*) malloc(sizeof(double)*(m*n));
				double* new_C = (double*) malloc(sizeof(double)*(m*2*n));
				if (temp == 10) printf("!!!!%d: %d,%d\n",rank,m,n);
				MPI_Recv(C_right, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
				for(int j = 0; j < m; j++) {
					copy(*C + n*j, *C + (j+1)*n, new_C + j*(2*n));
					copy(C_right + n*j, C_right + (j+1)*n, new_C + j*(2*n)+n);
				}
				free(C_right);
				free(*C);
				*C = new_C;
				n *= 2;
		printf("rank %d roger %d: %d\n",rank, temp, colors[i]);
				continue;
			}

			if (colors[i] == 1) {
				double* C_bot = (double*) malloc(sizeof(double)*(m*n));
				double* new_C = (double*) malloc(sizeof(double)*(m*2*n));
				if (temp == 10) printf("!!!!%d: %d,%d\n",rank,m,n);
		
				MPI_Recv(C_bot, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
				copy(*C, *C + m*n, new_C);
				copy(C_bot, C_bot + m*n, new_C + m*n);
				free(C_bot);
				free(*C);
				*C = new_C;
				m *= 2;
		printf("rank %d roger %d: %d\n",rank, temp, colors[i]);
				continue;
			}

			if (colors[i] == 2) {

				double* new_C = (double*) malloc(sizeof(double)*(m*n));
				if (temp == 10) printf("!!!!%d: %d,%d\n",rank,m,n);

				//MPI_Reduce(*C, new_C, m*n, MPI_DOUBLE, MPI_SUM, rank, comm01);
				MPI_Recv(new_C, m*n, MPI_DOUBLE, temp, 0, comm, MPI_STATUS_IGNORE);
/*
				for (int j = 0; j < m*n; j++) {
					(*C)[j] += new_C[j];
				}
*/
				vdAdd(m*n, *C, new_C, *C);
				free(new_C);

				k *= 2;
		printf("rank %d roger %d: %d\n",rank, temp, colors[i]);
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
