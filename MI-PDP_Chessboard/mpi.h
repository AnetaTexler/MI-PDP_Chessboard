#pragma once

/* Just fake header */

const int MPI_COMM_WORLD = 0;
const int MPI_INT = 0;
const int MPI_CHAR = 0;
const int MPI_ANY_SOURCE = 0;
const int MPI_ANY_TAG = 0;

struct MPI_Status
{
	int MPI_TAG;
	int MPI_SOURCE;
};

void MPI_Init(void *argc, char ***argv) {}
void MPI_Comm_rank(void *x, void *my_rank) {}
void MPI_Comm_size(void *y, void *p) {}

void MPI_Recv(void *message, int a, int * b, int *c, int d, int *q, MPI_Status *status) {}
void MPI_Send(void * a, int b, int * c, int d, int e, int * g) {}
void MPI_Iprobe(int * a, int b, int * c, int * d, MPI_Status *status) {}
void MPI_Finalize() {}
void MPI_Barrier(void * a) {}