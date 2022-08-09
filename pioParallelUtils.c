#ifndef _PIO_PARALLEL_UTILS_H_
#define _PIO_PARALLEL_UTILS_H_

#include <stdlib.h>

#ifndef ENABLE_MPI

// No MPI
int pio_nprocs() { return 1; }
int pio_myrank() { return 0; }
void pio_init_comm(int id) {}
void pio_exit_comm() {}
void pio_Allgather_i64(int64_t myValue, int64_t *result) { *result = myValue; }

#else
#include "mpi.h"
#include <stdlib.h>

static int nprocs = 1;
static int myrank = 0;

static MPI_Comm myworld = NULL;

int pio_nprocs() { return nprocs; }

int pio_myrank() { return myrank; }

void pio_init_comm(int id) {
  myworld = MPI_Comm_f2c(id);
  (void) MPI_Comm_size(myworld, &nprocs);
  (void) MPI_Comm_rank(myworld, &myrank);
}

void pio_exit_comm() { myworld = NULL; }

int pio_Allgather_i64(int64_t myValue, int64_t *result) {
  return MPI_Allgather(&myValue, 1, MPI_INTEGER8, result, 1, MPI_INTEGER8,
                       myworld);
}
#endif

#endif
