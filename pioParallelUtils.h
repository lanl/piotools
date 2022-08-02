#ifndef _PIO_PARALLEL_UTILS_H_
#define _PIO_PARALLEL_UTILS_H_

extern int pio_nprocs();
extern int pio_myrank();
extern void pio_init_comm(int id);
extern void pio_exit_comm();
extern void pio_Allgather_i64(int64_t myValue, int64_t *result);

#endif



