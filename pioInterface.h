//========================================================================================
// (C) (or copyright) 2022. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S. Department of Energy/National Nuclear
// Security Administration. All rights in the program are reserved by Triad
// National Security, LLC, and the U.S. Department of Energy/National Nuclear
// Security Administration. The Government is granted for itself and others
// acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
// in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do
// so.
//========================================================================================

#ifndef _PIOINTERFACE_H_
#define _PIOINTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>

/**
 * C functions defined in pioInterface()
 * Version: 220509: initial definition
 **/

extern void pio_init(const int ID, const char *fname, const int verbose);
extern void pio_release(const int ID);

extern int pio_nCell(const int ID);
extern int pio_nDim(const int ID);
extern int pio_nMat(const int ID);

extern int64_t pio_width(const int ID, const char *var);
extern int64_t pio_length(const int ID, const char *var);

extern int64_t *pio_daughter(const int ID);
extern double *pio_center(const int ID, int index);

extern bool pio_exists(const int ID, const char *var, const int index);

extern double *pio_get_d(const int ID, const char *var_name, const int index);
extern double *pio_get_range_d(int ID, const char *var_name, int index,
                               int64_t start, int64_t nCount);
extern void pio_release_d(double *ptr);

extern int64_t *pio_get_i64(const int ID, const char *var_name,
                            const int index);
extern int64_t *pio_get_range_i64(int ID, const char *var_name, int index,
                                  int64_t start, int64_t nCount);
extern void pio_release_i64(int64_t *ptr);

extern double **pio_get_matvar_d(int ID, const char *field);
extern double **pio_get_range_matvar_d(int ID, const char *field,
                                       int64_t iStart, int64_t nCount);
extern void pio_release_2d_d(double **ptr, int n);

extern double *pio_get_matvar_index_d(int ID, const char *field, int index);
extern double *pio_get_matvar_index_range_d(int ID, const char *field,
                                            int index, int64_t iStart,
                                            int64_t nCount);
#ifdef __cplusplus
}
#endif

#endif
