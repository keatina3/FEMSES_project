#ifndef _UTILS_H_
#define _UTILS_H_

void output_csv(char *fname, Mesh &M, float *u, int order);
void assign_ptrs(float*** arr_ptr, float** arr, int n, int m);
void assign_ptrs(int*** arr_ptr, int** arr, int n, int m);

#endif
