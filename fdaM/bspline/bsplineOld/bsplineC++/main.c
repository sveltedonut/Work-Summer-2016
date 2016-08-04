//
//  main.c
//  
//
//  Created by Harry Jiang on 2016-06-13.
//
//

#include <stdio.h>
#include <stdlib.h>
#include "Bspline.h"

int main() {
    int i,j;
    double *x = (double*)malloc(sizeof(double)*21);
    for (i = 0; i < 21; i++) {
        x[i] = i*0.1;
    }
    double breaks[] = {0, 2/3, 4/3, 2};
    double *basismat = (double*)calloc(21*4,sizeof(double));
    
    bspline(21, x, 4, breaks, 2, 1, basismat);
    for (i = 0; i < 21; i++) {
        for (j = 0; j < 4; j++) {
            printf("%f ", basismat[j*21 + i]);
        }
        printf("\n");
    }
    printf("\n");
    
    return 0;
}