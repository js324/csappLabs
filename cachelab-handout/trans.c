/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    //32 sets, 1 lne/direct mapped, 32 byte block (8 elements)
    //gonna have blocks of B, look at A one sliver/row at a time
    if (M == 32) {
        for (int BCol = 0; BCol < 32; BCol+=8) { //big incrementer over B (goes over whole array)
            for (int BRow = 0; BRow < 32; BRow += 8) {
                for (int i = BCol; i < BCol+8; i++) { //incrementer to go row by row of block aka down the column 
                    int e1, e2, e3, e4, e5, e6, e7, e8;
                    e1 = A[i][BRow];
                    e2 = A[i][BRow+1];
                    e3 = A[i][BRow+2];
                    e4 = A[i][BRow+3];
                    e5 = A[i][BRow+4];
                    e6 = A[i][BRow+5];
                    e7 = A[i][BRow+6];
                    e8 = A[i][BRow+7];
                    B[BRow][i] = e1;
                    B[BRow+1][i] = e2;
                    B[BRow+2][i] = e3;
                    B[BRow+3][i] = e4;
                    B[BRow+4][i] = e5;
                    B[BRow+5][i] = e6;
                    B[BRow+6][i] = e7;
                    B[BRow+7][i] = e8;
                    //loop unroll of 8 because one block can store 8 elements so we load A twice as much if we increment by 4 intsead
                    //we can store 32 "blocks"/subarrays of length 4 of B at once
                        //but we because we are limited on variables (12 max) and block size, we should update 8 blocks at a time
                        //though doesn't really matter (just means we have to load A and reload one of B's lines 4x more if +8 vs +32)
                }
            }
        }
    }
    else if (M == 64) {
        // int cnter;
        int e1, e2, e3, e4, e5, e6, e7, e8;
        for (int BRow = 0; BRow < 64; BRow += 8) {    
            for (int BCol = 0; BCol < 64; BCol+=8) { //big incrementer over B (goes over whole array)
                for (int i = BCol; i < BCol+4; i++) { //incrementer to go row by row of block aka down the column 
                    
                    //takes sliver of 8 elements from row in A
                    //why above ^ approach won't work again
                        //Due to the expanded size specifically at 64 elements, there is now thrasing between elements 4 spaces apart 
                            //The big difference is that at 64 elements (unlike 56, 48, 32, etc.), there is a large enough gap between the addrs of rows that it affects the tag bits
                            //So essentially difference in address between row 0 and row 3 is 1024 bytes or 0x400 which affects the tag bit and will map to same set
                                //In previous element sizes, all the set bits of the 8 columns/rows we looked at in B were unique, hence why it worked
                    e1 = A[i][BRow];
                    e2 = A[i][BRow+1];
                    e3 = A[i][BRow+2];
                    e4 = A[i][BRow+3]; 
                    e5 = A[i][BRow+4];
                    e6 = A[i][BRow+5];
                    e7 = A[i][BRow+6];
                    e8 = A[i][BRow+7];
                    B[BRow][i] = e1;
                    B[BRow+1][i] = e2;
                    B[BRow+2][i] = e3;
                    B[BRow+3][i] = e4;
                    B[BRow][i+4] = e5;
                    B[BRow+1][i+4] = e6;
                    B[BRow+2][i+4] = e7;
                    B[BRow+3][i+4] = e8;
                }
                for (int i = BRow; i < BRow+4; i++) { 
                    //Think of B now as a 8x8 block with 4x4 quadrants
                    //First upper left quadra of B is already set
                    
                    //save line from upper right quad
                    e1 = B[i][BCol+4];
                    e2 = B[i][BCol+5];
                    e3 = B[i][BCol+6];
                    e4 = B[i][BCol+7];
                    
                    //2 things: Sets same upper right quad line and also loads next 4 rows of A for other loop, 
                    B[i][BCol+4] = A[BCol+4][i];
                    B[i][BCol+5] = A[BCol+5][i];
                    B[i][BCol+6] = A[BCol+6][i];
                    B[i][BCol+7] = A[BCol+7][i];

                    //Set line of bottom left quad 
                    B[i+4][BCol] = e1;
                    B[i+4][BCol+1] = e2;
                    B[i+4][BCol+2] = e3;
                    B[i+4][BCol+3] = e4;
                    
                    //Finish with bottom right quad by taking from second block of A that we pulled from above ^
                    B[i+4][BCol+4] = A[BCol+4][i+4];
                    B[i+4][BCol+5] = A[BCol+5][i+4];
                    B[i+4][BCol+6] = A[BCol+6][i+4];
                    B[i+4][BCol+7] = A[BCol+7][i+4];
                    
                }
                // printf("Curr Bcol: %d BRow: %d\n", BCol, BRow);
                // printf("Curr B: %d\n", B[2][0]);
                
            }
        }
    }
    else if (M == 61) { //61 x 67 or other
        //5x7            -> 7 x 5
        //x x x x x x x     x x x Y x 
        //x x x x x x x     x x x Y x 
        //x x x x x x x     x x x Y x 4x4 is done
        //x x x x x x x     Y Y Y Y x 
        //x x x x x x x     x x x x x 
        //                  x x x x x
        //                  x x x x x

        //INPUT ARRAY is 67 x 61
        for (int BCol = 0; BCol < 56; BCol+=8) { //big incrementer over B (goes over whole array)
            for (int BRow = 0; BRow < 56; BRow += 8) {
                for (int i = BCol; i < BCol+8; i++) { //incrementer to go row by row of block aka down the column of A 
                    int e1, e2, e3, e4, e5, e6, e7, e8;
                    e1 = A[i][BRow];
                    e2 = A[i][BRow+1];
                    e3 = A[i][BRow+2];
                    e4 = A[i][BRow+3];
                    e5 = A[i][BRow+4];
                    e6 = A[i][BRow+5];
                    e7 = A[i][BRow+6];
                    e8 = A[i][BRow+7];
                    B[BRow][i] = e1;
                    B[BRow+1][i] = e2;
                    B[BRow+2][i] = e3;
                    B[BRow+3][i] = e4;
                    B[BRow+4][i] = e5;
                    B[BRow+5][i] = e6;
                    B[BRow+6][i] = e7;
                    B[BRow+7][i] = e8;
                    
                }
            }
        }
        //1-56x56-61 needs to be done
        //56-67x61 needs to be done

        //this does rows 1-56 x columns 56-60
        for (int i = 0; i < 56; i++) {
            for (int j = 56; j < 61; j+=5) {
                int e1, e2, e3, e4, e5;
                // e5, e6, e7, e8;
                e1 = A[i][j];
                e2 = A[i][j+1];
                e3 = A[i][j+2];
                e4 = A[i][j+3];
                e5 = A[i][j+4];
                // e6 = A[i][j+5];
                // e7 = A[i][j+6];
                // e8 = A[i][j+7];
                //56,57-63
                B[j][i] = e1;
                B[j+1][i] = e2;
                B[j+2][i] = e3;
                B[j+3][i] = e4;
                B[j+4][i] = e5;
                // B[j+5][i] = e6;
                // B[j+6][i] = e7;
                // B[j+7][i] = e8;
            }
        }

        //this gets rows 56-64 and its first 56 elements
        for (int j = 0; j < 56; j+=8) {
            for (int i = 56; i < 64; i++) {
                int e1, e2, e3, e4, e5, e6, e7, e8;
                e1 = A[i][j];
                e2 = A[i][j+1];
                e3 = A[i][j+2];
                e4 = A[i][j+3];
                e5 = A[i][j+4];
                e6 = A[i][j+5];
                e7 = A[i][j+6];
                e8 = A[i][j+7];
                //56,57-63
                B[j][i] = e1;
                B[j+1][i] = e2;
                B[j+2][i] = e3;
                B[j+3][i] = e4;
                B[j+4][i] = e5;
                B[j+5][i] = e6;
                B[j+6][i] = e7;
                B[j+7][i] = e8;
            }
        }
        //finishes off rows 56-64 and does cols 56-61
        for (int i = 56; i < 64; i++) {
            int j = 56;
            //this doesnt seem to affect misses at all, should double check
            int e1, e2, e3, e4, e5;
            e1 = A[i][j];
            e2 = A[i][j+1];
            e3 = A[i][j+2];
            e4 = A[i][j+3];
            e5 = A[i][j+4];

            B[j][i] = e1;
            B[j+1][i] = e2;
            B[j+2][i] = e3;
            B[j+3][i] = e4;
            B[j+4][i] = e5;
            
        }
        //does last rows 65-67 and first 56 columns
        for (int j = 0; j < 56; j+= 8) {
            for (int i = 64; i < 67; i++) {
                int e1, e2, e3, e4, e5, e6, e7, e8;
                e1 = A[i][j];
                e2 = A[i][j+1];
                e3 = A[i][j+2];
                e4 = A[i][j+3];
                e5 = A[i][j+4];
                e6 = A[i][j+5];
                e7 = A[i][j+6];
                e8 = A[i][j+7];
                //56,57-63
                B[j][i] = e1;
                B[j+1][i] = e2;
                B[j+2][i] = e3;
                B[j+3][i] = e4;
                B[j+4][i] = e5;
                B[j+5][i] = e6;
                B[j+6][i] = e7;
                B[j+7][i] = e8;
            }
        }
        for (int i = 64; i < 67; i++) {
            for (int j = 56; j < 61; j++) {
                B[j][i] = A[i][j];
            }
        }
    }
}
//FINAL NOTES ABOUT ABOVE:
    //Definitely not the cleanest, final result is 1997 misses which is full points
        //Possibly missed optimizations, need to look at comments
    //Doesn't really break the more than 12 int rule, in any given scope there is only 12 ints (though I assumed int parameters didnt coun't)
        //Loop unrolling not necessary, could have just added extra loop with no effect on caching (spatial locality)



/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int v1,v0;
    int blockForCol, blockForRow, row, col;
    for (blockForCol = 0; blockForCol < M; blockForCol += 16) {
            for (blockForRow = 0; blockForRow < N; blockForRow += 16) {
                for (row = blockForRow; (row < N) && (row < blockForRow + 16); row++) {
                    for (col = blockForCol; (col < M) && (col < blockForCol + 16); col++) {
                        if (row != col) {
                            B[col][row] = A[row][col];
                        } else {
                            v0 = A[row][col];
                            v1 = row;
                        }
                    }
                    if (blockForRow == blockForCol) {
                        B[v1][v1] = v0;
                    }
                }
            }
        }
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}