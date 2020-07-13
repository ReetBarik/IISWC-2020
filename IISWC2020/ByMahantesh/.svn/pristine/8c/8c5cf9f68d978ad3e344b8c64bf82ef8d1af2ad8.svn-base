// ***********************************************************************
//
//            Grappolo: A C++ library for graph clustering
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory
//
// ***********************************************************************
//
//       Copyright (2014) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#include "defs.h"

//Print out the linear arrangement (gap) statistics:
void computeMinLAScores(graph *G) {
    printf("Within computeMinLAScores()\n");
    long    NV        = G->numVertices;
    long    NE        = G->numEdges;
    long    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    long    minGap = NV, maxGap = 0;
    double  sum = 0, sum_sq = 0;
    
    for (long v = 0; v < NV; v++) {
        long adj1 = vtxPtr[v];
        long adj2 = vtxPtr[v+1];
        for (long j = adj1; j < adj2; j++) {
            long w = vtxInd[j].tail;
            if(w > v) { //Process each edge only once: upper triangle of adj matrix
                long gap = abs(w - v);
                //printf("Gap: %ld\n", gap);
                sum_sq += (double) (gap*gap);
                sum    += (double) (gap);
                if (gap > maxGap)
                    maxGap = gap;
                if ( gap < minGap )
                    minGap = gap;
            } //End of if()
        } //End of for (w)
    } //End of for(v)
    double average  = (double) sum / NE;
    double avg_sq   = (double) sum_sq / NE;
    double variance = fabs(avg_sq - (average*average));
    double std_dev  = sqrt(variance);
    
    printf("*******************************************\n");
    printf("Linear Arrangement (edge gap) Statistics :\n");
    printf("*******************************************\n");
    printf("Number of vertices   :  %ld\n", NV);
    printf("Number of edges      :  %ld\n", NE);
    printf("*******************************************\n");
    printf("Minimum Gap          :  %ld\n", minGap);
    printf("Average Gap          :  %lf\n", average);
    printf("Maximum Gap          :  %ld\n", maxGap);
    printf("Total   Gap          :  %g \n", sum);
    printf("Expected value of X^2:  %g\n", avg_sq);
    printf("Variance is          :  %g\n", variance);
    printf("Standard deviation   :  %g\n", std_dev);
    printf("*******************************************\n");
    
} //End of computeMinLAScores()

//Build a reordering scheme for vertices based on clustering
//Group all communities together and number the vertices contiguously
//N = Number of vertices
//C = Community assignments for each vertex stored in an order
//old2NewMap = Stores the output of this routine
void buildOld2NewMap(long N, long *C, long *commIndex) {
    printf("Within buildOld2NewMap(%ld) function...\n", N);
    printf("WARNING: Assumes that communities are numbered contiguously\n");
    assert(N > 0);
    //Compute number of communities:
    //Assume zero is a valid community id
    long nC=-1;
    bool isZero = false;
    bool isNegative = false;
    for(long i = 0; i < N; i++) {
        if(C[i] == 0)
            isZero = true; //Check if zero is a valid community
        if(C[i] < 0)
            isNegative = true; //Check if zero is a valid community
        if (C[i] > nC) {
            nC = C[i];
        }
    }
    printf("Largest community id observed: %ld\n", nC);
    if(isZero) {
        printf("Zero is a valid community id\n");
        nC++;
    }
    if(isNegative) {
        printf("Some vertices have not been assigned communities\n");
        nC++; //Place to store all the unassigned vertices
    }
    assert(nC>0);
    printf("Number of unique communities in C= %d\n", nC);
    
    //////////STEP 1: Create a CSR-like datastructure for communities in C
    long * commPtr = (long *) malloc ((nC+1) * sizeof(long)); assert(commPtr != 0);
    //long * commIndex = (long *) malloc (N * sizeof(long)); assert(commIndex != 0);
    long * commAdded = (long *) malloc (nC * sizeof(long)); assert(commAdded != 0);
    
    // Initialization
#pragma omp parallel for
    for(long i = 0; i < nC; i++) {
        commPtr[i] = 0;
        commAdded[i] = 0;
    }
    commPtr[nC] = 0;
    // Count the size of each community
#pragma omp parallel for
    for(long i = 0; i < N; i++) {
        if(C[i] < 0) { //A negative value
            __sync_fetch_and_add(&commPtr[nC],1); //Unassigned vertices
        } else { //A positive value
            if(isZero)
                __sync_fetch_and_add(&commPtr[C[i]+1],1); //Zero-based indexing
            else
                __sync_fetch_and_add(&commPtr[C[i]],1); //One-based indexing
        }
    }//End of for(i)
    //Prefix sum:
    for(long i=0; i<nC; i++) {
        commPtr[i+1] += commPtr[i];
    }
    //Group vertices with the same community in an order
#pragma omp parallel for
    for (long i=0; i<N; i++) {
        long tc = (long)C[i];
        if(tc < 0) { //A negative value
            tc = nC-1;
            long Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < N);
            commIndex[Where] = i; //The vertex id
        } else {
            if(!isZero)
                tc--; //Convert to zero based index
            long Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < N);
            commIndex[Where] = i; //The vertex id
        }
    }
    printf("Done building structure for C...\n");
    
    //////////STEP 2: Create the old2New map:
    //This step will now be handled outside the routine
    /*
     #pragma omp parallel for
     for (long i=0; i<N; i++) {
     old2NewMap[commIndex[i]] = i;
     //printf("(%ld) --> (%ld)\n", commIndex[i]+1, i+1);
     }
     */
    //Cleanup:
    free(commPtr); free(commAdded);
    //free(commIndex);
}//End of buildOld2NewMap()
