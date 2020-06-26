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

#ifndef _DEFS_H
#define _DEFS_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <unistd.h> //For getopts()

#define MilanRealMax HUGE_VAL       // +INFINITY
#define MilanRealMin -MilanRealMax  // -INFINITY

#define PRINT_DETAILED_STATS_
//#define PRINT_TERSE_STATS_

typedef struct comm
{
  long size;
  double degree;
}Comm;

typedef struct
{
    long cid;       //community ID
    double Counter; //Weight relative to that community
} mapElement;

typedef struct /* the edge data structure */
{
  long head;
  long tail;
  //long weight;
  double weight;
} edge;

typedef struct /* the graph data structure */
{
  long numVertices;        /* Number of columns                                */
  long sVertices;          /* Number of rows: Bipartite graph: number of S vertices; T = N - S */
  long numEdges;           /* Each edge stored twice, but counted once        */
  long * edgeListPtrs;     /* start vertex of edge, sorted, primary key        */
  edge * edgeList;         /* end   vertex of edge, sorted, secondary key      */
} graph;

struct clustering_parameters 
{
  const char *inFile; //Input file
  int ftype;  //File type

  bool strongScaling; //Enable strong scaling
  bool output; //Printout the clustering data
  bool VF; //Vertex following turned on
  bool coloring; //If coloring is turned on
  bool replaceMap; //If map data structure is replaced with a vector
  bool threadsOpt;
  double C_thresh; //Threshold with coloring on
  long minGraphSize; //Min |V| to enable coloring
  double threshold; //Value of threshold
       
  clustering_parameters();
  void usage();    
  
  //Define in parseInputParameter.cpp
  bool parse(int argc, char *argv[]);
};


/////////////////// FUNCTION CALLS ////////////////////
void displayGraphCharacteristics(graph *G);
void displayGraph(graph *G);
void displayGraphEdgeList(graph *G);
void displayGraphEdgeList(graph *G, FILE* out);
//Graph Clustering (Community detection)
void runMultiPhaseLouvainAlgorithm(graph *G, long *C_orig, int coloring, int replaceMap, long minGraphSize,
                                   double threshold, double C_threshold, int numThreads, int threadsOpt);

double parallelLouvianMethod(graph *G, long *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr);

double parallelLouvianMethodScale(graph *G, long *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr);

double parallelLouvianMethodNoMap(graph *G, long *C, int nThreads, double Lower,
                                  double thresh, double *totTime, int *numItr);

double algoLouvainWithDistOneColoring(graph* G, long *C, int nThreads, int* color,
                                      int numColor, double Lower, double thresh, double *totTime, int *numItr);

double algoLouvainWithDistOneColoringNoMap(graph* G, long *C, int nThreads, int* color,
                                           int numColor, double Lower, double thresh, double *totTime, int *numItr);

//***  Clustering Utility Functions ***//
//Distance-1 Coloring
int algoDistanceOneVertexColoring(graph *G, int *vtxColor, int nThreads, double *totTime);
int algoDistanceOneVertexColoringOpt(graph *G, int *vtxColor, int nThreads, double *totTime);


void buildColorSize(long NVer, int *vtxColor, int numColors, long *colorSize);
void computeVariance(long NVer, int numColors, long *colorSize);
void equitableDistanceOneColorBased(graph *G, int *vtxColor, int numColors, long *colorSize,
                                    int nThreads, double *totTime, int type);

//Other 
inline void Visit(long v, long myCommunity, short *Visited, long *Volts, 
				  long* vtxPtr, edge* vtxInd, long *C);
long buildCommunityBasedOnVoltages(graph *G, long *Volts, long *C, long *Cvolts);
void buildNextLevelGraph(graph *Gin, graph *Gout, long *C, long numUniqueClusters);
long renumberClustersContiguously(long *C, long size);
double buildNextLevelGraphOpt(graph *Gin, graph *Gout, long *C, long numUniqueClusters, int nThreads);
//Vertex following functions:
long vertexFollowing(graph *G, long *C);
double buildNewGraphVF(graph *Gin, graph *Gout, long *C, long numUniqueClusters);

//***  Utility Functions ***//
void duplicateGivenGraph(graph *Gin, graph *Gout);
void writeEdgeListToFile(graph *G, FILE* out);

//Compute Gini coefficient
double computeGiniCoefficient(long *colorSize, int numColors); 

//Random Number Generation:
void generateRandomNumbers(double *RandVec, long size);

void displayGraph(graph *G);
void displayGraphCharacteristics(graph *G);
graph* convertDirected2Undirected(graph *G);

void segregateEdgesBasedOnVoltages(graph *G, long *Volts);
void writeGraphPajekFormat(graph* G, char * filename);
void writeGraphMatrixMarketFormatSymmetric(graph* G, char *filename);
void writeGraphPajekFormatWithNodeVolts(graph* G, long *Cvolts, char * filename);
void writeGraphBinaryFormat(graph* G, char * filename); //Binary (each edge once)
void writeGraphMetisSimpleFormat(graph* G, char *filename); //Metis format; no weights

//File parsers:
void parse_Dimacs9FormatDirectedNewD(graph * G, char *fileName);
long removeEdges(long NV, long NE, edge *edgeList);
void SortNodeEdgesByIndex2(long NV, edge *list1, edge *list2, long *ptrs);
void SortEdgesUndirected2(long NV, long NE, edge *list1, edge *list2, long *ptrs);

void loadMetisFileFormat(graph *G, const char* filename); //Metis (DIMACS#10)
void parse_MatrixMarket(graph * G, char *fileName);       //Matrix-Market
void parse_MatrixMarket_Sym_AsGraph(graph * G, char *fileName);

void parse_Dimacs1Format(graph * G, char *fileName);      //DIMACS#1 Challenge format
void parse_Dimacs9FormatDirectedNewD(graph * G, char *fileName); //DIMACS#9 Challenge format
void parse_PajekFormat(graph * G, char *fileName);        //Pajek format (each edge stored only once
void parse_PajekFormatUndirected(graph * G, char *fileName);
void parse_DoulbedEdgeList(graph * G, char *fileName);

void parse_EdgeListBinary(graph * G, char *fileName); //Binary: Each edge stored only once 
void parse_SNAP(graph * G, char *fileName);

void loadADJMatrixFormat(graph *G, const char* filename); //Metis (DIMACS#10)

//For reading power grid data
long* parse_MultiKvPowerGridGraph(graph * G, char *fileName); //Four-column format


//Graph partitioning with Metis:
void MetisGraphPartitioner( graph *G, long *VertexPartitioning, int numParts );


void writeGraphBinaryFormatNew(graph* G, char * filename, long weighted); 
void parse_EdgeListBinaryNew(graph * G, char *fileName);
#endif
