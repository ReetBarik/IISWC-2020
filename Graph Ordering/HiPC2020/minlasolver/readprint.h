/*
#include "cmpfuncs.h"
#include <time.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <iostream>

using namespace std;
*/
void dot_write(TGraphC& G, leda::string f);
void dot_format_print(TGraphC& G, leda::string f);
void load_graph(TGraphC & G, leda::string instance_file, leda::string grtype);
void save_graph(TGraphC & G, leda::string output_file, leda::string grtype);
void load_metis_solution( TGraphC & G, leda::string title);
void hmetis_format_print( TGraphC & G, leda::string title);
void metis_format_print( TGraphC & G, leda::string title);
void graph_format_print( TGraphC & G, leda::string title);
void rmf_format_print( TGraphC & G, leda::string title);
void elist_format_print( TGraphC & G, leda::string title);
void row_graph_print( TGraphC & G, leda::string title);
void print_part_solution( TGraphC & G, leda::string title);
void matrix_print( TGraphC & G, leda::string title);
void graph_print( TGraphC & G, leda::string title);
void graph_print_4netresp( TGraphC & G, leda::string title);
double print_graph_info(TGraphC & G);
void gml_format_print(TGraphC & G, leda::string f);
void gml_write(leda::string f, TGraphC & G);
void gml_write_with_iw(leda::string f, TGraphC & G);
void gml_write_with_pw(leda::string f, TGraphC & G);
void gml_write_with_w(leda::string f, TGraphC & G);
void gml_weights_write(leda::string f, TGraphC & G,edge_array<double> & new_edge_costs);
void gml_ordered_write(leda::string f, TGraphC & G);
void read_part_solution(TGraphC & G, leda::string filename);
void read_graph(TGraphC & G, leda::string filename);
void read_snap(TGraphC & G, leda::string filename);
void read_elist(TGraphC & G, leda::string filename);
void read_banjo(TGraphC & G, leda::string filename);
void read_lemon(TGraphC & G, leda::string filename);
void print_lemon(TGraphC & G, leda::string filename);
void store_graph_arrangement(TGraphC & G, leda::string filename);
void read_arrangement(TGraphC & G, leda::string filename);
void print_diff_save_ArrId_and_ArrId(TGraphC & G);
void read_mtx(TGraphC & G, leda::string filename);
void read_rmfext(TGraphC & G, leda::string filename);
void read_rmf(TGraphC & G, leda::string filename);
void read_rmf2(TGraphC & G, leda::string filename);
void read_ine(TGraphC & G, leda::string filename);
void matlab_sparse_laplacian(TGraphC & G, leda::string title);
void scilab_sparse_laplacian(TGraphC & G, leda::string title);
void read_snap(TGraphC & G, leda::string filename);
