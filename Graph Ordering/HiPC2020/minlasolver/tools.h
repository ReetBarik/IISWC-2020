#include "cmpfuncs.h"
#include <stdarg.h>

extern leda::string i2s(int a);

extern leda::string d2s(double a);

extern leda::string getTime ();

extern bool char_search(leda::string s, char c);

extern double s2d(leda::string s);

extern int s2i(leda::string s);

extern edge is_edge(node v, node w,  TGraphC & G);

extern node second_adj_for_edge(edge e, node v, TGraphC & G);

extern void MY_error(char *fmt, ...);

extern double gen_prob();

extern void densify_graph(TGraphC & G);

extern void sparsify_graph(TGraphC & G);