//#include "cmpfuncs.h"
//#include <stdarg.h>
#include "tools.h"

leda::string i2s(int a)
{
  char s[20];
  sprintf(s, "%d", a);
  return (leda::string)s;
}

leda::string d2s(double a)
{
  char s[20];
  sprintf(s, "%f", a);
  return (leda::string)s;
}

leda::string getTime ()
{
  char buffer[SIZE];
  time_t curtime;
  struct tm *loctime;

  curtime = time (NULL);

  loctime = localtime (&curtime);

  return (leda::string)asctime (loctime);
}

bool char_search(leda::string s, char c)
{
  for(int i=0; i<s.length(); i++)
    if(s[i]==c)
      return true;
  return false;
}

double s2d(leda::string s)
{

  int i;
  char* r;
  r = (char*)malloc(s.length());
  //  cerr << "len="<< s.length() <<endl;
  for(i=0; i<s.length(); i++)
    r[i]=s[i];
  //  cerr << "-- " <<r << endl;

  return atof(r);
}

int s2i(leda::string s)
{

  int i;
  char r[s.length()+1];
  //  r = (char*)malloc(sizeof(char)*s.length());
  //  cerr << "len="<< s.length() <<", strlen(r)=" << strlen(r) << endl;
  for(i=0; i<s.length(); i++)
    r[i]=s[i];
  r[i] = '\0';
  //  cerr << "-- " <<r << endl;

  return atoi(r);
}

edge is_edge(node v, node w,  TGraphC & G)
{
  edge e;
  if(G.degree(v)<G.degree(w))
    {
      forall_adj_edges(e, v)
        if((G.source(e)==w)||(G.target(e)==w))
          return e;
    }
  else
    {
      forall_adj_edges(e, w)
        if((G.source(e)==v)||(G.target(e)==v))
          return e;
    }
  return nil;
}

node second_adj_for_edge(edge e, node v, TGraphC & G)
{
  node u = (G.source(e)==v)?G.target(e):G.source(e);
  return u;
}

void MY_error(char *fmt, ...)
{
  va_list argp;
  fprintf(stderr, "error: ");
  va_start(argp, fmt);
  fprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");
}

double gen_prob()
{
  int rnd1 = ((int)rand())%1000;
  int rnd2 = ((int)rand())%1000;
  double r = (rnd1<=rnd2)?((double)rnd1/(double)rnd2):((double)rnd2/(double)rnd1);
  return r;
}

void densify_graph(TGraphC & G) {
  edge e, f;
  node v, w, u;
  list<node> l1, l2;
  list<edge> l3, l4;
  
  TMP_CMP_GRAPHC = &G;
  G.sort_edges(&cmp_erw);
  
  double av_rw = 0;
  forall_edges(e, G) {
	av_rw+=G[e].rw;
  }
  av_rw = av_rw / (double)G.number_of_edges();
  
  int i=0;
  forall_edges(e, G) {
	if(G[e].imaginary==false) {
	  v = G.source(e);
	  w = G.target(e);
	  
	  forall_adj_edges(f, v) {			
		u = second_adj_for_edge(f, v, G);
		edge g = is_edge(u, w, G);
		if((g==nil)&&(G[f].imaginary==false)&&(G[e].rw > 1.5*av_rw)&&(G[f].rw > 1.5*av_rw)&&(u!=w)) {
		  CEdge ne; ne.w = (G[e].w + G[f].w)/2.0; ne.rw = (G[e].rw + G[f].rw)/2.0;
		  ne.imaginary = true;
		  G.new_edge(w, u, ne);			
		  i++;	
		}
	  }
	  
	  forall_adj_edges(f, w) {			
		u = second_adj_for_edge(f, w, G);
		edge g = is_edge(u, v, G);
		if((g==nil)&&(G[f].imaginary==false)&&(G[e].rw > 1.5*av_rw)&&(G[f].rw > 1.5*av_rw)&&(u!=v)) {
		  CEdge ne; ne.w = (G[e].w + G[f].w)/2.0; ne.rw = (G[e].rw + G[f].rw)/2.0;
		  ne.imaginary = true;
		  G.new_edge(v, u, ne);			
		  i++;	
		}
	  }
	  
	}
  }
  cout << "densified edges "  << i << endl;
}

void sparsify_graph(TGraphC & G) {
	if (Params.sparsification_coeff_type==abs_coeff)
		cerr << "sparsification: absolute" << endl;
	else if (Params.sparsification_coeff_type==rel_coeff)
		cerr << "sparsification: relative" << endl;
	else
	{
		cerr << "sparsification: bad type" << endl;
		exit(-1);
	}
  edge e, f;
  node v, w, u;
  list<node> l1, l2;
  list<edge> l3, l4;
  
  TMP_CMP_GRAPHC = &G;
  G.sort_edges(&cmp_erw);
  
  double av_rw = 0;
  forall_edges(e, G) {
	 // node v, w;
	 // v = G.source(e); w = G.target(e);
	 // cerr << G[e].rw << " "  << G.degree(v) << " " << G.degree(w) << endl;
	av_rw+=G[e].rw;
  }
  av_rw = av_rw / (double)G.number_of_edges();
  
  int i=0;
  forall_edges(e, G) {
	if(((Params.sparsification_coeff_type==rel_coeff)&&(G[e].rw < Params.sparsification_coeff * av_rw)) ||
	   ((Params.sparsification_coeff_type==abs_coeff)&&(G[e].rw < Params.sparsification_coeff)))
	 {
	  G.hide_edge(e);
	  i++;
	}
  }
  cout << "sparsified edges "  << i << endl;
}

void squize_nodes(TGraphC & G) {
	cerr << "before squizing " << G.number_of_nodes() << endl;
	node v;
	forall_nodes(v, G) {
		if(G.degree(v)==0) {
			G.del_node(v);
		}
	}
	int i=1;
	forall_nodes(v, G) {
		G[v].initial_id = i;
		G[v].ArrId = i;
		G[v].name = i2s(i);
		i++;
	}
	cerr << "after squizing " << G.number_of_nodes() << endl;
}
