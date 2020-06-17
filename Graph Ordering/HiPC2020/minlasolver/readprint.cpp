#include "cmpfuncs.h"
#include <time.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <iostream>

using namespace std;

#define str10K 100000

StatInfo Params;

void dot_write(TGraphC& G, leda::string f) {
  file_ostream o(f);
  o << "strict graph \"disjoint_union( ,  )\" {" << endl;

  node v, w;
  edge e;
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	o << G[v].initial_id << " -- " <<G[w].initial_id << ";\n";
  }
  o << "}\n";
  o.close();
}

void dot_format_print(TGraphC& G, leda::string f) {
  file_ostream o(f);
  o << "strict graph \"disjoint_union( ,  )\" {" << endl;

  node v, w;
  edge e;
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	o << G[v].initial_id << " -- " <<G[w].initial_id << ";\n";
  }
  o << "}\n";
  o.close();
}

void load_metis_solution( TGraphC & G, leda::string title) {
  cerr << "Reading solution vector file" << endl;
  FILE * inFile;

  node v;
  
  inFile = fopen(title, "r");
  int sol;
  forall_nodes(v, G) {
  	fscanf(inFile, "%d\n", &sol);
	G[v].ArrId = sol + 1;
  }
  fclose(inFile);
}

void hmetis_format_print( TGraphC & G, leda::string title) {
  file_ostream o(title);
  
  node v,w; edge e;
  o << G.number_of_edges() << " " << G.number_of_nodes() << " " << "11" << endl;
  
  forall_edges(e, G) {
  	v = G.source(e); w = G.target(e);
	o << 1+G[e].w << " " << G[v].metis_initial_id << " " << G[w].metis_initial_id << endl;
  }
  
  
  forall_nodes(v, G)
    {		
		o << (int)G[v].w  << endl;
    }
 
  o.close();

}


void metis_format_print( TGraphC & G, leda::string title) {
  file_ostream o(title);
  
  node v,w; edge e;
  o << G.number_of_nodes() << " " << G.number_of_edges() << " " << "11" << endl;
  /*
  forall_edges(e, G) {
  	v = G.source(e); w = G.target(e);
	cout  << G[v].initial_id << " " << G[w].initial_id << endl;
  }
  */
  forall_nodes(v, G)
    {		
		o << (int)G[v].w  << " ";
		forall_adj_edges(e, v) {
			w = second_adj_for_edge(e, v, G);
			o << G[w].metis_initial_id << " " << 1+(int)G[e].w << " ";
		}
      o << endl;
    }
 
  o.close();

}

void graph_format_print( TGraphC & G, leda::string title)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w; edge e;
  o << G.number_of_nodes() << " " << G.number_of_edges() << endl;
  forall_nodes(v, G)
    {
		forall_adj_edges(e, v) {
			w = second_adj_for_edge(e, v, G);
			o << G[w].initial_id << " ";
		}
      o << endl;
    }
 
  o.close();
}

void snap_format_print( TGraphC & G, leda::string title)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w; edge e;
  o << "# Nodes: " << G.number_of_nodes() << "\tEdges: " << G.number_of_edges() << endl;
  o << "# SrcNId\tDstNId" << endl;
  forall_edges(e, G)
    {
      v = G.source(e);
      w = G.target(e);
      o << G[v].initial_id << "\t" << G[w].initial_id << endl;
      o << G[w].initial_id << "\t" << G[v].initial_id << endl;
    }
 
  o.close();
}

void rmf_format_print( TGraphC & G, leda::string title)
{
//  TMP_CMP_GRAPHC = &G;
//  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w; edge e;
  
  o << "p ghct " << G.number_of_nodes() << " " << G.number_of_edges() << endl;
  forall_edges(e, G)
    {
      v = G.source(e); w = G.target(e);
      o << "e " << G[v].initial_id << " " << G[w].initial_id << " " << G[e].w << endl;
    }
 
  o.close();
}

void urmf_format_print( TGraphC & G, leda::string title)
{
//  TMP_CMP_GRAPHC = &G;
//  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w; edge e;
  
  o << "p ghct " << G.number_of_nodes() << " " << G.number_of_edges() << endl;
  forall_edges(e, G)
    {
      v = G.source(e); w = G.target(e);
      o << "e " << G[v].initial_id << " " << G[w].initial_id << " " << G[e].w << endl;
    }
 
  o.close();
}

void elist_format_print( TGraphC & G, leda::string title)
{
//  TMP_CMP_GRAPHC = &G;
//  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w; edge e;
  
  forall_edges(e, G)
    {
      v = G.source(e); w = G.target(e);
      o << G[v].initial_id << " " << G[w].initial_id << endl;
    }
 
  o.close();
}

void matlab_sparse_laplacian(TGraphC & G, leda::string title) {
  file_ostream o(title);
  
  o << "L=sparse(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  
  node v,w;
  edge e;
  
  node_array<int> A(G);
  int i=0;
  forall_nodes(v, G) {
    i++;
    A[v] = i;
  }
  
  forall_edges(e, G) {
    v = G.source(e); w = G.target(e);
    o << "L(" << A[v] << ", " << A[w] << ") = " << -G[e].w << ";" << endl;
    o << "L(" << A[w] << ", " << A[v] << ") = " << -G[e].w << ";" << endl;
  }
  o << "L = L - diag(sum(L));" << endl;
 
  o.close();  
}

void scilab_sparse_laplacian(TGraphC & G, leda::string title) {
  file_ostream o(title);
  
  o << "L=speye(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  
  node v,w;
  edge e;
  
  node_array<int> A(G);
  int i=0;
  forall_nodes(v, G) {
    i++;
    A[v] = i;
  }
  
  forall_edges(e, G) {
    v = G.source(e); w = G.target(e);
    o << "L(" << A[v] << ", " << A[w] << ") = " << -G[e].w << ";" << endl;
    o << "L(" << A[w] << ", " << A[v] << ") = " << -G[e].w << ";" << endl;
  }
  o << "L = L - speye(" << G.number_of_nodes() << ", " << G.number_of_nodes() <<");" << endl;
  o << "L = L - diag(sum(L, \"r\"));" << endl;
 
  o.close();  
}

void row_graph_print( TGraphC & G, leda::string title)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w;
  forall_nodes(v, G)
    {
      o << G[v].ArrId << endl;
    }
 
  o.close();
}
void print_part_solution( TGraphC & G, leda::string title)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_initial_id);

  file_ostream o(title);

  node v,w;
  forall_nodes(v, G)
    {
      o << G[v].ArrId << endl;
	  cout << G[v].ArrId << " ";
    }
 cout << endl;
  o.close();
}
void matrix_print( TGraphC & G, leda::string title)
{
file_ostream o(title);
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
	node v, w;
	edge e;
	  int i=1;
	  forall_nodes(v, G) {
	  	G[v].initial_id = i;
		i++;
	  }
	  forall_edges(e, G) {
	  	v = G.source(e); w = G.target(e);
		o << G[v].initial_id << " " << G[w].initial_id << " " << 1 << endl;
		o << G[w].initial_id << " " << G[v].initial_id << " " << 1 << endl;
	  }
}

void graph_print( TGraphC & G, leda::string fname)
{
  cerr << "Saving output arrangment in " << fname << endl;
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  
  file_ostream o(fname);
  
  node v;

  forall_nodes(v, G)
    {
      o << G[v].initial_id << endl;
    }
  o.close();  
}
void graph_print_4netresp( TGraphC & G, leda::string title)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  
  node v; edge e;
  forall_nodes(v, G) {
  	G[v].initial_id = G[v].ArrId;
  }
	file_ostream o(title);
	o << G.number_of_nodes() << " " << G.number_of_edges() << endl;
	
	forall_nodes(v, G)
    {
		forall_adj_edges(e, v) {
		node w = second_adj_for_edge(e, v, G);
			o << G[w].initial_id << " ";
		}
		o << endl;
	}
}
double print_graph_info(TGraphC & G)
{
  cerr << "************* Level " << Params.CURRENT_LEVEL << ", nodes=" << G.number_of_nodes() << ", #edges=" << G.number_of_edges();
  node v;
  double f=0;
  int mind=G.number_of_nodes();
  int maxd=-1;
  forall_nodes(v, G)
    {
      f = f+G.degree(v);
      if(maxd<G.degree(v))
        maxd=G.degree(v);
      if(mind>G.degree(v))
        mind=G.degree(v);
    }
  cerr << ", AvDeg=" << f/(double)G.number_of_nodes() << ", MaxDeg=" << maxd << ", MinDeg=" << mind << " ****************" << endl;
  //exit(1);
  return f/(double)G.number_of_nodes();
}

void gml_format_print(TGraphC & G, leda::string f)
{
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [" << endl; //\nversion 2\nCreator \"safro\"\ndirected 0\n";
      //o << "node_style [\n name \"default_node_style\" \n style [\ngraphics [\nw 16.0\n h 16.0\n\n]\n]\n]\n" ;
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].initial_id << "\n";
          o <<"label \"" << G[v].initial_id << "\"\n";
	  /*
          if(G[v].status==seed)
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\nfill \"#cccccc\"\n]\n";
            }
          else
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\n]\n";
            }
*/
  //        o <<"LabelGraphics [\ntype \"text\"\n]\n]\n";
	  o <<"]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"]\n";
        }
      o << "]\n";

}


void gml_write(leda::string f, TGraphC & G)
{
  if(Params.gml_out==true)
    {
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [\nversion 2\nCreator \"safro\"\ndirected 0\n";
      o << "node_style [\n name \"default_node_style\" \n style [\ngraphics [\nw 16.0\n h 16.0\n\n]\n]\n]\n" ;
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].ArrId << "\n";
          o <<"label \"" << G[v].initial_id << "(w=" << G[v].w<< ")\"\n";
          if(G[v].status==seed)
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\nfill \"#cccccc\"\n]\n";
            }
          else
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\n]\n";
            }

          o <<"LabelGraphics [\ntype \"text\"\n]\n]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"]\n";
        }
      o << "]\n";
    }

}

void gml_write_with_iw(leda::string f, TGraphC & G)
{
  if(Params.gml_out==true)
    {
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [\nversion 2\nCreator \"safro\"\ndirected 0\n";
      o << "node_style [\n name \"default_node_style\" \n style [\ngraphics [\nw 16.0\n h 16.0\n\n]\n]\n]\n" ;
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].ArrId << "\n";
          o <<"label \"" << G[v].initial_id << "\"\n";
          if(G[v].status==seed)
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\nfill \"#ff0000\"\n]\n";
            }
          else
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\n]\n";
            }

          o <<"LabelGraphics [\ntype \"text\"\n]\n]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"label \"" << /*e_costs[e] << "," <<*/ G[e].iw << "\"\n";
          o <<"]\n";
        }
      o << "]\n";
    }
}
void gml_write_with_pw(leda::string f, TGraphC & G)
{
cerr << "no pw" << endl;
exit(1);
/*
  if(Params.gml_out==true)
    {
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [\nversion 2\nCreator \"safro\"\ndirected 0\n";
      o << "node_style [\n name \"default_node_style\" \n style [\ngraphics [\nw 16.0\n h 16.0\n\n]\n]\n]\n" ;
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].ArrId << "\n";
          o <<"label \"" << G[v].initial_id << "\"\n";
          if(G[v].status==seed)
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\nfill \"#ff0000\"\n]\n";
            }
          else
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\n]\n";
            }

          o <<"LabelGraphics [\ntype \"text\"\n]\n]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"label \"" <<  G[e].pw << "\"\n";
          o <<"]\n";
        }
      o << "]\n";
    }
	*/
}
void gml_write_with_w(leda::string f, TGraphC & G)
{
  if(Params.gml_out==true)
    {
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [\nversion 2\nCreator \"safro\"\ndirected 0\n";
      o << "node_style [\n name \"default_node_style\" \n style [\ngraphics [\nw 16.0\n h 16.0\n\n]\n]\n]\n" ;
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].ArrId << "\n";
          o <<"label \"" << G[v].initial_id << "\"\n";
          if(G[v].status==seed)
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\nfill \"#ff0000\"\n]\n";
            }
          else
            {
              o << "fill \"#ffff00\"" <<"\n";
              o << "graphics [\n]\n";
            }

          o <<"LabelGraphics [\ntype \"text\"\n]\n]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"label \"" << /*e_costs[e] << "," <<*/ G[e].w << "\"\n";
          o <<"]\n";
        }
      o << "]\n";
    }
}

void gml_weights_write(leda::string f, TGraphC & G,edge_array<double> & new_edge_costs)
{
  if(Params.gml_out==true)
    {
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [\ndirected 0\n";
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].ArrId << "\n";
          o <<"label \"" << G[v].ArrId << "; " << ((G[v].status==seed)?("C-"):(" ")) <<G[v].initial_id << "(" << G[v].S_value<< ")" << "\"\n]\n";

          //      o <<"LabelGraphics [\ntype ""text""\n]\n]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"label \""<< G[e].w <<"\"\n";
          o <<"]\n";
        }
      o << "]\n";
    }
}

void gml_ordered_write(leda::string f, TGraphC & G)
{
  if(Params.gml_out==true)
    {
      file_ostream o(f);
      node v;
      edge e;


      o << "graph [\ndirected 0\n";
      forall_nodes(v, G)
        {
          o <<"node [\n";
          o <<"id " << G[v].ArrId << "\n";
          o <<"label \"" << G[v].ArrId << ", init=" << G[v].initial_id << "\"\n]\n";

          //      o <<"LabelGraphics [\ntype ""text""\n]\n]\n";
        }

      forall_edges(e, G)
        {
          int s = G[G.source(e)].ArrId;
          int t = G[G.target(e)].ArrId;
          o <<"edge [\n";
          o <<"source " << s << "\n";
          o <<"target " << t << "\n";
          o <<"]\n";
        }
      o << "]\n";
    }
}

void read_part_solution(TGraphC & G, leda::string filename)
{
  cerr << "Reading solution vector file" << endl;
  FILE * inFile;

  node v;
  
  inFile = fopen(filename, "r");
  int sol;
  int i=1;
  forall_nodes(v, G) {
  	fscanf(inFile, "%d\n", &sol);
	G[v].part_sol = sol+1;
	
	if(G[v].initial_id!=i) {
		cerr << "bad initial id" << endl;
		exit(-1);
	}
	i++;
  }
  fclose(inFile);
  
  node w; edge e;
  double c = 0;
  forall_edges(e, G) {
  	v = G.source(e);
	w = G.target(e);
	
	if(G[v].part_sol!=G[w].part_sol)
		c++;
  }
  cerr << "solution " << c << endl;
}

void read_graph(TGraphC & G, leda::string filename)
{
  cerr << "Reading .graph file" << endl;
  int Nv, Ne;
  int i;
  
  G.clear();
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
  getline (inFile, line);
  while(line.find('%')!=std::string::npos) {
	getline (inFile, line);
  }
  
  istringstream input(line);
  int toknum = 0;
  int edgeformat = 0;
  
  while (!input.eof()) {
	input >> token;
	if(toknum==0) {
	  Nv = atoi(token.c_str());
	} else if(toknum==1) {
	  Ne = atoi(token.c_str());
	} else if(toknum == 2) {
	  if((token.length()==0)||(token=="0")) {
	    edgeformat = 0; // unweighted everything
	  }
	  else if(token=="11") {
	    edgeformat = 2; // weighted nodes, weighted edges
	  }
	  else if(atoi(token.c_str())!=0) {
		cerr  << "strange header in graph file -" << token << "- " << endl;
		exit(-1);
	  }
	}
	toknum++;
  }
  
  cerr << "first line: nodes " << Nv << " edges " << Ne << endl;
  
  for(i=0; i<Nv; i++)
  {
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	G.new_node(n);
  }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  
  for(int n1=1; n1<=Nv; n1++)
  {
	getline (inFile, line);
	
	if(line.find('%')==std::string::npos)
	{
	  if(edgeformat==0) {
		v1 = AN[n1];
		istringstream tokens(line);		
		while (!tokens.eof()) {
		  tokens >> token;		  
		  v2 = AN[atoi(token.c_str())];
		  if(is_edge(v1,v2,G)==nil)
		  {
			CEdge ne; ne.w = 1.0; ne.real_w = 1.0;
			G.new_edge(v1, v2, ne);             			
		  }
		}
	  }
	  else if(edgeformat==2) {
	    v1 = AN[n1];
	    istringstream tokens(line);		
	    tokens >> token;
	    G[v1].w = atof(token.c_str());
	    while (!tokens.eof()) {
	      tokens >> token; // reading node number  
	      v2 = AN[atoi(token.c_str())];
	      tokens >> token; // reading edge weights
	      double weight = atof(token.c_str());
	      if(is_edge(v1,v2,G)==nil) {
		CEdge ne; ne.w = weight; ne.real_w = weight;
		G.new_edge(v1, v2, ne);             			
	      }
	    }
	  }
	}
  }
  cerr << "Edges added" << endl;
  inFile.close();
  
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
	G[v].captured=false;
	G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
  
  cerr << "min M " << minM << endl;
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}

void read_mtx(TGraphC & G, leda::string filename)
{
  cerr << "Reading .mtx file" << endl;
  int Nv, Ne;
  int i;
  
  G.clear();
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
  getline (inFile, line);
  while(line.find('%')!=std::string::npos) {
	getline (inFile, line);
  }
  
  istringstream input(line);
  int toknum = 0;
  int edgeformat = 0;
  
  while (!input.eof()) {
	input >> token;
	if(toknum==0) {
	  Nv = atoi(token.c_str());
	} else if(toknum==1) {
	  if(Nv!=atoi(token.c_str())) {
	    cerr << "rectangular matrix?" << endl;
	    exit(-1);
	  }
	} else if(toknum==2) {
	  Ne = atoi(token.c_str());
	}	
	else if(toknum == 3) {
	  if(atoi(token.c_str())!=0) {
		cerr  << "strange header in graph file" << endl;
		exit(-1);
	  }
	}
	toknum++;
  }
  
  cerr << "first line: nodes " << Nv << " edges " << Ne << endl;
  
  for(i=0; i<Nv; i++)
  {
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	G.new_node(n);
  }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  edge e;
    while(!inFile.eof())//(!feof(inFile))
    {
      int n1, n2;
      double w;
      node v, v1, v2;

      getline (inFile, line);
      istringstream input(line);
      
      input >> token; n1 = atoi(token.c_str());
      input >> token; n2 = atoi(token.c_str()); 
      input >> token; w = atof(token.c_str()); 
      
      //fscanf(inFile, "%d %d %lf\n", &n1, &n2, &w);
      v1 = AN[n1];//get_node_by_num(G, n1);
      v2 = AN[n2];//get_node_by_num(G, n2);

      CEdge ne; ne.w = w;
      ne.real_w = w;

      if((v1!=v2)&&(w!=0)&&((e=is_edge(v1, v2, G))==nil)) {
	  G.new_edge(v1, v2, ne); 
      }
    }

  cerr << "Edges added" << endl;
  inFile.close();
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
	G[v].captured=false;
	G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
  
  cerr << "min M " << minM << endl;
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}


void read_banjo(TGraphC & G, leda::string filename)
{
  cerr << "Reading banjo file" << endl;
  int Nv, Ne;
  int i;
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
  getline (inFile, line);
  while(line.find("Network #1, score")!=std::string::npos) {
	getline (inFile, line);
  }
  
  istringstream input(line);
  int toknum = 0;
  int edgeformat = 0;
  
  while (!input.eof()) {
	input >> token;
	if(toknum==0) {
	  Nv = atoi(token.c_str());
	} else if(toknum==1) {
	  Ne = atoi(token.c_str());
	} else if(toknum == 2) {
	  if(atoi(token.c_str())!=0) {
		cerr  << "strange header in graph file" << endl;
		exit(-1);
	  }
	}
	toknum++;
  }
  
  cerr << "first line: nodes " << Nv << " edges " << Ne << endl;
  
  for(i=0; i<Nv; i++)
  {
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	G.new_node(n);
  }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  
  for(int n1=1; n1<=Nv; n1++)
  {
	getline (inFile, line);
	
	if(line.find('%')==std::string::npos)
	{
	  if(edgeformat==0) {
		v1 = AN[n1];
		istringstream tokens(line);		
		while (!tokens.eof()) {
		  tokens >> token;		  
		  v2 = AN[atoi(token.c_str())];
		  if(is_edge(v1,v2,G)==nil)
		  {
			CEdge ne; ne.w = 1.0; ne.real_w = 1.0;
			G.new_edge(v1, v2, ne);             			
		  }
		}
	  }
	}
  }
  cerr << "Edges added" << endl;
  inFile.close();
  
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
	G[v].captured=false;
	G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
  
  cerr << "min M " << minM << endl;
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}


void read_lemon(TGraphC & G, leda::string filename)
{
  cerr << "Reading lemon file" << endl;
  int Nv, Ne;
  int i;
  
  G.clear();
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
  getline (inFile, line);
  while(line.find('%')!=std::string::npos) {
	getline (inFile, line);
  }
  
  istringstream input(line);
  int toknum = 0;
  int edgeformat = 0;
  
  while (!input.eof()) {
	input >> token;
	if(toknum==1) {
	  Nv = atoi(token.c_str());
	} 
	toknum++;
  }
  
  cerr << "first line: nodes " << Nv << endl;
  
  for(i=0; i<Nv; i++)
  {
    getline (inFile, line);    
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	n.status = fine;
	n.contr_edges = 0;
	node u = G.new_node(n);
	
    istringstream input(line);
    int toknum = 0;
    while (!input.eof()) {
	input >> token;
	if(toknum==2) {
	  G[u].w = atof(token.c_str());
	} 
	else if(toknum==3) {
	  G[u].T_max = atof(token.c_str());
	}
	toknum++;
    }	
   // cerr << G[u].initial_id << " " << G[u].w << " " << G[u].T_max << endl;
  }
  
  getline (inFile, line);    
  
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  
  while(!inFile.eof()) {
	getline (inFile, line);
	//cerr << line << endl;
	if(line.find('%')==std::string::npos) {		
		istringstream tokens(line);		
		int toknum = 0;
		while (!tokens.eof()) {
		  tokens >> token;		  
		  if(toknum==0) {
		    v1 = AN[atoi(token.c_str())];
		  }
		  else if(toknum==1) {
		    v2 = AN[atoi(token.c_str())];
		  }
		  else if(toknum==2) {
		    double w = atof(token.c_str());
		    if(is_edge(v1,v2,G)==nil) {
			CEdge ne; ne.w = w; ne.real_w = w;
			//cerr << w << endl;
			G.new_edge(v1, v2, ne);             			
		    }		    
		  }
		  toknum++;
		}
	}
  }
  
  cerr << "Edges added" << endl;
  inFile.close();
  
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
      G[v].contr_edges = 0;
	G[v].captured=false;
	//G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
  
  cerr << "min M " << minM << endl;
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}

void print_lemon(TGraphC & G, leda::string filename)
{
  
  int Nv, Ne;
  int i;
  
  file_ostream o(filename);
  
  o << "*Vertices " << G.number_of_nodes() << endl;
  node v, w;
  forall_nodes(v, G) {
    o << G[v].initial_id << " \"" << G[v].initial_id << "\" " << G[v].w << " " << G[v].T_max << endl;
  }
  edge e;
  
  o << "*Edges" << endl;
  forall_edges(e, G) {
    v = G.source(e);
    w = G.target(e);
    
    o << G[v].initial_id << " " << G[w].initial_id << " " << G[e].w << endl;
  }
  
  o.close();
}


void store_graph_arrangement(TGraphC & G, leda::string filename)
{
  FILE * inFile;

  inFile = fopen(filename, "aw");

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_initial_id);

  node v;
  forall_nodes(v, G)
    fprintf(inFile, "%d ", G[v].ArrId);
  fprintf(inFile, "\n");

  fclose(inFile);
}

void read_arrangement(TGraphC & G, leda::string filename)
{
  FILE * inFile;
  char s1[STR_LEN_1K];

  inFile = fopen(filename, "r");

  int pos;
  node v = G.first_node();
  while(!feof(inFile))
    {
      fscanf(inFile, "%d\n", &pos);
      G[v].ArrId = pos;
      G[v].save_ArrId = pos;
      if(G.last_node()!=v)
        v = G.succ_node(v);
    }

  fclose(inFile);
}

void print_diff_save_ArrId_and_ArrId(TGraphC & G)
{
  double avd = 0;
  int mxd = 0;
  int mnd = 0;
  array<int> vdiff(G.number_of_nodes()/100);
  node v;
  forall_nodes(v, G)
    {
      int d = abs(G[v].ArrId - G[v].save_ArrId);
      avd+=d;
      if(d>mxd) mxd=d;
      if(d<mnd) mnd=d;
      vdiff[d/100]++;
    }
  avd = avd/(double)G.number_of_nodes();
  cerr << "MIN DIST = " << mnd << endl;
  cerr << "MAX DIST = " << mxd << endl;
  cerr << "AVG DIST = " << avd << endl;
  for(int i=0; i<20; i++)
    cerr << "DIST FROM " << 100*i << " TO " << 100*(i+1) << " # of NODES = " << vdiff[i] << endl;
}
/*
  void check_same_coordinates(TGraphC & G)
  {
  node v, w;
  int count  = 0;
  for(int i=0; i<1000; i++)
  {
  w = G.choose_node();
  forall_nodes(v, G)
  {
  if(fabs(G[v].interpol_coord-G[w].interpol_coord)<DOUBLE_DIFF)
  count++;
  }
  }
  cerr << "REAL COORDS CHECKING, SAMPLING=1000 : " << (double)count/1000.0 << endl;
  }
*/
void read_mtx_old(TGraphC & G, leda::string filename)
{
  ifstream  inFile;
  char s1[STR_LEN_1K];
  char s2[STR_LEN_1K];
  int Nv, Ne;
  int i;

  cerr << "reading mtx file" << endl;
  G.clear();

  std::string str;
  inFile.open(filename);
  //inFile.getline(s1, 1024);
  inFile >> Nv >> Ne;
  cerr << "Nodes: " << Nv << " Edges: " << Ne << endl;

  for(i=0; i<Nv; i++)
    {
      Cnode n;
      n.ArrId = i+1;
      n.initial_id = i+1;
      G.new_node(n);
    }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
    {
      AN[i] = x;
      x = G.succ_node(x);
    }

  cerr << "Nodes defined" << endl;

  for(int en = 0; en < Ne; en++) {
    int n1, n2;
    double w;
    node v, v1, v2;

    inFile >> n1 >> n2 >> w;
    //cout << n1 << " " << n2 << " " << w << endl;
    if((n1!=n2)&&(w!=0)) {
      v1 = AN[n1];//get_node_by_num(G, n1);
      v2 = AN[n2];//get_node_by_num(G, n2);

      CEdge ne; ne.w = fabs(w);
      G.new_edge(v1, v2, ne); // !!!!!!!!!!! should be w here
    }
  }
  cerr << "Edges added" << endl;
  inFile.close();


  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
    {
      G[v].captured=false;
      G[v].w = 1;
      G[v].status = fine;
      int a = rand(); int b = rand();
      G[v].M = ((double)random())/((double)RAND_MAX);//(dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b));
      if(G[v].M<minM) minM = G[v].M;
    }
  cerr << "min M " << minM << endl;
  //    {
  //      G[v].ArrId = T[i];
  //      i++;
  //    }
  SI += "nodes=" + i2s(G.number_of_nodes()) + ", edges=" + i2s(G.number_of_edges()) + "\n";
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
/*
  cout << G.number_of_nodes() << " " << G.number_of_edges() << endl;
  edge e;
  node w;
  forall_edges(e, G) {
    v = G.source(e);
    w = G.target(e);
    cout << G[v].initial_id << " " << G[w].initial_id << " " << G[e].w << endl;
  }
  exit(1);
  */
}
void read_rmfext(TGraphC & G, leda::string filename) {
  FILE * inFile;
  //char s1[STR_LEN_1K];
  //char s2[STR_LEN_1K];
  char * s1 = (char*) malloc(STR_LEN_1K);
  char * s2 = (char*) malloc(STR_LEN_1K);
  
  int Nv, Ne;
  int i;

  G.clear();

  inFile = fopen(filename, "r");
  fscanf(inFile, "%s %s %d %d\n", s1, s2, &Nv, &Ne);

  //  array<int> T(Nv);
  //  for(i=0; i<Nv; i++)
  //    T[i]=i+1;
  //  T.permute();

  for (i = 0; i < Nv; i++) {
    Cnode n;
    n.ArrId = i + 1;
    n.initial_id = i + 1;
    //		n.name = i2s(i+1);
    G.new_node(n);
  }
  array<node> AN(Nv + 1);
  node x = G.first_node();
  for (i = 1; i <= Nv; i++) {
    AN[i] = x;
    x = G.succ_node(x);
  }

  cerr << "Nodes defined" << endl;
  int ww = 0;
  int ed_ex = 0;
  while (!feof(inFile)) {
    int n1, n2, name1, name2;
    double w;
    node v, v1, v2;

    fscanf(inFile, "%s %d %d %lf\n", s1, &n1, &n2, &w);

    v1 = AN[n1];//get_node_by_num(G, n1);
    v2 = AN[n2];//get_node_by_num(G, n2);

    G[v1].name = i2s(name1);
    G[v2].name = i2s(name2);

    if((v1!=v2)&&(is_edge(v1, v2, G) == nil)) {
      CEdge ne;
      ne.w = fabs(w);
      ne.real_w = fabs(w);
      G.new_edge(v1, v2, ne); // !!!!!!!!!!! should be w here
    } else {
      ww++;
    }
  }
  cerr << "Edges added" << endl;
  cerr << "existed edges " << ed_ex << endl;
  fclose(inFile);

  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
    {
      G[v].captured=false;
      G[v].w = 1;
      G[v].status = fine;
      int a = rand(); int b = rand();
      G[v].M = ((double)random())/((double)RAND_MAX);//(dblmin(a, b)/dmax(a,b))*(dblmin(a, b)/dmax(a,b));
      if(G[v].M<minM) minM = G[v].M;
    }
  cerr << "min M " << minM << endl;
  //    {
  //      G[v].ArrId = T[i];
  //      i++;
  //    }
  SI += "nodes=" + i2s(G.number_of_nodes()) + ", edges=" + i2s(
							       G.number_of_edges()) + "\n";
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes()
       << ", # of edges = " << G.number_of_edges() << endl;

  if (Ne != G.number_of_edges() + ww) {
    cerr << "Error in # of edges : Ne=" << Ne << "; G(E)="
	 << G.number_of_edges() + ww << endl;
    exit(1);
  }
  Nv = G.number_of_nodes();
  Ne = G.number_of_edges();
  
  free(s1);
  free(s2);
}

void read_rmf(TGraphC & G, leda::string filename)
{
  FILE * inFile;
  //char s1[STR_LEN_1K];
  //char s2[STR_LEN_1K];
  char * s1 = (char*) malloc(STR_LEN_1K);
  char * s2 = (char*) malloc(STR_LEN_1K);
  int Nv, Ne;
  int i;

  G.clear();

  inFile = fopen(filename, "r");
  fscanf(inFile, "%s %s %d %d\n", s1, s2, &Nv, &Ne);
  bool directed = false;
  if((leda::string)s1=="d") {
    directed = true;
  }

  //  array<int> T(Nv);
  //  for(i=0; i<Nv; i++)
  //    T[i]=i+1;
  //  T.permute();

  for(i=0; i<Nv; i++)
    {
      Cnode n;
      n.ArrId = i+1;
      n.initial_id = i+1;
      n.name = i2s(i+1);
      G.new_node(n);
    }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
    {
      AN[i] = x;
      x = G.succ_node(x);
    }
  edge e;
  int dir_e = 0;
  cerr << "Nodes defined" << endl;
  //  int ww=0;
  int self_loops = 0;
  while(!feof(inFile))
    {
      //      if(ww%1000==0)
      //        cerr << ww << endl;

      int n1, n2;
      double w;
      node v, v1, v2;

      fscanf(inFile, "%s %d %d %lf\n", s1, &n1, &n2, &w);
      /*
	forall_nodes(v, G)
        {
	if(G[v].ArrId == n1)
	v1 = v;
	if(G[v].ArrId == n2)
	v2 = v;
        }
      */
      v1 = AN[n1];//get_node_by_num(G, n1);
      v2 = AN[n2];//get_node_by_num(G, n2);

      CEdge ne; ne.w = w;
      ne.real_w = w;

      if(v1==v2) {
		  self_loops++;
      }
      else if(directed==true) {
		if((e=is_edge(v1, v2, G))!=nil) {
			G[e].w+=w;
			G[e].real_w+=w;
			dir_e++;
			//cerr << n1 << " " << n2 << endl;
			//exit(1);
		}
		else
			G.new_edge(v1, v2, ne); 
      }
      else  if(directed==false) {
		  if((e=is_edge(v1, v2, G))==nil) 
			G.new_edge(v1, v2, ne); // !!!!!!!!!!! should be w here
	  }
	  else {
		  cerr << "strange format" << endl;
		  exit(1);
	  }
      //      ww++;
    }
  cerr << "Edges added" << endl;
  cerr << "self loops " << self_loops << endl;
  fclose(inFile);


  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
    {
      G[v].captured=false;
      G[v].w = 1;
      G[v].status = fine;
      int a = rand(); int b = rand();
      G[v].M = ((double)random())/((double)RAND_MAX);//(dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b));
      if(G[v].M<minM) minM = G[v].M;
    }
  cerr << "min M " << minM << endl;
  //    {
  //      G[v].ArrId = T[i];
  //      i++;
  //    }
  SI += "nodes=" + i2s(G.number_of_nodes()) + ", edges=" + i2s(G.number_of_edges()) + "\n";
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
  if(Ne!=G.number_of_edges())
    {
      cerr << "Error in # of edges : Ne=" << Ne <<"; G(E)=" << G.number_of_edges() << endl;
      //exit(1);
    }
  cerr << "directed double edges " << dir_e << endl;
/*
  forall_nodes(v, G)
    if((G.degree(v)==0))//||(G[v].initial_id==2961))
      G.del_node(v);
	  */
  /*
  node_array<int> compnum(G);
  int scc = STRONG_COMPONENTS(G, compnum);
  cerr << scc << " connected components" << endl;
  if(scc>1) {
    int cc_num;
    int cc_size = 0;
    int curr_size = 0;
    for(int cc=0; cc<scc-1; cc++) {
      forall_nodes(v, G) {
	if(compnum[v]==cc)
	  curr_size++;
      }
      if(curr_size>cc_size) {
	cc_num = cc;
	cc_size = curr_size;
      }
    }
    forall_nodes(v, G) {
      if(compnum[v]!=cc_num)
	G.del_node(v);
    }
    cerr << "small cc removed" << endl;
  }
  */
  /*
  int r=1;
  forall_nodes(v, G) {
    G[v].initial_id = r;
    r++;
  }
  */
  free(s1);
  free(s2);

}
void read_rmf2(TGraphC & G, leda::string filename)
{
  FILE * inFile;
  char s1[STR_LEN_1K];
  char s2[STR_LEN_1K];
  int Nv, Ne;
  int i;

  G.clear();

  inFile = fopen(filename, "r");
  fscanf(inFile, "%d %d\n", &Ne, &Nv);

  //  array<int> T(Nv);
  //  for(i=0; i<Nv; i++)
  //    T[i]=i+1;
  //  T.permute();
  cerr << Nv << "\t" << Ne << endl;
  for(i=0; i<Nv; i++)
    {
      Cnode n;
      n.ArrId = i+1;
      n.initial_id = i+1;
      G.new_node(n);
    }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
    {
      AN[i] = x;
      x = G.succ_node(x);
    }

  cerr << "Nodes defined" << endl;
  //  int ww=0;
  while(!feof(inFile))
    {
      //      if(ww%1000==0)
      //        cerr << ww << endl;

      int n1, n2;
      double w;
      node v, v1, v2;

      fscanf(inFile, " %lf %d %d\n", &w, &n1, &n2);
      cerr << n1 << "\t" << n2 << endl;
      /*
	forall_nodes(v, G)
        {
	if(G[v].ArrId == n1)
	v1 = v;
	if(G[v].ArrId == n2)
	v2 = v;
        }
      */
      v1 = AN[n1];//get_node_by_num(G, n1);
      v2 = AN[n2];//get_node_by_num(G, n2);

      CEdge ne; ne.w = w;
      G.new_edge(v1, v2, ne); // !!!!!!!!!!! should be w here
      //      ww++;
    }
  cerr << "Edges added" << endl;
  fclose(inFile);


  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
    {
      G[v].captured=false;
      G[v].w = 1;
      G[v].status = fine;
      int a = rand(); int b = rand();
      G[v].M = ((double)random())/((double)RAND_MAX);//(dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b));
      if(G[v].M<minM) minM = G[v].M;
    }
  cerr << "min M " << minM << endl;
  //    {
  //      G[v].ArrId = T[i];
  //      i++;
  //    }
  SI += "nodes=" + i2s(G.number_of_nodes()) + ", edges=" + i2s(G.number_of_edges()) + "\n";
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
  if(Ne!=G.number_of_edges())
    {
      cerr << "Error in # of edges : Ne=" << Ne <<"; G(E)=" << G.number_of_edges() << endl;
      exit(1);
    }
}

void read_ine(TGraphC & G, leda::string filename)
{
  FILE * inFile;
  //char s1[STR_LEN_1K];
  //char s2[STR_LEN_1K];
  char * s1 = (char*) malloc(STR_LEN_1K);
  char * s2 = (char*) malloc(STR_LEN_1K);
  int Nv, Ne;
  int i;

  G.clear();

  inFile = fopen(filename, "r");
  fscanf(inFile, "%d %d\n", &Nv, &Ne);

  bool directed = false;

  for(i=0; i<Nv; i++)
    {
      Cnode n;
      n.ArrId = i+1;
      n.initial_id = i+1;
      n.name = i2s(i+1);
      G.new_node(n);
    }
  array<node> AN(Nv+1);
  node x = G.first_node();
  
  int a1, a2, a3;
  for(i=1; i<=Nv; i++)
    {
      AN[i] = x;
      x = G.succ_node(x);
      fscanf(inFile, "%d\t%d\t%d\n", &a1, &a2, &a3);
    }
  edge e;
  int dir_e = 0;
  cerr << "Nodes defined" << endl;
  //  int ww=0;
  while(!feof(inFile))
    {
      //      if(ww%1000==0)
      //        cerr << ww << endl;

      int n1, n2;
      double w;
      node v, v1, v2;

      fscanf(inFile, "%d\t%d\t%lf\n", &n1, &n2, &w);
      n1++; n2++;
      //cerr << n1 << " " << n2 << " " << w << endl;
      /*
	forall_nodes(v, G)
        {
	if(G[v].ArrId == n1)
	v1 = v;
	if(G[v].ArrId == n2)
	v2 = v;
        }
      */
      v1 = AN[n1];//get_node_by_num(G, n1);
      v2 = AN[n2];//get_node_by_num(G, n2);

      CEdge ne; ne.w = w;
      ne.real_w = w;

      if(v1==v2) {
      }
      else if(directed==true) {
	if((e=is_edge(v1, v2, G))!=nil) {
	  G[e].w+=w;
	  G[e].real_w+=w;
	  dir_e++;
	}
	else
	  G.new_edge(v1, v2, ne); 
      }
      else
	G.new_edge(v1, v2, ne); // !!!!!!!!!!! should be w here
      //      ww++;
    }
  cerr << "Edges added" << endl;
  fclose(inFile);


  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
    {
      G[v].captured=false;
      G[v].w = 1;
      G[v].status = fine;
      int a = rand(); int b = rand();
      G[v].M = ((double)random())/((double)RAND_MAX);//(dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b));
      if(G[v].M<minM) minM = G[v].M;
    }
  cerr << "min M " << minM << endl;
  //    {
  //      G[v].ArrId = T[i];
  //      i++;
  //    }
  SI += "nodes=" + i2s(G.number_of_nodes()) + ", edges=" + i2s(G.number_of_edges()) + "\n";
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
  if(Ne!=G.number_of_edges())
    {
      cerr << "Error in # of edges : Ne=" << Ne <<"; G(E)=" << G.number_of_edges() << endl;
      //exit(1);
    }
  cerr << "directed double edges " << dir_e << endl;
/*
  forall_nodes(v, G)
    if((G.degree(v)==0))//||(G[v].initial_id==2961))
      G.del_node(v);
	  */
  /*
  node_array<int> compnum(G);
  int scc = STRONG_COMPONENTS(G, compnum);
  cerr << scc << " connected components" << endl;
  if(scc>1) {
    int cc_num;
    int cc_size = 0;
    int curr_size = 0;
    for(int cc=0; cc<scc-1; cc++) {
      forall_nodes(v, G) {
	if(compnum[v]==cc)
	  curr_size++;
      }
      if(curr_size>cc_size) {
	cc_num = cc;
	cc_size = curr_size;
      }
    }
    forall_nodes(v, G) {
      if(compnum[v]!=cc_num)
	G.del_node(v);
    }
    cerr << "small cc removed" << endl;
  }
  */
  /*
  int r=1;
  forall_nodes(v, G) {
    G[v].initial_id = r;
    r++;
  }
  */
  free(s1);
  free(s2);

}

void read_snap(TGraphC & G, leda::string filename)
{
  cerr << "Reading snap file" << endl;
  int Nv, Ne;
  int i;
  
  G.clear();
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
    
  leda::list<int> v1list;
  leda::list<int> v2list;

  Nv = 0;
  
  while(!inFile.eof()) {
    getline (inFile, line);
    if(line.find('#')==std::string::npos) {      
      istringstream input(line);
      while (!input.eof()) {
	input >> token;
	int v1 = atoi(token.c_str());
	input >> token;
	int v2 = atoi(token.c_str());
	v1list.push_back(v1+1);
	v2list.push_back(v2+1);
	if(v1+1>Nv) Nv = v1 + 1;
	if(v2+1>Nv) Nv = v2 + 1;
      }
    }
  }
  cerr << "total nodes: " << Nv << endl;
  
  for(i=0; i<Nv; i++)
  {
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	G.new_node(n);
  }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  
  list_item it1, it2; it2 = v2list.first();
  forall_items(it1, v1list) {
    node n1 = AN[v1list[it1]];
    node n2 = AN[v2list[it2]];
    if((is_edge(n1,n2,G)==nil)&&(n1!=n2)) {
	CEdge ne; ne.w = 1; ne.real_w = 1;
	G.new_edge(n1, n2, ne);             			
    }
    it2 = v2list.succ(it2);
  }
  
  cerr << "Edges added" << endl;

  inFile.close();
  
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
	G[v].captured=false;
	G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
    
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}


void read_elist(TGraphC & G, leda::string filename)
{
  cerr << "Reading snap file" << endl;
  int Nv, Ne;
  int i;
  
  G.clear();
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
    
  leda::list<int> v1list;
  leda::list<int> v2list;

  Nv = 0;
  
  while(!inFile.eof()) {
    getline (inFile, line);
    if(line.find('#')==std::string::npos) {      
      istringstream input(line);
      while (!input.eof()) {
	input >> token;
	int v1 = atoi(token.c_str());
	input >> token;
	int v2 = atoi(token.c_str());
	v1list.push_back(v1+1);
	v2list.push_back(v2+1);
	if(v1+1>Nv) Nv = v1+1;
	if(v2+1>Nv) Nv = v2+1;
      }
    }
  }
  cerr << "total nodes: " << Nv << endl;
  
  for(i=0; i<Nv; i++)
  {
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	G.new_node(n);
  }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  
  list_item it1, it2; it2 = v2list.first();
  forall_items(it1, v1list) {
    node n1 = AN[v1list[it1]];
    node n2 = AN[v2list[it2]];
    if((is_edge(n1,n2,G)==nil)&&(n1!=n2)) {
	CEdge ne; ne.w = 1; ne.real_w = 1;
	G.new_edge(n1, n2, ne);             			
    }
    it2 = v2list.succ(it2);
  }
  
  cerr << "Edges added" << endl;

  inFile.close();
  
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
	G[v].captured=false;
	G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
    
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}

void read_welist(TGraphC & G, leda::string filename)
{
  cerr << "Reading snap file" << endl;
  int Nv, Ne;
  int i;
  
  G.clear();
  
  ifstream inFile (filename);
  if (!inFile.is_open()) {
	cerr << "Unable to open file " << filename << endl;
	exit(-1);
  }
  
  std::string line, token;
    
  leda::list<int> v1list;
  leda::list<int> v2list;
  leda::list<double> wlist;

  Nv = 0;
  
  while(!inFile.eof()) {
    getline (inFile, line);
    if(line.find('#')==std::string::npos) {      
      istringstream input(line);
      while (!input.eof()) {
	input >> token;
	int v1 = atoi(token.c_str());
	input >> token;
	int v2 = atoi(token.c_str());
	input >> token;
	double w = atof(token.c_str());
	v1list.push_back(v1);
	v2list.push_back(v2);
	wlist.push_back(w);
	if(v1+1>Nv) Nv = v1+1;
	if(v2+1>Nv) Nv = v2+1;
      }
    }
  }
  cerr << "total nodes: " << Nv << endl;
  
  for(i=0; i<Nv; i++)
  {
	Cnode n;
	n.ArrId = i+1;
	n.initial_id = i+1;
	G.new_node(n);
  }
  array<node> AN(Nv+1);
  node x = G.first_node();
  for(i=1; i<=Nv; i++)
  {
	AN[i] = x;
	x = G.succ_node(x);
  }
  
  cerr << "Nodes defined" << endl;
  node v1, v2;
  
  list_item it1, it2, it3; it2 = v2list.first(); it3 = wlist.first();
  forall_items(it1, v1list) {
    node n1 = AN[v1list[it1]];
    node n2 = AN[v2list[it2]];
    double w = wlist[it3];
    edge e = is_edge(n1,n2,G);
    if((e==nil)&&(n1!=n2)) {
	CEdge ne; ne.w = w; ne.real_w = w;
	G.new_edge(n1, n2, ne);             			
    }
    else if(e!=nil) {
      G[e].w = xdmax(G[e].w, w);
    }
    it2 = v2list.succ(it2);
    it3 = wlist.succ(it3);
  }
  
  cerr << "Edges added" << endl;

  inFile.close();
  
  
  node v;
  //  i=0;
  double minM = 10;
  forall_nodes(v, G)
  {
	G[v].captured=false;
	G[v].w = 1;
	G[v].status = fine;
	int a = rand(); int b = rand();
	G[v].M = ((double)random())/((double)RAND_MAX);//pow((dmin(a, b)/dmax(a,b))*(dmin(a, b)/dmax(a,b)), 3);
	if(G[v].M<minM) minM = G[v].M;
  }
    
  cerr << "Graph readed, # of nodes = " << G.number_of_nodes() << ", # of edges = " << G.number_of_edges() << endl;
	
}


void load_graph(TGraphC & G, leda::string instance_file, leda::string grtype) {
  if(grtype(0,3)=="-rmf") {
    read_rmf(G, instance_file);
  }
  else if(grtype(0,5)=="-elist") {
    read_elist(G, instance_file);
  }  
  else if(grtype(0,6)=="-welist") {
    read_welist(G, instance_file);
  }  
  else if(grtype(0,5)=="-graph") {
    read_graph(G, instance_file);
  }
  else if(grtype(0,3)=="-mtx") {
    read_mtx(G, instance_file);
  }
  else if(grtype(0,5)=="-lemon") {
    read_lemon(G, instance_file);
  }
  else if(grtype(0,4)=="-snap") {
    read_snap(G, instance_file);
  }
  else {
    cerr << "bad input format" << endl;
    exit(-1);
  }
}

void save_graph(TGraphC & G, leda::string output_file, leda::string grtype) {
  if(grtype(0,3)=="-rmf") {
    rmf_format_print(G, output_file);
  }
  else if(grtype(0,4)=="-urmf") {
    urmf_format_print(G, output_file);
  }
  else if(grtype(0,5)=="-graph") {
    graph_format_print(G, output_file);
  }
  else if(grtype(0,5)=="-elist") {
    elist_format_print(G, output_file);
  }
  else if(grtype(0,5)=="-snap") {
    snap_format_print(G, output_file);
  }
  else if(grtype(0,3)=="-gml") {
    gml_format_print(G, output_file);
  }
  else if(grtype(0,3)=="-dot") {
    dot_format_print(G, output_file);
  }
  else if(grtype(0,4)=="-lmat") {
    matlab_sparse_laplacian(G, output_file);
  }
  else if(grtype(0,4)=="-lsci") {
    scilab_sparse_laplacian(G, output_file);
  }
  else if(grtype(0,5)=="-lemon") {
    cerr << "bad lemon output format" << endl;
  }
  else {
    cerr << "bad output format" << endl;
    exit(-1);
  }
}

