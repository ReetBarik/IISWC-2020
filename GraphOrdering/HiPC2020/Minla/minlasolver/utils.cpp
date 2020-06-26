#include "cmpfuncs.h"
//#include "LEDA/basic_graph_alg.h"
#include <time.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <iostream>
using namespace std;
#define str10K 100000



void save_arrangement(TGraphC & G)
{
  node v;
  forall_nodes(v, G)
	G[v].save_ArrId=G[v].ArrId;
}

void local_save_arrangement(TGraphC & G, node u, node v)
{
  G[v].save_ArrId=G[v].ArrId;
  G[u].save_ArrId=G[u].ArrId;
}


void restore_arrangement(TGraphC & G/*, array<int> & MinArr*/)
{
  /*
  node v;
  int i=0;
  forall_nodes(v, G)
  {
	G[v].ArrId=MinArr[i];
	i++;
}
*/
  node v;
  forall_nodes(v, G)
	G[v].ArrId=G[v].save_ArrId;
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  //define_S_values(G);
}

void save_best_lcc_arrangement(TGraphC & G)
{
  cerr << "Bad call : save_best_lcc_arrangement" << endl;
  exit(1);
  /*
  node v;
  forall_nodes(v, G)
	G[v].bestLCC_ArrId=G[v].save_lcc_ArrId;
  */
}

void restore_best_lcc_arrangement(TGraphC & G)
{
  //cerr << "Bad call : restore_best_lcc_arrangement" << endl;
  //exit(1);
  
  node v;
  forall_nodes(v, G)
	G[v].ArrId=G[v].lcc_ArrId;
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  //define_S_values(G);
  
  
}


void save_and_clear_graph(TGraphC & G)
{
  /*
  leda::string fname = "/tmp/minla_" + i2s(S) + "_" + i2s(Params.prog_id);
  file_ostream o(fname);
  node v;
  forall_nodes(v, G)
	G[v].save(o);
  edge e;
  forall_edges(e, G)
  {
	G[e].save(o);
}
*/
}



node get_node_by_num(TGraphC & G, int n)
{
  if(n<G.number_of_nodes()/2)
  {
	node x = G.first_node();
	
	int i=1;
	
	while(i!=n)
	{
	  x = G.succ_node(x);
	  i++;
	}
	
	return x;
  }
  else
  {
	node x = G.last_node();
	
	int i=G.number_of_nodes();
	
	while(i!=n)
	{
	  x = G.pred_node(x);
	  i--;
	}
	
	return x;
  }
}

		  
		  void put_arrangement_into_G(TGraphC & G, list<int> & A)
		  {
			list_item it = A.first();
			node v;
			
			forall_nodes(v, G)
			{
			  cerr << "node id:" << G[v].initial_id << " will be " << A[it] << endl;
			  G[v].ArrId = A[it];
			  it = A.succ(it);;
			}
			
			TMP_CMP_GRAPHC = &G;
			G.sort_nodes(&cmp_ArrId);
			//define_S_values(G);
			//  forall_nodes(v, G)
			//    cerr << G[v].initial_id << endl;
		  }
		  