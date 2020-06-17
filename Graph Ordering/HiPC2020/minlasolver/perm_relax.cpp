#include "cmpfuncs.h"
#include <algorithm>

int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}

double minimize_dist_iter(TGraphC& G, int dist, double lin_arr, node & ret_first_node)
{
  int i,j;

  node V[dist];
  int A[dist];
  int Amin[dist];

  node q = ret_first_node;
  for(i=0; i<dist; i++)
    {
      A[i]=i;
      Amin[i]=i;
      V[i] = q;
      q = G.succ_node(q);
    }
  //  cerr << "nachali" << endl;
  node p = G.pred_node(ret_first_node);
  int N = factorial(dist);

  list<edge> el; edge e;
  for(i=0; i<dist; i++)
    forall_adj_edges(e, V[i])
      el.push_back(e);
  el.sort();
  el.unique();
  //  cerr << "nachali2" << endl;
  double new_lin_arr = lin_arr;
  double min_lin_arr = lin_arr;
  for(i=0; i<N; i++)
    {
      
      //      cerr << "---" << i << endl;
      list_item it;
      forall_items(it, el)
        {
          e = el[it]; node a = G.source(e); node b = G.target(e);
          new_lin_arr-=G[e].w*fabs(G[a].S_value - G[b].S_value);
        }
      //      cerr << "- edges" << endl;  
      node pred = p; int Id = (pred==nil)?0:G[pred].ArrId;
      for(j=0; j<dist; j++)
        {
          G[V[A[j]]].ArrId = Id+1; Id++;
          G.set_node_position(V[A[j]], pred);
          pred = V[A[j]];
        }
      //      cerr << "do nodes " << endl;
      //      cerr << G[V[A[0]]].initial_id << endl;
      //      cerr << G[V[A[j-1]]].initial_id << endl;        
      define_S_values_after_swap(V[A[0]], V[A[j-1]], G);
      //      cerr << "posle" << endl;
      forall_items(it, el)
        {
          e = el[it]; node a = G.source(e); node b = G.target(e);
          new_lin_arr+=G[e].w*fabs(G[a].S_value - G[b].S_value);
        }
      //      cerr << "+ edges" << endl;  
      if(new_lin_arr <= min_lin_arr)
        {
          min_lin_arr = new_lin_arr;
          for(j=0; j<dist; j++)
            Amin[j] = A[j];
        }
      
      next_permutation(A, A+dist);
      /*      
      node x;
      forall_nodes(x, G)
        cerr << G[x].initial_id << ", " ;
      cerr << endl;
      */
    }
  //  cerr << "minimiziruem" << endl;
  list_item it;
  forall_items(it, el)
    {
      e = el[it]; node a = G.source(e); node b = G.target(e);
      new_lin_arr-=G[e].w*fabs(G[a].S_value - G[b].S_value);
    }
  //  cerr << "minimiziruem2" << endl;
  node pred = p; int Id = (pred==nil)?0:G[pred].ArrId;
  //   cerr << "0" << endl;
  for(j=0; j<dist; j++)
    {
      //      cerr << "0.5-\n" ;cerr << Amin[j] << ", \n" ; cerr << V[Amin[j]] << ", \n"; cerr << G[V[Amin[j]]].ArrId<< endl;
      G[V[Amin[j]]].ArrId = Id+1; Id++;
      //      cerr << "1" << endl;
      G.set_node_position(V[Amin[j]], pred);//cerr << "2" << endl;
      pred = V[Amin[j]];//cerr << "3" << endl;
    }
  //  cerr << "minimiziruem3" << endl;
  define_S_values_after_swap(V[Amin[0]], V[Amin[j-1]], G);
  //  cerr << "minimiziruem4" << endl;
  forall_items(it, el)
    {
      e = el[it]; node a = G.source(e); node b = G.target(e);
      new_lin_arr+=G[e].w*fabs(G[a].S_value - G[b].S_value);
    }
  //  cerr << "minimiziruem5" << endl;
  ret_first_node = V[Amin[0]];
  //  cerr << "BEFORE MIN : " << lin_arr << ", AFTER MIN : " << new_lin_arr << endl;
  return new_lin_arr;
}

double minimize_dist(TGraphC & G, int dist)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);
  
  double lin_arr = calc_laC(G);
  double start_arr = lin_arr;
  
  
  double lin_arr_begin_of_cycle = lin_arr+1;
  int number_of_cycles = 0;
  while ((fabs(lin_arr_begin_of_cycle-lin_arr)>0.0000001)&&(number_of_cycles < Params.relax_sweeps))
    {
      node ret_first_node = G.first_node();
      number_of_cycles++;
      lin_arr_begin_of_cycle = lin_arr;
      int u=0;
      while(G[ret_first_node].ArrId+dist<=G.number_of_nodes())
        {
          u++;
          //      cerr << "Start from " << G[ret_first_node].ArrId << endl;
          lin_arr = minimize_dist_iter(G, dist, lin_arr, ret_first_node);
          ret_first_node = G.succ_node(ret_first_node);
        }
      //      cerr << u << endl;
    }

  cerr << "DISTANCE " << dist << " MINIMIZATION : " << start_arr << " -> " << lin_arr << endl;
  return lin_arr;
}
