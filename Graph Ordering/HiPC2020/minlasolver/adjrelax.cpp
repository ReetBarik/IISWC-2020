#include "cmpfuncs.h"

double adjacent_minimization(TGraphC & G, double & lin_arr, int upto_dist)
{
  node u, v;
  int flips_done = 0;
          
  //  array<node> NodesArr(G.number_of_nodes());
  cerr << "START  ADJ MIN : " << lin_arr << endl;
  /*
  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }
  */        
  //  min_lin_arr = lin_arr;
  
  double old_lin_arr = lin_arr;
  for(int i=0; i<5; i++)
    {
  edge e;
  forall_edges(e, G)
    {
      u = G.source(e);
      v = G.target(e);

      if(fabs(G[u].ArrId - G[v].ArrId) < upto_dist)
        {
          
          double new_lin_arr = flip_two_nodesC(u, v, G, old_lin_arr);
          //          cerr << "FLIP DONE = " << old_lin_arr << " " <<  new_lin_arr << endl;
          if(new_lin_arr > old_lin_arr)            
              old_lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
          else
            {
              flips_done++;
              /*
              if(new_lin_arr <= min_lin_arr)
                {
                  min_lin_arr = new_lin_arr;
                  
                  save_arrangement(G);
                }
              */
              old_lin_arr = new_lin_arr;
              //              NodesArr[i] = v;
              //              NodesArr[i+dist] = u;
            }
        }      
    }
    }  
  cerr << "FINISH ADJ MIN : " << old_lin_arr << ", flips : " << flips_done << endl;
  lin_arr = old_lin_arr;
  return old_lin_arr;
}
