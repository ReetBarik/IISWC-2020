#include "cmpfuncs.h"
double quad(double x)
{
  //cerr << V_Cycle_Number << endl;
  double ret = x;
  //for(int i=0; i<V_Cycle_Number+1; i++)
    ret=ret*x;
  return ret;
}
double calc_energy_by_interpol_coords(TGraphC & G)
{
  double ret = 0;
  edge e;
 
  forall_edges(e, G)
    {
      if(Params.improved_la_function==true)
        {
          
          node u = G.source(e);
          node v = G.target(e);          
          ret +=  G[e].w*quad(G[u].interpol_coord - G[v].interpol_coord);
        }
      else
        {
          cerr << "Error : Old Objective" << endl;
          exit(1);
          ret += quad(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
        }
    }
  return ret;  
}

double calc_real_laC(TGraphC & G)
{
  double ret = 0;
  edge e;

  forall_edges(e, G)
    {
      //      cerr << G[e].w  << "; " << G[G.source(e)].ArrId << ", " <<  G[G.target(e)].ArrId << endl;
      ret += fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
    }
  
  return ret;  
}

double calc_laC(TGraphC & G)
{
  double ret = 0;
  edge e;
  //  Params.improved_la_function=false;
  //  cerr << "--------------------\n";
  forall_edges(e, G)
    {
      if(Params.improved_la_function==true)
        {
          
          node u = G.source(e);
          node v = G.target(e);
          ret +=  G[e].w*fabs(G[u].S_value - G[v].S_value);
          //          ret +=  G[e].w*fabs(G[u].ArrId - G[v].ArrId)+fabs(G[u].S_value - G[v].S_value);
        }
      else
        {        
          ret += fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
        }
    }
  //  cerr << "--------------------\n";
  return ret;  
}

void define_S_values(TGraphC & G)
{
  if(Params.improved_la_function==true)
    {
      node v = G.first_node();

      G[v].S_value = G[v].w/2.0;
      v = G.succ_node(v);
      while(v!=nil)
        {
          node u = G.pred_node(v);
          G[v].S_value = G[u].S_value + G[u].w/2.0+G[v].w/2.0;
          v = G.succ_node(v);
        }
    }
  //  else
  //    cerr << "define_S_values : ERROR" << endl;
}

void define_S_values_after_swap(node s, node t, TGraphC & G)
{
  //  cout << G[s].ArrId << ", " << G[t].ArrId << endl;
  if(Params.improved_la_function==true)
    {
      node v = s;
      
      if(G.first_node()==v)
        G[v].S_value = G[v].w/2.0;
      else
        {
          node w = G.pred_node(v);
          G[v].S_value = G[w].S_value + G[w].w/2.0+G[v].w/2.0;
        }
      
      v = G.succ_node(v);
      if(v==nil)
        {
          cerr << "define_S_values_after_swap() : ERROR" << endl;
          exit(1);
        }
        
      while((v!=nil)&&(v!=G.succ_node(t)))
        {
          node u = G.pred_node(v);
          G[v].S_value = G[u].S_value + G[u].w/2.0+G[v].w/2.0;
          v = G.succ_node(v);
        }
    }
  //  else
  //      cerr << "define_S_values_after_swap : ERROR" << endl;

  //  if(abs(G[s].ArrId-G[t].ArrId)!=1)
  //    cerr << G[s].ArrId << "\t" << G[t].ArrId << endl;
}

