#include "cmpfuncs.h"

void decrease_lin_arr_by_nlist(TGraphC & G, list<node> & gap, double & lin_arr)
{
  edge e;

  list_item p;
  forall_items(p,gap)
    forall_adj_edges(e, gap[p])
    {
      if(G[e].InFlipFlag == false)
        {
          node x = G.source(e); node y = G.target(e);
          lin_arr -= G[e].w*fabs(G[x].S_value-G[y].S_value);
          G[e].InFlipFlag = true;
        }
    }
}

void increase_lin_arr_by_nlist(TGraphC & G, list<node> & gap, double & lin_arr)
{
  edge e;

  list_item p;
  forall_items(p,gap)
    forall_adj_edges(e, gap[p])
    {
      if(G[e].InFlipFlag == true)
        {
          node x = G.source(e); node y = G.target(e);
          lin_arr += G[e].w*fabs(G[x].S_value-G[y].S_value);
          G[e].InFlipFlag = false;
        }
    }
}

void decrease_lin_arr_by_gap(TGraphC & G, array<node> & gap, double & lin_arr)
{
  edge e;
  list<edge> ee;

  for(int i=0; i<gap.size(); i++)
    forall_adj_edges(e, gap[i])
      ee.push_back(e);

  ee.sort();
  ee.unique();

  list_item it;
  forall_items(it,ee)
    {
      e=ee[it];
      node x = G.source(e); node y = G.target(e);
      lin_arr -= G[e].w*fabs(G[x].S_value-G[y].S_value);
    }
}

void increase_lin_arr_by_gap(TGraphC & G, array<node> & gap, double & lin_arr)
{
  edge e;
  list<edge> ee;

  for(int i=0; i<gap.size(); i++)
    forall_adj_edges(e, gap[i])
      ee.push_back(e);

  ee.sort();
  ee.unique();

  list_item it;
  forall_items(it,ee)
    {
      e=ee[it];
      node x = G.source(e); node y = G.target(e);
      lin_arr += G[e].w*fabs(G[x].S_value-G[y].S_value);
    }
}



void dummy_flip_two_nodesC(node v, node u, TGraphC & G)
{  
  list_item it;

  int tmp = G[u].ArrId; G[u].ArrId = G[v].ArrId; G[v].ArrId = tmp;
         
  node xp = G.pred_node(u);
  node yp = G.pred_node(v);

  G.set_node_position(u, yp);
  G.set_node_position(v, xp);
      
  if(G[u].ArrId < G[v].ArrId)
    define_S_values_after_swap(u, v, G);
  else
    define_S_values_after_swap(v, u, G);
}

double flip_two_nodesC(node v, node u, TGraphC & G, double current_la)
{
  
  if(Params.improved_la_function==true)
    {
      node p;
      edge e;
      /*
      cerr << "perevernem : "<< G[u].initial_id << " and " << G[v].initial_id << endl;
      cerr << "nachalniy poryadok : " ;
      forall_nodes(p, G)
        cerr << G[p].initial_id << "(" << G[p].S_value << "), ";
      cerr << endl;
      */
      double lin_arr = current_la;
      
      list<node> node_to_recalc;
      node q = (G[u].ArrId<G[v].ArrId)?u:v;
      node last = (G[u].ArrId<G[v].ArrId)?v:u;
      
      while((q!=nil)&&(q!=G.succ_node(last)))
        {
          node_to_recalc.push_back(q);
          q = G.succ_node(q);
        }
      
      list_item it;
      /*
      forall_items(it, node_to_recalc)
        {
          p = node_to_recalc[it];
          forall_adj_edges(e, p)
          lin_arr = lin_arr - fabs(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
          } bilo */
      decrease_lin_arr_by_nlist(G, node_to_recalc, lin_arr);
      
      /*            
      forall_adj_edges(e, v)
        { 
          //          if((G.source(e) != u)&&(G.target(e) != u))
            lin_arr = lin_arr - fabs(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
        }
      
      forall_adj_edges(e, u)
        //        if((G.source(e) != v)&&(G.target(e) != v)) //{ cout << "+";
          lin_arr -= fabs(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
      */
      int tmp = G[u].ArrId; G[u].ArrId = G[v].ArrId; G[v].ArrId = tmp;
      //      TMP_CMP_GRAPHC = &G;  
      //            G.sort_nodes(&cmp_nC);
      //            define_S_values(G);

      
      node xp = G.pred_node(u);
      node yp = G.pred_node(v);

      G.set_node_position(u, yp);
      G.set_node_position(v, xp);
      
      if(G[u].ArrId < G[v].ArrId)
        define_S_values_after_swap(u, v, G);
      else
        define_S_values_after_swap(v, u, G);
      
      increase_lin_arr_by_nlist(G, node_to_recalc, lin_arr);
      /*
      forall_items(it, node_to_recalc)
        {
          p = node_to_recalc[it];
          forall_adj_edges(e, p)
            lin_arr = lin_arr + fabs(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
        }
        bilo
      */
      /*
      forall_adj_edges(e, v)
        //        if((G.source(e) != u)&&(G.target(e) != u))  
          lin_arr += fabs(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
      
      forall_adj_edges(e, u)
        //        if((G.source(e) != v)&&(G.target(e) != v))
        lin_arr += fabs(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
      
      */


      
      //      cerr << "konchili : " << lin_arr << endl;
      /*
      cerr << "konechniy poryadok : " ;
      forall_nodes(p, G)
        cerr << G[p].initial_id << "(" << G[p].S_value << "), ";
      cerr << endl;
      */
      return lin_arr;
    }
  else
    {
      //       cerr << "Error : The minla energy function should be improved" << endl;
      int tmp;
      edge e;
      
      double lin_arr = current_la;//calc_laC(G);
      
      forall_inout_edges(e, v)
        { 
          if((G.source(e) != u)&&(G.target(e) != u))
            lin_arr = lin_arr - fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
        }
      
      forall_inout_edges(e, u)
        if((G.source(e) != v)&&(G.target(e) != v)) //{ cout << "+";
          lin_arr -= fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
      
      tmp = G[u].ArrId; G[u].ArrId = G[v].ArrId; G[v].ArrId = tmp;
      
      forall_inout_edges(e, v)
        if((G.source(e) != u)&&(G.target(e) != u))  
          lin_arr += fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
      
      forall_inout_edges(e, u)
        if((G.source(e) != v)&&(G.target(e) != v))
          lin_arr += fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;

      return lin_arr;
    }
  
}
void flip_two_nodes_no_cost_calc(node v, node u, TGraphC & G)
{ 
  node p;
  edge e;

  list<node> node_to_recalc;
  node q = (G[u].ArrId<G[v].ArrId)?u:v;
  node last = (G[u].ArrId<G[v].ArrId)?v:u;
      
  while((q!=nil)&&(q!=G.succ_node(last)))
    {
      node_to_recalc.push_back(q);
      q = G.succ_node(q);
    }
  
  list_item it;

  //  decrease_lin_arr_by_nlist(G, node_to_recalc, lin_arr);
      
  int tmp = G[u].ArrId; G[u].ArrId = G[v].ArrId; G[v].ArrId = tmp;
        
  node xp = G.pred_node(u);
  node yp = G.pred_node(v);
  
  G.set_node_position(u, yp);
  G.set_node_position(v, xp);
  
  if(G[u].ArrId < G[v].ArrId)
    define_S_values_after_swap(u, v, G);
  else
    define_S_values_after_swap(v, u, G);
  
  //  increase_lin_arr_by_nlist(G, node_to_recalc, lin_arr);

}

double flip_two_lcc_nodesC(node v, node u, TGraphC & G, double current_la)
{
  cerr << "Bad call : flip_two_lcc_nodesC" << endl;
  exit(1);
  /*  
  if(Params.improved_la_function==true)
    {
      node p;
      edge e;
      double lin_arr = current_la;
      
      list<node> node_to_recalc;
      
      if(G[u].bestLCC_ArrId<G[v].bestLCC_ArrId)
        {
          node q=u;
          while((q!=nil)&&(q!=G.succ_node(v)))
            {
              node_to_recalc.push_back(q);
              q = G.succ_node(q);
            }
        }
      else
        {
          node q=v;
          while((q!=nil)&&(q!=G.succ_node(u)))
            {
              node_to_recalc.push_back(q);
              q = G.succ_node(q);
            }
        }
      
      list_item it;

      decrease_lin_arr_by_nlist(G, node_to_recalc, lin_arr);
      
      int tmp = G[u].bestLCC_ArrId; G[u].bestLCC_ArrId = G[v].bestLCC_ArrId; G[v].bestLCC_ArrId = tmp;
      
      node xp = G.pred_node(u);
      node yp = G.pred_node(v);

      G.set_node_position(u, yp);
      G.set_node_position(v, xp);
      
      if(G[u].bestLCC_ArrId < G[v].bestLCC_ArrId)
        define_S_values_after_swap(u, v, G);
      else
        define_S_values_after_swap(v, u, G);
      
      increase_lin_arr_by_nlist(G, node_to_recalc, lin_arr);

      return lin_arr;
    }
  else
    {
      //       cerr << "Error : The minla energy function should be improved" << endl;
      int tmp;
      edge e;
      
      double lin_arr = current_la;//calc_laC(G);
      
      forall_inout_edges(e, v)
        { 
          if((G.source(e) != u)&&(G.target(e) != u))
            lin_arr = lin_arr - fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
        }
      
      forall_inout_edges(e, u)
        if((G.source(e) != v)&&(G.target(e) != v)) //{ cout << "+";
          lin_arr -= fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
      
      tmp = G[u].ArrId; G[u].ArrId = G[v].ArrId; G[v].ArrId = tmp;
      
      forall_inout_edges(e, v)
        if((G.source(e) != u)&&(G.target(e) != u))  
          lin_arr += fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;
      
      forall_inout_edges(e, u)
        if((G.source(e) != v)&&(G.target(e) != v))
          lin_arr += fabs(G[G.source(e)].ArrId - G[G.target(e)].ArrId) * G[e].w;

      return lin_arr;
    }
  */  
}
