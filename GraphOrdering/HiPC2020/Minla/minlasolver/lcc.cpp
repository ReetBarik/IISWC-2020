#include "cmpfuncs.h"
#include <time.h>
#define LCC_DOUBLE_DIFF 0.00000000000000000001
node lcc_first;

list<list<node> > equal_A_vals;
list<two_tuple<node ,node> > lcc_tuples;

double calc_lcc_cost(TGraphC & G)
{
  double ret = 0;
  edge e;
  forall_edges(e, G)
    {
      node u = G.source(e);
      node v = G.target(e);
      ret +=  G[e].w*fabs(G[u].lcc_S_value - G[v].lcc_S_value);
    }
  return ret;  
}

void lcc_print(TGraphC & G)
{
  cerr << "LCC : ";
 
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_lcc_ArrId);
  node v;
  forall_nodes(v, G)
    cerr << G[v].initial_id << ", ";
  cerr << endl;

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
}

void lcc_init(TGraphC & G)
{
  node v;
  node prev = nil;
  
  //  double SM=0;
  forall_nodes(v, G)
    {
      G[v].lcc_ArrId = G[v].ArrId;
      //      G[v].lcc_SM = SM+G[v].M;
      //      SM=G[v].lcc_SM;

      if(prev!=nil)
        G[prev].lcc_right_adj = v;
      G[v].lcc_left_adj = prev;
      prev = v;
    }
  
  G[G.last_node()].lcc_right_adj = nil;

  lcc_first = G.first_node();
  
  //  lcc_print(G);
}

void recalc_lcc_SM(TGraphC & G)
{
  /*
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_bestLCCArrID);
  define_S_values(G);
  */
  node v = lcc_first;
  //  cerr << "first lcc node " << G[v].initial_id << endl;
  
  double SM=0;
  while(v!=nil)
    {
      //      cerr << G[v].initial_id << ", lcc=" << G[v].lcc_ArrId << endl;
      G[v].lcc_SM = SM+G[v].M;
      SM=G[v].lcc_SM;
      v = G[v].lcc_right_adj;
    }
  
}

void calculate_SM(TGraphC & G)
{
  node v;
  double SM = 0;
  //  cerr << "SM calc : ";
  forall_nodes(v, G)
    {
      //      cerr << G[v].initial_id << ", ";
      G[v].SM = SM+G[v].M;
      SM=G[v].SM;
    }
  //  cerr << endl;
}


void define_bound_lcc_points(TGraphC & G)
{
  node v, u;
  forall_nodes(v, G)
    {
      G[v].lcc_flag = false;
      G[v].lcc_A = G[v].SM - G[v].lcc_SM;
    }
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_lcc_A);

  /*  
  cerr << "AVALS : ";
  forall_nodes(v, G)
    cerr << G[v].lcc_A << "(" << G[v].initial_id << "), ";
  cerr << endl;
  */
  
  v = G.succ_node(G.first_node());
  list<node> A_val_list; A_val_list.clear();
  A_val_list.push_back(G.first_node());
  
  while(v!=G.last_node())
    {
      u = A_val_list[A_val_list.last_item()];
      if(fabs(G[u].lcc_A - G[v].lcc_A)<LCC_DOUBLE_DIFF)
        {
          A_val_list.push_back(v);
        }
      else
        {
          if(A_val_list.size()>1)
            equal_A_vals.push_back(A_val_list);
          A_val_list.clear(); 
          A_val_list.push_back(v);         
        }
      v = G.succ_node(v);  
    }
  
  if(A_val_list.size()>1)
    equal_A_vals.push_back(A_val_list);
}

void define_lcc_S_values(TGraphC & G)
{
  node v = G.first_node();
  
  G[v].lcc_S_value = G[v].w/2.0;
  v = G.succ_node(v);
  while(v!=nil)
    {
      node u = G.pred_node(v);
      G[v].lcc_S_value = G[u].lcc_S_value + G[u].w/2.0+G[v].w/2.0;
      v = G.succ_node(v);
    }
}

void print_lcc_error(TGraphC & G, node v1, node v2)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);

  node o = v1;
  node l = v1;

  while (o!=G.succ_node(v2))
    {
      cerr << "ArrId=" << G[o].ArrId << ", SM=" << G[o].SM << ", A=" << G[o].lcc_A << endl;
      o = G.succ_node(o);
    }
  cerr << "--------------------\n";
  while (l!=G[G.succ_node(v2)].lcc_right_adj)
    {
      cerr << "lcc_ArrId=" << G[l].lcc_ArrId << ", SM=" << G[l].SM << ", A=" << G[l].lcc_A << endl;
      l = G[l].lcc_right_adj;
    }
}

double sort_and_clear_by_lcc_ArrId_inside_A_vals(TGraphC & G)
{
  //  cerr << "sort : " << equal_A_vals.size() << endl;
  list_item it;
  
  TMP_CMP_GRAPHC = &G;
  forall_items(it, equal_A_vals)
    {      
      equal_A_vals[it].sort(&cmp_lcc_ArrId);
      list<node> tmp_node_list; tmp_node_list.clear();
      list_item it2;
      //      cerr << "sublist len = " << equal_A_vals[it].size() << endl;
      int open_flag = 1; // 1 - closed; 2 - opened
      forall_items(it2, equal_A_vals[it])
        {
          //          if(S==0)
          //          cerr << "id=" << G[equal_A_vals[it][it2]].initial_id << ", lcc_ArrId=" << G[equal_A_vals[it][it2]].ArrId << ", ArrId=" << G[equal_A_vals[it][it2]].lcc_ArrId << ", A=" << G[equal_A_vals[it][it2]].lcc_A << ", lccSM=" <<  G[equal_A_vals[it][it2]].lcc_SM << ", SM=" <<  G[equal_A_vals[it][it2]].SM <<endl;
          list_item prev_it2 = equal_A_vals[it].pred(it2);
          list_item succ_it2 = equal_A_vals[it].succ(it2);
          //          if(S==0)
          //          cerr << prev_it2 << ", " << succ_it2 << "(" << G[equal_A_vals[it][succ_it2]].ArrId << ")" << endl;
          
          if(
             //             ((it2==equal_A_vals[it].first_item())&&(G[succ_it2].ArrId!=G[it2].ArrId+1))||
             ((it2==equal_A_vals[it].last_item())&&(G[equal_A_vals[it][prev_it2]].ArrId!=G[equal_A_vals[it][it2]].ArrId-1))||
             ((succ_it2!=nil)&&(G[equal_A_vals[it][succ_it2]].ArrId!=G[equal_A_vals[it][it2]].ArrId+1)&&(open_flag==1))||
             ((succ_it2!=nil)&&(G[equal_A_vals[it][succ_it2]].ArrId==G[equal_A_vals[it][it2]].ArrId+1)&&(open_flag==2))
             )
            {                            
              tmp_node_list.push_back(equal_A_vals[it][it2]);
              //              cerr << "done " << open_flag << endl;
              open_flag = 3 - open_flag;                          
            }
        }
      //      equal_A_vals[it] = tmp_node_list;
      //      cerr << "tmp list constructed" << endl;    
      forall_items(it2, tmp_node_list)
        {        
          if(it2!=tmp_node_list.last_item())
            {              
              node v1 = tmp_node_list[it2]; node v2 = tmp_node_list[tmp_node_list.succ(it2)] ;
              // peresmotret' !!!
              if((fabs(G[v1].ArrId - G[v2].ArrId)>2)&&(G[v1].ArrId - G[v2].ArrId==G[v1].lcc_ArrId - G[v2].lcc_ArrId))
                {
                  two_tuple<node ,node> scsp(v1, v2);
                  lcc_tuples.push_back(scsp);
                  /*
                  if(G[v1].ArrId - G[v2].ArrId!=G[v1].lcc_ArrId - G[v2].lcc_ArrId)
                    {
                      cerr << "figovo" << endl;
                      cerr << "arr : " << G[v1].ArrId << "(" << G[v1].initial_id << "), " << G[v2].ArrId << "(" << G[v2].initial_id << ")" << endl;
                      cerr << "lcc : " << G[v1].lcc_ArrId << ", " << G[v2].lcc_ArrId << endl;
                      print_lcc_error(G, v1, v2);
                      //                      lcc_print(G);
                      //                      graph_print(G, "order : ");
                      
                      exit(1);
                    }
                  */
                }
            }
        }
      
    }

  lcc_tuples.sort(&cmp_lcc_tuples_by_length);

  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_lcc_ArrId);
  define_lcc_S_values(G);

  double lcc_lin_arr =  calc_lcc_cost(G);

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

  return lcc_lin_arr;
}

void check_lcc_cost(TGraphC & G)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_lcc_ArrId);
  define_lcc_S_values(G);

  double lcc_lin_arr =  calc_lcc_cost(G);
  cerr << "Check LCC cost : " << lcc_lin_arr << endl;

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G); 
}

double calc_subpermutation_cost(node v1, node v2, TGraphC & G)
{
  
  edge e;
  node v = G.succ_node(v1);

  double cost = 0;
  double to_remove_from_cost = 0;
  
  while (v!=v2)
    {
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          double e_cost = G[e].w*fabs(G[w].S_value - G[v].S_value);
          cost+=e_cost;
          if((G[w].ArrId>G[v1].ArrId)&&(G[w].ArrId<G[v2].ArrId))
            to_remove_from_cost+=e_cost;
        }
      v = G.succ_node(v);
    }

  cost = cost - to_remove_from_cost / 2.0;

  return cost;
}

double calc_subpermutation_cost_in_lcc(node v1, node v2, TGraphC & G)
{
  edge e;
  node v = G[v1].lcc_right_adj;

  double cost = 0;
  double to_remove_from_cost = 0;
  
  while (v!=v2)
    {
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          double e_cost = G[e].w*fabs(G[w].lcc_S_value - G[v].lcc_S_value);
          cost+=e_cost;
          if((G[w].lcc_ArrId>G[v1].lcc_ArrId)&&(G[w].lcc_ArrId<G[v2].lcc_ArrId))
            to_remove_from_cost+=e_cost;
        }
      v = G[v].lcc_right_adj;
    }

  cost = cost - to_remove_from_cost / 2.0;

  return cost;
}

double calc_sublcc_cost_in_order(node v1, node v2, TGraphC & G)
{
  edge e;
  node v = G[v1].lcc_right_adj;
  
  // redefine S values inside lcc given an order
  while (v!=v2)
    {
      G[v].save_S_value = G[v].S_value;
      G[v].S_value = G[G[v].lcc_left_adj].S_value + G[G[v].lcc_left_adj].w/2.0 + G[v].w / 2.0;
      v = G[v].lcc_right_adj;
    }

  double cost = 0;
  double to_remove_from_cost = 0;
  
  v = G[v1].lcc_right_adj;
  while (v!=v2)
    {
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          double e_cost = G[e].w*fabs(G[w].S_value - G[v].S_value);
          cost+=e_cost;
          if((G[w].ArrId>G[v1].ArrId)&&(G[w].ArrId<G[v2].ArrId))
            to_remove_from_cost+=e_cost;
        }
      v = G[v].lcc_right_adj;
    }

  cost = cost - to_remove_from_cost / 2.0;

  // restore S values inside lcc given an order
  v = G[v1].lcc_right_adj;
  while (v!=v2)
    {
      G[v].S_value = G[v].save_S_value;
      v = G[v].lcc_right_adj;
    }
      
  return cost;
}
double calc_subperm_cost_in_lcc(node v1, node v2, TGraphC & G)
{
  edge e;
  node v = G.succ_node(v1);
  
  // redefine S values inside lcc given an order
  while (v!=v2)
    {
      G[v].save_S_value = G[v].lcc_S_value;
      G[v].lcc_S_value = G[G.pred_node(v)].lcc_S_value + G[G.pred_node(v)].w/2.0 + G[v].w / 2.0;
      v = G.succ_node(v);
    }

  double cost = 0;
  double to_remove_from_cost = 0;
  
  v = G.succ_node(v1);
  while (v!=v2)
    {
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          double e_cost = G[e].w*fabs(G[w].lcc_S_value - G[v].lcc_S_value);
          cost+=e_cost;
          if((G[w].lcc_ArrId>G[v1].lcc_ArrId)&&(G[w].lcc_ArrId<G[v2].lcc_ArrId))
            to_remove_from_cost+=e_cost;
        }
      v = G.succ_node(v);
    }

  cost = cost - to_remove_from_cost / 2.0;

  // restore S values inside lcc given an order
  v = G.succ_node(v1);
  while (v!=v2)
    {
      G[v].lcc_S_value = G[v].save_S_value;
      v = G.succ_node(v);
    }
      
  return cost;
}

void insert_subperm_into_lcc(node v1, node v2, TGraphC & G)
{
  node v = G.succ_node(v1);

  //  cerr << "start insert : " << G[v].initial_id << endl;

  G[v1].lcc_right_adj = v;
  while (v!=v2)
    {
      //      cerr << "insert : " << G[v].initial_id << " lcc_arrid=" << G[G.pred_node(v)].lcc_ArrId + 1 << endl;
      G[v].lcc_S_value = G[G.pred_node(v)].lcc_S_value + G[G.pred_node(v)].w/2.0 + G[v].w / 2.0;
      G[v].lcc_ArrId = G[G.pred_node(v)].lcc_ArrId + 1;
      G[v].lcc_left_adj = G.pred_node(v);
      G[v].lcc_right_adj = G.succ_node(v);
      G[G.pred_node(v)].lcc_right_adj = v;
      G[G.succ_node(v)].lcc_left_adj = v;
      v = G.succ_node(v);
    }
  G[v2].lcc_left_adj = G.pred_node(v2);
}

void insert_sublcc_into_order(node v1, node v2, TGraphC & G)
{
  node v = G[v1].lcc_right_adj;

  int start_ArrId = G[v1].ArrId;
  int i = 1;
  while (v!=v2)
    {
      G[v].S_value = G[G[v].lcc_left_adj].S_value + G[G[v].lcc_left_adj].w/2.0 + G[v].w / 2.0;
      G[v].ArrId = start_ArrId + i;
      G.set_node_position(v, G[v].lcc_left_adj);
      v = G[v].lcc_right_adj;
      i++;
    } 
}
void check_lcc_depth(node v1, node v2, TGraphC & G)
{
  node v = v1;

  while(v!=v2)
    {
      G[v].lcc_depth++;
      v = G.succ_node(v);
    }
  

}

bool check_similarity(node v1, node v2, TGraphC & G)
{
  node v = v1;

  while(v!=v2)
    {
      if(G.succ_node(v)!=G[v].lcc_right_adj)
        return false;
      v = G.succ_node(v);
    }
  
  return true;
}

double all_lcc;
double all_avg_len;
double good_lcc;
double good_avg_len;
double good_lcc_avg_impr;

void try_to_fit(node v1, node v2, TGraphC & G, double & lcc_cost, double & order_cost)
{
  all_lcc++;
  all_avg_len+=fabs(G[v1].ArrId - G[v2].ArrId);

  check_lcc_depth(v1, v2, G);

  if(check_similarity(v1, v2, G)==true)
    return;
  
  //   double subperm_in_order = calc_subpermutation_cost(v1, v2, G);
  //  cerr << "sp in o" << endl;
  double    sublcc_in_lcc = calc_subpermutation_cost_in_lcc(v1, v2, G);
  //  cerr << "sp in lcc" << endl;
  //  double  sublcc_in_order = calc_sublcc_cost_in_order(v1, v2, G);
  //  cerr << "sl in o" << endl;
  double   subperm_in_lcc = calc_subperm_cost_in_lcc(v1, v2, G);
  //  cerr << "sp in lcc" << endl;
 
  if(subperm_in_lcc < sublcc_in_lcc)
    {
      good_lcc++;
      good_avg_len+=fabs(G[v1].ArrId - G[v2].ArrId);
      good_lcc_avg_impr+=fabs(subperm_in_lcc - sublcc_in_lcc);

      //      cerr << "subperm_in_lcc (" <<   subperm_in_lcc << ") < sublcc_in_lcc (" << sublcc_in_lcc << ")" << endl;
      insert_subperm_into_lcc(v1, v2, G);
      lcc_cost = lcc_cost - sublcc_in_lcc + subperm_in_lcc;
    }
  /*    
  if(sublcc_in_order < subperm_in_order)
    {
      //      cerr << "sublcc_in_order (" << sublcc_in_order <<") < subperm_in_order (" << subperm_in_order << ")" << endl;
      insert_sublcc_into_order(v1, v2, G);
      order_cost = order_cost - subperm_in_order + sublcc_in_order;
    }
  */
}

void lcc_update(TGraphC & G, double & lin_arr, double & lcc_cost)
{
   all_lcc=0;
   all_avg_len=0;
   good_lcc=0;
   good_avg_len=0;
   good_lcc_avg_impr=0;

  long t;

  //  cerr << "start lcc update" << endl;

  equal_A_vals.clear();
  lcc_tuples.clear();
  
  double start_lin_arr = lin_arr;
  node v, u;
  edge e;
  double new_lcc_lin_arr;

  calculate_SM(G);

  //  cerr << "SM calculated" << endl;

  recalc_lcc_SM(G); // no sort by lcc_ArrId, recalculate its SM

  //  cerr << "lcc_SM calculated" << endl;
 
  define_bound_lcc_points(G);

  //  cerr << "bound points defined" << endl;

 
  new_lcc_lin_arr = sort_and_clear_by_lcc_ArrId_inside_A_vals(G); // put into lcc_tuples the possible pairs
  //  cerr << "Current LCC cost = " << new_lcc_lin_arr << endl;
  double old_lcc = new_lcc_lin_arr;

  //  cerr << "A vals sorted and clear" << endl;
  
  list_item it;
  //  cerr << "# tuples " << lcc_tuples.size() << endl;
  forall_items(it, lcc_tuples)
    {
      //     lcc_print(G); 
      //      graph_print( G,"Before Order : ");
      //      cerr << "start trying to fit " << G[lcc_tuples[it].first()].initial_id << " and " << G[lcc_tuples[it].second()].initial_id << endl;
      try_to_fit(lcc_tuples[it].first(), lcc_tuples[it].second(), G, new_lcc_lin_arr, start_lin_arr);
      //      lcc_print(G); 
      //      graph_print( G,"After Order : ");
    }

  // TMP_CMP_GRAPHC = &G;
  // G.sort_nodes(&cmp_ArrId);
  lin_arr = start_lin_arr;
  lcc_cost = new_lcc_lin_arr;
  
  all_avg_len=all_avg_len/all_lcc;   
  good_avg_len=good_avg_len/good_lcc;
  good_lcc_avg_impr=good_lcc_avg_impr/good_lcc;

  //  cerr << "#LCC = " << all_lcc << "; AVG LEN = " << all_avg_len<<endl;
  //  cerr << "#GOOD LCC = " << good_lcc << "; AVG LEN = " << good_avg_len << "; AVG IMPR = " << good_lcc_avg_impr << endl;

  double all_n = 0;
  double all_d = 0;
  double max_d = 0;
  forall_nodes(v, G)
    {
      if(G[v].lcc_depth>0)
        {
          all_n++;
          all_d+=G[v].lcc_depth;
          if(G[v].lcc_depth>max_d)
            max_d = G[v].lcc_depth;
        }
      G[v].lcc_depth = 0;
    }
  //  cerr << "AVG LCC DEPTH = "<<  all_d/all_n << "; MAX LCC DEPTH = " << max_d << endl;
  

  //  cerr << "LCC improvment " << old_lcc-lcc_cost << endl;
  //  check_lcc_cost(G);
}
