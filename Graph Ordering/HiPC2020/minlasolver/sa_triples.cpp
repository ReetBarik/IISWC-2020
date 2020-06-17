#include "cmpfuncs.h"

double calc_T_start2(double avg_change, double real_part_of_bad_changes)
{ return (-1.0)* (avg_change) / log(real_part_of_bad_changes/100.0);
}

double get_second_triple_order(node u, node v, node w, TGraphC & G, double lin_arr)
{
  int i = G[u].ArrId - 1;
  
  array<triple_res> all_res(6);

  triple_res res;

  res.cost = lin_arr; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId;
  all_res[0] = res;

  double l1 = flip_two_nodesC(u, v, G, lin_arr);
  res.cost = l1; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[1] = res;

  double l2 = flip_two_nodesC(u, w, G, l1);
  res.cost = l2; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId;  
  all_res[2] = res;                                                                      
  double l3 = flip_two_nodesC(v, w, G, l2);
  res.cost = l3; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[3] = res;
  
  double l4 = flip_two_nodesC(v, u, G, l3);
  res.cost = l4; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[4] = res;                                                                      
  double l5 = flip_two_nodesC(w, u, G, l4);
  res.cost = l5; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[5] = res;
  
  all_res.sort(&cmp_triple);

  double diff=(all_res[5].cost - all_res[0].cost)/5.0;
  
  int j;
  for(j=1; j<6; j++)
    {
      if(all_res[j].cost-all_res[0].cost>1)
        break;
    }
  if(j<6)
    res = all_res[j];
  else
    res = all_res[0];

  //  res = all_res[5];

  G[u].ArrId = res.u_pos;
  G[v].ArrId = res.v_pos;
  G[w].ArrId = res.w_pos;

  node min_node, mid_node, max_node;
  
  if(G[u].ArrId==i+1)
    min_node = u;
  else if(G[v].ArrId==i+1)
    min_node = v;
  else if(G[w].ArrId==i+1)
    min_node = w;
  
  if(G[u].ArrId==i+2)
    mid_node = u;
  else if(G[v].ArrId==i+2)
    mid_node = v;
  else if(G[w].ArrId==i+2)
    mid_node = w;


  if(G[u].ArrId==i+3)
    max_node = u;
  else if(G[v].ArrId==i+3)
    max_node = v;
  else if(G[w].ArrId==i+3)
    max_node = w;

  node tmp_start = (G.first_node()==u)?nil:(G.pred_node(u));
  G.set_node_position(min_node, tmp_start);
  G.set_node_position(mid_node, min_node);  
  G.set_node_position(max_node, mid_node);  

  define_S_values_after_swap(min_node, max_node, G);

  //  cerr << "---------------------------------------------- " << res.cost << endl;
  if(res.cost<0)
    exit(1);
  
  return res.cost;
}

double get_two_best_triple_orders(node u, node v, node w, TGraphC & G, double lin_arr,triple_res & first, triple_res & second)
{
  int i = G[u].ArrId - 1;
  
  array<triple_res> all_res(6);

  triple_res res;

  res.cost = lin_arr; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId;
  all_res[0] = res;
  //  cerr << "cost 0 " << all_res[0].cost << endl;

  double l1 = flip_two_nodesC(u, v, G, lin_arr);
  res.cost = l1; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[1] = res;
  //  cerr << "cost 1 " << all_res[1].cost << endl;

  double l2 = flip_two_nodesC(u, w, G, l1);
  res.cost = l2; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId;  
  all_res[2] = res;
  //  cerr << "cost 2 " << all_res[2].cost << endl;
    
  double l3 = flip_two_nodesC(v, w, G, l2);
  res.cost = l3; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[3] = res;
  //  cerr << "cost 3 " << all_res[3].cost << endl;
  
  double l4 = flip_two_nodesC(v, u, G, l3);
  res.cost = l4; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[4] = res;
  //  cerr << "cost 4 " << all_res[4].cost << endl;
  
  double l5 = flip_two_nodesC(w, u, G, l4);
  res.cost = l5; res.u_pos = G[u].ArrId; res.v_pos = G[v].ArrId; res.w_pos = G[w].ArrId; 
  all_res[5] = res;
  //  cerr << "cost 5 " << all_res[5].cost << endl;
  
  all_res.sort(&cmp_triple);

  double diff=(all_res[5].cost - all_res[0].cost)/5.0;
 
  int j;
  for(j=1; j<6; j++)
    {
      if(all_res[j].cost-all_res[0].cost>1)
        break;
    }
  if(j<6)
    second = all_res[j];
  else
    second = all_res[0];
  
  //  second = all_res[5];
  
  res = all_res[0];
  first = all_res[0];
  
  //  cerr << "first cost " << first.cost << endl;
  //  cerr << "second cost " << second.cost << endl;
  
  G[u].ArrId = res.u_pos;
  G[v].ArrId = res.v_pos;
  G[w].ArrId = res.w_pos;

  node min_node, mid_node, max_node;
  
  if(G[u].ArrId==i+1)
    min_node = u;
  else if(G[v].ArrId==i+1)
    min_node = v;
  else if(G[w].ArrId==i+1)
    min_node = w;
  
  if(G[u].ArrId==i+2)
    mid_node = u;
  else if(G[v].ArrId==i+2)
    mid_node = v;
  else if(G[w].ArrId==i+2)
    mid_node = w;


  if(G[u].ArrId==i+3)
    max_node = u;
  else if(G[v].ArrId==i+3)
    max_node = v;
  else if(G[w].ArrId==i+3)
    max_node = w;

  node tmp_start = (G[u].ArrId==1)?nil:(G.pred_node(u));
  G.set_node_position(min_node, tmp_start);
  G.set_node_position(mid_node, min_node);  
  G.set_node_position(max_node, mid_node);  

  define_S_values_after_swap(min_node, max_node, G);

  //  cerr << "---------------------------------------------- " << res.cost << endl;
  if(res.cost<0)
    exit(1);
  
  return res.cost;
}

double part_of_changes_on_triples(TGraphC & G, double & avg_change, double start_lin_arr)
{
  /*
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_nC);
  define_S_values(G);
  */
  node v, u, w;
  
  u = G.first_node(); v = G.succ_node(u); w = G.succ_node(v);
  
  int counter = 0;
  double change_sum = 0;
  //  double lin_arr;
  while(w!=nil)
    {
      //      graph_print(G, "order : ");
      int u_pos = G[u].ArrId; int v_pos = G[v].ArrId; int w_pos = G[w].ArrId;
      node tmp_start = (G.first_node()==u)?nil:(G.pred_node(u));
      
      double new_lin_arr = get_second_triple_order(u, v, w, G, start_lin_arr);
      //      graph_print(G, "order : ");
      if(new_lin_arr - start_lin_arr > DOUBLE_DIFF)
        {
          counter++;
          change_sum+=fabs(start_lin_arr - new_lin_arr);
        }
      
      G[u].ArrId = u_pos; G[v].ArrId = v_pos; G[w].ArrId = w_pos;
      
      G.set_node_position(u, tmp_start);
      G.set_node_position(v, u);  
      G.set_node_position(w, v);  

      define_S_values_after_swap(u, w, G);
      
      //      lin_arr = start_lin_arr;
      
      v = G.succ_node(v);
      u = G.succ_node(u);
      w = G.succ_node(w);
    }

  if(counter == 0)
    avg_change = 0;
  else
    avg_change = change_sum / (double)counter;

  return 100.0*(double)counter/(double)(G.number_of_nodes()-1);
}

void accept_tripleflip(node & u, node & v, node & w, TGraphC & G, triple_res & first, array<node>&NodesArr, double new_lin_arr, double & min_lin_arr, double & lin_arr, int i)
{
  if(new_lin_arr <= min_lin_arr)
    {
      min_lin_arr = new_lin_arr;
      
      save_arrangement(G);
    }
  lin_arr = new_lin_arr;
  
  if(first.u_pos==i+1) NodesArr[i] = u;
  else if(first.v_pos==i+1) NodesArr[i] = v;
  else if(first.w_pos==i+1) NodesArr[i] = w;
  
  if(first.u_pos==i+2) NodesArr[i+1] = u;
  else if(first.v_pos==i+2) NodesArr[i+1] = v;
  else if(first.w_pos==i+2) NodesArr[i+1] = w;
  
  if(first.u_pos==i+3) NodesArr[i+2] = u;
  else if(first.v_pos==i+3) NodesArr[i+2] = v;
  else if(first.w_pos==i+3) NodesArr[i+2] = w;
}

int sa_one_cycle_triples(TGraphC & G, double T, double & lin_arr, double & min_lin_arr)
{
  //  cerr << "Nachali one cycle(T=" << T << ") : " << lin_arr << "\t" << min_lin_arr << endl;
  node u, v, w;
  int flips_done = 0;
          
  array<node> NodesArr(G.number_of_nodes());
  
  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }
            
  for(i=0; i<G.number_of_nodes()-3; i++)
    {
      //      cerr << "position : " << i << ", (one cycle)" << endl;
      //      graph_print(G, "one cycle order : ");
      //      graph_print(G, ">>> ");
      u = NodesArr[i];
      v = NodesArr[i+1];
      w = NodesArr[i+2];

      node pred = (i==0)?nil:(G.pred_node(u));
      /*
      if(pred==nil)
        cerr << "pos. predaccessor nil" << endl;
      else
        cerr << "pos. predaccessor " << G[pred].initial_id << endl; 
      */
      //      double old_lin_arr = lin_arr;
      triple_res first,second;
      //      cerr << "take triple : " << G[u].initial_id << ", " << G[v].initial_id << ", " << G[w].initial_id << endl;
      double new_lin_arr = get_two_best_triple_orders(u, v, w, G, lin_arr, first, second);
      //      cerr << "got new_lin_arr = " << new_lin_arr << ", linarr=" << lin_arr << ", sec=" << second.cost << ", first=" << first.cost << endl;
      if(lin_arr - new_lin_arr > DOUBLE_DIFF)
        {
          accept_tripleflip(u, v, w, G, first, NodesArr, new_lin_arr, min_lin_arr, lin_arr, i);
          flips_done++;

          //          cerr << "accepted " << lin_arr << endl;
        }
      else
        {
          //      cerr << "new = " << new_lin_arr << ", linarr=" << lin_arr << ", sec=" << second.cost << endl;
          if(move_with_min_prob(1, pow(2.718, -(second.cost - lin_arr)/T))==false)
            {
              //              cerr << "Make second cost move, arr="<< second.cost << ", P=" <<pow(2.718, -(second.cost - lin_arr)/T) <<endl;
              lin_arr = second.cost;
              G[u].ArrId = second.u_pos;
              G[v].ArrId = second.v_pos;
              G[w].ArrId = second.w_pos;
              
              node min_node, mid_node, max_node;
              
              if(G[u].ArrId==i+1)
                min_node = u;
              else if(G[v].ArrId==i+1)
                min_node = v;
              else if(G[w].ArrId==i+1)
                min_node = w;
              
              if(G[u].ArrId==i+2)
                mid_node = u;
              else if(G[v].ArrId==i+2)
                mid_node = v;
              else if(G[w].ArrId==i+2)
                mid_node = w;
              
              
              if(G[u].ArrId==i+3)
                max_node = u;
              else if(G[v].ArrId==i+3)
                max_node = v;
              else if(G[w].ArrId==i+3)
                max_node = w;
              
              G.set_node_position(min_node, pred);
              G.set_node_position(mid_node, min_node);  
              G.set_node_position(max_node, mid_node);
              
              NodesArr[i] = min_node;
              NodesArr[i+1] = mid_node;
              NodesArr[i+2] = max_node;
              
              define_S_values_after_swap(min_node, max_node, G);
              /*
                if(fabs(old_lin_arr - lin_arr)>DOUBLE_DIFF)
                {
                cerr << "ERRRROR!!!\n";
                exit(1);
                }
              */
            }
          else
            {
              //          cerr << "Accept first cost move\n";
              accept_tripleflip(u, v, w, G, first, NodesArr, new_lin_arr, min_lin_arr, lin_arr, i);
              flips_done++;

              //              cerr << "Accepted(Make second cost move), arr="<< lin_arr << ", P=" <<pow(2.718, -(second.cost - lin_arr)/T) <<endl;
            }      
        }
    } 
  //  exit(1);
  return flips_done;
}



double sa_on_triples(TGraphC & G)
{
}/*
  node u,v;
  double lin_arr = calc_laC(G);
  cerr << "SA_TRIPLES: Before = " << lin_arr << endl;
  
  double min_lin_arr = lin_arr;
  
  save_arrangement(G);  
  
  if((Params.use_lcc==true)&&(S<=1))
    lcc_init(G);
  
  bool last_minla_changed = false;
  double permit_perc; double end_perc;
  
  permit_perc = Params.start_sa_triples_perc;
  end_perc    = Params.finish_sa_triples_perc;


  int HC_STEPS;
  if(Params.const_number_of_hc ==true)
    HC_STEPS = Params.hot_cold_steps;
  else
    {
      HC_STEPS =  (int)((double)Params.hot_cold_steps*(double)G.number_of_nodes()/(double)sourceG_NofNodes)+1;
      if(HC_STEPS < 5)
        HC_STEPS = 5;
    }

  //  cerr << "Enter to HC steps" << endl;
  
  for(int T_stage = 0; T_stage<HC_STEPS; T_stage++)
    {
      //      cerr << "T STAGE " << T_stage << endl;
      double avg_change = 0;
      double p = part_of_changes_on_triples(G, avg_change,lin_arr);
      double T_start, T_end;
      //      cerr << "part of changes poschitali" << endl;
      //      cerr << avg_change << ", " << p << ", " << lin_arr << endl;
      
      T_start = calc_T_start2(avg_change, permit_perc);
      T_end = calc_T_start2(avg_change, end_perc);

      //      cerr << "TStart=" << T_start << ", TEnd=" << T_end << ", avgchange=" << avg_change << ", %ch=" << p << endl;
      
      for(double T=T_start; T>T_end; T=T*Params.sa_alpha)
        int flips_done = sa_one_cycle_triples(G, T, lin_arr, min_lin_arr);
      
      if(Params.sa_triples_improve==true)
        lin_arr = triples_minimization(G);
      else
        lin_arr = only_neighbors_relaxation(G);
      
      cerr << "ENERGY AFTER HC" << T_stage << " = " << lin_arr << endl;
      
      if(lin_arr <= min_lin_arr)
        {
          last_minla_changed = true;
          min_lin_arr = lin_arr;
          save_arrangement(G);
        }
      else
        last_minla_changed = false;
      
      if((Params.use_lcc==true)&&(S<=1))
        lcc_update(G, lin_arr);
      
    }
  
  restore_arrangement(G);
  lin_arr = calc_laC(G);
  
  if((Params.use_lcc==true)&&(S<=1))
    {
      restore_best_lcc_arrangement(G);
      double lcc_lin_arr = calc_laC(G);
      
      if(lcc_lin_arr-lin_arr > DOUBLE_DIFF)
        { 
          restore_arrangement(G);
          lin_arr = calc_laC(G);
        }
      else
        {
          lin_arr=lcc_lin_arr;
        }
    }
  
  cerr << "SA_TRIPLES: After = " << lin_arr << endl;
  
  return lin_arr;
}
*/
