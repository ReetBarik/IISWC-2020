#include "cmpfuncs.h"
#include <time.h>
bool move_with_min_prob(double q, double p)
{
  double min = (q<p)?q:p;
  //  cerr << "minimum " << min << endl;
  int rnd1 = ((int)rand())%1000;
  int rnd2 = ((int)rand())%1000;
  double r = (rnd1<=rnd2)?((double)rnd1/(double)rnd2):((double)rnd2/(double)rnd1);

  //  cerr << q << ", " << p << ", " << min << ", " << r << endl;
  if(min<r)
    return false;
  else
    return true;
}

void get_permit_perc(TGraphC & G, int dist, double & permit_perc, double & end_perc)
{
  if(dist==1)
    {
      if(Params.const_number_of_sa1pr==true)
        permit_perc = Params.start_sa1perc;
      else
        {
          permit_perc = ((double)Params.start_sa1perc*(double)G.number_of_nodes()/(double)sourceG_NofNodes)+1;
          if(permit_perc<20.0)
            permit_perc = 20.0;
        }
      end_perc = Params.finish_sa1perc;
    }
  else
    {
      permit_perc = Params.start_sa2perc;
      end_perc = Params.finish_sa2perc;
    }
  //    cerr << "PERMIT. PERC. : " << permit_perc << endl;
}

int calc_HC_steps(TGraphC & G)
{
  int HC_STEPS;
  if(Params.reverse_sa==false)
    {
      if(Params.const_number_of_hc ==true)
        HC_STEPS = Params.hot_cold_steps;
      else
        {
          HC_STEPS =  (int)((double)Params.hot_cold_steps*(double)G.number_of_nodes()/(double)sourceG_NofNodes)+1;
          if(HC_STEPS < 5)
            HC_STEPS = 5;
        }
    }
  else
    {      
      HC_STEPS = Params.CURRENT_LEVEL+Params.hot_cold_steps;
      //      if(S==0)
      //        HC_STEPS=Params.hot_cold_steps;
      cerr << "REVERSE SA : #Heating-Cooling steps =" << HC_STEPS << endl;
    }

  return HC_STEPS;
}

double light_minimization(TGraphC & G)
{
  double lin_arr;
  
  if(Params.strict_minimization==1)
    {
      lin_arr = triples_minimization(G);
    }
  else if(Params.strict_minimization==0)
    {
      lin_arr = only_neighbors_relaxation(G);
    }
  else if(Params.strict_minimization==2)
    {
      lin_arr = wmm_minimization(G);
    }
   else if(Params.strict_minimization==3)
    {
      lin_arr = insert_minimization(G, false);
    }
  
 

  return lin_arr;
}
double minimization(TGraphC & G)
{
  double lin_arr;
  
  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);
  
  if(Params.strict_minimization==1)
    {
      lin_arr = triples_minimization(G);
    }
  else if(Params.strict_minimization==0)
    {
      lin_arr = only_neighbors_relaxation(G);
    }
  else if(Params.strict_minimization==2)
    {
      lin_arr = wmm_minimization(G);
    }
  else if(Params.strict_minimization==3)
    {
      lin_arr = insert_minimization(G, false);
    }
  
  

  return lin_arr;
}

void check_save_ArrId(TGraphC & G)
{
  exit(0);
  cerr << "DEBUG CHECKING" << endl;
  node x,y;
  forall_nodes(x, G)
    forall_nodes(y, G)
    if((G[x].initial_id!=G[y].initial_id)&&(G[x].save_ArrId==G[y].save_ArrId))
      {
        cerr <<G[x].save_ArrId<<endl;
        graph_print(G, "!!! : ");
        exit(1);
      }
}

int new_sa_one_cycle_dist(TGraphC & G, double T, double & lin_arr, double & min_lin_arr, int dist, bool on_fines_only)
{
  //  graph_print(G, "before : " );
  node u, v;
  int flips_done = 0;
          
  array<node> NodesArr(G.number_of_nodes());
  
  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }
            
  for(i=0; i<G.number_of_nodes()-dist; i++)
    {
      u = NodesArr[i];
      v = NodesArr[i+dist];

      if((G[u].status==seed)&&(G[v].status==seed)&&(on_fines_only==true))
        {
          
        }
      else
        {
          ////////diff in coordinates/////
          double coord_diff = 1;
          if(Params.use_coord_in_sa == true)
            {
              if(dist==1)
                {
                  coord_diff = G[v].w+G[u].w;
                }
              else if(dist==2)
                {
                  node w = NodesArr[i+1];
                  coord_diff = 2*G[w].w+G[v].w+G[u].w+fabs(G[u].w-G[v].w);
                }
              else
                {
                  coord_diff = 1;
                }
            }
          else
            Params.sa_coor = 1;
          ////////////////////////////////
          double old_lin_arr = lin_arr;
          double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
          //           check_save_ArrId(G);
          double delta = new_lin_arr - lin_arr;
          if(move_with_min_prob(1, pow(2.718, -delta/(coord_diff*Params.sa_coor*T)))==false)
            {
              lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
              if(fabs(old_lin_arr - lin_arr)>DOUBLE_DIFF)
                {
                  cerr << "ERRRROR!!!\n";
                  exit(1);
                }
            }
          else
            {
              //              cerr << "arr=" << new_lin_arr << "; min=" << min_lin_arr << "; accepted flip : u=" << G[u].ArrId << ", v=" << G[v].ArrId << endl;
              //              check_save_ArrId(G);
              flips_done++;
              if(new_lin_arr <= min_lin_arr)
                {
                  min_lin_arr = new_lin_arr;              
                  //                  local_save_arrangement(G, u, v);
                  save_arrangement(G);
                  //                  cerr << "saving u=" << G[u].ArrId << ", v=" << G[v].ArrId << endl;
                }
              lin_arr = new_lin_arr;
              NodesArr[i] = v;
              NodesArr[i+dist] = u;
            }
        }
    }
  //  graph_print(G, "after : " );
  //  check_save_ArrId(G);
  
  return flips_done;
}
double calc_averages(TGraphC & G, double & avg_change, double & avg_coord_diff, int dist)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);
  
  node v, u;
  double lin_arr = calc_laC(G);

  // If the graph is less than dist
  if(G.number_of_nodes()<dist+1)
    {
      node v1 = G.first_node();
      node v2 = G.last_node();
      double new_la = flip_two_nodesC(v1, v2, G, lin_arr);
      lin_arr = flip_two_nodesC(v1, v2, G, new_la);
      avg_change = fabs(new_la-lin_arr);
      return 95.0;
    }
  
  double start_lin_arr = lin_arr;
  
  u = G.first_node(); v = G.first_node();
  
  for(int i=0; i<dist; i++)
        v = G.succ_node(v);
  
  int counter = 0;
  double change_sum = 0;
  double change_coord_diff = 0;
  node w;
  while(v!=nil)
    {
      if(dist==2)
        w = G.succ_node(u);
      
      double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
      
      if(new_lin_arr > lin_arr)
        {
          counter++;
          change_sum+=fabs(lin_arr - new_lin_arr);
          ////////diff in coordinates/////
          double coord_diff = 0;
          if(dist==1)
            coord_diff = G[v].w+G[u].w;
          else if(dist==2)                    
            coord_diff = 2*G[w].w+G[v].w+G[u].w+fabs(G[u].w-G[v].w);          
          else
            coord_diff = 1;
          ////////////////////////////////
          change_coord_diff+=coord_diff;
        }
      //      change_sum+=fabs(lin_arr - new_lin_arr);
      
      lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
      v = G.succ_node(v);
      u = G.succ_node(u);
    }

  if(counter == 0)
    avg_change = 0;
  else
    {
      avg_change = change_sum / (double)counter;
      avg_coord_diff = change_coord_diff / (double)counter;
    }

  return 100.0*(double)counter/(double)(G.number_of_nodes()-1);
}

double calc_T(double avg_change, double avg_coord_diff, double real_part_of_bad_changes)
{
  if(Params.use_coord_in_sa == true)
    return (-1.0)* (avg_change/avg_coord_diff) / log(real_part_of_bad_changes/100.0);
  else
   return (-1.0)* (avg_change) / log(real_part_of_bad_changes/100.0); 
}

double new_sa_on_distance(TGraphC & G, int dist, bool on_fines_only)
{

  //  graph_print(G, "aga : ");
  double lcc_cost = 0;
  
  node u,v;
  //  cerr << "start" << endl;
  double lin_arr = calc_laC(G);
  cerr << "SA_DIST: Before = " << lin_arr << endl;

  double min_lin_arr = lin_arr;
  save_arrangement(G);  

  if(Params.use_lcc==true)
    lcc_init(G);
 
  bool last_minla_changed = false;
  double permit_perc; double end_perc;
  //  cerr << "before get_permitP" << endl;
  get_permit_perc(G, dist, permit_perc, end_perc);

  int HC_STEPS = (on_fines_only==true)?5:calc_HC_steps(G);
  //  cerr << "calc HC" << endl;
  
  int number_of_sweeps = 0;
  double avg_minimization_time = 0;
  double avg_lcc_time = 0;
  double avg_sweep_time = 0;
  double savg_minimization_time = 0;
  double savg_lcc_time = 0;
  double savg_sweep_time = 0;
  
  for(int T_stage = 0; T_stage<HC_STEPS; T_stage++)
    {
      //      cerr << "Start new HC sweep VVVVVVVVVV" << endl;
      clock_t start = clock();
      long t = time(0);
      
      double avg_change = 0;
      double avg_coord_diff = 0;
      double p = calc_averages(G, avg_change, avg_coord_diff, dist);
      double T_start, T_end;
      //      cerr << "STEP " << T_stage << endl;
      T_start = calc_T(avg_change, avg_coord_diff, permit_perc);
      T_end = calc_T(avg_change, avg_coord_diff, end_perc);

      int c=0;
      for(double T=T_start; T>T_end; T=T*Params.sa_alpha)
        {
          c++;
          int flips_done = new_sa_one_cycle_dist(G, T, lin_arr, min_lin_arr, dist, on_fines_only);
        }
      avg_sweep_time+= ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
      savg_sweep_time+=time(0)-t;
      
      //      cerr << "SA : Number of passes in one HC sweep = " << c << endl;;
      number_of_sweeps+=c;

      start = clock();
      t = time(0);
      lin_arr = light_minimization(G);
      avg_minimization_time+= ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
       savg_minimization_time+=time(0)-t;
      
      if(lin_arr <= min_lin_arr)
        {
          last_minla_changed = true;
          min_lin_arr = lin_arr;
          save_arrangement(G);
        }
      else
        last_minla_changed = false;

      start = clock();
      t = time(0);
      if(Params.use_lcc==true)
        {
          lcc_update(G, lin_arr, lcc_cost);
          cerr << "ORDER COST = " << lin_arr << "\t LCC COST = " << lcc_cost << endl;
        }
      avg_lcc_time+=((double)clock() - (double)start) / (double)CLOCKS_PER_SEC;
      savg_lcc_time+=time(0)-t;

      //      graph_print(G, "jopa ");
    }
  
  cerr << "In each heating-cooling step : " << (float)number_of_sweeps / (float) HC_STEPS << " sweeps\n";
  cerr << "In each heating-cooling step avg. lcc time: " << (float)savg_lcc_time / (float) HC_STEPS << "\n";
  cerr << "In each heating-cooling step avg. minimization time: " << (float)savg_minimization_time / (float) HC_STEPS << "\n";
   cerr << "In each heating-cooling step avg. sweep time: " << (float)savg_sweep_time / (float) HC_STEPS << "\n";
 
   

  restore_arrangement(G);

  //  graph_print(G, "rest ");

  lin_arr = calc_laC(G);
    
  if((Params.use_lcc==true)&&(lcc_cost<lin_arr))
    {
      restore_best_lcc_arrangement(G);
      lin_arr = calc_laC(G);
    }
  
  cerr << "SA_DIST: After = " << lin_arr << endl;

  
  return lin_arr;
}
