#include "cmpfuncs.h"
#include <time.h>

//#define COLD_STEPS 30

class medge {
public:
  int group;
  double w;
  node s; // source
  node t; // target
  
  friend ostream& operator<<(ostream& o, const medge& s)
    { 
      return o;
    }
  
  friend istream& operator>>(istream& i, medge& s)
    {
      return i;
    }
  
};

typedef UGRAPH<Cnode, medge> TGraphWMM;

node second_adj_for_edge(edge e, node v, TGraphWMM & G)
{
  node u = (G.source(e)==v)?G.target(e):G.source(e);
  return u;
}

double one_pass_matching_minimization(TGraphC & G, double lin_arr)
{
  //  cerr << "BEFORE ONE PASS : " << lin_arr << endl;
  node u, v, w;

  TGraphWMM MG;
  forall_nodes(v, G)
    {
      u = MG.new_node(G[v]);
      MG[u].G_prev_ptr = v;
    }
  double new_lin_arr;
  u = MG.first_node();
  v = MG.succ_node(u);
  //  cerr << "1" << endl;
  while(v!=nil)
    {
      new_lin_arr = flip_two_nodesC(MG[u].G_prev_ptr, MG[v].G_prev_ptr, G, lin_arr);

      if(new_lin_arr < lin_arr)
        {
          medge e; e.w =  lin_arr - new_lin_arr;
          e.s = u; e.t = v;
          MG.new_edge(u, v, e);
        }
      
      new_lin_arr = flip_two_nodesC(MG[u].G_prev_ptr, MG[v].G_prev_ptr, G, new_lin_arr);
      v = G.succ_node(v);
      u = G.succ_node(u);
    }

  // distance 2
  u = MG.first_node();
  v = MG.succ_node(u);
  w = MG.succ_node(v);
  //  cerr << "1" << endl;
  while(w!=nil)
    {
      new_lin_arr = flip_two_nodesC(MG[u].G_prev_ptr, MG[w].G_prev_ptr, G, lin_arr);

      if(new_lin_arr < lin_arr)
        {
          medge e; e.w =  lin_arr - new_lin_arr;
          e.s = u; e.t = w;
          MG.new_edge(u, w, e);
        }
      
      new_lin_arr = flip_two_nodesC(MG[u].G_prev_ptr, MG[w].G_prev_ptr, G, new_lin_arr);
      w = G.succ_node(w);
      u = G.succ_node(u);
    }
  /*
  // distance 3
  u = MG.first_node();
  v = MG.succ_node(u); v = MG.succ_node(v);
  w = MG.succ_node(v);
  //  cerr << "1" << endl;
  while(w!=nil)
    {
      new_lin_arr = flip_two_nodesC(MG[u].G_prev_ptr, MG[w].G_prev_ptr, G, lin_arr);

      if(new_lin_arr < lin_arr)
        {
          medge e; e.w =  lin_arr - new_lin_arr;
          e.s = u; e.t = w;
          MG.new_edge(u, w, e);
        }
      
      new_lin_arr = flip_two_nodesC(MG[u].G_prev_ptr, MG[w].G_prev_ptr, G, new_lin_arr);
      w = G.succ_node(w);
      u = G.succ_node(u);
    }
  */
  //////////////////////////////////
  list<medge> M1; double M1w = 0;
  list<medge> M2; double M2w = 0;
  //  cerr << "2" << endl;
  int g = 1;

  edge e;
  while(MG.number_of_edges() > 0)
    {
      /*
      node x = MG.choose_node();
      while(MG.degree(x)==0)
        x = MG.choose_node();
      */
      node x = MG.first_node();
      while(MG.degree(x)==0)
        x = MG.succ_node(x);
      


      
      //      cerr << "a" << endl;
      while(MG.degree(x)>0)
        {
          double max_w = -1;
          edge max_e;
          forall_adj_edges(e, x)
            {
              if(MG[e].w > max_w)
                {
                  max_e = e;
                  max_w = MG[e].w;
                }
            }
          //          cerr << "x" << endl;  
          node y = second_adj_for_edge(max_e, x, MG);
          //          cerr << "x2" << endl;  
          if(g==1)
            {
              medge me; me = MG[max_e];
              M1.push_back(me); M1w+=MG[max_e].w; }
          else
            {
              medge me; me = MG[max_e];
              M2.push_back(me); M2w+=MG[max_e].w; }
          //          cerr << "y" << endl;
          g = 3 - g;
          MG.del_edge(max_e);
          x = y;
        }
    }
  //  cerr << "3" << endl;
  list_item it;

  if(M2w > M1w)
    {
      M1 = M2;
      M1w = M2w;
    }

  new_lin_arr = lin_arr;
  medge me;
  forall_items(it, M1)
    {
      me = M1[it];
      new_lin_arr = flip_two_nodesC(MG[me.s].G_prev_ptr, MG[me.t].G_prev_ptr, G, new_lin_arr);
    }

  //  cerr << "AFTER ONE PASS : " << new_lin_arr << endl;
  return new_lin_arr;
}

double part_of_changes_only_neghbors(TGraphC & G, double & avg_change)
{
  //  clock_t start = clock();
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);
  
  node v, u;
  double lin_arr = calc_laC(G);
  
  double start_lin_arr = lin_arr;
  
  u = G.first_node();
  v = G.succ_node(u);
  int counter = 0;
  double change_sum = 0;
  while(v!=nil)
    {
      double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
      
      if(new_lin_arr > lin_arr)
        {
          counter++;
          change_sum+=fabs(lin_arr - new_lin_arr);
        }
      //      change_sum+=fabs(lin_arr - new_lin_arr);
      
      lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
      //      cerr << start_lin_arr << "\t" << lin_arr << endl;
      //      if((lin_arr!=start_lin_arr)&&(fabs(lin_arr-start_lin_arr)>0.000001))
      //        {
      //          cerr << "SA relax : ERROR in flipping neighbors." << fabs(lin_arr-start_lin_arr) << endl;
      //          exit(1);
      //        }
      v = G.succ_node(v);
      u = G.succ_node(u);
    }

  if(counter == 0)
    avg_change = 0;
  else
    avg_change = change_sum / (double)counter;

  //  cerr << "\tPart of changes : " << (double)counter/(double)(G.number_of_nodes()-1) << endl;
  //  cerr << "Coarse nodes time = " << ((double) (clock() - start)) / CLOCKS_PER_SEC << endl;
  return 100.0*(double)counter/(double)(G.number_of_nodes()-1);
}

double part_of_changes_dist(TGraphC & G, double & avg_change, double & avg_coord_diff, int dist)
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
/*
double part_of_changes_fines_dist(TGraphC & G, double & avg_change, int dist)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_nC);
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
  while(v!=nil)
    {
      if((G[u].color!=green)||(G[v].color!=green))
        {
          double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
      
          if(new_lin_arr > lin_arr)
            {
              counter++;
              change_sum+=fabs(lin_arr - new_lin_arr);
            }
          //      change_sum+=fabs(lin_arr - new_lin_arr);
      
          lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
        }
      v = G.succ_node(v);
      u = G.succ_node(u);
    }

  if(counter == 0)
    avg_change = 0;
  else
    avg_change = change_sum / (double)counter;

  //  cerr << "\tPart of changes : " << (double)counter/(double)(G.number_of_nodes()-1) << endl;
  //  cerr << "Coarse nodes time = " << ((double) (clock() - start)) / CLOCKS_PER_SEC << endl;
  return 100.0*(double)counter/(double)(G.number_of_nodes()-1);
}
*/
// at the beginning of this function the nodes should be sorted by ArrId
int sa_one_cycle(TGraphC & G, double T, double & lin_arr, double & min_lin_arr)
{
  node u, v;
  int flips_done = 0;
          
  array<node> NodesArr(G.number_of_nodes());
  
  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }
            
  for(i=0; i<G.number_of_nodes()-1; i++)
    {
      //      int pos = rand()%(G.number_of_nodes()-1);
      
      u = NodesArr[i];
      v = NodesArr[i+1];
      
      
      double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
      //      cerr << "T=" << T << endl;
      if(move_with_min_prob(1, pow(2.718, -(new_lin_arr-lin_arr)/(T)))==false)  
        lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
      else
        {
          flips_done++;
          if(new_lin_arr <= min_lin_arr)
            {
              //              cerr << "SA1 saving : " << new_lin_arr << endl;
              min_lin_arr = new_lin_arr;
              save_arrangement(G/*, MinArr*/);
            }
          lin_arr = new_lin_arr;
          NodesArr[i] = v;
          NodesArr[i+1] = u;
        }      
    }
  return flips_done;
}
int sa_one_cycle_dist(TGraphC & G, double T, double & lin_arr, double & min_lin_arr, int dist)
{
  //  cerr << "Nachali one cycle : " << lin_arr << "\t" << min_lin_arr << endl;
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
      //      int pos = rand()%(G.number_of_nodes()-1);
      
      u = NodesArr[i];
      v = NodesArr[i+dist];
      
      ////////diff in coordinates/////
      double coord_diff = 0;
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
      ////////////////////////////////
      
      double old_lin_arr = lin_arr;
      double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
      //      cerr <<  lin_arr << "\tnovoe : " <<new_lin_arr<< endl;
      //      cerr << pow(2.718, -(new_lin_arr-lin_arr)/T) << endl;
      //      cerr << Params.sa_coor << endl;
      if(move_with_min_prob(1, pow(2.718, -(new_lin_arr-lin_arr)/(Params.sa_coor*coord_diff*T)))==false)
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
          flips_done++;
          if(new_lin_arr <= min_lin_arr)
            {
              //              cerr << "SA2-one cycle saving : " << new_lin_arr << endl;
              min_lin_arr = new_lin_arr;
              
              save_arrangement(G/*, MinArr*/);
              /*
              node p;
              forall_nodes(p, G)
                cerr << G[p].initial_id << ", ";
              cerr << endl;
              TMP_CMP_GRAPHC = &G;
              G.sort_nodes(&cmp_nC);
              define_S_values(G);
              double a_arr = calc_laC(G);
              cerr << "proverochka : " << a_arr << endl;
              forall_nodes(p, G)
                cerr << G[p].initial_id << ", ";
              cerr << endl;
              */
            }
          lin_arr = new_lin_arr;
          NodesArr[i] = v;
          NodesArr[i+dist] = u;
        }      
    }
  //  cerr << flips_done << endl;
  //  cerr << lin_arr << endl;
  //  exit(1);
  return flips_done;
}
/*
int sa_fines_one_cycle_dist(TGraphC & G, double T, double & lin_arr, double & min_lin_arr, int dist)
{
  //  cerr << "Nachali one cycle : " << lin_arr << "\t" << min_lin_arr << endl;
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
      //      int pos = rand()%(G.number_of_nodes()-1);
      
      u = NodesArr[i];
      v = NodesArr[i+dist];

      if((G[u].color!=green)||(G[v].color!=green))
        {
          double old_lin_arr = lin_arr;
          double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
          //      cerr <<  lin_arr << "\tnovoe : " <<new_lin_arr<< endl;
          //      cerr << pow(2.718, -(new_lin_arr-lin_arr)/T) << endl;
          if(move_with_min_prob(1, pow(2.718, -(new_lin_arr-lin_arr)/T))==false)
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
              flips_done++;
              if(new_lin_arr <= min_lin_arr)
                {
                  //              cerr << "SA2-one cycle saving : " << new_lin_arr << endl;
                  min_lin_arr = new_lin_arr;
                  
                  save_arrangement(G);
                }
              lin_arr = new_lin_arr;
              NodesArr[i] = v;
              NodesArr[i+dist] = u;
            }
        }      
    }
  //  cerr << flips_done << endl;
  //  cerr << lin_arr << endl;
  //  exit(1);
  return flips_done;
}
*/
double calc_T_start(double avg_change, double avg_coord_diff, double real_part_of_bad_changes)
{
  return (-1.0)* (avg_change/avg_coord_diff) / log(real_part_of_bad_changes/100.0);
}
double calc_T_start_fines(double avg_change, double real_part_of_bad_changes)
{
  return (-1.0)* (avg_change) / log(real_part_of_bad_changes/100.0);
}

double sa_on_distance(TGraphC & G, int dist)
{
  /*
  node u,v;
  double lin_arr = calc_laC(G);
  cerr << "SA_DIST: Before = " << lin_arr << endl;

  double min_lin_arr = lin_arr;
  
  save_arrangement(G);  

  if((Params.use_lcc==true)&&(G.number_of_nodes()>100))
    lcc_init(G);
 
  bool last_minla_changed = false;
  double permit_perc; double end_perc;
  
  if(dist==1)
    {
      if(Params.const_number_of_sa1pr==true)
        {
          permit_perc = Params.start_sa1perc;//+65.0*(S+0.1)/6.0;
        }
      else
        {
          //          permit_perc = ((double)Params.start_sa1perc*(double)G.number_of_nodes()/(double)sourceG_NofNodes)+1;
          permit_perc = ((double)Params.start_sa1perc-S*10.0)+1;
          if(permit_perc<20.0)
            permit_perc = 20.0;
        }
      end_perc = Params.finish_sa1perc;//+10.0/(S+0.1)/6.0;
    }
  else
    {
      permit_perc = Params.start_sa2perc;//+65.0*(S+0.1)/6.0;
      end_perc = Params.finish_sa2perc;//+10.0/(S+0.1)/6.0;
    }
  cerr << "PERMIT. PERC. : " << permit_perc << endl;

  int HC_STEPS;
  if(Params.const_number_of_hc ==true)
    HC_STEPS = Params.hot_cold_steps;
  else
    {
      HC_STEPS =  (int)((double)Params.hot_cold_steps*(double)G.number_of_nodes()/(double)sourceG_NofNodes)+1;
      if(HC_STEPS < 5)
        HC_STEPS = 5;
    }
  
   for(int T_stage = 0; T_stage<HC_STEPS; T_stage++)
    {
      double avg_change = 0;
      double avg_coord_diff = 0;
      double p = part_of_changes_dist(G, avg_change, avg_coord_diff, dist);
      double T_start, T_end;

      T_start = calc_T_start(avg_change, avg_coord_diff, permit_perc);
      T_end = calc_T_start(avg_change, avg_coord_diff, end_perc);
      
      for(double T=T_start; T>T_end; T=T*Params.sa_alpha)
          int flips_done = sa_one_cycle_dist(G, T, lin_arr, min_lin_arr, dist);

      minimization(G);

      cerr << "ENERGY AFTER HC" << T_stage << "\t";
      
      if(lin_arr <= min_lin_arr)
        {
          last_minla_changed = true;
          min_lin_arr = lin_arr;
          save_arrangement(G);
        }
      else
        last_minla_changed = false;
                  
      if((Params.use_lcc==true)&&(G.number_of_nodes()>100))
        lcc_update(G, lin_arr);
    }

   restore_arrangement(G);
   lin_arr = calc_laC(G);
   
   if((Params.use_lcc==true)&&(G.number_of_nodes()>100))
     {
       restore_best_lcc_arrangement(G);
       double lcc_lin_arr = calc_laC(G);
       
       if(lcc_lin_arr-lin_arr>0.0000001)
         {
           restore_arrangement(G);
           lin_arr = calc_laC(G);
         }
       else
         lin_arr=lcc_lin_arr;       
    }
  
   cerr << "SA_DIST: After = " << lin_arr << endl;
  
  
  return lin_arr;
  */
}

double sa_on_distance(TGraphC & G, int dist, int HC)
{
  /*
  node u,v;
  double lin_arr = calc_laC(G);
  cerr << "SA_DIST: Before = " << lin_arr << endl;
  double min_lin_arr = lin_arr;
  
  //  cerr << "SA_DIST: Before = " << lin_arr << endl;
  save_arrangement(G);
  
  if((Params.use_lcc==true)&&(S<=1))
    lcc_init(G);
 
  bool last_minla_changed = false;
  double permit_perc; double end_perc;
  if(dist==1)
    {
      permit_perc = Params.start_sa1perc;//+65.0*(S+0.1)/6.0;
      end_perc = Params.finish_sa1perc;//+10.0/(S+0.1)/6.0;
    }
  else
    {
      permit_perc = Params.start_sa2perc;//+65.0*(S+0.1)/6.0;
      end_perc = Params.finish_sa2perc;//+10.0/(S+0.1)/6.0;
    }
  cerr << "PERMIT. PERC. : " << permit_perc << endl;
  
  //   cerr << (double)G.number_of_nodes() << ", " << (double)sourceG_NofNodes<< ",  " <<COLD_STEPS1 << endl;
  for(int T_stage = 0; T_stage<HC; T_stage++)
    {
      //      cerr << "Level : " << S << ", cs=" << T_stage << ", minla=" << min_lin_arr << ", linarr=" << lin_arr << endl;
      double avg_change = 0;
      double avg_coord_diff = 0;
      double p = part_of_changes_dist(G, avg_change, avg_coord_diff, dist);
      double T_start, T_end;

      //      if(last_minla_changed==false)
      //        permit_perc = permit_perc * 1.2;
      //      else
      //        permit_perc = permit_perc / 1.2;

      //      cerr << permit_perc << endl;
      T_start = calc_T_start(avg_change, avg_coord_diff, permit_perc);
      T_end = calc_T_start(avg_change, avg_coord_diff, end_perc); 
      
      for(double T=T_start; T>T_end; T=T*Params.sa_alpha)
          int flips_done = sa_one_cycle_dist(G, T, lin_arr, min_lin_arr, dist);
      
      lin_arr = minimization(G);

      cerr << "ENERGY AFTER HC" << T_stage << " = " << lin_arr << endl;
      //      lin_arr = triples_minimization(G);  
      
      if(lin_arr <= min_lin_arr)
        {
          last_minla_changed = true;
          min_lin_arr = lin_arr;
          save_arrangement(G);
        }
      else
        last_minla_changed = false;
                  
      if((Params.use_lcc==true)&&(S<=1))
        {
          //          if(T_stage==0)
          //            lcc_init(G);
          //          else
            lcc_update(G, lin_arr);
        }
      
      //      cerr << "LevelA : " << S << ", cs=" << -1 << ", minla=" << min_lin_arr << ", linarr=" << lin_arr << endl;
    }
  //  cerr << "Level : " << S << ", cs=" << -1 << ", minla=" << min_lin_arr << ", linarr=" << lin_arr << endl;
  restore_arrangement(G);
  lin_arr = calc_laC(G);
  //  cerr << "Best after SA : " << lin_arr << endl;
  
  if((Params.use_lcc==true)&&(S<=1))
    {
      restore_best_lcc_arrangement(G);
      double lcc_lin_arr = calc_laC(G);
      //  cerr << "Best LCC after SA : " << lcc_lin_arr << endl;
      
      if(lcc_lin_arr-lin_arr>0.0000001)
        {
          restore_arrangement(G);
          lin_arr = calc_laC(G);
        }
      else
        {
          //      cerr << "---------------------- GOT LCC BEST ----------------------------" << endl;
          lin_arr=lcc_lin_arr;
        }
    }
  
  cerr << "SA_DIST: After = " << lin_arr << endl;
  
  
  return lin_arr;
  */
}
/*
double sa_on_fines_distance(TGraphC & G, int dist)
{
  node u,v;
  double lin_arr = calc_laC(G);
  cerr << "SA_ON_FINES_DIST: Before = " << lin_arr << endl;
  double min_lin_arr = lin_arr;
  
  //  cerr << "SA_DIST: Before = " << lin_arr << endl;
  save_arrangement(G);
  
  if((Params.use_lcc==true)&&(S<=1))
    lcc_init(G);
 
  bool last_minla_changed = false;
  double permit_perc; double end_perc;
  if(dist==1)
    { permit_perc = 20; end_perc = 1.0; }
  else
    { permit_perc = 20; end_perc = 1.0; }
  
  int HC_STEPS;
  if(Params.const_number_of_hc ==true)
    HC_STEPS = Params.hot_cold_steps;
  else
    {
      HC_STEPS =  (int)((double)Params.hot_cold_steps*(double)G.number_of_nodes()/(double)sourceG_NofNodes)+1;
      if(HC_STEPS < 5)
        HC_STEPS = 5;
    }
  HC_STEPS = 5;
  
  cerr << "|V|="<< G.number_of_nodes() << ", |V|source=" << (double)sourceG_NofNodes<< ",  HC = " << HC_STEPS << endl;
  for(int T_stage = 0; T_stage<HC_STEPS; T_stage++)
    {
      //      cerr << "Level : " << S << ", cs=" << T_stage << ", minla=" << min_lin_arr << ", linarr=" << lin_arr << endl;
      double avg_change = 0;
      
      double p = part_of_changes_fines_dist(G, avg_change, dist);
      double T_start, T_end;

      //      if(last_minla_changed==false)
      //        permit_perc = permit_perc * 1.2;
      //      else
      //        permit_perc = permit_perc / 1.2;

      //      cerr << permit_perc << endl;
      T_start = calc_T_start_fines(avg_change,permit_perc);
      T_end = calc_T_start_fines(avg_change, end_perc); 
      
      for(double T=T_start; T>T_end; T=T*Params.sa_alpha)
          int flips_done = sa_fines_one_cycle_dist(G, T, lin_arr, min_lin_arr, dist);
      
      minimization(G);

      cerr << "ENERGY AFTER HC" << T_stage << " = " << lin_arr << endl;
      //      lin_arr = triples_minimization(G);  
      
      if(lin_arr <= min_lin_arr)
        {
          last_minla_changed = true;
          min_lin_arr = lin_arr;
          save_arrangement(G);
        }
      else
        last_minla_changed = false;
                  
      if((Params.use_lcc==true)&&(S<=1))
        {
          //          if(T_stage==0)
          //            lcc_init(G);
          //          else
            lcc_update(G, lin_arr);
        }
      
      //      cerr << "LevelA : " << S << ", cs=" << -1 << ", minla=" << min_lin_arr << ", linarr=" << lin_arr << endl;
    }
  //  cerr << "Level : " << S << ", cs=" << -1 << ", minla=" << min_lin_arr << ", linarr=" << lin_arr << endl;
  restore_arrangement(G);
  lin_arr = calc_laC(G);
  //  cerr << "Best after SA : " << lin_arr << endl;
  
  if((Params.use_lcc==true)&&(S<=1))
    {
      restore_best_lcc_arrangement(G);
      double lcc_lin_arr = calc_laC(G);
      //  cerr << "Best LCC after SA : " << lcc_lin_arr << endl;
      
      if(lcc_lin_arr-lin_arr>0.0000001)
        {
          restore_arrangement(G);
          lin_arr = calc_laC(G);
        }
      else
        {
          //      cerr << "---------------------- GOT LCC BEST ----------------------------" << endl;
          lin_arr=lcc_lin_arr;
        }
    }
  
  cerr << "SA_ON_FINES_DIST: After = " << lin_arr << endl;
  
  
  return lin_arr;
}
*/
