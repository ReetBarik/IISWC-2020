// COARSENING NODES : double_dummy=weighted_deg

#include "cmpfuncs.h"
//#include "LEDA/basic_graph_alg.h"
#include <time.h>
#include <math.h>

double A0,V1,E0;
extern double mwf_minimization(TGraphC & G, int dist);
extern double calc_reset_mwf(TGraphC & G);
extern double wb_raise_power(TGraphC & G, int T, double curr_bnd, int pow_type);extern double insert_minimization_wb(TGraphC & G, bool for_fines_only);

#define COLOR(X)  ((G[X].status==seed) ? 'c' : 'f')
double calculate_wb_sval(TGraphC & G)
{
  node v;
  edge e;
  double total_wb = 0;
  forall_nodes(v, G)
    {
      double max_edge = 0;
      forall_adj_edges(e, v)
        {
	  node p = G.source(e);
	  node q = G.target(e);
	  if(((p==v)&&(G[q].S_value < G[p].S_value))||((q==v)&&(G[p].S_value < G[q].S_value)))
	    {
	      double ew = (G[G.source(e)].S_value - G[G.target(e)].S_value) *(G[G.source(e)].S_value - G[G.target(e)].S_value) * G[e].w;
	      if(ew > max_edge)
		max_edge = ew;
	    }
        }
      total_wb += max_edge;
    }
  return total_wb;
}

double calculate_wf_sval(TGraphC & G)
{
  list_item it, it2;
  node v, w;
  edge e;
  double total_wf = 0;
  list<node> ln; ln.clear();
  list<node> subg; subg.clear();
  forall_nodes(v, G)
    {
      subg.push_back(v);
      forall_adj_nodes(w, v)
	ln.push_back(w);
      list<node> ltmp = ln;
      ltmp.push_back(v);
      ltmp.sort();
      ltmp.unique();
      forall_items(it, subg)
	{
	  forall_items(it2, ltmp)
	    {
	      if(ltmp[it2]==subg[it])
		ltmp.erase(it2);
	    }
	}
      ltmp.push_back(v);
      ltmp.sort();
      ltmp.unique();
      total_wf += ltmp.size() * ltmp.size();

    }

  return sqrt(total_wf/(double)(G.number_of_nodes()));
}
double calculate_wf_sval2(TGraphC & G)
{
  list_item it, it2;
  node v, w;
  edge e;
  double total_wf = 0;


  array<list <node> > A(G.number_of_nodes());
  forall_edges(e, G)
    {
      node L = G.source(e);
      node R = G.target(e);
      if(G[L].ArrId > G[R].ArrId) {node TTT = L; L=R; R=TTT;}
      int l = G[L].ArrId;
      int r = G[R].ArrId;

      for(int i=l; i<r; i++)
	A[i].push_back(R);
    }
  int i=0;
  forall_nodes(v, G)
    {
      A[i].push_back(v);
      A[i].sort();
      A[i].unique();
      total_wf+=A[i].size() * A[i].size();
      i++;
    }
  return sqrt(total_wf/(double)(G.number_of_nodes()));
}


/*
void recalculate_edges(TGraphC & G)
{
  edge e;
  double s=0;
  forall_edges(e, G)
    {
      node v = G.source(e);
      node w = G.target(e);
      s+=fabs(G[v].S_value - G[w].S_value);
    }
  s = s/(double)G.number_of_edges();

  forall_edges(e, G)
    {  node v = G.source(e);
      node w = G.target(e);

      if(fabs(G[v].S_value - G[w].S_value)<s*0.1)
        G[e].w = G[e].w*2.0;
      else if (fabs(G[v].S_value - G[w].S_value)>s*1.9)
        G[e].w = G[e].w/2.0;

    }
}
double fix_edge_relative_weight(TGraphC & G, edge e)
{
  //  return G[e].rw;
  //  cerr << "fix_edge_relative_weight : error" << endl;
  //  exit(1);
  //  return G[e].w;
  if((V_Cycle_Number>0)&&(S==0))
    {
      node u = G.source(e);
      node v = G.target(e);

      double alpha= (double)V_Cycle_Number/(double)(Params.V_Cycle_Iter_Num-1);

      //return G[e].w;
      //      double alpha = (V_Cycle_Number);

      return (G[e].w*pow(fabs(G[u].S_value-G[v].S_value), alpha));//+G[e].prev_w;//G[e].w;
    }
  else
    return G[e].w;
}

void check_placement(TGraphC & G)
{
  coords_by_linear_order(G);

  node v; edge e;

  forall_nodes(v, G)
    {
      double lsum = 0;
      double rsum = 0;
      double l=0; double r=0;
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          if(G[w].ArrId<G[v].ArrId)
            {
              lsum+=G[e].w;
              l++;
            }
          else
            {
              rsum+=G[e].w;
              r++;
            }
        }
      if((lsum/(lsum+rsum) < 0.4)||(lsum/(lsum+rsum) > 0.6))
        {
          cerr << "Left edges of v ("<<l <<","<< r << ") : " << lsum/(lsum+rsum) << endl;
          minimize_vertex_coordinate(v, G);

        }
    }
  double lin_arr = linear_order_by_coord(G);
  cerr << "Res = " << lin_arr << endl;

}
*/
//#define PART_OF_GREATEST_WEIGHT 0.01
double cpu_time;
clock_t prev_clock;

extern double segment_minimization(TGraphC & G, double la);
extern double contraction_minimization(TGraphC & G, double la);
extern double segment_min(TGraphC & G, double la, double thr);

//double A0,V1,E0;

void update_time()
{
  clock_t now = clock();
  cpu_time+=(double)(now - prev_clock)/((double)CLOCKS_PER_SEC)/60.0;
  prev_clock = now;
}

int MAX_LEVEL;

int basis_solutions_passed;

array<int> solutions_on_level;

list<TGraphC> AllLevelsGraphs;
list<edge_array<double> > AllLevelsEdgeCosts;
double fix_edge_relative_weight(TGraphC & G, edge e)
{
  //  return G[e].w;
  if((V_Cycle_Number>0)&&(Params.CURRENT_LEVEL==0))
    {
      node u = G.source(e);
      node v = G.target(e);

      double alpha= (double)V_Cycle_Number/(double)(Params.V_Cycle_Iter_Num-1);
      //            if(V_Cycle_Number==1)
      //              alpha = 1;
      //            else if(V_Cycle_Number==2)
      //              alpha = 2;
     return G[e].w/pow(fabs(G[u].S_value-G[v].S_value), alpha);
    }
  else
    return G[e].w;
}



void find_all_edges(node s, TGraphC & G, TGraphC & H, node_array<double> & edges_accum)
{


  list<node> edges_to;
  node x,y,z;
  edge e,f,g;

  node sH = G[s].ptr_to_coarse;
  //  cerr << "Taking v=" << G[s].initial_id << endl;
  list_item e_it, g_it;
  forall_items(e_it, G[s].iw_links)
    {

      double e_iw = G[s].iw_links[e_it].second();
      x = G[s].iw_links[e_it].first();
      //      cerr << "v's iw link to " << G[x].initial_id << "with iw=" << e_iw << endl;

      forall_adj_edges(f, x)
        {
          double f_w = G[f].w;
          y = second_adj_for_edge(f, x, G);

          //          cerr <<  "Taking y=" << G[y].initial_id << endl;
          //           cerr << "iw link to y" <<  "with w=" << f_w << endl;
          if(y!=s)
            {
              forall_items(g_it, G[y].iw_links)
                {
                  double g_iw = G[y].iw_links[g_it].second();
                  z = G[y].iw_links[g_it].first();

                  if((z!=s)&&(z!=x)&&(G[z].status==seed))
                    {

                      //                                    cerr <<  "Taking" << endl;
                  //                  cerr << "iw link to z" <<  "with iw=" << g_iw << endl;
                      //                      if(((G[s].initial_id==828)&&(G[z].initial_id==337))||
                      //                         ((G[z].initial_id==828)&&(G[s].initial_id==337)))
                      //                                            cerr << G[s].initial_id << "---(" << e_iw << ")---" << G[x].initial_id << "^" << COLOR(x) << "---(" << f_w << ")---" << G[y].initial_id << "^" << COLOR(y) << "---(" << g_iw << ")---" << G[z].initial_id << "^" << COLOR(z) << endl;

                      if(G[s].initial_id<G[z].initial_id)
                        {
                          edges_to.push_back(G[z].ptr_to_coarse);
                          //                      cerr << "da" << endl;
                          edges_accum[G[z].ptr_to_coarse]+=e_iw*f_w*g_iw;
                          //                      cerr << "opsa" << endl;
                        }
                    }
                }
            }
        }
    }

  list_item d_it;
  //  cerr << "adding" << endl;
  //  edges_to.sort();
  //  edges_to.unique();

  forall_items(d_it, edges_to)
    {
      //      cout << H[sH].initial_id << "---" << H[edges_to[d_it]].initial_id << endl;
      //      if(H[sH].initial_id<H[edges_to[d_it]].initial_id)
      //        {
      if(edges_accum[edges_to[d_it]]>0)
        {
          CEdge ne; ne.w = edges_accum[edges_to[d_it]];
          H.new_edge(edges_to[d_it], sH, ne);
          //          cerr << "added " << H[edges_to[d_it]].initial_id << "---" << H[sH].initial_id << " with w=" << ne.w << endl;
          //        }
          edges_accum[edges_to[d_it]] = 0;
        }
    }
  //  cerr << "added\n";
}

void filter_edges(TGraphC & H)
  // double_dummy = adj_edges_w
{
  // Recalculation of filtering thershold
  double filtering_threshold = (double)Params.threshold_edge_weight*pow(0.9,log(E0/(double)H.number_of_edges()));
  if(filtering_threshold > Params.threshold_edge_weight)
    filtering_threshold = Params.threshold_edge_weight;
  if(filtering_threshold < 0.00001)
    filtering_threshold = 0.00001;
  if(Params.CURRENT_LEVEL==0)
    filtering_threshold = Params.threshold_edge_weight;
  cerr << "Recalculated filtering threshold = " << filtering_threshold << endl;
  //////////////////////////////////////////

  edge e,f;
  node u, v;
  forall_edges(e, H)
    {
      u = H.source(e);
      v = H.target(e);
      /*
            if(H[u].double_dummy < H[e].w)
              H[u].double_dummy=H[e].w;
            if(H[v].double_dummy < H[e].w)
              H[v].double_dummy=H[e].w;
      */
      H[u].double_dummy+=H[e].w;
      H[v].double_dummy+=H[e].w;
    }

  long N = H.number_of_nodes();
  if(Params.use_edges_abs_filter == false)
    {
      //      int c=0;
      forall_edges(e, H)
        {
          u = H.source(e);
          v = H.target(e);
          double u_sum = H[u].double_dummy;
          //          double u_avg = u_sum / (double)H.degree(u);
          if(u_sum*filtering_threshold /*Params.rel_filter(N)*/ /*threshold_edge_weight*/>H[e].w)//u_sum/100
            {
              double v_sum = H[v].double_dummy;
              //              double v_avg = v_sum / (double)H.degree(v);
              if(v_sum*filtering_threshold /*Params.rel_filter(N)*//*threshold_edge_weight*/>H[e].w)//v_sum/100
                {
                  //                  if((H[e].w<u_avg/2.0) && (H[e].w<v_avg/2.0))
                    H[e].be_deleted = true;
                  //              cerr << H[e].w << " (" << u_sum << ", " << v_sum << ", " << H.degree(v) << ", " << H.degree(u) << "), ";
                }
            }
        }
      //  cerr << endl;
      //      cerr << "Using relative value filtering : " << c << " edges\n";
    }
  else if(Params.use_edges_abs_filter == true)
    {
      forall_edges(e, H)
        H[e].w = H[e].w/2.0;

      int c=0;
      forall_edges(e, H)
        if(H[e].w<Params.threshold_edge_weight)
          {
            H[e].be_deleted = true;
            c++;
          }
      cerr << "Using absolute value filtering : " << c << " edges\n";
      forall_edges(e, H)
        if(H[e].be_deleted == true)
          H.del_edge(e);
      forall_edges(e, H)
        {
          if(H[e].be_deleted==true)
            {
              cerr << "BE DELETED" << endl;
              exit(1);
            }
        }
      return;
    }

  forall_edges(e, H)
    {
      if(H[e].be_deleted == true)
        {
          H.del_edge(e);
          //          cerr << H[e].w << ",";
        }

      //      else
      //        H[e].w = H[e].w/2.0;
    }
  //  cerr << endl;
  //  forall_edges(e, H)
  //    {
      /*
      if(H[e].be_deleted==true)
        {
          cerr << "BE DELETED" << endl;
          exit(1);
        }
      */
  //       H[e].w = H[e].w/2.0;
  //    }

}
/*
void old_by_avg_filter_edges(TGraphC & G, TGraphC & H)
{
  edge e,f;
  node u, v;
  forall_edges(e, H)
    {
      u = H.source(e);
      v = H.target(e);

      H[u].adj_edges_w+=H[e].w;
      H[v].adj_edges_w+=H[e].w;
    }

  if(Params.use_edges_abs_filter == false)
    {
      //      int c=0;
      forall_edges(e, H)
        {
          u = H.source(e);
          v = H.target(e);
          double u_sum = H[u].adj_edges_w;
          //          forall_adj_edges(f, u)
          //            u_sum+=H[f].w;
          double u_avg = u_sum / (double)G.degree(u);
          if(u_avg/100.0*Params.threshold_edge_weight>H[e].w)
            {
              double v_sum = H[v].adj_edges_w;
              //              forall_adj_edges(f, v)
              //                v_sum+=H[f].w;
              double v_avg = v_sum / (double)G.degree(v);
              if(v_avg/100.0*Params.threshold_edge_weight>H[e].w)
                {
                  //                  if((H[e].w<u_avg/2.0) && (H[e].w<v_avg/2.0))
                    H[e].be_deleted = true;
                  //              cerr << H[e].w << " (" << u_sum << ", " << v_sum << ", " << H.degree(v) << ", " << H.degree(u) << "), ";
                }
            }
        }
      //  cerr << endl;
      //      cerr << "Using relative value filtering : " << c << " edges\n";
    }
  else if(Params.use_edges_abs_filter == true)
    {
      forall_edges(e, H)
        H[e].w = H[e].w/2.0;

      int c=0;
      forall_edges(e, H)
        if(H[e].w<Params.threshold_edge_weight)
          {
            H[e].be_deleted = true;
            c++;
          }
      cerr << "Using absolute value filtering : " << c << " edges\n";
      forall_edges(e, H)
        if(H[e].be_deleted == true)
          H.del_edge(e);
      forall_edges(e, H)
        {
          if(H[e].be_deleted==true)
            {
              cerr << "BE DELETED" << endl;
              exit(1);
            }
        }
      return;
    }

  forall_edges(e, H)
    {
      if(H[e].be_deleted == true)
        {
          H.del_edge(e);
          //          cerr << H[e].w << ",";
        }
      //      else
      //        H[e].w = H[e].w/2.0;
    }
  //  cerr << endl;
  forall_edges(e, H)
    {

//      if(H[e].be_deleted==true)
//        {
//          cerr << "BE DELETED" << endl;
//          exit(1);
//        }

       H[e].w = H[e].w/2.0;
    }

}
*/
void construct_new_graph(TGraphC & G, TGraphC & H)
{
  node v;
  node u;

  //  Params.threshold_coarse_edges = Params.threshold_edge_weight;

  /*
  forall_nodes(v, H)
    {
      u = H[v].G_prev_ptr;
      find_all_edges(u, G, H);


      H[v].double_dummy = 0;
    }
*/
  edge e;
node_array<double> edges_accumulator(H);

  // Adding edges from paths of length 1
/*
  forall_edges(e, G)
    {
      u = G.source(e); v = G.target(e);
      if((G[u].status==seed)&&(G[v].status==seed))
        {
          CEdge h; h.w = 2*G[e].w;
          H.new_edge(G[u].ptr_to_coarse, G[v].ptr_to_coarse, h);
          edges_accumulator[G[u].ptr_to_coarse] = 0;
        }
    }
*/
  forall_nodes(v, H)
    {
      u = H[v].G_prev_ptr;
      find_all_edges/*_from_paths_of_lenght_2_3*/(u, G, H, edges_accumulator);

      H[v].double_dummy = 0;

    }

  //  cerr << "Edges defined" << endl;
  //  cerr << "Removing double edges" << endl;
  //  forall_edges(e, H)
  //    if(H[H.source(e)].initial_id<H[H.target(e)].initial_id)
  //      H.hide_edge(e);

  cerr << "Start filtering with " << H.number_of_edges() << endl;
  clock_t filter_start = clock();
  filter_edges(H);
  cerr << "Only Filtering time : " << (double)(clock()-filter_start)/(double)CLOCKS_PER_SEC << endl;

  //  G.restore_all_edges();

  // edge e;
  //  forall_edges(e, H)
  //    H[e].w = H[e].w/2.0;


  /*
  forall_nodes(v, H)
    {
      H[v].total_edges_weight = 0;
      forall_adj_edges(e, v)
        H[v].total_edges_weight+=H[e].w;
    }
  */

  /*
  double nseg = 50;
  double max=H[H.first_edge()].w;
  double min=H[H.first_edge()].w;
  double sum=0;
  forall_edges(e, H)
    {
      if(H[e].w > max) max = H[e].w;
      if(H[e].w < min) min = H[e].w;
      sum+=H[e].w;
    }
  double avg = sum / (double)H.number_of_edges();
  double seg = (max-min)/nseg;

  leda::string filename = i2s(Params.CURRENT_LEVEL);
  filename = "stat_e"+filename+".txt";
  FILE * inFile = fopen(filename, "w");
  fprintf(inFile, "max\t%f\n", max);
  fprintf(inFile, "min\t%f\n", min);
  fprintf(inFile, "seg\t%f\n", seg);

  int distr[(int)nseg];
  for(int j=0; j<nseg; j++) distr[j]=0;

  forall_edges(e, H)
    {
      distr[(int)(H[e].w/seg)]++;
      if((int)(H[e].w/seg)==0)
        {
          //          cerr << H[e].p1 << ", " << H[e].p2 << ", " << H[e].p3 << ", " << endl;
        }
    }

  for(int i=0; i<nseg; i++)
    fprintf(inFile, "seg%d\t%d\n", i,distr[i]);

  fclose(inFile);
  //  forall_edges(e, H)
  //    H[e].w = H[e].w/2.0;
  */
}


void calculate_vertex_weights(TGraphC & G, TGraphC & H)
{
  node v, w, u;
  edge e;

  forall_nodes(v, H)
    {
      w = H[v].G_prev_ptr;
      if (H[v].initial_id != G[w].initial_id)
        {
          cerr << "jopa prejopa" << endl;
          exit(1);
        }
      if(G[w].status != seed)
        {
          cerr << "Vertex weight : ERROR in G_prev_ptr" << endl;
          exit(1);
        }

      double edges_to_fine_sum = 0;
      forall_adj_edges(e, w)
        {
          u = second_adj_for_edge(e, w, G);
          if(G[u].status!=seed)
            {
              edges_to_fine_sum+=G[e].iw*G[u].w;
            }
        }
      H[v].w = G[w].w + edges_to_fine_sum;
      //G[w].w = H[v].w;
    }
}

void check_node_weights(TGraphC & G)
{
  //return ;
  double aa=0;
  node v;
  forall_nodes(v, G)
    aa+=G[v].w;
  cerr << "NODES WEIGHT = " << aa << endl;
  return ;
}

void basis_solution_iterations(TGraphC & G)
{
  //  cerr << "a" << endl;
  MAX_LEVEL = 0;
  node v;
  //  cerr << "b" << endl;
  basis_solutions_passed = 0;

//  list<LA> results;
  //  cerr << "c" << endl;
  double res;
  // cycle over the permitted basis solutions
  cerr << basis_solutions_passed << "\t" <<  Params.basis_solutions_number << endl;
//  while(basis_solutions_passed < Params.basis_solutions_number)
    {
      cerr << "starting basis " << basis_solutions_passed << endl;
      Params.CURRENT_LEVEL = -1;

      if(basis_solutions_passed==0)
        {
          //          cerr << "jopsa" << endl;
          res = method_s3(G, true);
        }
      else
        res = method_s3(G, false);
      
      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_ArrId);
      
      cerr << "final solution " << calc_laC(G) << endl;
      /*
      cerr << "b=[";
      forall_nodes(v, G) {
	cerr << G[v].initial_id << " ";
      }
      cerr << "];"<<endl;
       print_Laplacian(G, "lap.m");
      */
      //exit(1);

      Params.CURRENT_LEVEL = 0;
      int was_function = Params.relax_sweeps;//Params.improved_la_function;
      Params.relax_sweeps = 100; //Params.improved_la_function = false;
            res = triples_minimization(G);
      Params.relax_sweeps =  was_function;//Params.improved_la_function = was_function;

      res = calc_real_laC(G);

     // LA current_la;
     // current_la.weight = res;
     // current_la.V.init(G);
      //forall_nodes(v, G) {
	//cerr << G[v].ArrId << " ";
        //        G[v].ArrIds[basis_solutions_passed-1] = G[v].ArrId;
        //  current_la.V[v] = G[v].ArrId;
      //}
     // cerr << endl;
      
   //   cerr << "Storing solution : " << current_la.weight << endl;
 //     results.push_back(current_la);

      //      cerr << "Final result for BASIS Solution " << basis_solutions_passed << " = " << current_la.weight << endl;
    }
  // choose the best basis solution
//  results.sort(&cmp_LA);
//  list_item la_it = results.first();

  //cerr << "After sorting : " << results[la_it].weight << endl;

  // set the arrangement of G to the best
//  forall_nodes(v, G)
  //  {

    //  G[v].ArrId = results[la_it].V[v];
     // if(G[v].ArrId==0) cerr << "0 " ;
   // }
   // cerr << endl;

  // resort the nodes of G according to the best choosen solution
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);
}

void reset_all_vertices(TGraphC & G)
{
  node v;

  forall_nodes(v, G)
    {
      G[v].status = fine;
    }
}
extern leda::string instance_file;

void mwf_graph_print(TGraphC & G, leda::string filename)
{
  exit(1);
  FILE * outFile;
  outFile = fopen("tempedges.dat", "w");

  cout << filename << endl;
  FILE * allFile;
  allFile = fopen(filename, "w");

  //  TMP_CMP_GRAPHC = &G;
  //  G.sort_nodes(&cmp_ArrId);

    edge e;
    int c=0;
  forall_edges(e, G)
    {
      node v = G.source(e);
      node w = G.target(e);

      if(G[v].ArrId>G[w].ArrId)
	{
	  node t=v; v=w; w=t;
	}

      node t = G.succ_node(v);
      while(t!=w)
	{
	  fprintf(outFile, "%d %d 1.0\n", G[v].initial_id, G[t].initial_id);
	  c++;
	  t = G.succ_node(t);
	}
    }

  fclose(outFile);
  cout << "aga" << endl;
  fprintf(allFile, "p ghct %d %d\n", G.number_of_nodes(), c);
  fclose(allFile);
  //  cout << filename << endl;
  system("cat tempedges.dat  >>"+filename);
  //  system("\\rm -rf tempedges.dat");

}

double V_Cycle_Iterations(TGraphC & G)
{
  node v;
//  list<LA> results;

  for(int i=0; i<Params.V_Cycle_Iter_Num; i++)
    {
      cerr << "Starting " << i << "-th V-Cycle Iteration" << endl;
      Params.amount_of_coarse_nodes = start_threshold_on_nodes;

      AllLevelsGraphs.clear();
      AllLevelsEdgeCosts.clear();

      reset_all_vertices(G);

      V_Cycle_Number = i;
      //      cerr << "aga" << endl;
      basis_solution_iterations(G);

      //cerr << "After all basis solutions, |V|=" << G.number_of_nodes() << endl;

      double res = calc_real_laC(G);

      cerr << "After all basis solutions, cost = " << res << endl;
/*
      LA current_la;
      current_la.weight = res;
      current_la.V.init(G);
      forall_nodes(v, G) 
        current_la.V[v] = G[v].ArrId;

      results.push_back(current_la);
*/
//      best_previous_V_Cycle_LA = current_la;      
      cerr << "END OF " << i << " V-CYCLE\n" << endl;
    }

//  results.sort(&cmp_LA);
//  list_item la_it = results.first();

  // set the arrangement of G to the best
//  forall_nodes(v, G) {
//    G[v].ArrId = results[la_it].V[v];
    //cerr <<G[v].ArrId << endl;  
//  }

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

  double res = calc_real_laC(G);
  cerr << ">>>>>>>>> BEST COST = " << res << " <<<<<<<<<\n";

  //double mwfstart = calc_reset_mwf(G);
  //cerr << "mwfstart = " << mwfstart << endl;

  //int i=0;
  /*
  forall_nodes(v, G)
    {
      G[v].ArrId = G.number_of_nodes()-i;
      i++;
    }
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  double mwfstart_rev = calc_reset_mwf(G);
  cerr << "mwfstart = " << mwfstart_rev << endl;
  */




//cerr << "Workbound before min = " << calculate_wb_sval(G) << endl;
  //   double lll = wb_raise_power(G, 3, -1, 2);//insert_minimization_wb(G, ALL);
     //cerr << "Workbound after windows = " << lll << endl;
    // lll = insert_minimization_wb(G, ALL);

     // lll = wb_raise_power(G, 10, -1, 2);
     // lll = insert_minimization_wb(G, ALL);


     // lll = wb_raise_power(G, 10, -1, 2);
     // lll = insert_minimization_wb(G, ALL);



     //cerr << "Workbound after min = " << calculate_wb_sval(G) << endl;

     //     mwfstart=insert_minimization_wb(G, false);

//i=0;
  //   forall_nodes(v, G)
    //   {
	// G[v].ArrId = G.number_of_nodes()-i;
	 //i++;
       //}
  //TMP_CMP_GRAPHC = &G;
  //G.sort_nodes(&cmp_ArrId);
  //define_S_values(G);


       //mwfstart = calc_reset_mwf(G);
     //cerr << "mwfend = " << mwfstart << endl;


     //leda::string ttt="mwf";
     //     mwf_graph_print(G, ttt+instance_file);

  //  cerr << "ogo : " << results[la_it].weight << endl;

  //return results[la_it].weight;
  return res;
}

void calc_oldpw(TGraphC & G) {
	//  edge_array<double> ret(G);
	edge e;
	node v, u;
	/*
	 forall_edges(e, G)
	 if(((G[G.source(e)].color != green)&&(G[G.target(e)].color != green))||
	 ((G[G.source(e)].color == green)&&(G[G.target(e)].color == green)))
	 G[e].pw = fix_edge_relative_weight(G, e);//G[e].w;
	 */
	forall_edges(e, G)
	G[e].pw = 0;

	forall_nodes(v, G) {
		if (G[v].status!=seed) {
			double edges_to_green_sum = 0;
			forall_adj_edges(e, v) {
				u = second_adj_for_edge(e, v, G);
				if (G[u].status==seed)
					edges_to_green_sum+=fix_edge_relative_weight(G, e);
			}
			forall_adj_edges(e, v) {
				u = second_adj_for_edge(e, v, G);
				if (G[u].status==seed)
					G[e].pw = fix_edge_relative_weight(G, e)/edges_to_green_sum;
			}
		}
	}
}


void coarsening_preparations(TGraphC & G, TGraphC & H)
{
	if((Params.alg_dist==wag_and_gs)||(Params.alg_dist==gauss_seidel)||(Params.iw_calc_type==iw_alg))
	{
		//nodes_algdist_relaxation2d(G);
		//algdist_Jacobi_underrelax(G);
	  /*
	  if(Params.sym_nonsym==symmetric)
	    algdist_symJacobi_underrelax(G);
	  else
	    algdist_nonsymJacobi_underrelax(G);
	  */

	  //algdist_minla_relax(G);
	  algdist_symJacobi_underrelax(G);

	}

	double sum_fines_wd = sort_fine_nodes_by_wdegree(G);

	add_heavy_nodes(G, H, sum_fines_wd);

	if (Params.alg_dist==gauss_seidel)
	{
		coarse_nodes_alg_sum(G, H);
		//wag_coarse_nodes(G, H);
		//coarse_nodes_alg_max(G, H);
		//nodes_algdist_compatible_relaxation2d(G);

	}
	else if(Params.alg_dist==wag_and_gs)
		coarse_nodes_alg_and_wag_sum(G, H);
	else if(Params.alg_dist==jacobian)
	{
		cerr << "Jacobi is not supported" << endl;
		exit(1);
	}
	else if(Params.alg_dist==noalgdist)
	{
		//nodes_algdist_relaxation2d(G);
		//coarse_nodes_alg_sum(G, H);
		wag_coarse_nodes(G, H);
	}

	calc_oldpw(G);
	clear_iw_links(G);

	//if((Params.alg_dist==wag_and_gs)||(Params.alg_dist==gauss_seidel))
	//	nodes_algdist_compatible_relaxation2d(G);

	if(Params.iw_calc_type==iwbamg)
	{
		exit(1);
		if(interpol_order(G)>1)
		{
			//wag_calc_iw(G);
			new_bamg_alg_calc_iw_by_rw(G);
			//bamg_calc_iw(G);
		}
		else
			minla_alg_calc_iw_by_rw(G);
			//alg_calc_iw_by_rw_and_w(G);
			//new_bamg_alg_calc_iw_by_rw(G);
	}
	else if(Params.iw_calc_type==iw_alg)
	{
		//nodes_algdist_compatible_relaxation2d(G);
		alg_calc_iw_by_rw(G, H);
		//alg_calc_iw_by_rw_and_w(G);
	}
	else if(Params.iw_calc_type==iw_alg_and_wag)
		alg_calc_iw_by_rw_and_w(G);
	else if(Params.iw_calc_type==iw_wag)
	{
		//Params.alg_dist=gauss_seidel;
		//nodes_algdist_relaxation2d(G);
		//new_bamg_alg_calc_iw_by_rw(G);
		//alg_calc_iw_by_rw(G);
		//Params.alg_dist=noalgdist;
		wag_calc_iw(G);
	}
}

double method_s3( TGraphC & G,
                  bool first_time_down /* for multiple basis solutions */)

{
  clock_t start, end;
  node a;

  Params.CURRENT_LEVEL++;
  double avg_deg = print_graph_info(G);
  //    if(S==1)
  //    {
  //      graph_print(G, "order :");

  //    }
  //  if(first_time_down==true)
  //    AllLevelsGraphs.push_back(G);



  //  cerr << "Entering into the level " << S << ", nodes="<<G.number_of_nodes()<<", edges=" << G.number_of_edges() << endl;
  node v, u;
  edge e;

  //  check_node_weights(G);

  int pos, i, j;
  list_item it;
  double lin_arr;

  leda::string fn = i2s(Params.CURRENT_LEVEL);
  fn = fn+".gml";

  int level = Params.CURRENT_LEVEL;

  // Construct smaller graph H
  if(G.number_of_nodes() > Params.min_graph_size)
    {
      TGraphC H;
      //      edge_array<double> new_edge_costs(G);

      if(first_time_down==true)
        {

          clock_t start = clock();
	coarsening_preparations(G, H);

          construct_new_graph(G, H);
          cerr << "Coarse graph construction. time (sec) : " << ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC << endl;


          if(Params.CURRENT_LEVEL==0)
            {
              cerr << "resaving : " << G.number_of_edges() << endl;
              E0=G.number_of_edges();
              V1=H.number_of_nodes();
              A0=avg_deg;
            }


          start = clock();
          //          cerr << "Vertex weight : Start" << endl;
          calculate_vertex_weights(G, H);
          cerr << "Coarse nodes weights def. time (sec) : " << ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC << endl;
          update_time();

          cerr << "TOTAL CPU TIME = " << cpu_time << " M." << endl;

/*
          if(G.number_of_nodes()<1000)
            {
              cerr << "SCC : Start" << endl;
              node_array<int> compnum(H);
              int scc = STRONG_COMPONENTS(H, compnum);
              if(scc!=1)
                cerr << ">>>>>>> Error : SCC=" << scc << endl;
            }
*/
          if((Params.basis_solutions_number>1)&&(Params.V_Cycle_Iter_Num>1))
            //          if(Params.basis_solutions_number>1)
            {
              cerr << "Store a fine graph : Start" << endl;
              AllLevelsGraphs.push_back(G);
              //              G.clear();
            }
          else
            {
              //  exit(1);
              TGraphC Dummy;
              Dummy.clear();
              AllLevelsGraphs.push_back(Dummy);
            }

          //          gml_write_with_iw("opsa_iw_"+i2s(S)+".gml", G);
          //          gml_write_with_pw("opsa_pw_"+i2s(S)+".gml", G);
          //          gml_write_with_w("opsa_w_"+i2s(S)+".gml", G);
        }
      else
        {
          cerr << "Current S=" << Params.CURRENT_LEVEL << endl;
          H = AllLevelsGraphs[AllLevelsGraphs.get_item(Params.CURRENT_LEVEL+1)];
        }


      double H_lin_arr = method_s3(H, first_time_down);

      Params.CURRENT_LEVEL--;
      //      if(Params.basis_solutions_number>1)
        //        {
      //  G = AllLevelsGraphs[AllLevelsGraphs.get_item(S)];
      //  cerr << "Restored G : " << G.number_of_nodes() << endl;
      //        }

      cerr << "S=" << Params.CURRENT_LEVEL << "; Level " << level << " interpolation and relaxation" << endl;
      //      graph_print(H,"COARSE BEFORE INTERPOLATION : " );

      start = clock();
      interpolate_s2(G, H);
      //      cerr << "Interpolation time (sec) : " << ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC << endl;

      //      check_same_coordinates(G);

      lin_arr = linear_order_by_coord(G);

      //    recalculate_edges(G);

       cerr << "Interpolation time + Nodes pumping time (sec) : " << ((double)(clock() - start)) / (double)CLOCKS_PER_SEC << endl;
      cerr << "cost after interpol = " << lin_arr << endl;
update_time();

      dummyNodesArr.resize(G.number_of_nodes());
      int p=0;
      forall_nodes(v, G)
        {
          dummyNodesArr[p] = v;
          p++;
        }

      //insert_minimization(G, true);
update_time();
      //      forall_nodes(v, G)
      //        G[v].after_interpol_ArrId = G[v].ArrId;

      //      print_diff_save_ArrId_and_ArrId(G);

      coords_by_linear_order(G);


      if(Params.check_best_relax == true)
        {
        //lin_arr = minimization(G);
          /*
          relaxation_of_fines(G);
          lin_arr = linear_order_by_coord(G);
          lin_arr = minimization(G);
          //          cerr << "1" << endl;
          if(Params.use_sa_in_relax==true)
            lin_arr = new_sa_on_distance(G, 1, true);
          //            lin_arr = sa_on_fines_distance(G, 1);
          lin_arr = minimization(G);
          //          cerr << "2" << endl;

          node_array<int> A(G); forall_nodes(v, G) A[v]=G[v].ArrId;
          double after_comp_rel = lin_arr;

          cerr << "ENERGY AFTER COMP.REL. " << lin_arr << endl;

          relaxation_gauss_seidel(G);

          lin_arr = linear_order_by_coord(G);
          lin_arr = minimization(G);
          //          cerr << "3" << endl;
          if(Params.use_sa_in_relax==true)
          lin_arr = new_sa_on_distance(G, 1, true);
          //lin_arr = sa_on_fines_distance(G, 1);
          lin_arr = minimization(G);
          //          cerr << "4" << endl;

          double after_gs_rel = lin_arr;
          cerr << "ENERGY AFTER G-S REL. " << lin_arr << endl;

          if(after_gs_rel>after_comp_rel)
            forall_nodes(v, G) G[v].ArrId=A[v];

          lin_arr = minimization(G);
         */
        }
      else
        {
          if(Params.use_compatible_relaxation==true)
            relaxation_of_fines(G);
          if(Params.use_gauss_seidel_relaxation==true)
            //            relaxation_gauss_seidel(G);
            relaxation_of_gs(G);
          //          lin_arr = linear_order_by_coord(G);
          //if(S==0)
          //lin_arr = segment_minimization(G, lin_arr);
          cerr << "After segment min : " << lin_arr << endl;
          lin_arr = minimization(G);
        }

   update_time();

   //   if(S==0)
   //     {
       cout << "ENERGY AFTER (COMP/G-S)-RELAXATIONS " << lin_arr << endl;
       //       cerr << "Random arr = " << G.number_of_nodes() * G.number_of_edges() / 3.0 << endl;
       //       exit(1);
       //     }

      //      lin_arr = insert_minimization(G);

      //      check_same_coordinates(G);
      //      print_diff_save_ArrId_and_ArrId(G);
      //      exit(1);
      //      graph_print(G, "AFTER RELAXATION :");

      // start = clock();

      //      for(int y=0; y<10; y++)
      //      segment_relaxation(G);
       cerr << "Before segment minimization = " << lin_arr << endl;
       //        lin_arr = segment_minimization(G, lin_arr);
        cerr << "After segment minimization = " << lin_arr << endl;

       if((Params.insert_sa_dist>0)&&(((Params.use_sa_only_in_last_two_levels == true) && (Params.CURRENT_LEVEL<2))||(Params.use_sa_only_in_last_two_levels == false)))
        {
          long t = time(0);

          //          if(G.number_of_nodes()>11)
            lin_arr = insert_sa(G, lin_arr);

          //if(Params.do_sa_1==true)
          //            lin_arr = new_sa_on_distance(G, 1, false);

           //if(Params.do_sa_2==true)
           // lin_arr = new_sa_on_distance(G, 2, false);
          cerr << "Time for entire SA (sec) : " << ((double)time(0) - (double)t) << endl;
        }

      update_time();
      cerr << "CURRENT TOTAL CPU TIME = " << cpu_time << " M." << endl;

      //  relaxation_of_gs(G);

      //      lin_arr = triples_minimization(G);
      //      lin_arr = wmm_minimization(G);

      //      cerr << "opsa" << endl;
      //      lin_arr = adjacent_minimization(G, lin_arr, 100);

      //      if (level==0)
      //        G.sort_nodes(&cmp_nC);
      /*
      graph_print(G, "ogogo : ");
      cerr << "cost = " << lin_arr << endl;
      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_initial_id);
      int t=1;
      forall_nodes(v, G)
        {
          G[v].ArrId = t;
          t++;
        }
      define_S_values(G);
      lin_arr = calc_laC(G);
      */
      //      graph_print(G, "ogogo : ");
      //      cerr << "cost = " << lin_arr << endl;

      //      check_placement(G);
      //      check_placement(G);
      //      check_placement(G);
      /*
      if(S==0)
        {
          edge f;
          forall_edges(f, G)
            G[f].w = 1;
          lin_arr = calc_laC(G);

        }
      */

      //  if(S==0) graph_print(G, "not final : ");

      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
      lin_arr = segment_minimization(G, lin_arr);

      relaxation_of_gs(G);
      lin_arr = minimization(G);
      */


              //  lin_arr = contraction_minimization(G, lin_arr);




      //if(S==0) graph_print(G, "final : ");


                  //if(S==0)

      //  {
      /*
          double t=0.2;
          for(int i=1; i<51; i++)
            {
              t+=0.3*(double)i/50.0;

              lin_arr = segment_min(G, lin_arr, t);
              minimization(G);

              // relaxation_of_gs(G);minimization(G);
            }
      */
      //      graph_print(G, "final : ");
          /*
          if(S==0)
                 {
       lin_arr = segment_min(G, lin_arr);
          minimization(G);
          //relaxation_of_gs(G);minimization(G);
        }
          */
          //          lin_arr = segment_minimization(G, 8);
          //          minimization(G);
          cerr << "AFTER SEG MIN = " << lin_arr << endl;
          //        }

    }
  else if(basis_solutions_passed==0)
    {
      basis_solutions.resize(Params.CURRENT_LEVEL+1);
      solutions_on_level.resize(Params.CURRENT_LEVEL+1);
      for(int i=solutions_on_level.size()-1; i>0; i--)
        if(i==solutions_on_level.size()-1)
          solutions_on_level[i] = Params.basis_solutions_number;
        else if(i==0)
          solutions_on_level[i] = Params.basis_solutions_number;//1;
        else
          solutions_on_level[i] = Params.basis_solutions_number;//(int)((float)solutions_on_level[i+1]/2.0);

      AllLevelsGraphs.push_back(G);

      MAX_LEVEL = level;

      node v;
      lin_arr = solve_exact_minla_C(G);

      list_item it = basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].L.first();

      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_initial_id);

      forall_nodes(v, G)
        {
          G[v].ArrId = basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].L[it];
          it = basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].L.succ(it);
        }

      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_ArrId);
      define_S_values(G);

      //      cerr << "Takes basis solution with C=" << calc_laC(G) << endl;


      cerr << "order of vertices " ;
      forall_nodes(v, G)
        cerr << G[v].initial_id << "(" << G[v].w << "), ";
      cerr << endl;

      cerr << "Takes basis solution with C=" << calc_laC(G) << endl;

      basis_solutions[Params.CURRENT_LEVEL].erase(basis_solutions[Params.CURRENT_LEVEL].first());
      basis_solutions_passed++;
    }
  else
    {
      cerr << "Ready basis solution" << endl;
      list_item it = basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].L.first();

      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_initial_id);

      forall_nodes(v, G)
        {
          G[v].ArrId = basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].L[it];
          it = basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].L.succ(it);
        }

      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_ArrId);

      define_S_values(G);

      forall_nodes(v, G)
        cerr << G[v].initial_id << " (" << G[v].w << ")\n";
      cerr << endl;

      cerr << "Takes basis solution with C=" << calc_laC(G) << endl;


      basis_solutions[Params.CURRENT_LEVEL].erase(basis_solutions[Params.CURRENT_LEVEL].first());
      basis_solutions_passed++;
    }
  return lin_arr;
}

