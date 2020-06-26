#include <time.h>
#include "cmpfuncs.h"
#include "random_variate.h"
//#include "LEDA/dictionary.h"
//#include "LEDA/matrix.h"
#define SEGMENT_SIZE 10
#define ONLY_FINE true
#define ALL false
//#define MY_MAX_INT 100000000
/*
double Insert_SA_Avg_Change[21];
double Insert_SA_Counter[21];
double Current_T[21];
*/

extern double E0;

array<list<double> > EChanges;
array<node> dummyNodesArr;

void volumeorder2interpolcoord(TGraphC & G)
{
  linear_order_by_coord(G);
  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

  node v;
  forall_nodes(v, G)
    G[v].interpol_coord = G[v].S_value;
}
/*
void get_vertex_median_bounds(node & v, TGraphC & G, int & left, int & right)
{
  node w;
  edge e;
  list<two_tuple<edge,node> > l;

  forall_adj_edges(e, v)
    {
      w = second_adj_for_edge(e, v, G);
      two_tuple<edge, node> two(e, w);
      l.push_back(two);
    }
  if(l.length()==1)
    {
      
    }
    G[v].interpol_coord = G[l[l.first_item()].second()].interpol_coord;
}
*/

// Returns previous to median element
list_item list_median(list<two_tuple<edge,node> > & lst, TGraphC & G, double esum, double & min_l, double & min_r)
{
  list_item it, min_ptr = nil;

  double l = 0;
  double r = esum;
  double min_dist = esum;
  forall_items(it, lst)
    {
      //      cerr << "x" << endl;
      edge e = lst[it].first();
      double ecost = G[e].w;// * fabs(G[G.source(e)].interpol_coord - G[G.target(e)].interpol_coord);
      l+=ecost;
      r-=ecost;

      //      cerr << l << "\t" << r << "\t" << min_dist << endl;
      if(fabs(l-r)<=min_dist)
        {
          min_dist = fabs(l-r);
          min_ptr = it;
          min_l = l;
          min_r = r;
        }
      else
        {
          return min_ptr;
        }
    }

  //  cerr << "LIST MEDIAN : return median" << endl;
  min_l = 0;
  min_r = esum;
  return lst.get_item(floor(lst.size()/2.0));
}

void minimize_seg(node b1, node b2, TGraphC & G, double & thr)
{
  double p1 = G[b1].S_value; double p2 = G[b2].S_value;
  list<two_tuple<edge,node> > l;

  node v = b1; node w;
  edge e;
  double lsum=0;
  double rsum=0;
  double esum=0;
  while((v!=nil)&&(v!=G.succ_node(b2)))
    {
      //      cerr << "For " << G[v].initial_id << " (" << G[v].ArrId << ") deg=" << G.degree(v) << endl;
      forall_adj_edges(e, v)
        {
          w = second_adj_for_edge(e, v, G);
          if((G[w].S_value<p1)||(G[w].S_value>p2))
            {
              two_tuple<edge, node> two(e, w);
              l.push_back(two);

              if(G[w].S_value<p1)
                {
                  lsum+=fabs(G[v].S_value - G[w].S_value)*G[e].w;
                  esum+=G[e].w;
                }
              if(G[w].S_value>p2)
                {
                  rsum+=fabs(G[v].S_value - G[w].S_value)*G[e].w;
                  esum+=G[e].w;
                }
            }
          
        }
      v = G.succ_node(v);
    }
  //  cerr << "REL : " << lsum/(lsum+rsum) << "\t" << rsum/(lsum+rsum) << endl;
  if(((lsum/(lsum+rsum)>thr)&&(rsum/(lsum+rsum)>thr))||(lsum/(lsum+rsum)==0)||(rsum/(lsum+rsum)==0))
    return;
  else
    {
    }
    //    cerr << "l size : " << l.size() << endl;
  
  TMP_CMP_GRAPHC = &G;
  l.sort(&cmp_two_tuple);
  l.unique(&cmp_two_tuple);
  //  cerr << "l size : " << l.size() << endl;
  //  cerr << "SIZE : " << fabs(G[b1].ArrId - G[b2].ArrId) << endl;
  
  double right_edges, left_edges;
  list_item mid_ptr = list_median(l, G, esum, left_edges, right_edges);

  //  cerr << mid_ptr << endl;

  list_item it2;
  
  if(right_edges>left_edges)
    {
      if(mid_ptr==l.last())
        it2=l.last();
      else
        it2 = l.succ(mid_ptr);
      //      forall_items(it2, l)
      //        if(G[l[it2].second()].S_value > G[b2].S_value)
      //          break;
      
      //          cerr << "r:near to " << G[l[it2].second()].initial_id << endl;
      //      G[l[it2].second()].lh_pts++;
      double p=G[l[it2].second()].S_value - 1.0/100000.0;
      v = b1;
      while((v!=nil)&&(v!=G.succ_node(b2)))
        {                  
          G[v].interpol_coord = p+0.00000001;
          p+=0.00000001;
          v = G.succ_node(v);
        }
    }
  else if(right_edges<left_edges)
    {
      it2 = mid_ptr;
      //      forall_items(it2, l)
      //        if((G[l[it2].second()].S_value < G[b1].S_value)&&(G[l[l.succ(it2)].second()].S_value>G[b2].S_value))
      //          break;
      
      //            cerr << "l:near to " << G[l[it2].second()].initial_id << "(" << G[l[it2].second()].ArrId << ")" << endl;
      //      G[l[it2].second()].rh_pts++;
      double p = G[l[it2].second()].S_value + 1.0/100000.0 ;
      v = b1;
      while((v!=nil)&&(v!=G.succ_node(b2)))
        {                  
          G[v].interpol_coord = p+0.00000001;
          //                    cerr << "\t Node " << G[v].initial_id << " got " << G[v].interpol_coord << endl;
          p+=0.00000001;
          v = G.succ_node(v);
        }
      
    }
  //  cerr << "end\n";
}

node NEXT_NODE(node & v, TGraphC & G)
{
  if(v==nil)
    return nil;
  else
    return G.succ_node(v);
}

double segment_min(TGraphC & G, double la, double thr)
{
  double s = 0;
  edge e;
  forall_edges(e, G)
    s+=G[e].w;
  s = s/(double)G.number_of_edges();

  double EW = 0;
  forall_edges(e, G)
    EW+=(G[e].w - s)*(G[e].w - s);
  EW = sqrt(EW / (double)G.number_of_edges());
  //cerr << s << "\t" << EW << endl;
  EW = s-2*EW;

  //  double EW = s / 5.0;
  
  //coords_by_lin_order(G);
  node x;
  forall_nodes(x, G)
    {
      G[x].seg_min_save_ArrId = G[x].ArrId;
      //      G[x].lh_pts = 0;
      //      G[x].rh_pts = 0;
    }
  
  double save_lin_arr = la;
  
  coords_by_linear_order(G);
  node b1 = G.first_node(); node b2 = G.succ_node(b1);
 

  while(b2!=nil)
    {
      if(((e=is_edge(b2, G.pred_node(b2) , G))!=nil)&&(G[e].w >= EW))
        b2 = NEXT_NODE(b2,G);              
      else
        {
          b2 = G.pred_node(b2);
          if((G[b2].ArrId - G[b1].ArrId>1))
            {
              //                            cerr << "Segment len = " << G[b2].ArrId - G[b1].ArrId << " from " << G[b1].initial_id << " to " << G[b2].initial_id << endl;
                minimize_seg(b1, b2, G, thr);
              //              double lin_arr = linear_order_by_coord(G);
              //   lin_arr = coords_by_linear_order(G);
            }
          b1 = NEXT_NODE(b2,G);
          b2 = NEXT_NODE(b1,G);          
        }
    }
  
   double lin_arr = linear_order_by_coord(G);
   lin_arr = minimization(G);
   
   //   double lcc_cost = save_lin_arr;
   
   //   lcc_update(G, lin_arr, lcc_cost);

     
   if(lin_arr - save_lin_arr/400.0> save_lin_arr)
     {
       forall_nodes(x, G)
         G[x].ArrId = G[x].seg_min_save_ArrId;

       TMP_CMP_GRAPHC = &G;
       G.sort_nodes(&cmp_ArrId);
       define_S_values(G);
       lin_arr =  save_lin_arr;
     }
   //   else if((lcc_cost<lin_arr)&&(lcc_cost<save_lin_arr))
   //    {
   //      restore_best_lcc_arrangement(G);
   //      lin_arr = calc_laC(G);      
   //    } 
   //   else
   lin_arr = coords_by_linear_order(G);

   
   return lin_arr;
}

double segment_minimization(TGraphC & G, double lin_arr)
{
  if(Params.segment_minimization <1)
    return lin_arr;
  
  node x;
  forall_nodes(x, G)
    G[x].seg_min_mainsave_ArrId = G[x].ArrId;

  lcc_init(G);
  
  double t=0.1;
  int c = Params.segment_minimization;
  double upto = 0.5;

  double lacost = lin_arr;
  double lcc_cost = lin_arr;
  double save_lin_arr = lin_arr;
  
  for(int i=1; i<c+1; i++)
    {
      t+=upto*(double)i/(double)c;

      //      cerr << "\tEntering segment_min, lacost=" << lacost << ", t=" << t << endl;
      lacost = segment_min(G, lacost, t);
      minimization(G);

      lcc_update(G, lacost, lcc_cost);
    }

  if((lcc_cost<lacost+1)&&(lcc_cost<save_lin_arr+1))
    {
      restore_best_lcc_arrangement(G);
      lacost = calc_laC(G); 
    }
  else if((save_lin_arr < lacost)&&(save_lin_arr < lcc_cost))
    {      
      forall_nodes(x, G)
         G[x].ArrId = G[x].seg_min_mainsave_ArrId;

       TMP_CMP_GRAPHC = &G;
       G.sort_nodes(&cmp_ArrId);
       define_S_values(G);
       lacost =  save_lin_arr;
     
       lacost = coords_by_linear_order(G);

      //      cerr << "Error : LCC is bigger " << lacost << "\t" << save_lin_arr << "\t" << lcc_cost << endl;
      //      exit(1);
    }
  else
    lacost = coords_by_linear_order(G);
    
    

  return lacost;
}
/*
double contraction_minimization(TGraphC & G, double lin_arr)
{
  double s = 0;
  edge e, f, g; node y;
  forall_edges(e, G)
    s+=G[e].w;
  s = s/(double)G.number_of_edges();
  double EW = s;

  cerr << "MED EDGE : " << EW << endl;

  TGraphC H = G;
  forall_nodes(y, H)
    H[y].contr_nodes.push_back(y);

  for(int k=0; k<1; k++)
    {
      cerr << k << "-------------------\n";
      array<edge> E(H.number_of_edges());
      edge_array<bool> Eb(H);
      forall_edges(e, H)
        Eb[e] = true;

      int l=0;
      forall_edges(e, H)
        {
          E[l]=e;
          l++;
        }
      
      TMP_CMP_GRAPHC = &H;
      E.sort(&cmp_ew);
      
      for(int i=0; i<E.size(); i++)
        {
          //cerr << i << endl;
          e = E[i];
          if((Eb[e] == true)&&(H[e].is_new!=true))
            {
              //              cerr << "w="<<H[e].w<<endl;
              
              node v = H.source(e); node w = H.source(e);
              if((fabs(H[v].ArrId - H[w].ArrId) == 1)&&(H[e].w >= EW))
                {
                  Cnode wv;
                  wv.w = H[w].w + H[v].w;
                  wv.ArrId = H[w].ArrId;
                  wv.contr_nodes = H[w].contr_nodes;
                  list_item it;
                  forall_items(it, H[v].contr_nodes)
                    wv.contr_nodes.push_back(H[v].contr_nodes[it]);
                  
                  node x = H.new_node(wv);
                  forall_adj_edges(f, w)
                    if((y=second_adj_for_edge(f, w, H))!=v)
                      {
                        CEdge ne; ne.w = H[f].w; ne.is_new = true;
                        if((g=is_edge(x,y,H))==nil)                          
                            H.new_edge(x,y,ne);                          
                        else                          
                          H[g].w+=ne.w;
                        Eb[g] = false;
                      }
                  
                  forall_adj_edges(f, v)
                    if((y=second_adj_for_edge(f, v, H))!=w)
                      {
                        CEdge ne; ne.w = H[f].w; ne.is_new = true;
                        if((g=is_edge(x,y,H))==nil)                          
                          H.new_edge(x,y,ne);
                        else                          
                          H[g].w+=ne.w;
                        Eb[g] = false;
                      }
                  H.del_node(w);
                  H.del_node(v);
                }
              Eb[e] = false;
            }
          
          //int j=1;
          //forall_nodes(y, H)
           // {
            //  H[y].ArrId = j;
             // j++;
           // }
          
        }

      //      cerr << "a " << endl;
      TMP_CMP_GRAPHC = &H;
      H.sort_nodes(&cmp_ArrId);
      int j=1;
      node z;
      forall_nodes(z, H)
        {
          H[z].ArrId = j;
          j++;
        }
      //      define_S_values(H);
    }
  define_S_values(H);
  relaxation_of_gs(H);
  double laH = minimization(H);
  cerr << laH << endl;

  list_item it;
  node v;
  forall_nodes(v, H)
    {
      double p = H[v].ArrId;
      forall_items(it, H[v].contr_nodes)
        {
          G[H[v].contr_nodes[it]].dArrId = p + 0.000001;
          p+=0.000001;
        }
    }

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_dArrId);
  int j=1;
  forall_nodes(v, G)
    {
      G[v].ArrId = j;
      j++;
    }
  define_S_values(G);
  double la = minimization(G);

  return la;
}
*/
/*
double segment_minimization(TGraphC & G, int len)
{
  coords_by_linear_order(G);
  //return la;
  node b1 = G.succ_node(G.first_node()); node b2;

  node v = b1; node w = G.succ_node(b1);
  for(int i=0; i<len-1; i++)
    w = G.succ_node(w);
  
    //  w = G.succ_node(w);w = G.succ_node(w);w = G.succ_node(w);w = G.succ_node(w);w = G.succ_node(w);
    //  v=w; w = G.succ_node(v);
    //  w = G.succ_node(w);w = G.succ_node(w);w = G.succ_node(w);w = G.succ_node(w);w = G.succ_node(w);
  
  
  while((v!=nil)&&(w!=G.pred_node(G.last_node())))
    {
      
    //      if(is_edge(v, w , G)!=nil)
    //        {
    //          v = G.succ_node(v);
    //          w = G.succ_node(w);
    //        }
    //      else
      
        {
          //  if(b1!=b2)
            {
              //  b2 = v;
              //              cerr << "aaaa" << endl;
              //              if((fabs(G[b1].ArrId-G[b2].ArrId)<20)&&(fabs(G[b1].ArrId-G[b2].ArrId)>5))
              //                {
                            node n1 = G.pred_node(v); node n2 = G.succ_node(w);
                            if((is_edge(n1,v, G)==nil)&&(is_edge(n2,w,G)==nil))
                {
                  //                  cerr << "eof segment : " << G[w].initial_id << endl;
                  minimize_seg(v, w, G);
                  double lin_arr = linear_order_by_coord(G);
                  //                  cerr << "LO by C " <<  lin_arr << endl;
                  lin_arr = coords_by_linear_order(G);
                  //                  cerr << "C by LO " <<  lin_arr << endl;
                  //                  return lin_arr;
                } 
                  //                  cerr << lin_arr << endl;
                  //                }
              //              cerr << "bbbb" << endl;
                  //              v = w;
                  //              w = G.succ_node(v);
                  //              b1 = v;
            }
        }

        v = G.succ_node(v);
        w = G.succ_node(w);
 
    }

  double lin_arr = linear_order_by_coord(G);
  return lin_arr;
}
*/

void /*old_good_*/minimize_vertex_coordinate(node v, TGraphC & G)
{
  node w;
  G[v].interpol_coord = 0;
  double sum_of_edges = 0;
  list<two_tuple<edge,node> > l;
  edge e;
  forall_adj_edges(e, v)
    {
      w = second_adj_for_edge(e, v, G);
      //                  G[v].interpol_coord+=G[e].w*G[w].interpol_coord;
      //                  sum_of_edges+=G[e].w;
      two_tuple<edge, node> two(e, w);
      l.push_back(two);
    }
  if(l.length()==1)
    G[v].interpol_coord = G[l[l.first_item()].second()].interpol_coord;
  else
    {
      TMP_CMP_GRAPHC = &G;
      l.sort(&cmp_two_tuple);
      
      list_item it1 = l.first_item();
      double left_edges = 0;
      double right_edges = 0;
      list_item it2 = l.succ(it1); while(it2!=nil) {right_edges+=G[l[it2].first()].w; it2 = l.succ(it2);}
      double min_diff = fabs(right_edges - left_edges);
      it2 = l.succ(it1);
      
      G[v].interpol_coord = -1;
      while(it2!=nil)
        {
          left_edges+=G[l[it1].first()].w;
          if(fabs(right_edges - left_edges)<min_diff)
            min_diff = fabs(right_edges - left_edges);
          else
            {
              G[v].interpol_coord = G[l[it1].second()].interpol_coord+(G[l[it2].second()].interpol_coord-G[l[it1].second()].interpol_coord)/100.0;
              break;
            }
          right_edges-=G[l[it2].first()].w;
          if(fabs(right_edges - left_edges)<min_diff)
            min_diff = fabs(right_edges - left_edges);
          else
            {
              G[v].interpol_coord = (G[l[it1].second()].interpol_coord+G[l[it2].second()].interpol_coord)/2.0;
              break;
            }
          it2 = l.succ(it2);
          it1 = l.succ(it1);
        }
      if(G[v].interpol_coord == -1)
        G[v].interpol_coord = G[l[it1].second()].interpol_coord-fabs(G[l[it1].second()].interpol_coord-G[l[l.pred(it1)].second()].interpol_coord)/100.0;
      
    }
  
  return;
}
void new_minimize_vertex_coordinate(node v, TGraphC & G)
{
  node w;
  G[v].interpol_coord = 0;
  double sum_of_edges = 0;
  
  list<two_tuple<edge,node> > l;
  
  edge e;
  forall_adj_edges(e, v)
    {
      w = second_adj_for_edge(e, v, G);
      two_tuple<edge, node> two(e, w);
      l.push_back(two);
      sum_of_edges+=G[e].w;
    }
  
  if(l.length()==1)
    G[v].interpol_coord = G[l[l.first_item()].second()].interpol_coord;
  else
    {
      TMP_CMP_GRAPHC = &G;
      l.sort(&cmp_two_tuple);
      
      list_item it1 = l.first_item();
      double left_edges = G[l[it1].first()].w;
      double right_edges = sum_of_edges - left_edges;
      double min_diff = fabs(right_edges - left_edges);
      list_item it2 = l.succ(it1);
      
      G[v].interpol_coord = -1;
      while(it2!=nil)
        {
          //          cerr << left_edges << "\t" << right_edges << endl;
          if(fabs(right_edges - left_edges)<=min_diff)
            min_diff = fabs(right_edges - left_edges);
          else
            {
              left_edges-=G[l[it1].first()].w;
              right_edges+=G[l[it1].first()].w;
              it1 = l.pred(it1); it2 = l.pred(it2);
              
              //              G[v].interpol_coord = G[l[it1].second()].interpol_coord+(G[l[it2].second()].interpol_coord-G[l[it1].second()].interpol_coord)/100.0;
              break;
            }

          it2 = l.succ(it2);
          it1 = l.succ(it1);
          if(it2!=nil)
            {
              left_edges+=G[l[it1].first()].w;
              right_edges-=G[l[it1].first()].w;
            }
        }
      if(it2==nil)
        { it1 = l.pred(it1); it2 = l.last(); }
      //      cerr << "done " << left_edges << "\t" << right_edges << endl;
      if(right_edges>left_edges)
        G[v].interpol_coord = G[l[it2].second()].interpol_coord - (G[l[it2].second()].interpol_coord - G[l[it1].second()].interpol_coord)/100000.0;        
      else if(right_edges<left_edges)
        G[v].interpol_coord = G[l[it1].second()].interpol_coord + (G[l[it2].second()].interpol_coord - G[l[it1].second()].interpol_coord)/100000.0 ;       
      else
        G[v].interpol_coord = G[l[it1].second()].interpol_coord + (G[l[it2].second()].interpol_coord - G[l[it1].second()].interpol_coord)/2.0;
    }
      if(G[v].interpol_coord == -1)
        {
          cerr << "minimize vertex : G[v].interpol_coord = -1" << endl;
          exit(1);
          //          G[v].interpol_coord = G[l[it1].second()].interpol_coord-fabs(G[l[it1].second()].interpol_coord-G[l[l.pred(it1)].second()].interpol_coord)/100.0;
        }
    
      //      cerr << "vishli" << endl;
  return;
}


// END : ordered, interpol_coord are calculated
void relaxation_of_fines(TGraphC & G)
  // double_dummy = prev_interpol_coord
{
  edge e;
  node v, w;

  double prev_la = calc_laC(G);
  forall_nodes(v, G)
    G[v].double_dummy = G[v].interpol_coord;
  
  cerr << "before compatible relaxation : " << prev_la << endl;
  clock_t start =  clock();
  int i;
  for(i=0; i<Params.comp_relaxation_sweeps+Params.CURRENT_LEVEL; i++)
    {
      
      forall_nodes(v, G)
        {
          if(G[v].status!=seed)
            {
              minimize_vertex_coordinate(v, G);
            }
        }
      
      double lin_arr = linear_order_by_coord(G);
     
      
      //      lin_arr = minimization(G);
      cerr << "cost after compatible step " << i << ", arr= " << lin_arr << "\n";
      //      lin_arr = insert_minimization(G, ONLY_FINE);
      
      coords_by_linear_order(G);
      
      node q;
      if(/*(i>0)&&*/(prev_la<lin_arr))
        {
          forall_nodes(q, G)
            G[q].interpol_coord = G[q].double_dummy;
          double lin_arr = linear_order_by_coord(G);
          cerr << "cost after reverse compatible step " <<  ", arr= " << lin_arr << endl; 
          break;
        }
      prev_la = lin_arr;
      
      forall_nodes(q, G)
        G[q].double_dummy = G[q].interpol_coord;
      //      coords_by_linear_order(G);
      
      //      print_diff_save_ArrId_and_ArrId(G);
      
      //    cerr << "Compatible relax. step time (sec) : " << ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC << endl;
    }
  insert_minimization(G, ONLY_FINE);
  cerr << "Compatible relax. time (sec) for " << i << " sweeps : " << ((double)(clock() - start)) / (double)CLOCKS_PER_SEC << endl;
}
/*
double imaginary_minimize_vertex_coordinate(node v, TGraphC & G)
{
  double interpol_coord_before = G[v].interpol_coord;
  double interpol_coord_after  = 0;
  
  node w;
  //  G[v].interpol_coord = 0;
  double sum_of_edges = 0;
  list<two_tuple<edge,node> > l; // list of all neighbors+edges
  edge e;
  forall_adj_edges(e, v)
    {
      w = second_adj_for_edge(e, v, G);
      //                  G[v].interpol_coord+=G[e].w*G[w].interpol_coord;
      //                  sum_of_edges+=G[e].w;
      two_tuple<edge, node> two(e, w);
      l.push_back(two);
    }
  if(l.length()==1)
    interpol_coord_after = G[l[l.first_item()].second()].interpol_coord;
  else
    {
      TMP_CMP_GRAPHC = &G;
      l.sort(&cmp_two_tuple);
      
      list_item it1 = l.first_item();
      double left_edges = 0;
      double right_edges = 0;
      list_item it2 = l.succ(it1); while(it2!=nil) {right_edges+=G[l[it2].first()].w; it2 = l.succ(it2);}
      double min_diff = fabs(right_edges - left_edges);
      it2 = l.succ(it1);
      
      interpol_coord_after = -1;
      while(it2!=nil)
        {
          left_edges+=G[l[it1].first()].w;
          if(fabs(right_edges - left_edges)<min_diff)
            min_diff = fabs(right_edges - left_edges);
          else
            {
              interpol_coord_after = G[l[it1].second()].interpol_coord+(G[l[it2].second()].interpol_coord-G[l[it1].second()].interpol_coord)/100.0;
              break;
            }
          right_edges-=G[l[it2].first()].w;
          if(fabs(right_edges - left_edges)<min_diff)
            min_diff = fabs(right_edges - left_edges);
          else
            {
              interpol_coord_after = (G[l[it1].second()].interpol_coord+G[l[it2].second()].interpol_coord)/2.0;
              break;
            }
          it2 = l.succ(it2);
          it1 = l.succ(it1);
        }
      if(interpol_coord_after == -1)
        interpol_coord_after = G[l[it1].second()].interpol_coord-fabs(G[l[it1].second()].interpol_coord-G[l[l.pred(it1)].second()].interpol_coord)/100.0;
      
    }

  
  return fabs(interpol_coord_before - interpol_coord_after);
}
*/
void construct_vertex_queue(TGraphC & G, array<node> & VQ)
{
  node v,w;
  int i=0;
  forall_nodes(v, G)
    {
      G[v].adj_diff_in_coords = 0;
      
      double prev = G[v].interpol_coord;
      minimize_vertex_coordinate(v, G);
      G[v].diff_in_coords = fabs(prev - G[v].interpol_coord);
      G[v].interpol_coord = prev;
      VQ[i] = v;
      i++;
    }

  edge e;
  forall_edges(e, G)
    {
      v = G.source(e);
      w = G.target(e);
      G[v].adj_diff_in_coords+=G[w].diff_in_coords;
      G[w].adj_diff_in_coords+=G[v].diff_in_coords;      
    }
  /*
  forall_nodes(v, G)
    {
      G[v].adj_diff_in_coords = 0;
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          G[v].adj_diff_in_coords+=G[w].diff_in_coords;
        }        
    }
  */
 
}

void relaxation_of_gs(TGraphC & G)
  // double_dummy = prev_interpol_coord
{
  edge e;
  node v, w;
  //  TMP_CMP_GRAPHC = &G;  
  //  G.sort_nodes(&cmp_coord_C);                
  coords_by_linear_order(G);
  double prev_la = calc_laC(G);
  forall_nodes(v, G)
    G[v].double_dummy = G[v].interpol_coord;

  cerr << "before gs la = " << prev_la << endl;

  clock_t start =  clock();
  int i;

  //  array<node> VQ(G.number_of_nodes());
  for(i=0; i<Params.gs_relaxation_sweeps + Params.CURRENT_LEVEL*2; i++)
    {
      forall_nodes(v, G)
        {
          minimize_vertex_coordinate(v/*VQ[j]*/, G);
        }
      
      
      double lin_arr = linear_order_by_coord(G);
      
      
      cerr << "cost after gs step " << i << ", arr= " << lin_arr << "\n";

      coords_by_linear_order(G);

      
      node q;
      
      if(prev_la<lin_arr)
        {
          forall_nodes(q, G)
            G[q].interpol_coord = G[q].double_dummy;
          lin_arr = linear_order_by_coord(G);
          cerr << "cost after reverse gs step " <<  ", arr= " << lin_arr << endl; 
          break;
        }
      else
        {
          prev_la = lin_arr;
       
          forall_nodes(q, G)
            G[q].double_dummy = G[q].interpol_coord;
        }
       //      coords_by_linear_order(G);

       
      //      print_diff_save_ArrId_and_ArrId(G);

      
    }

  insert_minimization(G, ALL);
   
  cerr << "Gauss-Seidel relax. time (sec) for " << i << " sweeps : " << ((double)(clock() - start)) / (double)CLOCKS_PER_SEC << endl;
}

double get_best_triple_order(TGraphC & G, array<node> & NodesArr, int i, double lin_arr)
{
  node u = NodesArr[i];
  node v = NodesArr[i+1];
  node w = NodesArr[i+2];

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
  res = all_res[0];

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

  node tmp_start = G.pred_node(u);
  G.set_node_position(min_node, tmp_start);
  G.set_node_position(mid_node, min_node);  
  G.set_node_position(max_node, mid_node);  

  NodesArr[i] = min_node;
  NodesArr[i+1] = mid_node;
  NodesArr[i+2] = max_node;
  
  define_S_values_after_swap(min_node, max_node, G);

  //  cerr << "---------------------------------------------- " << res.cost << endl;
  if(res.cost<0)
    exit(1);
  
  return res.cost;
}

double wmm_minimization(TGraphC & G)
{
  /*
  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_nC);
  define_S_values(G);
  */

  double lin_arr = calc_laC(G);
  int k;
  double old=lin_arr;
  for(k=0; k<50; k++)
    {
      double old_lin_arr = lin_arr;
      lin_arr = one_pass_matching_minimization(G, lin_arr);
      if(fabs(lin_arr - old_lin_arr)<1)
        break;
    }
  //  cerr << "WMM : IMPR=" << old-lin_arr << ", sweeps=" << k << endl;
  return lin_arr;
}

double triples_minimization(TGraphC & G)
{
  node u, v, w;
  
  double lin_arr = calc_laC(G);
  
  double start_la = lin_arr;
  
  
  //  cerr << "1" << endl;
  
  array<node> NodesArr(G.number_of_nodes());
  array<bool> triples_ff(G.number_of_nodes());

  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      triples_ff[i] = false;
      i++;
    }
  //  cerr << "2" << endl;

  //  graph_print(G);
  
  double lin_arr_begin_of_cycle = lin_arr+1;
  int number_of_cycles = 0;

  int ff_count = 0;
  
  while ((fabs(lin_arr_begin_of_cycle-lin_arr)>0.0000001)&&(number_of_cycles < Params.relax_sweeps))
    {
      //      cerr << "TRIPLES RELAX SWEEP START" << endl;
      //      graph_print(G);
      number_of_cycles++;
      lin_arr_begin_of_cycle = lin_arr;
      
      for(i=0; i<G.number_of_nodes()-3; i++)
        {
          //          cerr << i << endl;
          double prev_cost = lin_arr;

          if(triples_ff[i]==false)
            lin_arr = get_best_triple_order(G, NodesArr, i, lin_arr);
          else
            ff_count++;
          //          cerr << lin_arr << endl;
          if(fabs(lin_arr-prev_cost)<=0.00001)
            {
              triples_ff[i] = true;
            }
          else
            {
              int min_pos_ff = mymax(i-2, 0);
              int max_pos_ff = mymin(i+2, G.number_of_nodes()-3);
              for(int j=min_pos_ff; j<max_pos_ff+1; j++)
                triples_ff[j] = false;
            }
        }
    }
  //  cerr << "3" << endl;

  //  lin_arr = calc_laC(G);
  
  //  cerr << "TRIPLES MINIMIZATION FINISH : ff_count " << ff_count << endl;
  
  return lin_arr;
}

double get_best_lcc_triple_order(TGraphC & G, array<node> & NodesArr, int i, double lin_arr)
{
  /*
  node u = NodesArr[i];
  node v = NodesArr[i+1];
  node w = NodesArr[i+2];

  array<triple_res> all_res(6);

  triple_res res;

  res.cost = lin_arr; res.u_pos = G[u].bestLCC_ArrId; res.v_pos = G[v].bestLCC_ArrId; res.w_pos = G[w].bestLCC_ArrId;
  all_res[0] = res;

  double l1 = flip_two_lcc_nodesC(u, v, G, lin_arr);
  res.cost = l1; res.u_pos = G[u].bestLCC_ArrId; res.v_pos = G[v].bestLCC_ArrId; res.w_pos = G[w].bestLCC_ArrId; 
  all_res[1] = res;

  double l2 = flip_two_lcc_nodesC(u, w, G, l1);                                               
  res.cost = l2; res.u_pos = G[u].bestLCC_ArrId; res.v_pos = G[v].bestLCC_ArrId; res.w_pos = G[w].bestLCC_ArrId;  
  all_res[2] = res;                                                                      
  double l3 = flip_two_lcc_nodesC(v, w, G, l2);                                               
  res.cost = l3; res.u_pos = G[u].bestLCC_ArrId; res.v_pos = G[v].bestLCC_ArrId; res.w_pos = G[w].bestLCC_ArrId; 
  all_res[3] = res;
  
  double l4 = flip_two_lcc_nodesC(v, u, G, l3);                                               
  res.cost = l4; res.u_pos = G[u].bestLCC_ArrId; res.v_pos = G[v].bestLCC_ArrId; res.w_pos = G[w].bestLCC_ArrId;
  
  all_res[4] = res;
  
  double l5 = flip_two_lcc_nodesC(w, u, G, l4);                                               
  res.cost = l5; res.u_pos = G[u].bestLCC_ArrId; res.v_pos = G[v].bestLCC_ArrId; res.w_pos = G[w].bestLCC_ArrId; 
  all_res[5] = res;
  
  all_res.sort(&cmp_triple);
  res = all_res[0];

  G[u].bestLCC_ArrId = res.u_pos;
  G[v].bestLCC_ArrId = res.v_pos;
  G[w].bestLCC_ArrId = res.w_pos;

  

  node min_node, mid_node, max_node;
  
  if(G[u].bestLCC_ArrId==i+1)
    min_node = u;
  else if(G[v].bestLCC_ArrId==i+1)
    min_node = v;
  else if(G[w].bestLCC_ArrId==i+1)
    min_node = w;
  
  if(G[u].bestLCC_ArrId==i+2)
    mid_node = u;
  else if(G[v].bestLCC_ArrId==i+2)
    mid_node = v;
  else if(G[w].bestLCC_ArrId==i+2)
    mid_node = w;


  if(G[u].bestLCC_ArrId==i+3)
    max_node = u;
  else if(G[v].bestLCC_ArrId==i+3)
    max_node = v;
  else if(G[w].bestLCC_ArrId==i+3)
    max_node = w;

  node tmp_start = G.pred_node(u);
  G.set_node_position(min_node, tmp_start);
  G.set_node_position(mid_node, min_node);  
  G.set_node_position(max_node, mid_node);  

  NodesArr[i] = min_node;
  NodesArr[i+1] = mid_node;
  NodesArr[i+2] = max_node;
  
  define_S_values_after_swap(min_node, max_node, G);

  if(res.cost<0)
    exit(1);

    
  return res.cost;
  */
}


// should be sorted by bestLCCArrId first
double lcc_triples_minimization(TGraphC & G)
{
  node u, v, w;
  double lin_arr = calc_laC(G);
  double start_la = lin_arr;
  //  cerr << "start tr. min lcc=" << lin_arr << endl;
  array<node> NodesArr(G.number_of_nodes());

  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }

  double lin_arr_begin_of_cycle = lin_arr+1;
  int number_of_cycles = 0;
  
  while ((fabs(lin_arr_begin_of_cycle-lin_arr)>0.0000001)&&(number_of_cycles < Params.relax_sweeps))
    {
      //      cerr << "opsa" << endl;
      number_of_cycles++;
      lin_arr_begin_of_cycle = lin_arr;
      
      for(i=0; i<G.number_of_nodes()-3; i++)
        lin_arr = get_best_lcc_triple_order(G, NodesArr, i, lin_arr);
    }

  return lin_arr;
}

double only_neighbors_relaxation(TGraphC & G)
{
  node u, v;
  double lin_arr = calc_laC(G);
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);
  
  array<node> NodesArr(G.number_of_nodes());

  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }

  double lin_arr_begin_of_cycle = lin_arr+1;
  int number_of_cycles = 0;
  while ((lin_arr_begin_of_cycle>lin_arr)&&(number_of_cycles < Params.relax_sweeps))
    {
      //cerr << "aga" << endl;
      number_of_cycles++;
      lin_arr_begin_of_cycle = lin_arr;
      
      for(i=0; i<G.number_of_nodes()-1; i++)
        {
          u = NodesArr[i];
          v = NodesArr[i+1];
          
          double new_lin_arr = flip_two_nodesC(u, v, G, lin_arr);
          if(new_lin_arr>=lin_arr)
            lin_arr = flip_two_nodesC(u, v, G, new_lin_arr);
          else
            {
              lin_arr = new_lin_arr;
              NodesArr[i] = v;
              NodesArr[i+1] = u;
            }        
        }
      
    }

  lin_arr = calc_laC(G);

  //  cerr << "Number of cycles in neighbors relaxation = " << number_of_cycles << endl;

  return lin_arr;
}

double linear_order_by_coord(TGraphC & G)
{
  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_coord_C);                

  int last_num = 1;
  node v;
  forall_nodes(v, G)
    {
      G[v].ArrId = last_num;
      last_num++;
    }

  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_ArrId);                
  define_S_values(G);
  
  return calc_laC(G);
}

double coords_by_linear_order(TGraphC & G)
{
  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

  node v;
  forall_nodes(v, G)
    {
      G[v].interpol_coord = G[v].S_value;
    }

  TMP_CMP_GRAPHC = &G;  
  G.sort_nodes(&cmp_coord_C);                
  
  return calc_laC(G);
}

//int DIST = 3;

double s_left_dist;
double s_right_dist;

double find_node_best_position(TGraphC & G, array<node> & NodesArr, int p, double lin_arr)
{
  node v = NodesArr[p];
  //  sum_of_bad_changes = 0;
  //  num_of_bad_changes = 0;

  double start_lin_arr = lin_arr;
  double new_lin_arr = lin_arr;
  double min_lin_arr = lin_arr;

  int left_dist = (p<Params.insert_minimization_dist)?p:Params.insert_minimization_dist;
  int right_dist = (G.number_of_nodes()-1-p<Params.insert_minimization_dist)?(G.number_of_nodes()-1-p):Params.insert_minimization_dist;
  
  int left_adj = -1;
  int right_adj = G.number_of_nodes()+1;
  int i,j,k;

  if((Params.individual_minimization == true)&&(Params.CURRENT_LEVEL==0))
    {
      edge e;
      forall_adj_edges(e, v)
        {
          node w = second_adj_for_edge(e, v, G);
          int w_ArrId = G[w].ArrId;
          int v_ArrId = G[v].ArrId;
          if((w_ArrId < v_ArrId)&&(w_ArrId>left_adj))
            left_adj = w_ArrId;
          else if((w_ArrId > v_ArrId)&&(w_ArrId<right_adj))
            right_adj = w_ArrId;
        }
    
      if(left_adj != -1)
        left_dist = G[v].ArrId - left_adj;
      if(right_adj != G.number_of_nodes()+1)
        right_dist = right_adj - G[v].ArrId;
      
      s_left_dist+=left_dist;
      s_right_dist+=right_dist;
    }

  //  cerr << "L=" << left_dist << "; R=" << right_dist << " ----------------\n" << endl;graph_print(G, "");
      
  //  int num_of_flips = get_the_number_of_flips_to_minimize(G, NodesArr, p, lin_arr);

  int saved_min_i = 0;
 
  // Move a node to the leftmost position
  for(i=0; i<=left_dist-1; i++)
    {
      //      cerr << "Flip node : " << G[NodesArr[p-i]].initial_id << " - " << G[NodesArr[p-i-1]].initial_id << endl;
      new_lin_arr = flip_two_nodesC(NodesArr[p-i], NodesArr[p-i-1], G, new_lin_arr);
      node dummy = NodesArr[p-i]; NodesArr[p-i] = NodesArr[p-i-1]; NodesArr[p-i-1] = dummy;
      if(new_lin_arr <= min_lin_arr)
        {
          //          cerr << "\t Got min : " << new_lin_arr << endl;
          min_lin_arr = new_lin_arr;
          saved_min_i = -(i+1);
        }  
    }
  

  //min_lin_arr = new_lin_arr;
  
    
  // Do left+right flips
  
  
  for(i=0; i<left_dist+right_dist; i++)
    {
      if(i>=left_dist)
        {
          //          cerr << "Flip node : " << G[NodesArr[p-left_dist+i]].initial_id << " - " << G[NodesArr[p-left_dist+i+1]].initial_id << endl;
          new_lin_arr = flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G, new_lin_arr);
          node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;
          /*
            if(new_lin_arr > lin_arr)
            {
            sum_of_bad_changes+=fabs(new_lin_arr - lin_arr);
            num_of_bad_changes++;
            }
          */
          if(new_lin_arr <= min_lin_arr)
            {
              //              cerr << "\t Got min : " << new_lin_arr << endl;
              min_lin_arr = new_lin_arr;
              saved_min_i = i+1;
            }
        }
      else
        {
          //          cerr << "dFlip node : " << G[NodesArr[p-left_dist+i]].initial_id << " - " << G[NodesArr[p-left_dist+i+1]].initial_id << endl;
          dummy_flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G);
          node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;
          new_lin_arr = start_lin_arr;
        }
    }

  if(saved_min_i<=0)
    saved_min_i = left_dist + saved_min_i ;
  
  //  cerr << "saved_min_i=" << saved_min_i << endl;
  // Restore minimum
  for(j=0; j<(left_dist+right_dist-saved_min_i); j++)
    {
      //      cerr << "dFlip node : " << G[NodesArr[p+right_dist-j]].initial_id << " - " << G[NodesArr[p+right_dist-1-j]].initial_id << endl;
      dummy_flip_two_nodesC(NodesArr[p+right_dist-j], NodesArr[p+right_dist-1-j], G);      
      
      node dummy = NodesArr[p+right_dist-j]; NodesArr[p+right_dist-j] = NodesArr[p+right_dist-1-j]; NodesArr[p+right_dist-1-j] = dummy;
    }

  //  best_pos = p+right_dist-j;
  return min_lin_arr;
  
}
double find_nodes_averages(TGraphC & G, array<node> & NodesArr, int p, double lin_arr)
{
  //  cerr << p << ", " << lin_arr << endl;
  
  double new_lin_arr = lin_arr;
  double start_lin_arr = lin_arr;

  int i,j,k;
  node v = NodesArr[p];
  //  int left_dist = Params.insert_minimization_dist;
  //  int right_dist = Params.insert_minimization_dist;
  int left_dist = (p<Params.insert_sa_dist)?p:Params.insert_sa_dist;
  int right_dist = (G.number_of_nodes()-1-p<Params.insert_sa_dist)?(G.number_of_nodes()-1-p):Params.insert_sa_dist;

  // Move a node to the leftmost position
  for(i=0; i<=left_dist-1; i++)
    {
      new_lin_arr = flip_two_nodesC(NodesArr[p-i], NodesArr[p-i-1], G, new_lin_arr);
      node dummy = NodesArr[p-i]; NodesArr[p-i] = NodesArr[p-i-1]; NodesArr[p-i-1] = dummy;
    }

  //  cerr << "a" << endl;

  // Do left+right flips
  int saved_min_i = 0;

  //int global_index;
  for(i=0; i<left_dist+right_dist; i++)
    {
      /*
      int global_index;
      if(left_dist<Params.insert_minimization_dist)
        global_index = Params.insert_minimization_dist*2 - left_dist-right_dist+i;
      else
        global_index = i;
      */
      if(new_lin_arr > start_lin_arr)
        {
          //          Insert_SA_Avg_Change[i]+=new_lin_arr-start_lin_arr;
          //          Insert_SA_Counter[i]++;
          EChanges[i].push_back(new_lin_arr-start_lin_arr);
        }

      //      cerr << "flip " << p-left_dist+i << " to " << p-left_dist+i+1 << endl;
      new_lin_arr = flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G, new_lin_arr);
      node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;

    }
  if(new_lin_arr > start_lin_arr)
        {
          //          Insert_SA_Avg_Change[i]+=new_lin_arr-start_lin_arr;
          //          Insert_SA_Counter[i]++;

          EChanges[i].push_back(new_lin_arr-start_lin_arr);
        }
  
  saved_min_i = right_dist;
  //  cerr << saved_min_i << endl;
  // Restore minimum
  for(j=0; j<=(left_dist+right_dist-saved_min_i-1); j++)
    {
      //      cerr << "flip " << p+right_dist-j <<" >>> " << p+right_dist-1-j << endl;
      new_lin_arr = flip_two_nodesC(NodesArr[p+right_dist-j], NodesArr[p+right_dist-1-j], G, new_lin_arr);      
      
      node dummy = NodesArr[p+right_dist-j]; NodesArr[p+right_dist-j] = NodesArr[p+right_dist-1-j]; NodesArr[p+right_dist-1-j] = dummy;
    }
 
  if(fabs(new_lin_arr-start_lin_arr)>0.1)
    {
      cerr << fabs(new_lin_arr-start_lin_arr) << endl;
      cerr << "find nodes averages :" << endl;
      cerr << j << ", " << new_lin_arr << ", " << start_lin_arr << endl;
      exit(1);
    }
  
  return new_lin_arr;
  
}

void check_averages_for_insert_sa(TGraphC & G, double start_lin_arr, array<double> & Insert_SA_Avg_Change, array<double> & Insert_SA_Counter, array<double> & Current_T)
{
  
  node u, v, w;
  
  double lin_arr = start_lin_arr;
  
  double start_la = lin_arr;

  for(int k=0; k<Params.insert_sa_dist*2+1; k++)
    {
      Insert_SA_Avg_Change[k] = 0;
      Insert_SA_Counter[k] = 0;
    }
    
  array<node> NodesArr(G.number_of_nodes());

  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      i++;
    }

  //  cerr << "1" << endl;
  for(i=0; i<G.number_of_nodes(); i++)
    {
      if((i>Params.insert_sa_dist)&&(i<G.number_of_nodes()-Params.insert_sa_dist))
        lin_arr = find_nodes_averages(G, NodesArr, i, lin_arr);
    }
  //  cerr << "2" << endl;
  //  for(int k=0; k<Params.insert_sa_dist*2+1; k++)
  //      Insert_SA_Avg_Change[k] = Insert_SA_Avg_Change[k]/Insert_SA_Counter[k];
  
}

double insert_minimization(TGraphC & G, bool for_fines_only)
{
  int save_dist = Params.insert_minimization_dist;
  // Recalculation of insert minimization dist
  int minimization_dist = (int)((double)Params.insert_minimization_dist+log(E0/(double)G.number_of_edges())/2.0+0.5);
  if(minimization_dist < Params.insert_minimization_dist)
    minimization_dist = Params.insert_minimization_dist;
  //  if(interpolation_order > 100)
  //    interpolation_order = 100;
  if(Params.CURRENT_LEVEL==0)
    minimization_dist =  Params.insert_minimization_dist;
  
  Params.insert_minimization_dist = minimization_dist;
  cerr << "Recalculated minimization distance = " << Params.insert_minimization_dist << endl;
  //////////////////////////////////////////

  
  clock_t start_time = clock();
  
  node u, v, w;
  
  double lin_arr = calc_laC(G);
  
  double start_la = lin_arr;
  
  //  Insert_SA_Avg_Change = 0;
    
  array<node> NodesArr(G.number_of_nodes());
  //  array<node> dummyNodesArr(G.number_of_nodes());

  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      //      dummyNodesArr[i] = v;
      i++;
    }

  double lin_arr_begin_of_cycle = lin_arr+1;
  int number_of_cycles = 0;

  //  double all_sums_of_bad_changes = 0;
  //  int    all_num_of_bad_changes = 0;

  s_left_dist = 0;
  s_right_dist = 0;

  //  int best_pos;
  cerr << "Start min after recalc cost" << endl;
  while ((number_of_cycles==0)||(lin_arr_begin_of_cycle/lin_arr>(1.0001))&&(number_of_cycles < Params.relax_sweeps))
    {
      //      cerr << "DISTANCE RELAX SWEEP START" << endl;
     
      number_of_cycles++;
      lin_arr_begin_of_cycle = lin_arr;

      double nodes_moved = 0;
      double nodes_moved2 = 0;
      for(i=0; i<G.number_of_nodes(); i++)
        {
          double prev_cost = lin_arr;
          double sum_of_bad_changes = 0;
          int num_of_bad_changes = 0;

          v = dummyNodesArr[i];
          int p = G[v].ArrId;
          if((for_fines_only == false)||((for_fines_only == true)&&(G[v].status!=seed)))
            {
              
              lin_arr = find_node_best_position(G, NodesArr, G[v].ArrId-1, lin_arr);
              //              all_sums_of_bad_changes+=sum_of_bad_changes;
              //              all_num_of_bad_changes+=num_of_bad_changes;
            }

          if(G[v].ArrId!=p)
            {
              nodes_moved+=fabs(G[v].ArrId-p);
              nodes_moved2++;
              //              cerr << fabs(best_pos-p) <<endl;
            }
                    
        }
      //      cerr << "INSERT MINIMIZATION : " << nodes_moved2/G.number_of_nodes()*100.0 << " of nodes were moved; total ch =  " << nodes_moved/G.number_of_nodes()*100.0 << " of nodes were moved in this cycle" << endl;
      if(Params.many_cycles_in_min == false)
        break;
      //      cerr << "dist : " << lin_arr << endl;
    }
  //  cerr << "INSERT MINIMIZATION # of cycles : " << number_of_cycles << endl;
  if(Params.individual_minimization==true)
    cerr << "INSERT MINIMIZATION : L_DIST=" << s_left_dist/(double)(G.number_of_nodes()*number_of_cycles) << "; R_DIST=" << s_right_dist/(double)(G.number_of_nodes()*number_of_cycles) << endl;
  
  /*
    double lin_arr = calc_laC(G);
  if(fabs(oolin_arr-lin_arr)>0.001)
    {cerr << "asasasa" << endl;
      exit(1);}
  */
  //  Insert_SA_Avg_Change = all_sums_of_bad_changes / (double)all_num_of_bad_changes;

  cerr << "INSERT MINIMIZATION TIME : " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << endl;
  cerr << "INSERT MINIMIZATION IMPROVEMENT : " << start_la - lin_arr << endl;

  // restore parameter
  Params.insert_minimization_dist = save_dist;

  return lin_arr;  
}

//double qqq,rrr,bbb,ddd;
/*
double old_find_node_best_position_with_sa(TGraphC & G, array<node> & NodesArr, int p, double lin_arr, array<double> & Current_T)
{
  //  graph_print(G, "Before : " );
  double new_lin_arr = lin_arr;
  double min_lin_arr = lin_arr;
  double start_lin_arr = lin_arr;
   
  int i,j,k;
  //  node v = NodesArr[p];
  //  cerr << "p=" << p  << "; |V|=" << G.number_of_nodes() << endl;
  int left_dist = (p<Params.insert_minimization_dist)?p:Params.insert_minimization_dist;
  int right_dist = (G.number_of_nodes()-1-p<Params.insert_minimization_dist)?(G.number_of_nodes()-1-p):Params.insert_minimization_dist;

  //    cerr << "flip node " << G[NodesArr[p]].initial_id << endl;
 array<double> local_probs(left_dist+right_dist+1);  
 array<double> local_costs(left_dist+right_dist+1);  

 double probs_sum = 0;

  double probs_min = 1000;
  // Move a node to the leftmost position
  for(i=0; i<=left_dist-1; i++)
    {
      int global_index = left_dist - i - 1;
            
      new_lin_arr = flip_two_nodesC(NodesArr[p-i], NodesArr[p-i-1], G, new_lin_arr);
      node dummy = NodesArr[p-i]; NodesArr[p-i] = NodesArr[p-i-1]; NodesArr[p-i-1] = dummy;

      //      cerr << "g.i. " << global_index << endl;
      local_costs[global_index] = new_lin_arr;
      local_probs[global_index] = dblmin(1.0,pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index] )));
      
      if((i==0)&&(G[NodesArr[p-i]].was_flipped_with==NodesArr[p-i-1])&&(G[NodesArr[p-i-1]].was_flipped_with==NodesArr[p-i]))
        local_probs[global_index] = 0;        
        
      if(probs_min>1-local_probs[global_index])
        probs_min = 1-local_probs[global_index];
                
      probs_sum+=local_probs[global_index];
    
    }

  //  cerr << "===============================\n";
  //  graph_print(G, "Flips left done : " );
    
  // Do left+right flips
  int saved_min_i = 0;
  min_lin_arr = new_lin_arr;
 
  // array<double> local_costs(left_dist+right_dist);
 
  //   cerr << "L=" << left_dist << "; R="<< right_dist<< endl;
  
  //    cerr << "local probs" << endl;
  //  double lr_count = left_dist+right_dist;

  local_costs[left_dist] = start_lin_arr;
  
  for(i=0; i<left_dist+right_dist+1; i++)
    {
      //      cerr << "step " << i << endl;
      if(i>=left_dist)
        {
          int global_index;
          if(left_dist<Params.insert_minimization_dist)
            global_index = Params.insert_minimization_dist*2 - left_dist-right_dist+i;
          else
            global_index = i;
          
          local_probs[i] = dblmin(1.0,pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index] )));  //  /lr_count;
          
          if((i==left_dist-1)&&(G[NodesArr[p-left_dist+i]].was_flipped_with==NodesArr[p-left_dist+i+1])&&(G[NodesArr[p-left_dist+i+1]].was_flipped_with==NodesArr[p-left_dist+i]))
            {
              local_probs[i] = 0;
              //          cerr << "zanuhali" << endl;
            }
          else if((i==left_dist+1)&&(G[NodesArr[p-left_dist+i]].was_flipped_with==NodesArr[p-left_dist+i-1])&&(G[NodesArr[p-left_dist+i-1]].was_flipped_with==NodesArr[p-left_dist+i]))
            {
              local_probs[i] = 0;
              //cerr << "zanuhali" << endl;
            }
          
          if((i!=left_dist)&&(probs_min>1-local_probs[i]))
            {
              probs_min = 1-local_probs[i];
              
            }
          if(i!=left_dist)
            probs_sum+=local_probs[i];
          
          //      if(i!=left_dist)
          //        probs_sum+=local_probs[i];
          
          //      cerr << "current cost = " << new_lin_arr << endl;
          //                 cerr << "local probs " << local_probs[i] << ", diffdist=" << new_lin_arr-lin_arr << endl;//"; T=" << Current_T[global_index] << "; pow=" <<pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index])) << "; global=" << global_index << endl;
          //      if(p-left_dist+i+1 != G.number_of_nodes())
          //        cerr <<  "do flip " << G[NodesArr[p-left_dist+i]].initial_id << " and " << G[NodesArr[p-left_dist+i+1]].initial_id << endl ;
          
          if(i<left_dist+right_dist)
            {
              new_lin_arr = flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G, new_lin_arr);
              node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;
              local_costs[i+1] = new_lin_arr;
              
              //                      cerr << "new " << new_lin_arr << endl;
              
              //                if(new_lin_arr <= min_lin_arr)
              //                {
              //                min_lin_arr = new_lin_arr;
              //                saved_min_i = i+1;
              //                //              cerr << "got min " << saved_min_i << "; with cost=" << min_lin_arr << endl;
              //                }
              
            }
        }
      else
        {
          dummy_flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G);
          node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;
          new_lin_arr = start_lin_arr;
        }
    }
  
  if(probs_min==1000) probs_min = 0;
  local_probs[left_dist] = probs_min;
  probs_sum += probs_min;

  //  cerr << "probs sum=" << probs_sum << endl;
  for(i=0; i<local_probs.size(); i++)
    {
      local_probs[i] = local_probs[i]/probs_sum;
    }
  //  cerr << "Prob to stay " << local_probs[left_dist] << endl;
  //  graph_print(G, "Flips right done : " );
  // Restore minimum
  
  //  cerr << "bilo : " << lin_arr << "; " << "stalo : " << min_lin_arr << "; " << min_lin_arr - lin_arr << "; " << saved_min_i << "; left_dist=" << left_dist << endl;
  //  bool xxx = false;
  //  ddd++;
  //  if(saved_min_i!=left_dist)     
  //    {
      //      cerr << "ogo" << endl;
      //            cerr << saved_min_i << endl;
      //      bbb++;
      //      xxx=true;
      //    }
    //  else
  // {
  //        for(j=0; j<local_probs.size(); j++)
  //          cerr << local_probs[j] << "; cost=" << local_costs[j] << endl;
                 //      cerr << "nachali R; " << local_probs.size() << endl;
      random_variate R(local_probs);
      //        cerr << "jopsa" << endl;
      saved_min_i = R.generate();

      //      if(saved_min_i==left_dist)
      //        {
      //          qqq++;//=local_probs[left_dist];
      //        }
      //      rrr++;

        
      //      cerr << "generated index " << saved_min_i << "\n ";

      //                     cerr << "final pos = " << saved_min_i << endl;
      //    }

  new_lin_arr = local_costs[local_costs.size()-1];
  for(j=0; j<(left_dist+right_dist-saved_min_i); j++)
    {
      //      cerr << j <<", " << left_dist+right_dist-saved_min_i <<endl;
      //      cerr << "do return flip " << G[NodesArr[p+right_dist-j]].initial_id << " to " << G[NodesArr[p+right_dist-1-j]].initial_id << endl;
      dummy_flip_two_nodesC(NodesArr[p+right_dist-j], NodesArr[p+right_dist-1-j], G);      

      new_lin_arr = local_costs[local_costs.size()-j-2];
      node dummy = NodesArr[p+right_dist-j]; NodesArr[p+right_dist-j] = NodesArr[p+right_dist-1-j]; NodesArr[p+right_dist-1-j] = dummy;
    }  
if(saved_min_i-right_dist==1)
    {
      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j+1];
      G[NodesArr[p+right_dist-j+1]].was_flipped_with = NodesArr[p+right_dist-j];
    }
  else if(saved_min_i-right_dist==-1)
    {
      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j-1];
      G[NodesArr[p+right_dist-j-1]].was_flipped_with = NodesArr[p+right_dist-j];
    }
  //  cerr << "start : " << p << endl;
  //    cerr << "finished with : " << new_lin_arr << endl;
 
  //  graph_print(G, "After : ");
  //  exit(1);
  return new_lin_arr;
}
*/
double find_node_best_position_with_sa(TGraphC & G, array<node> & NodesArr, int p, double lin_arr, array<double> & Current_T)
{
  double new_lin_arr = lin_arr;
  double min_lin_arr = lin_arr;
  double start_lin_arr = lin_arr;
  //   graph_print(G, "Before : " );
  int i,j,k;
  int left_dist = (p<Params.insert_sa_dist)?p:Params.insert_sa_dist;
  int right_dist = (G.number_of_nodes()-1-p<Params.insert_sa_dist)?(G.number_of_nodes()-1-p):Params.insert_sa_dist;

  array<double> local_probs(left_dist+right_dist+1);  
  array<double> local_costs(left_dist+right_dist+1);  

  double probs_sum = 0;
  double probs_min = 1000;
  //       cerr << "flip node " << G[NodesArr[p]].initial_id << endl;

  //       cerr << "L=" <<left_dist << "; R=" << right_dist << endl;
  // Move a node to the leftmost position
  for(i=0; i<=left_dist-1; i++)
    {
      int global_index = Params.insert_sa_dist-i-1;//left_dist - i - 1;
            
      new_lin_arr = flip_two_nodesC(NodesArr[p-i], NodesArr[p-i-1], G, new_lin_arr);
      node dummy = NodesArr[p-i]; NodesArr[p-i] = NodesArr[p-i-1]; NodesArr[p-i-1] = dummy;

      local_costs[left_dist - i - 1] = new_lin_arr;
      local_probs[left_dist - i - 1] = dblmin(1.0,pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index] )));
      //      cerr << "gi = " << global_index << "; T=" << Current_T[global_index] << "; lpr=" << local_probs[left_dist - i - 1] << endl;
      
      if((i==0)&&(G[NodesArr[p-i]].was_flipped_with==NodesArr[p-i-1])&&(G[NodesArr[p-i-1]].was_flipped_with==NodesArr[p-i]))
        local_probs[left_dist - i - 1] = 0;        
        
      if(probs_min>1-local_probs[left_dist - i - 1])
        probs_min = 1-local_probs[left_dist - i - 1];
                
      probs_sum+=local_probs[left_dist - i - 1];
    
    }
  //graph_print(G, "Flips left done : " );
  // Do left+right flips
  int saved_min_i = 0;
  min_lin_arr = new_lin_arr;
 
  local_costs[left_dist] = start_lin_arr;
  
  for(i=0; i<left_dist+right_dist+1; i++)
    {
      if(i>=left_dist)
        {
          int global_index;
          if(left_dist<Params.insert_sa_dist)
            global_index = Params.insert_sa_dist*2 - left_dist-right_dist+i;
          else
            global_index = i;
          
          local_probs[i] = dblmin(1.0,pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index] ))); 
          
          if((i==left_dist-1)&&(G[NodesArr[p-left_dist+i]].was_flipped_with==NodesArr[p-left_dist+i+1])&&(G[NodesArr[p-left_dist+i+1]].was_flipped_with==NodesArr[p-left_dist+i]))
            local_probs[i] = 0;
          else if((i==left_dist+1)&&(G[NodesArr[p-left_dist+i]].was_flipped_with==NodesArr[p-left_dist+i-1])&&(G[NodesArr[p-left_dist+i-1]].was_flipped_with==NodesArr[p-left_dist+i]))
            local_probs[i] = 0;
          
          if((i!=left_dist)&&(probs_min>1-local_probs[i]))
            {
              probs_min = 1-local_probs[i];
              
            }
          if(i!=left_dist)
            probs_sum+=local_probs[i];
          
          
          if(i<left_dist+right_dist)
            {
              new_lin_arr = flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G, new_lin_arr);
              node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;
              local_costs[i+1] = new_lin_arr;              
            }
        }
      else
        {
          dummy_flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G);
          node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;
          new_lin_arr = start_lin_arr;
        }
    }
  
  if(probs_min==1000) probs_min = 0;
  local_probs[left_dist] = probs_min;
  probs_sum += probs_min;
  //  cerr << "probs_sum=" << probs_sum << endl;
  for(i=0; i<local_probs.size(); i++)
    local_probs[i] = local_probs[i]/probs_sum;

  //  for(j=0; j<local_probs.size(); j++)
  //              cerr << local_probs[j] << endl;
      
  random_variate R(local_probs);
  saved_min_i = R.generate();
  // cerr << saved_min_i << endl;
  new_lin_arr = local_costs[local_costs.size()-1];
  for(j=0; j<(left_dist+right_dist-saved_min_i); j++)
    {
      dummy_flip_two_nodesC(NodesArr[p+right_dist-j], NodesArr[p+right_dist-1-j], G);
      new_lin_arr = local_costs[local_costs.size()-j-2];
      node dummy = NodesArr[p+right_dist-j]; NodesArr[p+right_dist-j] = NodesArr[p+right_dist-1-j]; NodesArr[p+right_dist-1-j] = dummy;
    }  
  
  if(saved_min_i-right_dist==1)
    {
      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j+1];
      G[NodesArr[p+right_dist-j+1]].was_flipped_with = NodesArr[p+right_dist-j];
    }
  else if(saved_min_i-right_dist==-1)
    {
      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j-1];
      G[NodesArr[p+right_dist-j-1]].was_flipped_with = NodesArr[p+right_dist-j];
    }
  //graph_print(G, "After : " );
// if((left_dist==5)&&(right_dist==5))
// exit(1);
  return new_lin_arr;
}
/*
double very_old_find_node_best_position_with_sa(TGraphC & G, array<node> & NodesArr, int p, double lin_arr, array<double> & Current_T)
{
    graph_print(G, "Before : " );
  double new_lin_arr = lin_arr;
  double min_lin_arr = lin_arr;
  
  int i,j,k;
  node v = NodesArr[p];
  //  cerr << "p=" << p  << "; |V|=" << G.number_of_nodes() << endl;
  int left_dist = (p<Params.insert_sa_dist)?p:Params.insert_sa_dist;
  int right_dist = (G.number_of_nodes()-1-p<Params.insert_sa_dist)?(G.number_of_nodes()-1-p):Params.insert_sa_dist;

  //      cerr << "flip node " << G[NodesArr[p]].initial_id << endl;

  // Move a node to the leftmost position
  for(i=0; i<=left_dist-1; i++)
    {
      new_lin_arr = flip_two_nodesC(NodesArr[p-i], NodesArr[p-i-1], G, new_lin_arr);
      node dummy = NodesArr[p-i]; NodesArr[p-i] = NodesArr[p-i-1]; NodesArr[p-i-1] = dummy;
    }

  //  cerr << "===============================\n";
  //    graph_print(G, "Flips left done : " );
    
  // Do left+right flips
  int saved_min_i = 0;
  min_lin_arr = new_lin_arr;
 
  // array<double> local_costs(left_dist+right_dist);
  array<double> local_probs(left_dist+right_dist+1);  

  //   cerr << "L=" << left_dist << "; R="<< right_dist<< endl;
  
  //    cerr << "local probs" << endl;
  //  double lr_count = left_dist+right_dist;
  double probs_sum = 0;

  double probs_min = 1000;
  for(i=0; i<left_dist+right_dist+1; i++)
    {
      int global_index;
      if(left_dist<Params.insert_sa_dist)
        global_index = Params.insert_sa_dist*2 - left_dist-right_dist+i;
      else
        global_index = i;
      
      local_probs[i] = dblmin(1.0,pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index] )));  //  /lr_count;

      if((i==left_dist-1)&&(G[NodesArr[p-left_dist+i]].was_flipped_with==NodesArr[p-left_dist+i+1])&&(G[NodesArr[p-left_dist+i+1]].was_flipped_with==NodesArr[p-left_dist+i]))
        {
          local_probs[i] = 0;
          //          cerr << "zanuhali" << endl;
        }
      else if((i==left_dist+1)&&(G[NodesArr[p-left_dist+i]].was_flipped_with==NodesArr[p-left_dist+i-1])&&(G[NodesArr[p-left_dist+i-1]].was_flipped_with==NodesArr[p-left_dist+i]))
        {
          local_probs[i] = 0;
          //cerr << "zanuhali" << endl;
        }
      
      if((i!=left_dist)&&(probs_min>1-local_probs[i]))
        {
          probs_min = 1-local_probs[i];
          
        }
       if(i!=left_dist)
         probs_sum+=local_probs[i];
         
        //      if(i!=left_dist)
        //        probs_sum+=local_probs[i];

      //      cerr << "current cost = " << new_lin_arr << endl;
       //       cerr << "local probs " << local_probs[i] << ", diffdist=" << new_lin_arr-lin_arr << endl;//"; T=" << Current_T[global_index] << "; pow=" <<pow(2.718, -(new_lin_arr-lin_arr)/(Current_T[global_index])) << "; global=" << global_index << endl;
      //      if(p-left_dist+i+1 != G.number_of_nodes())
      //        cerr <<  "do flip " << G[NodesArr[p-left_dist+i]].initial_id << " and " << G[NodesArr[p-left_dist+i+1]].initial_id << endl ;

      if(i<left_dist+right_dist)
        {
          new_lin_arr = flip_two_nodesC(NodesArr[p-left_dist+i], NodesArr[p-left_dist+i+1], G, new_lin_arr);
          node dummy = NodesArr[p-left_dist+i]; NodesArr[p-left_dist+i] = NodesArr[p-left_dist+i+1]; NodesArr[p-left_dist+i+1] = dummy;

          //                    cerr << "new " << new_lin_arr << ", min " << min_lin_arr << endl;
          
          //          if(new_lin_arr <= min_lin_arr)
          //            {
          //              min_lin_arr = new_lin_arr;
          //              saved_min_i = i+1;
              //              cerr << "got min " << saved_min_i << "; with cost=" << min_lin_arr << endl;
              //            }
         
        }
    }
  if(probs_min==1000) probs_min = 0;
  local_probs[left_dist] = probs_min;
  probs_sum += probs_min;

  //  cerr << "probs sum=" << probs_sum << endl;
  for(i=0; i<local_probs.size(); i++)
    {
      local_probs[i] = local_probs[i]/probs_sum;
    }
  //  cerr << "Prob to stay " << local_probs[left_dist] << endl;
  //  graph_print(G, "Flips right done : " );
  // Restore minimum
  
  //  cerr << "bilo : " << lin_arr << "; " << "stalo : " << min_lin_arr << "; " << min_lin_arr - lin_arr << "; " << saved_min_i << "; left_dist=" << left_dist << endl;
  //  bool xxx = false;
  //  ddd++;
  //  if(saved_min_i!=left_dist)     
  //    {
      //      cerr << "ogo" << endl;
      //            cerr << saved_min_i << endl;
      //      bbb++;
      //      xxx=true;
      //    }
    //  else
  //{
  //            for(j=0; j<local_probs.size(); j++)
  //              cerr << local_probs[j] << endl;
                 //      cerr << "nachali R; " << local_probs.size() << endl;
      random_variate R(local_probs);
      //        cerr << "jopsa" << endl;
      saved_min_i = R.generate();

  
        
      //      cerr << "generated index " << saved_min_i << "\n ";

      
      //      if(gen_prob()<local_probs[index])
      //         saved_min_i = index;
      //      else
      //        saved_min_i = right_dist;
      
      
      //      array<int> tmp_probs(2);
      //          tmp_probs[0] = MY_MAX_INT - local_probs[index]; tmp_probs[1] = local_probs[index];
      //        
      //          random_variate R1(tmp_probs);
      //          int index2 = R1.generate();
      //          if(index2==0)
      //            saved_min_i = right_dist;
      //          else
      //            saved_min_i = index;
      
      //                       cerr << "final pos = " << saved_min_i << endl;
      //    }
  
  for(j=0; j<(left_dist+right_dist-saved_min_i); j++)
    {
      //      cerr << j <<", " << left_dist+right_dist-saved_min_i <<endl;
      //      cerr << "do return flip " << G[NodesArr[p+right_dist-j]].initial_id << " to " << G[NodesArr[p+right_dist-1-j]].initial_id << endl;
      new_lin_arr = flip_two_nodesC(NodesArr[p+right_dist-j], NodesArr[p+right_dist-1-j], G, new_lin_arr);      
      
      node dummy = NodesArr[p+right_dist-j]; NodesArr[p+right_dist-j] = NodesArr[p+right_dist-1-j]; NodesArr[p+right_dist-1-j] = dummy;
    }
   
    //  if(abs(saved_min_i-right_dist)==1)
    //    {
    //      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j+1];
    //      G[NodesArr[p+right_dist-j+1]].was_flipped_with = NodesArr[p+right_dist-j];
    //    }
  
 if(saved_min_i-right_dist==1)
    {
      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j+1];
      G[NodesArr[p+right_dist-j+1]].was_flipped_with = NodesArr[p+right_dist-j];
    }
  else if(saved_min_i-right_dist==-1)
    {
      G[NodesArr[p+right_dist-j]].was_flipped_with = NodesArr[p+right_dist-j-1];
      G[NodesArr[p+right_dist-j-1]].was_flipped_with = NodesArr[p+right_dist-j];
    }

  //  cerr << "start : " << p << endl;
  //  cerr << "saved : " << saved_min_i << endl;
 
  //   graph_print(G, "After : ");
   //   if((left_dist==5)&&(right_dist==5))
   //   exit(1);
  return new_lin_arr;
}
*/

double insert_sa_one_cycle(TGraphC & G, double & lin_arr, array<double> & Current_T)
{
  
  node u, v, w;
  
  array<node> NodesArr(G.number_of_nodes());
  //  array<node> dummyNodesArr(G.number_of_nodes());
  
  int i=0;
  forall_nodes(v, G)
    {
      NodesArr[i] = v;
      //      dummyNodesArr[i] = v;
      i++;
    }

  //  for(i=Params.insert_minimization_dist; i<G.number_of_nodes()-Params.insert_minimization_dist; i++)
  forall_nodes(v, G)
    {
      G[v].was_flipped_with = nil;
      //      G[v].touched = 0; // to remove this checker
    }

  for(i=0; i<dummyNodesArr.size(); i++)
    //  forall_nodes(v, G)
    {
      v = dummyNodesArr[i];
      //      G[v].touched++;  // to remove this checker
      lin_arr = /*very_old_*/find_node_best_position_with_sa(G, NodesArr, G[v].ArrId-1, lin_arr, Current_T);
    }
  /////////to remove this checker/////////////
  /*
  double x=0;
   forall_nodes(v, G)
    {
      if(G[v].touched!=1)
        {cerr << "jopsa" << endl; exit(1);}
      x+=G[v].touched;
       
    }
  */
   //cerr << x / (double)G.number_of_nodes() << endl;
  //  cerr << "===================================================================" << endl;
   /////////////////////
  
  /*
   for(i=0; i<G.number_of_nodes(); i++)
    {
      //if(G.number_of_nodes()>Params.insert_minimization_dist*2+1)
      lin_arr = find_node_best_position_with_sa(G, NodesArr, i, lin_arr, Current_T);
      //      cerr << "Finished with LA=" << lin_arr << endl;
    }
  */
  return lin_arr;  
}
double calc_T_insert_sa(double avg_change, double real_part_of_bad_changes)
{
   return (-1.0)* (avg_change) / log(real_part_of_bad_changes/100.0); 
}

void calc_vector_T_insert_sa(double real_part_of_bad_changes, array<double> & Insert_SA_Avg_Change, array<double> & Current_T)
{
  for(int i=0; i<Params.insert_sa_dist*2+1; i++)
    {
      double perc_to_permit =  real_part_of_bad_changes;//*pow(1.1,fabs(i-Params.insert_minimization_dist));
      
      Current_T[i] = /*pow(0.9,fabs(i-Params.insert_minimization_dist)) * */(-1.0)* ( Insert_SA_Avg_Change[i]) / log(perc_to_permit/100.0);
      //      cerr << "dist=" << fabs(i-Params.insert_minimization_dist) << "; %=" << perc_to_permit << ";\t T=" << Current_T[i] << ";\t AVG=" << Insert_SA_Avg_Change[i] << endl;
    }
  Current_T[Params.insert_sa_dist] = 0;
}

void get_permit_perc_insert_sa(double & permit_perc)
{
      permit_perc = Params.sa_permition_percent;
      //      cerr << "SA Perm % " << permit_perc << endl;
}

void ReloadT(array<double> & Current_T)
{
  //  cerr << "Reload T" << endl;
  for(int i=0; i<EChanges.size(); i++)
    {
      EChanges[i].sort();
      //      x.sort();
      int ptr = (int)((double)EChanges[i].size() * (double)Params.sa_permition_percent / 100.0);
      //     cerr << "EChanges[i].size=" << EChanges[i].size() << "; ptr calculated " << ptr << endl;
      double x;
      list_item it = EChanges[i].first();
      for(int j=0; j<ptr; j++)
        {       x = EChanges[i][it]; it = EChanges[i].succ(it); }
      Current_T[i] = (-1.0)* x / log(95.0/100.0);
    }
  Current_T[Params.insert_sa_dist] = 0;
  //  exit(0);
}

int recalc_HC_steps(TGraphC & G)
{
  int steps = (int)((double)Params.hot_cold_steps*log(E0/(double)G.number_of_edges())/2.0+0.5);
  if(steps < Params.hot_cold_steps)
    steps = Params.hot_cold_steps;
  if(steps>100)
    steps = 100;
  return steps;
}

void reverse_graph_order(TGraphC & G)
{
  node v;
  forall_nodes(v, G)
    G[v].ArrId = G.number_of_nodes()+1-G[v].ArrId;

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

}

double old_good_insert_sa(TGraphC & G, double arr_cost)
{

  int save_dist = Params.insert_sa_dist;
  // Recalculation of insert sa dist
  int sa_dist = (int)((double)Params.insert_sa_dist+log(E0/(double)G.number_of_edges())/2.0+0.5);
  if(sa_dist < Params.insert_sa_dist)
    sa_dist = Params.insert_sa_dist;
  if(Params.CURRENT_LEVEL==0)
    sa_dist =  Params.insert_sa_dist;
  
  Params.insert_sa_dist = sa_dist;
  cerr << "Recalculated sa distance = " << Params.insert_sa_dist << endl;
  //////////////////////////////////////////

  
  int save_insert_sa_dist = Params.insert_sa_dist;
  while(G.number_of_nodes()-1<Params.insert_sa_dist*2+1)
    Params.insert_sa_dist--;

  EChanges.resize(Params.insert_sa_dist*2+1); 

  cerr << "Start SA on distance +/- " << Params.insert_sa_dist << endl;
  array<double> Insert_SA_Avg_Change(Params.insert_sa_dist*2+1);
  array<double> Insert_SA_Counter(Params.insert_sa_dist*2+1);
  array<double> Current_T(Params.insert_sa_dist*2+1);

  //  graph_print(G, "aga : ");
  double lcc_cost = 0;
  node u,v;
  //  cerr << "start" << endl;
  double lin_arr = arr_cost;
  cerr << "SA_DIST: Before = " << lin_arr << endl;

  double min_lin_arr = lin_arr;
  //  save_arrangement(G);  

  if(Params.use_lcc==true)
    lcc_init(G);
 
  bool last_minla_changed = false;

  double permit_perc;
  //  cerr << "before get_permitP" << endl;
  get_permit_perc_insert_sa(permit_perc);

  int HC_STEPS = recalc_HC_steps(G);//Params.hot_cold_steps;
  cerr << "# of HC STEPS = " << HC_STEPS << endl;
  
  int number_of_sweeps = 0;

  double avg_min=0;
  double avg_lcc=0;
  double avg_sweep=0;
  
  for(int T_stage = 0; T_stage<HC_STEPS; T_stage++)
    {
      update_time();
      for(int k=0; k<EChanges.size(); k++)
        EChanges[k].clear();
      check_averages_for_insert_sa(G, lin_arr, Insert_SA_Avg_Change, Insert_SA_Counter, Current_T);
      ReloadT(Current_T);

      //      calc_vector_T_insert_sa(permit_perc, Insert_SA_Avg_Change, Current_T);
      /*           
      cerr << "Recalc T : " << endl;
      for(int i=0; i<Params.insert_sa_dist*2+1; i++)
        {
          cerr << Current_T[i] << endl;
        }
      */
      double sweep = time(0);     

      //      qqq=0; rrr=0; bbb=0; ddd=0;
      
      for(int p=0; p<Params.number_of_sweeps_in_one_hc; p++)
        {
          lin_arr = insert_sa_one_cycle(G, lin_arr, Current_T);
          for(int i=0; i<Params.insert_sa_dist*2+1; i++)
            Current_T[i] = Current_T[i] * Params.sa_alpha;          
        }
      avg_sweep+=(double)time(0)-sweep;
      //      cerr << "#vert. remained = " << qqq/rrr << endl;
      //      cerr << "impr = " << bbb/ddd << endl;

      cerr << "ORDER COST BEFORE MIN. = " << lin_arr << "\n";
      double min=time(0);
      lin_arr = insert_minimization(G, ALL);
      avg_min+=(double)time(0)-min;
      
      double lcc=time(0);
      if(Params.use_lcc==true)
        {
          lcc_update(G, lin_arr, lcc_cost);
          //          reverse_graph_order(G);
          //          lcc_update(G, lin_arr, lcc_cost);
          //          reverse_graph_order(G);
          cerr << "ORDER COST = " << lin_arr << "\t LCC COST = " << lcc_cost << endl;
          cerr << "------------------\n";
        }
      avg_lcc+=(double)time(0)-lcc;
    }
  
  cerr << "Average for minimization in one HC step = " << avg_min/(double)HC_STEPS << endl;
  cerr << "Average for lcc in one HC step = " << avg_lcc/(double)HC_STEPS << endl;
  cerr << "Average for SA = " << avg_sweep/(double)HC_STEPS<< endl;
    
  if((Params.use_lcc==true)&&(lcc_cost<lin_arr))
    {
      restore_best_lcc_arrangement(G);
      lin_arr = calc_laC(G);      
    }
  
  cerr << "SA_INSERT_DIST: After = " << lin_arr << endl;

  Params.insert_sa_dist = save_insert_sa_dist;


  // restore sa prams
  Params.insert_sa_dist = save_dist;

  return lin_arr;
}

double insert_sa(TGraphC & G, double arr_cost)
{

  int save_dist = Params.insert_sa_dist;
  // Recalculation of insert sa dist
  int sa_dist = (int)((double)Params.insert_sa_dist+log(E0/(double)G.number_of_edges())/2.0+0.5);
  if(sa_dist < Params.insert_sa_dist)
    sa_dist = Params.insert_sa_dist;
  if(Params.CURRENT_LEVEL==0)
    sa_dist =  Params.insert_sa_dist;
  
  Params.insert_sa_dist = sa_dist;
  cerr << "Recalculated sa distance = " << Params.insert_sa_dist << endl;
  //////////////////////////////////////////

  
  int save_insert_sa_dist = Params.insert_sa_dist;
  while(G.number_of_nodes()-1<Params.insert_sa_dist*2+1)
    Params.insert_sa_dist--;

  EChanges.resize(Params.insert_sa_dist*2+1); 

  cerr << "Start SA on distance +/- " << Params.insert_sa_dist << endl;
  array<double> Insert_SA_Avg_Change(Params.insert_sa_dist*2+1);
  array<double> Insert_SA_Counter(Params.insert_sa_dist*2+1);
  array<double> Current_T(Params.insert_sa_dist*2+1);

  //  graph_print(G, "aga : ");
  double lcc_cost = 0;
  node u,v;
  //  cerr << "start" << endl;
  double lin_arr = arr_cost;
  cerr << "SA_DIST: Before = " << lin_arr << endl;

  double min_lin_arr = lin_arr;
  //  save_arrangement(G);  

  if(Params.use_lcc==true)
    lcc_init(G);
  else
    save_arrangement(G);
 
  bool last_minla_changed = false;

  double permit_perc;
  //  cerr << "before get_permitP" << endl;
  get_permit_perc_insert_sa(permit_perc);

  int HC_STEPS = recalc_HC_steps(G);//Params.hot_cold_steps;
  cerr << "# of HC STEPS = " << HC_STEPS << endl;
  
  int number_of_sweeps = 0;

  double avg_min=0;
  double avg_lcc=0;
  double avg_sweep=0;

  double min_sa_cost = lin_arr;
  //  node_array<int> min_sa_order(G);
  
  
  for(int T_stage = 0; T_stage<HC_STEPS; T_stage++)
    {
      update_time();
      for(int k=0; k<EChanges.size(); k++)
        EChanges[k].clear();
      check_averages_for_insert_sa(G, lin_arr, Insert_SA_Avg_Change, Insert_SA_Counter, Current_T);
      ReloadT(Current_T);

      //      calc_vector_T_insert_sa(permit_perc, Insert_SA_Avg_Change, Current_T);
      /*           
      cerr << "Recalc T : " << endl;
      for(int i=0; i<Params.insert_sa_dist*2+1; i++)
        {
          cerr << Current_T[i] << endl;
        }
      */
      double sweep = time(0);     

      //      qqq=0; rrr=0; bbb=0; ddd=0;
      
      for(int p=0; p<Params.number_of_sweeps_in_one_hc; p++)
        {
          lin_arr = insert_sa_one_cycle(G, lin_arr, Current_T);
          for(int i=0; i<Params.insert_sa_dist*2+1; i++)
            Current_T[i] = Current_T[i] * Params.sa_alpha;          
        }
      avg_sweep+=(double)time(0)-sweep;
      //      cerr << "#vert. remained = " << qqq/rrr << endl;
      //      cerr << "impr = " << bbb/ddd << endl;

      cerr << "ORDER COST BEFORE MIN. = " << lin_arr << "\n";
      double min=time(0);
      lin_arr = insert_minimization(G, ALL);
      avg_min+=(double)time(0)-min;
      
      double lcc=time(0);
      if(Params.use_lcc==true)
        {
          lcc_update(G, lin_arr, lcc_cost);
          //          reverse_graph_order(G);
          //          lcc_update(G, lin_arr, lcc_cost);
          //          reverse_graph_order(G);
          if(lcc_cost>lin_arr)
            {
              lcc_init(G);
              lcc_cost = lin_arr;
            }
          cerr << "ORDER COST = " << lin_arr << "\t LCC COST = " << lcc_cost << endl;
          cerr << "------------------\n";
        }
      
      if((lin_arr < min_sa_cost)&&(Params.use_lcc==false))
        {
          save_arrangement(G);
          min_sa_cost = lin_arr;
        }
      avg_lcc+=(double)time(0)-lcc;
    }
  
  //  cerr << "Average for minimization in one HC step = " << avg_min/(double)HC_STEPS << endl;
  //  cerr << "Average for lcc in one HC step = " << avg_lcc/(double)HC_STEPS << endl;
  //  cerr << "Average for SA = " << avg_sweep/(double)HC_STEPS<< endl;
    
  if((Params.use_lcc==true)&&(lcc_cost<lin_arr))
    {
      restore_best_lcc_arrangement(G);
      lin_arr = calc_laC(G);      
    }

  if((Params.use_lcc==false)&&(min_sa_cost<lin_arr))
    {
      restore_arrangement(G);
      lin_arr = min_sa_cost;
    }
  
  cerr << "SA_INSERT_DIST: After = " << lin_arr << endl;

  Params.insert_sa_dist = save_insert_sa_dist;


  // restore sa prams
  Params.insert_sa_dist = save_dist;

  return lin_arr;
}
