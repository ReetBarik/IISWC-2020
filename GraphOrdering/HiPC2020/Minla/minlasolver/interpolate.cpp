#include "cmpfuncs.h"
int calc_window_start_size(TGraphC & G)
{
  double N = G.number_of_nodes();
  int s;

  if(N<120)
    {
      s = (int)((double)N/8.0);
      if(s==1)
        s++;
     
    }
  else
    s = 5;

  if(s>=N*9.0/20.0)
    s = (int)(s/2.0);
  
  return s;
  
  /*
  if(N < 1000)
    return (int)(10.0*N/1000.0)+5;
  else if(N<10000)
    return (int)(10.0*N/10000.0)+10;
  else if(N<100000)
    return (int)(30.0*N/100000.0)+20;
  else if(N<1000000)
    return (int)(150.0*N/1000000.0)+50;
  else
    return 200;
  */
}

void update_real_neighbors(node v, node left, node right, TGraphC & G)
{
  return;
  node l = left;
  node r = right;

  //  cerr << l << ", " << v << ", " << r << endl;
  
  if(l==r)
    {
      //      cerr << "aga : " << G[left].interpol_coord << " - " << G[v].interpol_coord <<  " - " << G[right].interpol_coord << endl;
      if(G[v].interpol_coord < G[l].interpol_coord)
        {
          if(G[l].left_real_neighbor==nil)
            {
              G[l].left_real_neighbor = v;
              G[v].right_real_neighbor = l;
              return;
            }        
          l = G[l].left_real_neighbor;          
        }
      else
        {
          if(G[r].right_real_neighbor==nil)
            {
              G[r].right_real_neighbor = v;
              G[v].left_real_neighbor = r;
              return;
            }
          //          cerr << "a" << endl;cerr << G[r].right_real_neighbor << endl;
          r = G[r].right_real_neighbor;
          //          cerr << G[r].right_real_neighbor << endl;
        }
    }
  
  if((l==nil)||(r==nil))
    {
      cerr << "Error : two nil real neighbors" << endl;
      exit(1);
    }
  else
    {
      //      cerr << G[l].interpol_coord << " - " << G[v].interpol_coord <<  " - " << G[r].interpol_coord << endl;
    }

  G[v].left_real_neighbor = l;
  G[v].right_real_neighbor = r;
  if(l!=nil)
    G[l].right_real_neighbor = v;
  if(r!=nil)
    G[r].left_real_neighbor = v;
  //  cerr << l << ", " << v << ", " << r << endl;
  //  cerr << "v konze : " << G[l].interpol_coord << " - " << G[v].interpol_coord <<  " - " << G[r].interpol_coord << endl;
  
}

void get_nodes_to_interpolate(TGraphC & G, double thresh, int number_of_fines)
{
  node v, w;
  edge e;
  int were_interpolated = 0;
  
  while(were_interpolated!=number_of_fines)
    {
      //cerr << "were interpolated " << were_interpolated << ", thresh=" << thresh << endl;
      forall_nodes(v, G)
        {
          if(G[v].interpol_coord == -1)
            {
              double count = 0;
              double interpol_coord = 0;
              double weights_sum = 0;
              list<two_tuple<edge,node> > l;
              forall_adj_edges(e, v)
                {
                  w = second_adj_for_edge(e, v, G);                  
                  if(G[w].interpol_coord != -1)
                    {
                      count++;
                      two_tuple<edge, node> two(e, w);
                      l.push_back(two);
                    }
                }
              if(count/(double)G.degree(v)>thresh)
                {
                  
                  if(l.length()==1)
                    {
                      if(G[l[l.first_item()].second()].right_real_neighbor==nil)
                        G[v].interpol_coord = G[l[l.first_item()].second()].interpol_coord+0.01;
                      else
                        G[v].interpol_coord = G[l[l.first_item()].second()].interpol_coord+(G[G[l[l.first_item()].second()].right_real_neighbor].interpol_coord - G[l[l.first_item()].second()].interpol_coord)/100.0;
                      update_real_neighbors(v,l[l.first_item()].second(), l[l.first_item()].second(), G);
                    }
                  else
                    {
                      TMP_CMP_GRAPHC = &G;
                      l.sort(&cmp_two_tuple);

                      list_item it;
                      //                      forall_items(it, l)
                      //                        cerr << G[l[it].second()].interpol_coord << "(e=" << G[l[it].first()].w  << "), ";
                      //                      cerr << endl;
                      
                      list_item it1 = l.first_item();
                      double left_edges = 0;
                      double right_edges = 0;
                      list_item it2 = l.succ(it1); while(it2!=nil) {right_edges+=G[l[it2].first()].w; it2 = l.succ(it2);}
                      double min_diff = fabs(right_edges - left_edges);
                      it2 = l.succ(it1);

                      G[v].interpol_coord = -1;
                      while(it2!=nil)
                        {
                          //                          cerr << "1. min_diff = " << min_diff << endl;
                          left_edges+=G[l[it1].first()].w;
                          if(fabs(right_edges - left_edges)<min_diff)
                            min_diff = fabs(right_edges - left_edges);
                          else
                            {
                              G[v].interpol_coord = G[l[it1].second()].interpol_coord+(G[l[it2].second()].interpol_coord-G[l[it1].second()].interpol_coord)/100.0;
                              update_real_neighbors(v, l[it1].second(), l[it1].second(), G);
                              break;
                            }
                          //                          cerr << "2. min_diff = " << min_diff << endl;
                          right_edges-=G[l[it2].first()].w;
                          if(fabs(right_edges - left_edges)<min_diff)
                            min_diff = fabs(right_edges - left_edges);
                          else
                            {
                              G[v].interpol_coord = (G[l[it1].second()].interpol_coord+G[l[it2].second()].interpol_coord)/2.0;
                              update_real_neighbors(v,l[it1].second(), l[it2].second(), G);
                              break;
                            }
                          //                          cerr << "3. min_diff = " << min_diff << endl;
                          it2 = l.succ(it2);
                          it1 = l.succ(it1);
                        }
                      if(G[v].interpol_coord == -1)
                        {                          
                          G[v].interpol_coord = G[l[it1].second()].interpol_coord-fabs(G[l[it1].second()].interpol_coord-G[l[l.pred(it1)].second()].interpol_coord)/100.0;
                          //                          cerr << "was not interpolated in cycle" << endl;
                        }
                      //                      cerr << G[v].interpol_coord << endl;
                      
                    }
                  
                  were_interpolated++;
                }
            }                        
        }
      thresh-=0.05;
    }
  
}

void interpolate_s2(TGraphC & G, TGraphC & H)
{
  list_item it;
 
  node v, w;
  edge e;

  cerr << "INTERPOL : " << H.number_of_nodes() << " -> " << G.number_of_nodes() << endl;
  TMP_CMP_GRAPHC = &H;
  H.sort_nodes(&cmp_ArrId);

  forall_nodes(v, G)
    G[v].interpol_coord = -1;
  
  ////////// new coarse coordinates ///////////
  double last_pos = 0;
  node left_neighbor = nil;
  node G_prev;
  
  forall_nodes(v, H)
    {      
      G_prev = H[v].G_prev_ptr;
      if(G[G_prev].initial_id!=H[v].initial_id)
        {
          cerr << "jopsa jopsa" << endl;
	  exit(1);
        }
      /*
      if((Params.basis_solutions_number>1)||(Params.V_Cycle_Iter_Num>1))
        {
          forall_nodes(G_prev, G)
            if(G[G_prev].initial_id==H[v].initial_id) break;
        }
      */
      G[G_prev].interpol_coord = last_pos + H[v].w/2.0;
      
      // update the real neighbors
      G[G_prev].left_real_neighbor = left_neighbor;
      if(left_neighbor != nil)
        G[left_neighbor].right_real_neighbor = G_prev;
      left_neighbor = G_prev;
      /////////////////////////////
      
      last_pos+=H[v].w;
    }
  G[G_prev].right_real_neighbor = nil;
  
  forall_nodes(v, G)
    {
      if((G[v].status==seed)&&(G[v].interpol_coord==-1))
        {
          cerr << "ERROR : Green node with interpolation coordinate = -1" << endl;
          exit(1);
        }
      if((G[v].status==seed)&&(G[v].left_real_neighbor==nil)&&(G[v].right_real_neighbor==nil))
        {
          cerr << "ERROR : Green node with nil real neighbors" << endl;
          exit(1);
        }
    }
      
  /////////////////////////////

  cerr << "|V(G)|=" << G.number_of_nodes() << " <<<=== |V(H)|=" << H.number_of_nodes() << endl;
  get_nodes_to_interpolate(G, 0.9, G.number_of_nodes()-H.number_of_nodes());

 
  
  return;
  /////////////////////////////
    
  forall_nodes(v, G)
    {
      if(G[v].status != seed)
        {          
          G[v].interpol_coord = 0;
          forall_adj_edges(e, v)
            {
              w = second_adj_for_edge(e, v, G);
              
              if(G[w].status == seed)
                G[v].interpol_coord+=G[w].interpol_coord*G[e].pw;              
            }
        }
    }
}
