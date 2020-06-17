#include "cmpfuncs.h"
#include <time.h>
#include <math.h>

void initialize_f2c_connections(node & v, TGraphC & G, list<edge> & tmp_ec_array) {
  edge e;
  tmp_ec_array.clear();
  forall_adj_edges(e, v) {
    node u = second_adj_for_edge(e, v, G);
    if (G[u].status==seed)
      tmp_ec_array.push_back(e);
  }
}

void choose_best_pair_of_f2c_connections( list<edge> & tmp_ec_array, edge & best_e1, edge & best_e2, double & best_cost, node & v, TGraphC & G) {
   // cerr << Params.min_part_size << endl;
  list_item eit1, eit2;
  if(tmp_ec_array.size() >= 2) {
    forall_items(eit1, tmp_ec_array) {
      forall_items(eit2, tmp_ec_array) {
	if(eit1<eit2) {
	  edge e1 = tmp_ec_array[eit1];
	  edge e2 = tmp_ec_array[eit2];
	  double cost = G[e1].rw + G[e2].rw;
	  node u1 = second_adj_for_edge( e1, v, G);
	  node u2 = second_adj_for_edge( e2, v, G);
	  if((G[u1].wd2 + G[v].w*G[e1].rw/cost <= Params.min_part_size)&&(G[u2].wd2 + G[v].w*G[e2].rw/cost <= Params.min_part_size)&&(cost > best_cost)) {
	    best_e1 = e1;
	    best_e2 = e2;
	    best_cost = cost;
	  }
	}
      }
    }
  }
}

void choose_best_pair_of_f2c_edgeweight_connections( list<edge> & tmp_ec_array, edge & best_e1, edge & best_e2, double & best_cost, node & v, TGraphC & G) {
  list_item eit1, eit2;
  if(tmp_ec_array.size() >= 2) {
    forall_items(eit1, tmp_ec_array) {
      forall_items(eit2, tmp_ec_array) {
	if(eit1<eit2) {
	  edge e1 = tmp_ec_array[eit1];
	  edge e2 = tmp_ec_array[eit2];
	  double cost = G[e1].w + G[e2].w;
	  node u1 = second_adj_for_edge( e1, v, G);
	  node u2 = second_adj_for_edge( e2, v, G);
	  if((G[u1].wd2 + G[v].w*G[e1].w/cost <= Params.min_part_size)&&(G[u2].wd2 + G[v].w*G[e2].w/cost <= Params.min_part_size)&&(cost > best_cost)) {
	    best_e1 = e1;
	    best_e2 = e2;
	    best_cost = cost;
	  }
	}
      }
    }
  }
}


void balanced_coarsening_operator_by_algdist(TGraphC & G, TGraphC & H) {

  edge e;
  node v, u;

  cerr << "Interpolation operator: balanced coarsening, algebraic distance" << endl;

  //cerr  << "mps " << Params.min_part_size << endl;
  // Recalculation od interpolation order
  int interpolation_order = interpol_order(G);
  //cerr << interpolation_order << endl; exit(1);
  int cc=0;
  int N = G.number_of_nodes();

  int added_fine_nodes = 0;

  double AVG_AGG_SIZE = 0;
  double TOTAL_SIZE = 0;
  forall_nodes(v, G) {
    TOTAL_SIZE += G[v].w;
    if (G[v].status!=seed) {
      G[v].wd2 = 0;
    }
    else {
      two_tuple<node, double> t(v, 1.0);
      G[v].iw_links.push_back(t);
      G[v].wd2 = G[v].w;
    }
  }
  AVG_AGG_SIZE = ceil(TOTAL_SIZE / (double)H.number_of_nodes()) * (1+Params.imbalance);
  cerr << "interpolation weights: total=" << TOTAL_SIZE << " |V_H|=" << H.number_of_nodes() << " " << AVG_AGG_SIZE << endl;
//exit(1);
  forall_edges(e, G)
    G[e].iw = 0;

  forall_nodes(v, G) {
    if (G[v].status!=seed) {
      //cerr << "start iw-conn fine id " << G[v].initial_id << " deg "<< G.degree(v) << endl;
      list<edge> tmp_ec_array;
      initialize_f2c_connections(v, G, tmp_ec_array);

      list_item it;

      TMP_CMP_GRAPHC = &G;
      //tmp_ec_array.sort(&cmp_edges);
      tmp_ec_array.sort(&cmp_edges_revrw);

      while(tmp_ec_array.size()>interpolation_order) {
		tmp_ec_array.erase(tmp_ec_array.last());
      }

      bool was_attached = false;
      //cerr  << "list len before pairs " << tmp_ec_array.length() << endl;
      //if(G.degree(v)!=1) {
		if(Params.CURRENT_LEVEL<Params.io2uptolevel) {
			edge best_e1, best_e2;
			double best_cost = -1;
			choose_best_pair_of_f2c_connections(tmp_ec_array, best_e1, best_e2, best_cost, v, G);

			if(best_cost > 0) {
				was_attached = true;
				tmp_ec_array.clear();
				tmp_ec_array.push_back(best_e1);
				tmp_ec_array.push_back(best_e2);
				node u1 = second_adj_for_edge( best_e1, v, G);
				node u2 = second_adj_for_edge( best_e2, v, G);

				G[u1].wd2 += G[v].w * G[best_e1].rw / best_cost;
				G[u2].wd2 += G[v].w * G[best_e2].rw / best_cost;
			}
		}
	//cerr  << "list len " << tmp_ec_array.length() << endl;
		list<edge> tmp_ec_array2 = tmp_ec_array;
		while((was_attached==false)&&(tmp_ec_array.size()>0)) {
		u = second_adj_for_edge( tmp_ec_array[tmp_ec_array.first()], v, G);
		if(G[u].wd2 + G[v].w <= Params.min_part_size) {
			while(tmp_ec_array.size()>1)
				tmp_ec_array.erase(tmp_ec_array.last());
			was_attached = true;
			G[u].wd2 += G[v].w;
			//cerr  << "attached" << endl;
			}
		else {/*
	    if((tmp_ec_array.size()==1)&&(H.number_of_nodes() > 0.8 * G.number_of_nodes())) {
	      u = second_adj_for_edge( tmp_ec_array2[tmp_ec_array2.first()], v, G);
	      was_attached = true;
	      G[u].wd2 += G[v].w;
	      break;
	    }
	    */
	    //cerr  << "erased first "   << G[u].wd2 <<"+"<< G[v].w<<" > "  << Params.min_part_size << endl;
			tmp_ec_array.erase(tmp_ec_array.first());
		}
	}
	//}
	//else {
	  //  was_attached = true;
	  //}
	 // cerr  << was_attached << " " << added_fine_nodes << endl;
	  if(was_attached==false) {
	    //cerr << "add coarse" << endl;
	    Cnode v_tmp;
	    v_tmp.ArrId = G[v].ArrId;
	    v_tmp.G_prev_ptr= v;
	    v_tmp.status = fine;
	    v_tmp.M = G[v].M;
	    v_tmp.initial_id = G[v].initial_id;
	    v_tmp.part_sol = G[v].part_sol;
	    node n = H.new_node(v_tmp);
	    G[v].status=seed;
	    G[v].ptr_to_coarse = n;
	    //		cerr << "added " << endl;
	    two_tuple<node, double> t(v, 1.0);
	    G[v].iw_links.push_back(t);

	    G[v].wd2 = G[v].w;
	    AVG_AGG_SIZE = ceil(TOTAL_SIZE / (double)H.number_of_nodes())* (1+Params.imbalance);
	    // cerr << H.number_of_nodes() << " " << AVG_AGG_SIZE << endl;
	    added_fine_nodes++;
	  }
	  else {
	    //cerr << "def fine" << endl;
	    double edges_to_seed_sum = 0;
	    node s, t;

	    forall_items(it, tmp_ec_array)
	    {
	      edges_to_seed_sum+=G[tmp_ec_array[it]].rw;
	    }
	    forall_items(it, tmp_ec_array) {
	      //G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;
	      G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].rw/edges_to_seed_sum;

	      s = G.source(tmp_ec_array[it]);
	      t = G.target(tmp_ec_array[it]);

	      two_tuple<node, double> t1(s, G[tmp_ec_array[it]].iw);
	      two_tuple<node, double> t2(t, G[tmp_ec_array[it]].iw);

	      G[s].iw_links.push_back(t2);
	      G[t].iw_links.push_back(t1);

	      if ((G[s].status==seed)&&(G[t].status==seed)) {
			cerr << "jopa" << endl;
			exit(1);
	      }
	      if ((G[s].status!=seed)&&(G[t].status!=seed)) {
			cerr << "jopa" << endl;
			exit(1);
	      }
	    }
	  }
    }
  }

  forall_nodes(v, G) {
    //if(G[v].status==seed) cerr << G[v].initial_id << " " << G[v].wd2 << endl;
    G[v].wd2 = 0;
  }
  cerr << "Interpolation weights: added " << added_fine_nodes << " fine nodes; AVG_AGG_SIZE " << AVG_AGG_SIZE  << endl;
 // exit(1);
}

void balanced_coarsening_operator_by_edgeweights(TGraphC & G, TGraphC & H) {

  edge e;
  node v, u;

  cerr << "Interpolation operator: balanced coarsening, edge weights" << endl;

  //cerr  << "mps " << Params.min_part_size << endl;
  // Recalculation od interpolation order
  int interpolation_order = interpol_order(G);
  int cc=0;
  int N = G.number_of_nodes();

  int added_fine_nodes = 0;

  double AVG_AGG_SIZE = 0;
  double TOTAL_SIZE = 0;
  forall_nodes(v, G) {
    TOTAL_SIZE += G[v].w;
    if (G[v].status!=seed) {
      G[v].wd2 = 0;
    }
    else {
      two_tuple<node, double> t(v, 1.0);
      G[v].iw_links.push_back(t);
      G[v].wd2 = G[v].w;
    }
  }
  AVG_AGG_SIZE = ceil(TOTAL_SIZE / (double)H.number_of_nodes()) * (1+Params.imbalance);
  cerr << "interpolation weights: total=" << TOTAL_SIZE << " |V_H|=" << H.number_of_nodes() << " " << AVG_AGG_SIZE << endl;

  forall_edges(e, G)
    G[e].iw = 0;

  forall_nodes(v, G) {
    if (G[v].status!=seed) {
      //cerr << "start iw-conn fine id " << G[v].initial_id << " deg "<< G.degree(v) << endl;
      list<edge> tmp_ec_array;
      initialize_f2c_connections(v, G, tmp_ec_array);

      list_item it;

      TMP_CMP_GRAPHC = &G;
      tmp_ec_array.sort(&cmp_edges);

      while(tmp_ec_array.size()>interpolation_order) {
	tmp_ec_array.erase(tmp_ec_array.last());
      }

      bool was_attached = false;
      //cerr  << "list len before pairs " << tmp_ec_array.length() << endl;
      //if(G.degree(v)!=1) {
	if(Params.CURRENT_LEVEL<Params.io2uptolevel) {
	  edge best_e1, best_e2;
	  double best_cost = -1;
	  choose_best_pair_of_f2c_edgeweight_connections(tmp_ec_array, best_e1, best_e2, best_cost, v, G);

	  if(best_cost > 0) {
	    was_attached = true;
	    tmp_ec_array.clear();
	    tmp_ec_array.push_back(best_e1);
	    tmp_ec_array.push_back(best_e2);
	    node u1 = second_adj_for_edge( best_e1, v, G);
	    node u2 = second_adj_for_edge( best_e2, v, G);

	    G[u1].wd2 += G[v].w * G[best_e1].w / best_cost;
	    G[u2].wd2 += G[v].w * G[best_e2].w / best_cost;
	  }
	}
	//cerr  << "list len " << tmp_ec_array.length() << endl;
	list<edge> tmp_ec_array2 = tmp_ec_array;
	while((was_attached==false)&&(tmp_ec_array.size()>0)) {
	  u = second_adj_for_edge( tmp_ec_array[tmp_ec_array.first()], v, G);
	  if(G[u].wd2 + G[v].w <= Params.min_part_size) {
	    while(tmp_ec_array.size()>1)
	      tmp_ec_array.erase(tmp_ec_array.last());
	    was_attached = true;
	    G[u].wd2 += G[v].w;
	    //cerr  << "attached" << endl;
	  }
	  else {/*
	    if((tmp_ec_array.size()==1)&&(H.number_of_nodes() > 0.8 * G.number_of_nodes())) {
	      u = second_adj_for_edge( tmp_ec_array2[tmp_ec_array2.first()], v, G);
	      was_attached = true;
	      G[u].wd2 += G[v].w;
	      break;
	    }
	    */
	    //cerr  << "erased first "   << G[u].wd2 <<"+"<< G[v].w<<" > "  << Params.min_part_size << endl;
	    tmp_ec_array.erase(tmp_ec_array.first());
	  }
	}
	//}
	//else {
	  //  was_attached = true;
	  //}
	  //cerr  << was_attached << " " << added_fine_nodes << endl;
	  if(was_attached==false) {
	    //cerr << "add coarse" << endl;
	    Cnode v_tmp;
	    v_tmp.ArrId = G[v].ArrId;
	    v_tmp.G_prev_ptr= v;
	    v_tmp.status = fine;
	    v_tmp.M = G[v].M;
	    v_tmp.initial_id = G[v].initial_id;
	    v_tmp.part_sol = G[v].part_sol;
	    node n = H.new_node(v_tmp);
	    G[v].status=seed;
	    G[v].ptr_to_coarse = n;
	    //		cerr << "added " << endl;
	    two_tuple<node, double> t(v, 1.0);
	    G[v].iw_links.push_back(t);

	    G[v].wd2 = G[v].w;
	    AVG_AGG_SIZE = ceil(TOTAL_SIZE / (double)H.number_of_nodes())* (1+Params.imbalance);
	    // cerr << H.number_of_nodes() << " " << AVG_AGG_SIZE << endl;
	    added_fine_nodes++;
	  }
	  else {
	    //cerr << "def fine" << endl;
	    double edges_to_seed_sum = 0;
	    node s, t;

	    forall_items(it, tmp_ec_array)
	    {
	      edges_to_seed_sum+=G[tmp_ec_array[it]].w;
	    }
	    forall_items(it, tmp_ec_array) {
	      //G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;
	      G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].w/edges_to_seed_sum;

	      s = G.source(tmp_ec_array[it]);
	      t = G.target(tmp_ec_array[it]);

	      two_tuple<node, double> t1(s, G[tmp_ec_array[it]].iw);
	      two_tuple<node, double> t2(t, G[tmp_ec_array[it]].iw);

	      G[s].iw_links.push_back(t2);
	      G[t].iw_links.push_back(t1);

	      if ((G[s].status==seed)&&(G[t].status==seed)) {
		cerr << "jopa" << endl;
		exit(1);
	      }

	      if ((G[s].status!=seed)&&(G[t].status!=seed)) {
		cerr << "jopa" << endl;
		exit(1);
	      }
	    }
	  }
    }
  }

  forall_nodes(v, G) {
    //if(G[v].status==seed) cerr << G[v].initial_id << " " << G[v].wd2 << endl;
    G[v].wd2 = 0;
  }
  cerr << "Interpolation weights: added " << added_fine_nodes << " fine nodes; AVG_AGG_SIZE " << AVG_AGG_SIZE  << endl;
}
