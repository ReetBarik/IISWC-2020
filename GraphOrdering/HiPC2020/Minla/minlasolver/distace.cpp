#include "cmpfuncs.h"
#include <time.h>
#include <math.h>

void calc_wdegree_futurevolume(TGraphC & G) {
  
  cerr << "compute weighted degree future volume" << endl;
  
  node v, u;
  edge e, f;
  
  //double io = interpol_order(G);
  
  forall_nodes(v, G) {
	G[v].future_edges_sum = 0;
	G[v].wd = 0;
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	
	G[v].future_edges_sum+=fix_edge_relative_weight(G, e);
	G[u].future_edges_sum+=fix_edge_relative_weight(G, e);
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	
	if (G[u].status!=seed)
	  G[v].wd+=G[u].w*(fix_edge_relative_weight(G, e)/G[u].future_edges_sum);
	
	if (G[v].status!=seed)
	  G[u].wd+=G[v].w*(fix_edge_relative_weight(G, e)/G[v].future_edges_sum);
  }
  
  double all_wd = 0;
  
  forall_nodes(v, G) {
	G[v].wd+=G[v].w;
	//	cout  << G[v].wd << " " << G.degree(v) << endl;
  }
  
}


void mloga_Jacobi_get_best_rw_value(int rit, TGraphC& G, node_array<array<double> > & rv){
	edge e, e1;
	node u, v, w, t;

	//normalization of a random vector
	double s = 0;
	forall_nodes(v, G) 
	  s+=rv[v][rit]*rv[v][rit];
	
	s = sqrt(s);
	forall_nodes(v, G)
	  rv[v][rit] = rv[v][rit] / s;

	forall_edges(e, G) {
		v = G.source(e);
		w = G.target(e);

		double rw = log(fabs(rv[v][rit]-rv[w][rit]))/log(2);

		//double rw = fabs(rv[v][rit]-rv[w][rit]);
		//double rw = xdmax(fabs(G[v].rv[rit]), fabs(G[w].rv[rit]));
		//double rw = pow(G[v].rv[rit]-G[w].rv[rit], 4);

		//G[e].rw = xdmax(G[e].rw, rw);

		G[e].rw += rw;
	}
}

void Jacobi_get_best_rw_value(int rit, TGraphC& G, node_array<array<double> > & rv){
  edge e, e1;
  node u, v, w, t;
  
  //normalization of a random vector
  if(Params.alg_dist_normalize==true) {
	double minv = 10;
	double maxv = -10;
	forall_nodes(v, G) {
	  minv = xdmin(rv[v][rit], minv);
	  maxv = xdmax(rv[v][rit], maxv);
	}
	if(minv!=maxv) {
	  forall_nodes(v, G)
		rv[v][rit] = (rv[v][rit] - minv)/(maxv-minv) - .5;
	}
  }

  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);

	double rw = fabs(rv[v][rit]-rv[w][rit]);
	
	//double rw = fabs(rv[v][rit]-rv[w][rit]);
	//double rw = xdmax(fabs(G[v].rv[rit]), fabs(G[w].rv[rit]));
	//double rw = pow(G[v].rv[rit]-G[w].rv[rit], 4);
	
	//G[e].rw = xdmax(G[e].rw, rw);
	
	G[e].rw += rw*rw;
  }
}

void mloga_Jacobi_correct_rw2(TGraphC & G) {
	edge e;
	node v, w;
	

	double min_rw = G[G.first_edge()].rw;
	forall_edges(e, G) {
	 
	  if((min_rw > G[e].rw)&&(G[e].rw!=-1.0/0.0))
	    min_rw = G[e].rw;
	}
		cerr << min_rw << endl;
	//	exit(1);
	if(min_rw > 0)
	  min_rw = 0;
	forall_edges(e, G) {
	 // G[e].rw = 1.0/*(double)Params.alg_distance_rvec*/ / G[e].rw;
		//G[e].rw = pow(2.17,G[e].rw);
	  if(G[e].rw == -1.0/-0.0)
	    G[e].rw = min_rw;
	  G[e].rw = 1.0 / (-min_rw + G[e].rw + 0.0000000001);

	}
}
void Jacobi_correct_rw2(TGraphC & G) {
  edge e;
  node v, w;

  double same_cl = 0;
	double diff_cl = 0;
	forall_edges(e, G) {
	  // G[e].rw = 1.0/*(double)Params.alg_distance_rvec*/ / G[e].rw;
	  //G[e].rw = pow(2.17,G[e].rw);
	  //  if(G[e].rw == -1.0/-0.0)
	  //  G[e].rw = min_rw;
	  //G[e].rw = 1.0 / (-min_rw + G[e].rw + 0.0000000001);
	  //cout << G[e].rw << " ";
	  v = G.source(e);
	  w = G.target(e);
	  
	  G[e].rw = 1.0 / (sqrt(G[e].rw)/(double)Params.alg_distance_rvec + 0.0000000001);

	}
	if(Params.alg_dist == asym_gauss_seidel) {
		forall_nodes(v, G)
			G[v].rw_sum = 0;
		forall_edges(e, G) {
			//G[G.source(e)].rw_sum += G[e].rw;
			//G[G.target(e)].rw_sum += G[e].rw;
			G[G.source(e)].rw_sum = xdmax(G[G.source(e)].rw_sum,G[e].rw);
			G[G.target(e)].rw_sum = xdmax(G[G.target(e)].rw_sum,G[e].rw);
		}	
		forall_edges(e, G) {
			G[e].rw = xdmax(G[e].rw / G[G.source(e)].rw_sum, G[e].rw / G[G.target(e)].rw_sum);
		}
	}
	
}


int algdist_symJacobi_underrelax(TGraphC & G) {
  int RVEC = Params.alg_distance_rvec;
  int GS_ITER = Params.alg_distance_iter;
  //srand(time(0));
  if ((Params.alg_dist==gauss_seidel)||(Params.alg_dist==asym_gauss_seidel)) {
	if(Params.problem==mloga)
	  cerr << "MLOGA algebraic distance" << endl;
	cerr << "Algebraic distance coarsening: Jacobi underrelaxation, RVEC="	<< RVEC << " ITER=" << GS_ITER << endl;
  }

node v, w;
  edge e;
  //	cerr << RVEC  << endl;
  node_array<array<double> > rv(G);
  
  forall_edges(e, G) {
	G[e].rw = 0;
  }
  
  //cerr << "node volumes " << endl;
  forall_nodes(v, G) {
	//if(S==2)
	//		cerr  << G[v].w << " ";
	rv[v].resize(RVEC);
	for (int i=0; i<RVEC; i++)  {
	  rv[v][i] = -.5+(double)rand()/(double)RAND_MAX;
	}
  }
  //	cerr << "a" << endl;
  //calc_alg_wdegree(G);
  //color_and_reorder_vertices(G);
  //reorder_vertices_boundary(G);
  node_array<double> jac_rv(G);
  
  if(Params.alg_dist_volumes==future_wdeg) 
	calc_wdegree_futurevolume(G);
  else {
	  cerr << "algdist: no future volume" << endl;
  }
  
  edge_array<double> norm_factor(G);
  if(Params.alg_dist_volumes==future_wdeg) {    
    forall_edges(e, G) {
      v = G.source(e); w = G.target(e);
      norm_factor[e] = sqrt(G[v].wd*G[w].wd);
    }
  }
  else if(Params.alg_dist_volumes==wdeg) {
    forall_edges(e, G) {
      v = G.source(e); w = G.target(e);
      norm_factor[e] = sqrt(G[v].w*G[w].w);
    }
  }
  else {
    forall_edges(e, G) {
      norm_factor[e] = 1;
    }    
  }
  
  cerr << "Start iterations" << endl;
  for (int rit=0; rit<RVEC; rit++) {		
	forall_nodes(v, G)
	  jac_rv[v] = rv[v][rit];
	for (int gsit=0; gsit<GS_ITER; gsit++) {
	  forall_nodes(v, G) {
		if(G.degree(v)!=0) {
		  double s=0;
		  double new_c = 0;
		  forall_adj_edges(e, v) {
			w = second_adj_for_edge(e, v, G);
						
			new_c+=rv[w][rit]*G[e].w/norm_factor[e];
			s+=G[e].w/norm_factor[e];
		  }
		  jac_rv[v] = new_c/s;
		}
		else
		  jac_rv[v] = 0;
	  }
	  forall_nodes(v, G)  {
		rv[v][rit] = .5 * jac_rv[v] + .5 * rv[v][rit];
	  }
	}
	//				legalize_rv(G, rit);
	if(Params.problem==mloga)
	  mloga_Jacobi_get_best_rw_value(rit, G, rv);
	else
	  Jacobi_get_best_rw_value(rit, G, rv);
  }
  
  if(Params.problem==mloga)
	mloga_Jacobi_correct_rw2(G);
  else
	Jacobi_correct_rw2(G);
}
