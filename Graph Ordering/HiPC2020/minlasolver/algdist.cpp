#include "cmpfuncs.h"
//#include "LEDA/basic_graph_alg.h"
//#include "LEDA/graph_misc.h"
#include <time.h>
#include <math.h>


void check_bad_rw(TGraphC& G) {
  cerr << "CHECKING BAD RW" << endl;
  int c = 0;
  int cs = 0;
  edge e;
  forall_edges(e, G) {
	if(G[e].rw==0) c++;
	if(G[e].rw<0.000000000000000001) cs++;
	if(std::isinf(G[e].rw))   { cerr << "infinite rw " << G[e].rw << endl; exit(1);}
	if(std::isnan(G[e].rw)) { cerr << "not a number rw " << G[e].rw << endl; exit(1);}
  }
  cerr << "Zero algdist edges: " << c << endl;
  cerr << "Less than 0.1e-17 algdist edges: " << cs << endl;
}

void print_Laplacian(TGraphC& G, leda::string f){
  node v, w;
  edge e;
  file_ostream o(f);
  o << "A = zeros(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  
  forall_edges(e, G)
  {
	v = G.source(e);
	w = G.target(e);
	
	o << "A(" << G[v].initial_id << ", " << G[w].initial_id << ")=" << -G[e].w << ";" << endl;
	o << "A(" << G[w].initial_id << ", " << G[v].initial_id << ")=" << -G[e].w << ";" << endl;
  }
  forall_nodes(v, G)
  {
	double s = 0;
	forall_adj_edges(e, v)
	{
	  s+=G[e].w;
	}
	o << "A(" << G[v].initial_id << ", " << G[v].initial_id << ")=" << s << ";" << endl;
  }
  
  o.close();
  /*
  forall_nodes(v, G)
  {
	forall_nodes(w, G)
	{
	  if(v!=w)
	  {
		if((e=is_edge(w, v, G))!=NULL)
		  o << G[e].w << " ";
		else
		  o << 0 << " ";
		
}
else
{
  double x=0;
  forall_adj_edges(e, v)
	x+=G[e].w;
  o << -x << " ";
}
}
o << endl;
}
cout << "Laplacian printed" << endl;
*/
}

void print_F2COperator(TGraphC& G, TGraphC& H, leda::string f){
  edge e;
  node v, w;
  file_ostream o(f);
  o << "M = zeros(" << G.number_of_nodes() << ", " << H.number_of_nodes() << ");" << endl;
  forall_edges(e, G)
  {
	v = G.source(e);
	w = G.target(e);
	
	if((G[v].status==seed)&&(G[w].status==fine)&&(G[e].iw!=0))
	{
	  o << "M(" << G[w].ArrId << ", " << H[G[v].ptr_to_coarse].ArrId << ")=" << G[e].iw << ";" << endl;
	}
	else if((G[w].status==seed)&&(G[v].status==fine)&&(G[e].iw!=0))
	{
	  o << "M(" << G[v].ArrId << ", " << H[G[w].ptr_to_coarse].ArrId << ")=" << G[e].iw << ";" << endl;
	}
  }
  forall_nodes(v, G)
  {
	if(G[v].status==seed)
	{
	  o << "M(" << G[v].ArrId << ", " << H[G[v].ptr_to_coarse].ArrId << ")=" << 1 << ";" << endl;
	}
  }
  o.close();
}
/*
void dot_write(TGraphC& G, leda::string f) {
  file_ostream o(f);
  o << "digraph AAA {\n";
  node v, w;
  edge e;
  forall_nodes(v, G)
	o << G[v].initial_id << "[label=\"" << G[v].initial_id << "\"];\n";
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	o << G[v].initial_id << "->" <<G[w].initial_id << "[label=\"" << G[e].w
	<< ", " << 1.0/G[e].rw << "\"];\n";
  }
  o << "}\n";
  o.close();
}
*/
double old_good_sort_fine_nodes_by_weighted_degree_s3(TGraphC & G) {
  if (Params.future_volume_sort==alg)
	cerr << "FVol sort by algebraic distance" << endl;
  else if (Params.future_volume_sort==nosort)
	cerr << "No FVol sort" << endl;
  else
	cerr << "FVol sort by weighted aggregation" << endl;
  
  //  return 0;
  node v, u;
  edge e, f;
  
  double io = interpol_order(G);
  
  forall_nodes(v, G) {
	G[v].future_edges_sum = 0;
	G[v].w_future_edges_sum = 0;
	G[v].wd = 0;
  }
  
  forall_edges(e, G) {
	
	//cout << G[e].rw << endl;
	v = G.source(e);
	u = G.target(e);
	//if((G[v].initial_id==16)||(G[u].initial_id==16))
	//{
	  //cout << G[e].rw << "\t" << G[e].w << endl;
	  //}
	  
	  double e_w;
	  
	  if (Params.future_volume_sort==alg)
		e_w = /*1.0/*/G[e].rw;
	  else
		e_w = fix_edge_relative_weight(G, e);
	  
	  G[v].future_edges_sum+=e_w;
	  G[u].future_edges_sum+=e_w;
	  G[v].w_future_edges_sum+=G[e].w;
	  G[u].w_future_edges_sum+=G[e].w;
  }
  //	exit(1);
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	double e_w;
	
	if (Params.future_volume_sort==alg)
	  e_w = /*1.0/*/G[e].rw;
	else
	  e_w = fix_edge_relative_weight(G, e);
	
	if (Params.future_volume_sort==alg) {
	  if (G[u].status!=seed) {
		//G[v].wd+=G[u].w * dblmin(1,(G.degree(u)*e_w/dblmin(1.0*io,ceil(G.degree(u)*Params.threshold_coarse_nodes_current)))/G[u].future_edges_sum);// G[u].w*e_w/G[u].future_edges_sum;
		G[v].wd+=G[u].w*(e_w/G[u].future_edges_sum + G[e].w/G[u].w_future_edges_sum);
		//        cerr << dblmin(1,(G.degree(u)*e_w/dblmin(1.0*io,ceil(G.degree(u)*Params.threshold_coarse_nodes_current)))/G[u].future_edges_sum) <<endl;
	  }
	  
	  if (G[v].status!=seed) {
		//G[u].wd+=G[v].w * dblmin(1,(G.degree(v)*e_w/dblmin(1.0*io,ceil(G.degree(v)*Params.threshold_coarse_nodes_current)))/G[v].future_edges_sum);//G[v].w*e_w/G[v].future_edges_sum;
		G[u].wd+=G[v].w*(e_w/G[v].future_edges_sum+G[e].w/G[v].w_future_edges_sum);
	  }
	} else {
	  if (G[u].status!=seed) {
		//G[v].wd+=G[u].w * dblmin(1,(G.degree(u)*e_w/dblmin(1.0*io,ceil(G.degree(u)*Params.threshold_coarse_nodes_current)))/G[u].future_edges_sum);// G[u].w*e_w/G[u].future_edges_sum;
		G[v].wd+=G[u].w*(e_w/G[u].future_edges_sum);
		//        cerr << dblmin(1,(G.degree(u)*e_w/dblmin(1.0*io,ceil(G.degree(u)*Params.threshold_coarse_nodes_current)))/G[u].future_edges_sum) <<endl;
	  }
	  
	  if (G[v].status!=seed) {
		//G[u].wd+=G[v].w * dblmin(1,(G.degree(v)*e_w/dblmin(1.0*io,ceil(G.degree(v)*Params.threshold_coarse_nodes_current)))/G[v].future_edges_sum);//G[v].w*e_w/G[v].future_edges_sum;
		G[u].wd+=G[v].w*(e_w/G[v].future_edges_sum);
	  }
	  
	}
  }
  
  double all_wd = 0;
  
  forall_nodes(v, G) {
	
	G[v].wd+=G[v].w;
	
	//      G[v].wd = 1; // TO REMOVE !!!!!!!!!!!!!!!!!!!!! G[v].wd+=G[v].w;
	
	if (G[v].status == seed)
	  G[v].wd = -1;
	else
	  all_wd+=G[v].wd;
  }
  
  if (/*(S!=0)&&*/(Params.random_reorder==true)) {
	//      cerr << "OK" << endl;
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  if ((Params.future_volume_sort!=nosort)) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_weighted_degree_C);
  }
  
  //  cerr << "NO FUTURE VOLUME SORT!" << endl;
  
  /*
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_nodes_degree);
  
  node w;
  forall_nodes(v, G)
  {
	cerr << G.degree(v) << " (" << G[v].wd << "); ";
	
	//      cerr << endl;
}
cerr << endl;
exit(1);
*/
  return all_wd;
}

double sort_fine_nodes_by_alg_and_wag_wdegree(TGraphC & G) {
  cerr << "FVol sort by mix of alg and wag" << endl;
  
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  node v, u;
  edge e, f;
  
  double io = interpol_order(G);
  
  forall_nodes(v, G) {
	G[v].future_edges_sum = 0;
	G[v].w_future_edges_sum = 0;
	G[v].wd = 0;
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	
	G[v].future_edges_sum+=G[e].rw;
	G[u].future_edges_sum+=G[e].rw;
	G[v].w_future_edges_sum+=fix_edge_relative_weight(G, e);
	G[u].w_future_edges_sum+=fix_edge_relative_weight(G, e);
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	
	if (G[u].status!=seed)
	  G[v].wd+=G[u].w*(.5*G[e].rw/G[u].future_edges_sum + .5*fix_edge_relative_weight(G, e)/G[u].w_future_edges_sum);
	
	if (G[v].status!=seed)
	  G[u].wd+=G[v].w*(.5*G[e].rw/G[v].future_edges_sum + .5*fix_edge_relative_weight(G, e)/G[v].w_future_edges_sum);
  }
  
  double all_wd = 0;
  
  forall_nodes(v, G) {
	
	G[v].wd+=G[v].w;
	
	if (G[v].status == seed)
	  G[v].wd = -1;
	else
	  all_wd+=G[v].wd;
  }
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_weighted_degree_C);
  
  return all_wd;
}

double calc_alg_wdegree(TGraphC & G) {
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  cerr << "FVol sort algebraic weighted degree" << endl;
  
  node v, u;
  edge e, f;
  
  double io = interpol_order(G);
  
  forall_nodes(v, G) {
	G[v].future_edges_sum = 0;
	G[v].wd = 0;
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	//cout << G[e].rw << endl;
	
	G[v].future_edges_sum+=G[e].rw;
	G[u].future_edges_sum+=G[e].rw;
  }
  forall_nodes(v, G) {
	if(G[v].future_edges_sum==0)
	  G[v].future_edges_sum = 1;
  }
  forall_edges(e, G) {
	
	v = G.source(e);
	u = G.target(e);
	//cout << G[u].future_edges_sum << endl;
	//cout << G[v].future_edges_sum << endl;
	if (G[u].status!=seed)
	{
	  G[v].wd+=G[u].w*(G[e].rw/G[u].future_edges_sum);
	  //if ((G[v].wd > 999999999)||(G[v].wd < .000000000000000001))
	  //cout << G[u].w << " " << G[e].rw << " " << G[u].future_edges_sum << endl;
	}
	if (G[v].status!=seed)
	{
	  G[u].wd+=G[v].w*(G[e].rw/G[v].future_edges_sum);
	  //if ((G[u].wd > 999999999)||(G[u].wd < .000000000000000001))
	  //cout << G[v].w << " " << G[e].rw << " " << G[v].future_edges_sum << endl;
	}
  }
  
  
  forall_nodes(v, G) {
	G[v].wd+=G[v].w;
  }
}

double sort_fine_nodes_by_alg_wdegree(TGraphC & G) {
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  cerr << "FVol sort algebraic weighted degree" << endl;
  
  node v, u;
  edge e, f;
  
  double io = interpol_order(G);
  
  forall_nodes(v, G) {
	G[v].future_edges_sum = 0;
	G[v].wd = 0;
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	u = G.target(e);
	//cout << G[e].rw << endl;
	
	G[v].future_edges_sum+=G[e].rw;
	G[u].future_edges_sum+=G[e].rw;
  }
  
  forall_edges(e, G) {
	
	v = G.source(e);
	u = G.target(e);
	//cout << G[u].future_edges_sum << endl;
	//cout << G[v].future_edges_sum << endl;
	if (G[u].status!=seed)
	{
	  G[v].wd+=G[u].w*(G[e].rw/G[u].future_edges_sum);
	  //if ((G[v].wd > 999999999)||(G[v].wd < .000000000000000001))
	  //cout << G[u].w << " " << G[e].rw << " " << G[u].future_edges_sum << endl;
	}
	if (G[v].status!=seed)
	{
	  G[u].wd+=G[v].w*(G[e].rw/G[v].future_edges_sum);
	  //if ((G[u].wd > 999999999)||(G[u].wd < .000000000000000001))
	  //cout << G[v].w << " " << G[e].rw << " " << G[v].future_edges_sum << endl;
	}
  }
  
  double all_wd = 0;
  
  forall_nodes(v, G) {
	G[v].wd+=G[v].w;
	
	if (G[v].status == seed)
	  G[v].wd = -1;
	else
	  all_wd+=G[v].wd;
  }
  
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_weighted_degree_C);
  //cout << "level " << S << endl;
  //forall_edges(e, G)
  //{
	//	v = G.source(e);
	//	node w = G.target(e);
	
	//cout << "sort: " << G[v].initial_id << "-" << G[w].initial_id << "(" //<< G.degree(w) << ")" << " w=" << G[e].w << " " << "rw=" << G[e].rw //<<  endl;
	//}
	//cout << "---------------" << endl;
	return all_wd;
}
int algdist_Jacobi_underrelax2(TGraphC & G) ;

double sort_fine_nodes_by_wag_wdegree(TGraphC & G) {
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  cerr << "FVol sort wag weighted degree" << endl;
  
  node v, u;
  edge e, f;
  
  double io = interpol_order(G);
  
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
	
	if (G[v].status == seed)
	  G[v].wd = -1;
	else
	  all_wd+=G[v].wd;
  }
  
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_weighted_degree_C);
  
  return all_wd;
}

double sort_fine_nodes_by_wdegree(TGraphC & G) {
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  double sum_fines_wd;
  if(Params.future_volume_sort==alg) {
	sum_fines_wd = sort_fine_nodes_by_alg_wdegree(G);
	//algdist_Jacobi_underrelax2(G);
  }
  else if(Params.future_volume_sort==mix_alg_wag)
	sum_fines_wd = sort_fine_nodes_by_alg_and_wag_wdegree(G);
  else if(Params.future_volume_sort==nosort)
	cerr << "sort_fine_nodes_by_wdegree: no future volume sort" << endl;
  else if(Params.future_volume_sort==wd)
	sum_fines_wd = sort_fine_nodes_by_wag_wdegree(G);
  
  return sum_fines_wd;
}



void add_heavy_nodes(TGraphC & G, TGraphC & H, double & sum_fines_wd) {
  cerr << "Adding global heavy nodes" << endl;
  cout << "sum_fines_wd = " << sum_fines_wd << endl;
  node v, u;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  double mid_wd = 0;
  //  forall_nodes(v, G)
  //    mid_wd+=G[v].wd;
  mid_wd = sum_fines_wd/*mid_wd*// (double)G.number_of_nodes();
  
  double sd_wd = 0;
  forall_nodes(v, G)
	sd_wd+=(G[v].wd - mid_wd)*(G[v].wd - mid_wd);
  sd_wd = sqrt(sd_wd / (double)G.number_of_nodes());
  
  forall_nodes(v, G) {
	//if (G[v].status != seed)
	//{
	  //cout << G[v].wd << " " <<Params.threshold_heavy_nodes*mid_wd << endl;
	  //}
	  if ((G[v].status != seed)&&(G[v].wd>=Params.threshold_heavy_nodes
		*/*sd_wd+*/mid_wd)) {
		//cerr << G[v].initial_id << "\t" << G.degree(v) << "\t" << G[v].wd << "\t" << mid_wd << endl;
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.M = G[v].M;
	  v_tmp.status = fine;
	  v_tmp.initial_id = G[v].initial_id;
	  
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  //           ret.push_back(v);
	  }
}
//cout << "Heavy node: |V(H)|=" << H.number_of_nodes() << endl;
//exit(1);
}
void mloga_add_locally_heavy_nodes(TGraphC & G, TGraphC & H) {
  if (Params.add_heavy==false) {
    cerr << "Don't add locally heavy vertices" << endl;
    return;
  }
  node v, u, w;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  forall_nodes(v, G) {
    //cout << G[v].wd << " " << G.degree(v) << endl;
    double avg_neigh_fvol = 0;
    forall_adj_nodes(w, v) {
      avg_neigh_fvol += G[w].wd;
    }
    avg_neigh_fvol = avg_neigh_fvol / (double)G.degree(v);
    if(G[v].wd > 5*avg_neigh_fvol) {
      Cnode v_tmp;
      v_tmp.ArrId = G[v].ArrId;
      v_tmp.G_prev_ptr = v;
      v_tmp.M = G[v].M;
      v_tmp.status = fine;
      v_tmp.initial_id = G[v].initial_id;
      
      node n = H.new_node(v_tmp);
      G[v].status=seed;
      G[v].ptr_to_coarse = n;		
      
      G[v].visited = true;
    }
  }
}
void old_add_locally_heavy_nodes(TGraphC & G, TGraphC & H) {
  cerr << "Add locally heavy vertices" << endl;
  node v, u, w;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  forall_nodes(v, G) {
	//cout << G[v].wd << " " << G.degree(v) << endl;
	double avg_neigh_fvol = 0;
	forall_adj_nodes(w, v) {
	  avg_neigh_fvol += G[w].wd;
	}
	avg_neigh_fvol = avg_neigh_fvol / (double)G.degree(v);
	if(G[v].wd > Params.threshold_heavy_nodes*avg_neigh_fvol) {
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.M = G[v].M;
	  v_tmp.status = fine;
	  v_tmp.initial_id = G[v].initial_id;
	  
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;		
	  
	  G[v].visited = true;
	}
  }
}

static int cmp_nddbl(const two_tuple<node, double> & p1,const two_tuple<node, double> & p2) {
  if (p1.second()<p2.second())
	return -1;
  else if (p1.second()>=p2.second())
	return 1;
  else
	return 0;
}
/*
void normalize_rvec_by_transform(int vec_num, TGraphC & G) {
  node v;
  double minv = 10;
  double maxv = -10;
  forall_nodes(v, G) {
	minv = xdmin(G[v].rv[vec_num], minv);
	maxv = xdmax(G[v].rv[vec_num], maxv);
	}
	if(minv==maxv)
	  return;
	forall_nodes(v, G)
	  G[v].rv[vec_num] = (G[v].rv[vec_num] - (minv+maxv)/2.0)/(maxv-minv);
	
	//minv = 10;
	//maxv = -10;
	//forall_nodes(v, G) {
	  //	minv = xdmin(G[v].rv2[vec_num], minv);
	  //	maxv = xdmax(G[v].rv2[vec_num], maxv);
	  //}
	  
	  //forall_nodes(v, G)
	  //	G[v].rv2[vec_num] = (G[v].rv2[vec_num] - (minv+maxv)/2.0)/(maxv-minv);
	  
	  }
	  */
void normalize_prev_rvec_by_transform(int vec_num, TGraphC & G) {
  cerr << "Old version: normalize_prev_rvec_by_transform" << endl;
  exit(-1);
  /*
  node v;
  double minv = 10;
  double maxv = -10;
  forall_nodes(v, G) {
	minv = xdmin(G[v].prev_rv[vec_num], minv);
	maxv = xdmax(G[v].prev_rv[vec_num], maxv);
}
if(minv==maxv)
  return;
forall_nodes(v, G)
  G[v].prev_rv[vec_num] = (G[v].prev_rv[vec_num] - (minv+maxv)/2.0)/(maxv-minv);
*/
}
/*
void new_normalize_rvec_by_transform(int vec_num, TGraphC & G) {
  //cerr << "no scaling" << endl;
  //return ;
  node v;
  double minv = 10;
  double maxv = -10;
  forall_nodes(v, G) {
	minv = xdmin(G[v].rv[vec_num], minv);
	maxv = xdmax(G[v].rv[vec_num], maxv);
	}
	
	if (fabs(minv-maxv)>0.000000001) {
	  double k=0.5/maxv;
	  double b=-.5-minv;
	  forall_nodes(v, G)
		G[v].rv[vec_num] = k*G[v].rv[vec_num]+b;
	  }
	  }
	  */
/*
void legalize_rv(TGraphC & G, int vec_num) {
  //  return;
  //  cerr << "Legalize random coordinates at vector " << vec_num << endl;
  
  array<two_tuple<node, double> > rv(G.number_of_nodes());
  node v;
  int i=0;
  double totalw = 0;
  forall_nodes(v, G) {
	two_tuple<node, double> p(v, G[v].rv[vec_num]);
	rv[i] = p;
	i++;
	totalw+=G[v].w;
	//  cerr << i << endl;
	}
	//cerr << "ogo" << endl;
	rv.sort(&cmp_nddbl);
	double N = G.number_of_nodes();
	double prev=-.5;
	for (i=0; i<N; i++) {
	  // cerr << i << endl;
	  v = rv[i].first();
	  // cerr << "a" << endl;
	  G[v].rv[vec_num] = prev+G[v].w/(2.0*totalw);//((double)(i+1))/N;
	  prev+=G[v].w/totalw;
	  }
	  //cerr << "done" << endl;
	  //rv.resize(0);
	  //cerr << "ogo" << endl;
	  }
	  */
/*
void normalize_rw_by_endpoints(TGraphC & G) {
  cerr << "Normalizing algebraic distances" << endl;
  node v, w;
  edge e;
  
  forall_nodes(v, G)
	G[v].all_edges_sum = 0;
  
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	
	G[v].all_edges_sum+=G[e].rw;
	G[w].all_edges_sum+=G[e].rw;
	//	G[v].w_all_edges_sum+=G[e].w;
	//	G[w].w_all_edges_sum+=G[e].w;
	
	}
	
	forall_edges(e, G) {
	  v = G.source(e);
	  w = G.target(e);
	  G[e].rw = (G[e].rw / sqrt(G[v].all_edges_sum*G[w].all_edges_sum));
	  //G[e].w = (G[e].w / sqrt(G[v].w_all_edges_sum*G[w].w_all_edges_sum));
	  //G[e].rw = G[e].rw*(G[v].all_edges_sum+G[w].all_edges_sum)/(G[v].all_edges_sum*G[w].all_edges_sum);
	  }
	  }
	  
	  void normalize_ew_by_endpoints(TGraphC & G) {
		
		cerr << "Normalizing edge weights" << endl;
		//exit(1);
		node v, w;
		edge e;
		
		
		forall_nodes(v, G)
		  G[v].all_edges_sum = 0;
		
		forall_edges(e, G) {
		  v = G.source(e);
		  w = G.target(e);
		  
		  G[v].all_edges_sum+=G[e].rw;
		  G[w].all_edges_sum+=G[e].rw;
		  
		  }
		  forall_edges(e, G) {
			v = G.source(e);
			w = G.target(e);
			G[e].rw = dblmin(G[e].rw / sqrt(G[v].all_edges_sum),G[e].rw / sqrt(G[w].all_edges_sum));//(G[e].rw / sqrt(G[v].all_edges_sum*G[w].all_edges_sum));
		  }
		  
		  }
		  */
void Jacobi_correct_rw(TGraphC & G) {
  edge e;
  node v, w;
  check_bad_rw(G);
  forall_edges(e, G)
	if(G[e].rw==0)
	  G[e].rw = 0.000000000001;
	
	forall_nodes(v, G) {
	  double w_v = G[v].w;
	  G[v].min_v = 100000000;
	  forall_adj_edges(e, v) {
		w = second_adj_for_edge(e, v, G);
		G[v].min_v = xdmin(G[v].min_v, G[e].rw * xdmin(w_v, G[w].w));
	  }
	}
	forall_edges(e, G) {
	  v = G.source(e);
	  w = G.target(e);
	  double ad1 = G[e].rw * xdmin(G[v].w, G[w].w)/G[v].min_v;
	  double ad2 = G[e].rw * xdmin(G[v].w, G[w].w)/G[w].min_v;
	  G[e].rw = 1.0/xdmax(ad1, ad2);
	}
	
	cerr << "Jacobi algebraic distances corrected" << endl;
}




/*
void normalize_and_save_edges(TGraphC & G) {
  node v, w;
  edge e;
  
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	
	G[v].all_edges_sum+=G[e].w;
	G[w].all_edges_sum+=G[e].w;
	
	G[e].old_W = G[e].w;
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	G[e].w = G[e].w / sqrt(G[v].all_edges_sum*G[w].all_edges_sum);
  }
}

void restore_after_normalization(TGraphC & G) {
  edge e;
  
  forall_edges(e, G) {
	G[e].w = G[e].old_W;
  }
}
*/
extern double current_e_entropy;
extern double current_e_sdv;

void rw_by_cossim(TGraphC& G, node_array<array<double> > & rv) {
  edge e;
  node v, w;
  int RVEC = Params.alg_distance_rvec;
  
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	
	G[e].rw = 0;
	
	double vs = 0;
	double ws = 0;
	for(int i=0; i<RVEC; i++) {
	  double vc = rv[v][i]+1;
	  double wc = rv[w][i]+1;
	  
	  G[e].rw+=vc * wc;
	  vs+=vc*vc;
	  ws+=wc*wc;
	}
	G[e].rw = 1.0/(1 - G[e].rw/(sqrt(vs)*sqrt(ws)));
  }
  
  /*
  TMP_CMP_GRAPHC = &G;
  G.sort_edges(&cmp_edges_rw);
  double same_cl = 0;
  double diff_cl = 0;
  int ed=0;
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	//ed++;
	if(G[v].part_sol==G[w].part_sol) {
	  same_cl++;
	  //cout << G[e].rw << " = ";// << same_cl/(double)G.number_of_edges() << endl;
}
else {
  diff_cl++;
  //	cout << G[e].rw << " ! ";// << diff_cl/(double)G.number_of_edges() << endl;
  
}
//cout <<  " " << diff_cl << endl;
//cout << same_cl/(double)G.number_of_edges() << " " << diff_cl/(double)G.number_of_edges() << endl;

}
//cout << endl;
//exit(1);
*/
}


/*
void calculate_edge_entropy(TGraphC & G) {
  cerr << "error: calculate_edge_entropy\n";
  exit(-1);
  
  int RVEC = Params.alg_distance_rvec;
  edge e;
  forall_edges(e, G)
  {
	G[e].H = 0;
	for(int h=0; h<5; h++)
	{
	  //if(G[e].rw_distr[h]>0)
	  //	G[e].H-=G[e].rw_distr[h]/(double)RVEC*log(G[e].rw_distr[h]/(double)RVEC);
	}
  }
}
*/
void color_and_reorder_vertices(TGraphC & G) {
  node v, w;
  edge e;
  forall_nodes(v, G)
	G[v].col = -1;
  forall_nodes(v, G)
  {
	list<int> adjcols; adjcols.clear();
	forall_adj_nodes(w, v)
	{
	  if(G[w].col!=-1)
		adjcols.push_back(G[w].col);
	}
	if(adjcols.size()==0)
	  G[v].col = 1;
	else
	{
	  adjcols.sort();
	  adjcols.unique();
	  int c = 1;
	  list_item it = adjcols.first();
	  while(G[v].col==-1)
	  {
		if(it==nil)
		  G[v].col = c;
		else if(c<adjcols[it])
		  G[v].col = c;
		else
		{
		  c++;
		  it = adjcols.succ(it);
		}
	  }
	}
	
  }
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_color);
}

void reorder_vertices_boundary(TGraphC & G) {
  node v, w;
  edge e;
  forall_nodes(v, G)
	G[v].wd = 0;
  forall_edges(e, G)
  {
	v = G.source(e); w = G.target(e);
	//G[v].wd+=G[e].w; G[w].wd+=G[e].w;
	G[v].wd+=1; G[w].wd+=1;
  }
  
  TMP_CMP_GRAPHC = &G;
  cerr << "return sort nodes?" << endl;
  exit(1);
  //G.sort_nodes(&cmp_weighted_degree_rev);
}
/*
void rescale_rvec(int vec_num, TGraphC & G) {
  node v;
  double minv = G.number_of_nodes()+100;
  double maxv = -100;
  forall_nodes(v, G) {
	minv = xdmin(G[v].rv[vec_num], minv);
	maxv = xdmax(G[v].rv[vec_num], maxv);
	}
	//cerr << "minv = " << minv << " maxv=" << maxv << endl;
	
	forall_nodes(v, G) {
	  G[v].rv[vec_num] = (double)G.number_of_nodes() * (G[v].rv[vec_num] - minv) / (double)(maxv - minv);
	  if(G[v].rv[vec_num] < 0) {
		cerr << maxv << " " << minv << endl;
		exit(1);
	  }
	  }
	  }
	  */

void random_walk(TGraphC & G) {
  node rv = G.choose_node();
  node v, u;
  edge e,f;
  /*
  forall_nodes(v, G) G[v].rmark = 0;
  for(int i=0; i<10*G.number_of_nodes(); i++) {
	G[v].rmark++;
	double s = 0;
	forall_adj_edges(e, v) {
	  s+=G[e].w;
}
double vv = gen_prob();
double p = 0;
bool flag = false;
forall_adj_edges(e, v) {
  p+=G[e].w;
  if((p > vv)&&(flag==false)) {
	f = e;
	flag = true;
}
}
v = second_adj_for_edge(e, v, G);
}
*/
}

int rvec_global_num;



int cossim_algdist_symJacobi_underrelax(TGraphC & G) {
  cerr << "cossim" << endl;
  exit(1); // correct alg_dist_volumes
  int RVEC = Params.alg_distance_rvec;
  int GS_ITER = Params.alg_distance_iter;
  //srand(time(0));
  if( (Params.alg_dist==gauss_seidel)||(Params.alg_dist==maximum_matching)) {
	cerr << "cosine Algebraic distance coarsening: Jacobi underrelaxation, RVEC="	<< RVEC << " ITER=" << GS_ITER << endl;
  }
  if(Params.sym_nonsym==symmetric) {
	//cerr << "Symmetric relaxation" << endl;
  }
  node v, w;
  edge e;
  //	cerr << RVEC  << endl;
  node_array<array<double> > rv(G);
  
  forall_edges(e, G) {
	G[e].rw = 0;
  }
  
  forall_nodes(v, G) {
	rv[v].resize(RVEC);
	for (int i=0; i<RVEC; i++)  {
	  rv[v][i] = -.5+(double)rand()/(double)RAND_MAX;
	}
  }
  node_array<double> jac_rv(G);
  
  //cout << "eee " << G.number_of_edges() << endl;
  for (int rit=0; rit<RVEC; rit++) {		
	//cout << "random vector " << rit << endl;
	forall_nodes(v, G)
	  jac_rv[v] = rv[v][rit];
	for (int gsit=0; gsit<GS_ITER; gsit++) {
	  forall_nodes(v, G) {
		if(G.degree(v)!=0) {
		  double s=0;
		  double new_c = 0;
		  forall_adj_edges(e, v) {
			w = second_adj_for_edge(e, v, G);
			
			
			if(Params.alg_dist_volumes==true) {
			  new_c+=rv[w][rit]*G[e].w/sqrt(G[v].w*G[w].w);
			  s+=G[e].w/sqrt(G[v].w*G[w].w);
			}
			else {
			  new_c+=rv[w][rit]*G[e].w;
			  s+=G[e].w;						
			}
		  }
		  jac_rv[v] = new_c/s;
		}
		else
		  jac_rv[v] = 0;
	  }
	  forall_nodes(v, G) 
		rv[v][rit] = .5 * jac_rv[v] + .5 * rv[v][rit];
	}
	//				legalize_rv(G, rit);
	Jacobi_get_best_rw_value(rit, G, rv);
  }
  
  rw_by_cossim(G, rv);	
  
}


int gauss_seidel_algdist_symJacobi_underrelax(TGraphC & G) {
  exit(-1); // correct algdist_volumes
  int RVEC = Params.alg_distance_rvec;
  int GS_ITER = Params.alg_distance_iter;
  //srand(time(0));
  if (Params.alg_dist==gauss_seidel) {
	cerr << "Algebraic distance coarsening: Jacobi underrelaxation, RVEC="	<< RVEC << " ITER=" << GS_ITER << endl;
  }
  if(Params.sym_nonsym==symmetric) {
	//cerr << "Symmetric relaxation" << endl;
  }
  node v, w;
  edge e;
  //	cerr << RVEC  << endl;
  node_array<double> rv(G);
  
  forall_edges(e, G) {
	G[e].rw = 0;
  }
  
  //	cerr << "a" << endl;
  //calc_alg_wdegree(G);
  //color_and_reorder_vertices(G);
  //reorder_vertices_boundary(G);
  node_array<double> jac_rv(G);
  
  //cout << "eee " << G.number_of_edges() << endl;
  for (int rit=0; rit<RVEC; rit++) {		
	forall_nodes(v, G) {
	  rv[v] = -.5+(double)rand()/(double)RAND_MAX;
	  jac_rv[v] = rv[v];
	}
	//cout << "random vector " << rit << endl;
	for (int gsit=0; gsit<GS_ITER; gsit++) {
	  forall_nodes(v, G) {
		if(G.degree(v)!=0) {
		  double s=0;
		  double new_c = 0;
		  forall_adj_edges(e, v) {
			w = second_adj_for_edge(e, v, G);
			
			if(Params.alg_dist_volumes==true) {
			  new_c+=rv[w]*G[e].w/sqrt(G[v].w*G[w].w);
			  s+=G[e].w/sqrt(G[v].w*G[w].w);
			}
			else {
			  new_c+=rv[w]*G[e].w;
			  s+=G[e].w;						
			}
		  }
		  rv[v] = new_c/s;
		}
		else
		  rv[v] = 0;
	  }
	  //forall_nodes(v, G) 
	  //  rv[v] = .5 * jac_rv[v] + .5 * rv[v];
	}
	
	if(Params.alg_dist_normalize==true) {
	  double minv = 10;
	  double maxv = -10;
	  forall_nodes(v, G) {
		minv = xdmin(rv[v], minv);
		maxv = xdmax(rv[v], maxv);
	  }
	  if(minv!=maxv) {
		forall_nodes(v, G)
		  rv[v] = (rv[v] - minv)/(maxv-minv) - .5;
	  }
	}
	forall_edges(e, G) {
	  v = G.source(e);
	  w = G.target(e);
	  double rw = fabs(rv[v]-rv[w]);
	  G[e].rw += rw*rw;
	}
  }
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	
	G[e].rw = 1.0 / (sqrt(G[e].rw)/(double)Params.alg_distance_rvec + 0.0000000001);
  }
}


int algdist_symJacobi_underrelax_stage2(TGraphC & G) {
  cerr << "stage2" << endl;
  exit(-1);
  int RVEC = Params.alg_distance_rvec;
  int GS_ITER = Params.alg_distance_iter;
  //srand(time(0));
  if (Params.alg_dist==gauss_seidel) {
	cerr << "Algebraic distance coarsening: Jacobi underrelaxation, RVEC="	<< RVEC << " ITER=" << GS_ITER << endl;
  }
  if(Params.sym_nonsym==symmetric) {
	//cerr << "Symmetric relaxation" << endl;
  }
  node v, w;
  edge e;
  //	cerr << RVEC  << endl;
  node_array<array<double> > rv(G);
  
  forall_edges(e, G) {
	G[e].rw = 0;
  }
  
  forall_nodes(v, G) {
	rv[v].resize(RVEC);
	for (int i=0; i<RVEC; i++)  {
	  if(G[v].status==seed)
		rv[v][i]=-.5;
	  else
		rv[v][i]=.5;
	}
  }
  
  node_array<double> jac_rv(G);
  
  //cout << "eee " << G.number_of_edges() << endl;
  for (int rit=0; rit<RVEC; rit++) {		
	//cout << "random vector " << rit << endl;
	for (int gsit=0; gsit<GS_ITER; gsit++) {
	  forall_nodes(v, G) {
		if(G[v].status!=seed) {
		  //cout << "a " << G.degree(v) << " ";
		  if(G.degree(v)!=0) {
			double s=0;
			double new_c = 0;
			forall_adj_edges(e, v) {
			  w = second_adj_for_edge(e, v, G);
			  if(Params.alg_dist_volumes==true) {
				new_c+=rv[w][rit]*G[e].w/sqrt(G[v].w*G[w].w);
				s+=G[e].w/sqrt(G[v].w*G[w].w);
			  }
			  else {
				new_c+=rv[w][rit]*G[e].w;
				s+=G[e].w;						
			  }
			}
			jac_rv[v] = new_c/s;
			//	cout << jac_rv[v] << " " << new_c << " " << s << "** ";	
		  }
		  else
			jac_rv[v] = 0;	
		  //cerr << jac_rv[v] << " " << new_c << " " << s << "** ";			
		}
	  }
	  //cout << endl;
	  forall_nodes(v, G) 
		rv[v][rit] = .5 * jac_rv[v] + .5 * rv[v][rit];
	}
	//				legalize_rv(G, rit);
	Jacobi_get_best_rw_value(rit, G, rv);
  }
  Jacobi_correct_rw2(G);	
  
  //forall_edges(e, G)
  //cout << G[e].rw << " ";
  //cout << endl;
  //normalize_ew_by_endpoints(G);
  //cout << "cc = " << cc << endl;
  //exit(1);
}


double HGdiv(TGraphC & G, TGraphC & H) {
  return (double)H.number_of_nodes()/(double)G.number_of_nodes();
}

double find_algdist_pivot(TGraphC & G, int p) {
  TMP_CMP_GRAPHC = &G;
  G.sort_edges(&cmp_edges_revrw);
  
  double ret_val;
  edge e;
  int ed=0;
  if(p==100) {
	return G[G.first_edge()].rw;
  } else {
	e = G.first_edge();
	double i=0;
	while(e!=nil) {
	  if(100.0*((double)G.number_of_edges()-i)/(double)G.number_of_edges() < p) {
		ret_val = G[e].rw;
		e = nil;				
	  }
	  else {
		i++;
		e = G.succ_edge(e);
	  }
	}
	return ret_val;
  }	
}


void localize_rw(TGraphC & G) {
  node u, v, w;
  edge e;
  forall_nodes(v, G) {
	G[v].wd = 0;
  }
  
  forall_edges(e, G)
  {
	u = G.source(e);
	v = G.target(e);
	
	G[u].wd+=G[e].rw;
	G[v].wd+=G[e].rw;
  }	
  
  forall_nodes(v, G) {
	G[v].wd = G[v].wd / (double)G.degree(v);
  }
  
  forall_edges(e, G)
  {
	u = G.source(e);
	v = G.target(e);
	
	G[e].rw = xdmax(G[e].rw / G[v].wd, G[e].rw / G[u].wd);		
  }	
  
}

bool bignodeprob(TGraphC & G, node & v, node & u, double & AV_NSIZE) {
  double p = gen_prob();
  
  
}

void algdist_matching(TGraphC & G, TGraphC & H){
  cerr << "Algebraic distance matching" << endl;
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  //localize_rw(G);
  
  TMP_CMP_GRAPHC = &G;
  //if(S<4)
  G.sort_edges(&cmp_erw/*_vols*/);
  //else
  //G.sort_edges(&cmp_edges);
  
  edge e, f;
  node u, v;
  
  double AV_NSIZE = 0;
  forall_nodes(v, G) {
	AV_NSIZE+=G[v].w;
  }
  AV_NSIZE = AV_NSIZE/(double)G.number_of_nodes();//Params.min_graph_size;	
  
  
  //double rw_pivot = find_algdist_pivot(G, 80);
  
  int gmatch = 0;
  int bmatch = 0;
  
  forall_edges(e, G)
  {
	u = G.source(e);
	v = G.target(e);
	
	if(/*(G[e].rw > rw_pivot)&&*/(G[u].matched!=true)&&(G[v].matched!=true)/*&&(G[v].w < 2*AV_NSIZE)&&(G[u].w < 2*AV_NSIZE)*/)
	{
	  //cerr << "merging " << G[u].w << " " << G[v].w << " " << G[e].rw << endl;
	  G[u].matched = true;
	  G[v].matched = true;
	  
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.matched = false;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  v_tmp.w = G[u].w + G[v].w;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  G[u].status=seed;
	  G[u].ptr_to_coarse = n;
	  
	  if((G[u].part_sol==G[v].part_sol) &&(G[v].part_sol!=-1)&&(G[u].part_sol!=-1)){
		gmatch++;
		H[n].part_sol = G[u].part_sol;
	  }
	  else {
		bmatch++;
		H[n].part_sol = -1;
	  }
	  
	}
  }
  cerr  << "good matches " << gmatch << " bad matches " << bmatch << endl;
  //exit(1);
  forall_nodes(v, G)
  {
	if(G[v].matched==false)
	{
	  G[v].matched = true;
	  
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.matched = false;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  v_tmp.w = G[v].w;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  
	  H[n].part_sol = G[v].part_sol;
	}
  }
  
  forall_edges(e, G) {
	u = G.source(e);
	v = G.target(e);
	
	node Hu = G[u].ptr_to_coarse;
	node Hv = G[v].ptr_to_coarse;
	
	if(Hu!=Hv) {
	  if((f=is_edge(Hu, Hv, H))!=nil) {
		H[f].w+=G[e].w;
		H[f].real_w+=G[e].w;
	  }
	  else {
		CEdge ne;
		ne.w = G[e].w;
		ne.real_w = G[e].w;
		H.new_edge(Hv, Hu, ne);
	  }
	}	
  }	
}

void maxmatching_coarsening(TGraphC & G, TGraphC & H){
  if (Params.random_reorder==true) {
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_rnd);
  }
  
  TMP_CMP_GRAPHC = &G;
  G.sort_edges(&cmp_erw);
  
  edge e;
  node u, v;
  forall_edges(e, G)
  {
	u = G.source(e);
	v = G.target(e);
	
	if((G[u].matched!=true)&&(G[v].matched!=true))
	{
	  G[u].matched = true;
	  G[v].matched = true;
	  
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  
	  // link on itself
	  two_tuple<node, double> t(v, 1.0);
	  G[v].iw_links.push_back(t);
	  
	  // links on fine and backward
	  two_tuple<node, double> t1(u, 1.0);
	  two_tuple<node, double> t2(v, 1.0);
	  
	  G[u].iw_links.push_back(t2);
	  G[v].iw_links.push_back(t1);
	}
  }
  
  forall_nodes(v, G)
  {
	if(G[v].matched==false)
	{
	  G[v].matched = true;
	  
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  
	  // link on itself
	  two_tuple<node, double> t(v, 1.0);
	  G[v].iw_links.push_back(t);
	}
  }
}

void wag_coarse_nodes(TGraphC & G, TGraphC & H/*, list<node> & ret*/) {
  
  cerr << "Weighted aggregation (classical)" << endl;
  
  node v, u;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  /*
  double sum_fines_wd = sort_fine_nodes_by_weighted_degree_s3(G, false);
  
  add_heavy_nodes(G, H, sum_fines_wd);
  */
  double sum_fines_wd = sort_fine_nodes_by_wdegree(G);
  
  cerr << "Heavy coarse nodes : " << H.number_of_nodes() << endl;
  
  list_item it;
  
  int H_size = H.number_of_nodes();
  //   int resort_done = false;
  //   cerr << T << endl;
  
  //  while(HGdiv(G, H)<0.6)
  {
	
	//      cerr << H.number_of_nodes() << " " << HGdiv(G, H) << " " << T << endl;
	forall_nodes(v, G) {
	  
	  if (G[v].status != seed)//&&(HGdiv(G, H)<0.6))
	  {
		//           cerr << "\t" << G[v].initial_id << endl;
		double edges_to_seed_sum = 0;
		double all_edges_sum = 0;
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  double e_w = fix_edge_relative_weight(G, e);
		  if (G[u].status==seed)
			edges_to_seed_sum+=e_w;
		  all_edges_sum+=e_w;
		}
		//      cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
		if (edges_to_seed_sum/all_edges_sum < T) {
		  //               cerr << G[v].initial_id << endl;
		  Cnode v_tmp;
		  v_tmp.ArrId = G[v].ArrId;
		  v_tmp.G_prev_ptr = v;
		  v_tmp.status = fine;
		  v_tmp.M = G[v].M;
		  v_tmp.initial_id = G[v].initial_id;
		  node n = H.new_node(v_tmp);
		  G[v].status=seed;
		  G[v].ptr_to_coarse = n;
		  H_size++;
		}
	  } //else
	  //break;
	}
	//      T = T+0.05;
  }
  /*
  forall_nodes(v, G)
  if(G[v].status==seed)
  cout << G[v].initial_id << endl;
  */
  //exit(1);
  /*
  forall_nodes(v, G)
  {
	
	if((G[v].status != seed)&&(HGdiv(G, H)<0.6))
	{
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  H_size++;
}
}

*/
  
  /*
  TMP_CMP_GRAPHC = &G;
  G.sort_edges(&cmp_edges);
  
  
  double kk=0;
  double cc = 0;
  node w;
  forall_edges(e, G)
  {
	v = G.source(e);
	w = G.target(e);
	leda::string v1s, v2s;
	if(G[v].status==seed)
	v1s = "c";
	else
	  v1s = "f";
	if(G[w].status==seed)
	v2s = "c";
	else
	  v2s = "f";
	if((G[w].status!=seed)&&(G[v].status!=seed))
	kk++;
	cc++;
	
	cout << G[e].w << " " << v1s << " " << v2s << " " << kk/cc <<endl;
}
cout << "------------------------------" << endl;
*/
  if (H.number_of_nodes()>Params.min_graph_size)
	return;
  else {
	if ((double)G.number_of_nodes()*Params.amount_of_coarse_nodes
	  <(double)Params.min_graph_size)
	  Params.amount_of_coarse_nodes = (double)Params.min_graph_size
	  /(double)G.number_of_nodes();
  }
  
  cerr << Params.amount_of_coarse_nodes << endl;
  bool enought_coarse_nodes=false;
  while (enought_coarse_nodes==false) {
	cerr << "T=" << T << "; ret len=" << H.number_of_nodes() << endl;
	
	if (fabs(H.number_of_nodes()-G.number_of_nodes()
	  *Params.amount_of_coarse_nodes)>0.0001) {
	  T = T+0.025;
	sort_fine_nodes_by_wdegree(G);
	
	forall_nodes(v, G) {
	  if ((G[v].status != seed)&&(fabs(H.number_of_nodes()
		-G.number_of_nodes()*Params.amount_of_coarse_nodes)
		>0.0001)) {
		double edges_to_seed_sum = 0;
	  double all_edges_sum = 0;
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);
		if (G[u].status==seed)
		  edges_to_seed_sum+=fix_edge_relative_weight(G, e);
		all_edges_sum+=fix_edge_relative_weight(G, e);
	  }
	  if (edges_to_seed_sum/all_edges_sum < T) {
		//                       cerr << G[v].initial_id << endl;
		Cnode v_tmp;
		v_tmp.ArrId = G[v].ArrId;
		v_tmp.G_prev_ptr = v;
		v_tmp.status = fine;
		v_tmp.M = G[v].M;
		v_tmp.initial_id = G[v].initial_id;
		node n = H.new_node(v_tmp);
		G[v].status=seed;
		G[v].ptr_to_coarse = n;
		//                       ret.push_back(v);
	  }
	  }
	}
	} else
	  enought_coarse_nodes = true;
  }
  
  return;
}

void coarse_nodes_alg_and_wag_sum(TGraphC & G, TGraphC & H/*, list<node> & ret*/) {
  
  cerr << "coarse nodes: algebraic distance and weighted aggretation sum"
  << endl;
  
  node v, u, w;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  /*
  double sum_fines_wd = sort_fine_nodes_by_weighted_degree_s3(G, false);
  
  add_heavy_nodes(G, H, sum_fines_wd);
  */
  
  double sum_fines_wd = sort_fine_nodes_by_wdegree(G);
  
  cerr << "Heavy coarse nodes : " << H.number_of_nodes() << endl;
  
  list_item it;
  
  int H_size = H.number_of_nodes();
  forall_nodes(v, G) {
	// cerr << "vdeg: " << G.degree(v) << endl;
	if (G[v].status != seed)//&&(HGdiv(G, H)<0.6))
	{
	  //           cerr << "\t" << G[v].initial_id << endl;
	  double edges_to_seed_sum = 0;
	  double all_edges_sum = 0;
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);
		double e_w = G[e].rw;//fix_edge_relative_weight(G, e);
		if (G[u].status==seed)
		  edges_to_seed_sum+=e_w;
		all_edges_sum+=e_w;
	  }
	  double w_edges_to_seed_sum = 0;
	  double w_all_edges_sum = 0;
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);
		double w_e_w = /*G[e].nw;*/fix_edge_relative_weight(G, e);
		if (G[u].status==seed)
		  w_edges_to_seed_sum+=w_e_w;
		w_all_edges_sum+=w_e_w;
	  }
	  
	  //      cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
	  
	  if ((.5*w_edges_to_seed_sum/w_all_edges_sum+.5*edges_to_seed_sum
		/all_edges_sum < T)&&(G.degree(v)>1)) {
		//               cerr << G[v].initial_id << endl;
	  //cerr << "added v" << endl;
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  H_size++;
	  }
	} else
	  break;
  }
  
  forall_edges(e, G) {
	v = G.source(e);
	w = G.target(e);
	
	if ((G.degree(v)==1)&&(G[w].status!=seed)) {
	  Cnode v_tmp;
	  v_tmp.ArrId = G[w].ArrId;
	  v_tmp.G_prev_ptr = w;
	  v_tmp.status = fine;
	  v_tmp.M = G[w].M;
	  v_tmp.initial_id = G[w].initial_id;
	  node n = H.new_node(v_tmp);
	  G[w].status=seed;
	  G[w].ptr_to_coarse = n;
	  H_size++;
	} else if ((G.degree(w)==1)&&(G[v].status!=seed)) {
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  H_size++;
	}
  }
  //return;
  if (H.number_of_nodes()>Params.min_graph_size)
	return;
  else {
	if ((double)G.number_of_nodes()*Params.amount_of_coarse_nodes
	  <(double)Params.min_graph_size)
	  Params.amount_of_coarse_nodes = (double)Params.min_graph_size
	  /(double)G.number_of_nodes();
  }
  
  cerr << Params.amount_of_coarse_nodes << endl;
  bool enought_coarse_nodes=false;
  while (enought_coarse_nodes==false) {
	cerr << "T=" << T << "; ret len=" << H.number_of_nodes() << endl;
	
	if (fabs(H.number_of_nodes()-G.number_of_nodes()
	  *Params.amount_of_coarse_nodes)>0.0001) {
	  T = T+0.025;
	sort_fine_nodes_by_wdegree(G);
	
	forall_nodes(v, G) {
	  if ((G[v].status != seed)&&(fabs(H.number_of_nodes()
		-G.number_of_nodes()*Params.amount_of_coarse_nodes)
		>0.0001)) {
		double edges_to_seed_sum = 0;
	  double all_edges_sum = 0;
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);
		if (G[u].status==seed)
		  edges_to_seed_sum+=G[e].rw;//fix_edge_relative_weight(G, e);
		all_edges_sum+=G[e].rw;//fix_edge_relative_weight(G, e);
	  }
	  double w_edges_to_seed_sum = 0;
	  double w_all_edges_sum = 0;
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);
		double w_e_w = fix_edge_relative_weight(G, e);
		if (G[u].status==seed)
		  w_edges_to_seed_sum+=w_e_w;
		w_all_edges_sum+=w_e_w;
	  }
	  
	  //      cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
	  
	  if (.5*w_edges_to_seed_sum/w_all_edges_sum+.5
		*edges_to_seed_sum/all_edges_sum < T) {
		
		//if (edges_to_seed_sum/all_edges_sum < T) {
  //                       cerr << G[v].initial_id << endl;
  Cnode v_tmp;
  v_tmp.ArrId = G[v].ArrId;
  v_tmp.G_prev_ptr = v;
  v_tmp.status = fine;
  v_tmp.M = G[v].M;
  v_tmp.initial_id = G[v].initial_id;
  node n = H.new_node(v_tmp);
  G[v].status=seed;
  G[v].ptr_to_coarse = n;
  //                       ret.push_back(v);
	  }
	  }
	}
	} else
	  enought_coarse_nodes = true;
  }
  
  return;
}

void gradual_coarse_nodes(TGraphC & G, TGraphC & H) {
  cerr << "coarse nodes: gradual coarsening" << endl;
  
  node v, u, w;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  
}


void add_locally_heavy_nodes(TGraphC & G, TGraphC & H) {
  
  cerr << "locally heavy coarse nodes: algebraic distances of edges" << endl;
  
  node v, u, w;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  cerr << "Heavy coarse nodes before adding locally heavy: " << H.number_of_nodes() << endl;
  cerr << "locally heavy threshold: " << T << endl;
  list_item it;
  
  double pivot_rw = find_algdist_pivot(G, 90);
  cout << "pivot rw " << pivot_rw  << endl;
  
  int H_size = H.number_of_nodes();
  forall_nodes(v, G) {
	// cerr << "vdeg: " << G.degree(v) << endl;
	if (G[v].status != seed)//&&(HGdiv(G, H)<0.6))
	{
	  //           cerr << "\t" << G[v].initial_id << endl;
	  double edges_to_seed_sum = 0;
	  double all_edges_sum = 0;
	  
	  //double edges_to_seed_gsum = 0;
	  //double all_edges_gsum = 0;
	  
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);				
		double e_w = G[e].rw;
		//double e_w = fix_edge_relative_weight(G, e);
		//double e_gw = fix_edge_relative_weight(G, e);
		if ((G[u].status!=seed)&&(e_w>pivot_rw))	{				
		  edges_to_seed_sum+=e_w;
		  //edges_to_seed_gsum+=e_gw;
		  
		}
		all_edges_sum+=e_w;
		//all_edges_gsum+=e_gw;
	  }
	  //cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
	  if ((edges_to_seed_sum/all_edges_sum> T/2.0)&&(H.number_of_nodes() < T*G.number_of_nodes()))//||(edges_to_seed_gsum/all_edges_gsum < T))/*&&(G.degree(v)>1)*/
	  {
		//               cerr << G[v].initial_id << endl;
		//cerr << "added v" << endl;
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
		H_size++;
	  }
	  //else
	  //	cout << G[v].initial_id << "\t" << G.degree(v) << "\t" << G[v].wd << endl;
	} //else
	//break;
  }
  cerr << "Heavy coarse nodes after adding locally heavy: " << H.number_of_nodes() << endl;
}

void deg1_neighbors_to_seeds(TGraphC & G, TGraphC & H) {
  node v, w;
  edge e;
  
  forall_nodes(v, G) {
	if((G.degree(v)==1)&&(G[v].status!=seed)) {
	  forall_adj_edges(e, v) {
		w = second_adj_for_edge(e, v, G);
		if(G[w].status!=seed) {
		  Cnode v_tmp;
		  v_tmp.ArrId = G[w].ArrId;
		  v_tmp.G_prev_ptr= w;
		  v_tmp.status = fine;
		  v_tmp.M = G[w].M;
		  v_tmp.initial_id = G[w].initial_id;
		  v_tmp.part_sol = G[w].part_sol;
		  node n = H.new_node(v_tmp);
		  G[w].status=seed;
		  G[w].ptr_to_coarse = n;
		  
		}
	  }
	}
  }
  		int deg1seeds=0;
		 forall_nodes(v, G) {
	   if((G[v].status==seed)&&(G.degree(v)==1)) {
	     
	     deg1seeds++;
	   }
	 }
	 cout << "deg1_neighbors_to_seeds: G seeds with degree 1 = " << deg1seeds << endl;
   
}


void coarse_nodes_alg_sum(TGraphC & G, TGraphC & H/*, list<node> & ret*/) {
  
  cerr << "coarse nodes: algebraic distance sum" << endl;
  
  node v, u, w;
  edge e;
  int i;
  double T = Params.threshold_coarse_nodes_current;
  
  
  //double sum_fines_wd = sort_fine_nodes_by_weighted_degree_s3(G, false);
  
  //add_heavy_nodes(G, H, sum_fines_wd);
  
  
  double sum_fines_wd = sort_fine_nodes_by_wdegree(G);
  
  cerr << "Heavy coarse nodes : " << H.number_of_nodes() << endl;
  
  list_item it;
  
  //double pivot_rw = find_algdist_pivot(G, 50);
    
  int H_size = H.number_of_nodes();
  forall_nodes(v, G) {
	//cerr << "=============== id " << G[v].initial_id << " vdeg: " << G.degree(v) << endl;
	if ((G[v].status != seed))//&&(G.degree(v)>1))
	{
	  //cerr << "\t" << G.degree(v) << endl;
	  double edges_to_seed_sum = 0;
	  double all_edges_sum = 0;
	  
	  double edges_to_seed_gsum = 0;
	  double all_edges_gsum = 0;
	  
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);				
		double e_w = G[e].rw;
		//cerr << e_w << " " ;
		//double e_w = fix_edge_relative_weight(G, e);
		double e_gw = fix_edge_relative_weight(G, e);
		if ((G[u].status==seed))//&&(e_w>pivot_rw))//&&(G[u].w < AV_CSIZE * IMB_FACTOR))
		{
		  //cout << "da" << endl;
		  edges_to_seed_sum+=e_w;
		  edges_to_seed_gsum+=e_gw;
		  
		}
		all_edges_sum+=e_w;
		all_edges_gsum+=e_gw;
	  }
	  //cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
	  
	  if ((edges_to_seed_sum/all_edges_sum < T))//||(edges_to_seed_gsum/all_edges_gsum < T))/*&&(G.degree(v)>1)*/
	  {
		//               cerr << G[v].initial_id << endl;
		//cerr << "added v" << endl;
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
		H_size++;
	  }
	  //else
	  //	cout << G[v].initial_id << "\t" << G.degree(v) << "\t" << G[v].wd << endl;
	} //else
	//break;
  }

		int deg1seeds=0;
		 forall_nodes(v, G) {
	   if((G[v].status==seed)&&(G.degree(v)==1)) {
	     
	     deg1seeds++;
	   }
	 }
	 cout << "G seeds with degree 1 = " << deg1seeds << endl;
   
  cerr << "before deg-1 neighbours : " << H.number_of_nodes() << endl;
  //deg1_neighbors_to_seeds(G, H);
  cerr << "after deg-1 neighbours : " << H.number_of_nodes() << endl;


  double max_rw = 0;
  double av_rw = 0;
  forall_edges(e, G)
  {
	av_rw += G[e].rw;
	if(max_rw < G[e].rw)
	  max_rw = G[e].rw;
  }
  av_rw = av_rw / (double)G.number_of_edges();
  cout << "av_rw=" << av_rw << ", max_rw=" << max_rw << endl;
  
  cerr << "coarse_nodes_alg_sum: coarse nodes " << H.number_of_nodes() << endl;
  cerr << "coarse_nodes_alg_sum: min size for external part "  << Params.external_part_size << endl;
  if (H.number_of_nodes()>Params.external_part_size)//Params.min_graph_size)
  return;
  
  
  //cerr << Params.amount_of_coarse_nodes << endl;
  bool enought_coarse_nodes=false;
  while (enought_coarse_nodes==false) {
	cerr << "T=" << T << "; |V_H|=" << H.number_of_nodes() << " |V_G|=" << G.number_of_nodes() << endl;
	
	//if (fabs(H.number_of_nodes()-G.number_of_nodes()*Params.amount_of_coarse_nodes)>0.0001) {
  if(H.number_of_nodes()<Params.external_part_size) {//Params.min_graph_size) {
  T = T+0.025;
  //sort_fine_nodes_by_wdegree(G);
  
  forall_nodes(v, G) {
	if ((G[v].status != seed)&&(H.number_of_nodes()<Params.external_part_size)) {//(fabs(H.number_of_nodes()-G.number_of_nodes()*Params.amount_of_coarse_nodes)>0.0001)) {
	  double edges_to_seed_sum = 0;
	  double all_edges_sum = 0;
	  forall_adj_edges(e, v) {
		u = second_adj_for_edge(e, v, G);
		if (G[u].status==seed)
		  edges_to_seed_sum+=G[e].rw;//fix_edge_relative_weight(G, e);
		all_edges_sum+=G[e].rw;//fix_edge_relative_weight(G, e);
	  }
	  if ((edges_to_seed_sum/all_edges_sum < T)&&(H.number_of_nodes()<Params.external_part_size)) {//Params.min_graph_size)) {
	  //                       cerr << G[v].initial_id << endl;
	  Cnode v_tmp;
	  v_tmp.ArrId = G[v].ArrId;
	  v_tmp.G_prev_ptr = v;
	  v_tmp.status = fine;
	  v_tmp.M = G[v].M;
	  v_tmp.initial_id = G[v].initial_id;
	  v_tmp.part_sol = G[v].part_sol;
	  node n = H.new_node(v_tmp);
	  G[v].status=seed;
	  G[v].ptr_to_coarse = n;
	  //                       ret.push_back(v);
	  }
	}
  }
  } else
	enought_coarse_nodes = true;
  }
  
  return;
  }
  
  void coarse_nodes_alg_max (TGraphC & G, TGraphC & H/*, list<node> & ret*/) {
	
	cerr << "coarse nodes: algebraic distance max" << endl;
	
	node v, u, w;
	edge e;
	int i;
	double T = Params.threshold_coarse_nodes_current;
	
	/*
	double sum_fines_wd = sort_fine_nodes_by_weighted_degree_s3(G, false);
	
	add_heavy_nodes(G, H, sum_fines_wd);
	*/
	
	double sum_fines_wd = sort_fine_nodes_by_wdegree(G);
	
	cerr << "Heavy coarse nodes : " << H.number_of_nodes() << endl;
	
	list_item it;
	
	int H_size = H.number_of_nodes();
	T = Params.threshold_coarse_nodes_current;;//0.8;
	forall_nodes(v, G) {
	  // cerr << "vdeg: " << G.degree(v) << endl;
	  if (G[v].status != seed)//&&(HGdiv(G, H)<0.6))
	  {
		//           cerr << "\t" << G[v].initial_id << endl;
		double edges_to_seed_sum = 0;
		//double all_edges_sum = 0;
		
		//double edges_to_seed_gsum = 0;
		//double all_edges_gsum = 0;
		double max_edge_to_seed = 0;
		
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  double e_w = G[e].rw;//fix_edge_relative_weight(G, e);
		  //double e_gw = fix_edge_relative_weight(G, e);
		  if (G[u].status==seed)
		  {
			edges_to_seed_sum+=e_w;
			//edges_to_seed_gsum+=e_gw;
			
			if(max_edge_to_seed < e_w)
			  max_edge_to_seed = e_w;
			
		  }
		  //all_edges_sum+=e_w;
		  //all_edges_gsum+=e_gw;
		}
		//      cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
		if ((max_edge_to_seed/edges_to_seed_sum > T)&&(G.degree(v)>1)) {
		  //               cerr << G[v].initial_id << endl;
		  //cerr << "added v" << endl;
		  
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
		  H_size++;
		  
		}
	  } //else
	  //break;
	}
	
	
	/*
	int added_counter = 0;
	forall_nodes(v, G) {
	if (G[v].status != seed)
	{
	  double edges_to_seed_sum = 0;
	  double edges_to_seed_c = 0;
	  double all_edges_sum = 0;
	  double max_to_seed = 0;
	  
	  forall_adj_edges(e, v) {
	u = second_adj_for_edge(e, v, G);
	double e_w = G[e].rw;
	if (G[u].status==seed)
	{
	  //cout << e_w << " ";
	  max_to_seed = xdmax(max_to_seed, e_w);
	  edges_to_seed_sum+=e_w;
	  edges_to_seed_c++;
  }
  }
  //cout << endl;
  //      cerr << edges_to_seed_sum << " " << all_edges_sum << endl;
  if ((edges_to_seed_sum/edges_to_seed_c*1.05 > max_to_seed)&&(G.degree(v)>1))
  {
	added_counter++;
	//cout << "adding " << edges_to_seed_sum/edges_to_seed_c << " " << max_to_seed << endl;
	Cnode v_tmp;
	v_tmp.ArrId = G[v].ArrId;
	v_tmp.G_prev_ptr= v;
	v_tmp.status = fine;
	v_tmp.M = G[v].M;
	v_tmp.initial_id = G[v].initial_id;
	node n = H.new_node(v_tmp);
	G[v].status=seed;
	G[v].ptr_to_coarse = n;
	H_size++;
  }
  }// else
  //break;
  }
  
  cerr << "Added " << added_counter << " vertices with weak algebraic connections" << endl;
  */
	
	
	/*
	forall_nodes(v, G)
	if(G[v].status==seed)
	cout << G[v].initial_id << endl;
	*/
	//exit(1);
	/*
	forall_nodes(v, G)
	{
	  
	  if((G[v].status != seed)&&(HGdiv(G, H)<0.6))
	  {
		Cnode v_tmp;
		v_tmp.ArrId = G[v].ArrId;
		v_tmp.G_prev_ptr = v;
		v_tmp.status = fine;
		v_tmp.M = G[v].M;
		v_tmp.initial_id = G[v].initial_id;
		node n = H.new_node(v_tmp);
		G[v].status=seed;
		G[v].ptr_to_coarse = n;
		H_size++;
  }
  }
  
  */
	
	if (H.number_of_nodes()>Params.min_graph_size)
	  return;
	else {
	  if ((double)G.number_of_nodes()*Params.amount_of_coarse_nodes
		<(double)Params.min_graph_size)
		Params.amount_of_coarse_nodes = (double)Params.min_graph_size
		/(double)G.number_of_nodes();
	}
	cerr << Params.amount_of_coarse_nodes << endl;
	bool enought_coarse_nodes=false;
	while (enought_coarse_nodes==false) {
	  cerr << "T=" << T << "; ret len=" << H.number_of_nodes() << " " << G.number_of_nodes() << endl;
	  
	  //if (fabs(H.number_of_nodes()-G.number_of_nodes()*Params.amount_of_coarse_nodes)>0.0001) {
	if(H.number_of_nodes()<Params.min_graph_size) {
	  T = T+0.025;
	  //sort_fine_nodes_by_wdegree(G);
	  
	  forall_nodes(v, G) {
		if ((G[v].status != seed)&&(H.number_of_nodes()<Params.min_graph_size)) {//(fabs(H.number_of_nodes()-G.number_of_nodes()*Params.amount_of_coarse_nodes)>0.0001)) {
		  double edges_to_seed_sum = 0;
		  double all_edges_sum = 0;
		  forall_adj_edges(e, v) {
			u = second_adj_for_edge(e, v, G);
			if (G[u].status==seed)
			  edges_to_seed_sum+=G[e].rw;//fix_edge_relative_weight(G, e);
			all_edges_sum+=G[e].rw;//fix_edge_relative_weight(G, e);
		  }
		  if ((edges_to_seed_sum/all_edges_sum < T)&&(H.number_of_nodes()<Params.min_graph_size)) {
			//                       cerr << G[v].initial_id << endl;
			Cnode v_tmp;
			v_tmp.ArrId = G[v].ArrId;
			v_tmp.G_prev_ptr = v;
			v_tmp.status = fine;
			v_tmp.M = G[v].M;
			v_tmp.initial_id = G[v].initial_id;
			node n = H.new_node(v_tmp);
			G[v].status=seed;
			G[v].ptr_to_coarse = n;
			//                       ret.push_back(v);
		  }
		}
	  }
	} else
	  enought_coarse_nodes = true;
	}
	
	return;
  }
  
  void calc_pw(TGraphC & G) {
	cerr << "Error: old calc_pw" << endl;
	exit(1);
	/*	
	edge e;
	node v, u;
	forall_edges(e, G)
		  G[e].pw = 0;
	
	forall_nodes(v, G) {
		  if (G[v].status!=seed) {
		  double edges_to_seed_sum = 0;
		  forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if (G[u].status==seed)
		  edges_to_seed_sum+=fix_edge_relative_weight(G, e);
  }
  forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if (G[u].status==seed)
		  G[e].pw = fix_edge_relative_weight(G, e)/edges_to_seed_sum;
  }
  }
  }
  */
  }
  
  void clear_iw_links(TGraphC & G) {
	node v;
	forall_nodes(v, G)
	  G[v].iw_links.clear();
	//cerr << "A" << endl;
  }
  
  int interpol_order(TGraphC & G) {
	/*
	int interpolation_order = (int)((double)Params.maximum_of_connections_between_fine_and_its_coarses+log(E0/(double)G.number_of_edges())+0.5);
	if(interpolation_order < Params.maximum_of_connections_between_fine_and_its_coarses)
		  interpolation_order = Params.maximum_of_connections_between_fine_and_its_coarses;
	if(interpolation_order > 500)
		  interpolation_order = 500;
	if(S==0)
		  interpolation_order = Params.maximum_of_connections_between_fine_and_its_coarses;
	cerr << "Recalculated interpolation order = " << interpolation_order << endl;
	
	return interpolation_order;
	*/
	
	//cerr << "Recalculated interpolation order = "	<< (int)xdmax(2,Params.maximum_of_connections_between_fine_and_its_coarses-S) << endl;
	//return (int)xdmax(2,Params.maximum_of_connections_between_fine_and_its_coarses-S);
	
	return Params.maximum_of_connections_between_fine_and_its_coarses;
	/*
	cerr << "Recalculated interpolation order = "	<< (int)xdmin(10,Params.maximum_of_connections_between_fine_and_its_coarses+S) << endl;
	return (int)xdmin(10,Params.maximum_of_connections_between_fine_and_its_coarses+S);
	*/
	//return Params.maximum_of_connections_between_fine_and_its_coarses;
  }
  
  void alg_calc_iw(TGraphC & G) {
	
	cerr << "Interpolation weights: sort by algebraic d istance" << endl;
	
	// Recalculationod interpolation order
	int interpolation_order = interpol_order(G);
	
	int cc=0;
	
	edge e;
	node v, u;
	int N = G.number_of_nodes();
	
	forall_edges(e, G)
	  G[e].iw = 0;
	
	forall_nodes(v, G) {
	  if (G[v].status!=seed) {
		
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if ((G[u].status==seed))//&&(G[u].initial_id<401))
		  tmp_ec_array.push_back(e);
		}
		
		TMP_CMP_GRAPHC = &G;
		tmp_ec_array.sort(&cmp_rnd_edges);
		tmp_ec_array.sort(&cmp_edges_revrw);
		list_item it;
		
		while (tmp_ec_array.size()>interpolation_order/*Params.get_maximum_of_connections_between_fine_and_its_coarses(N)*/)
		  tmp_ec_array.erase(tmp_ec_array.last());
		
		double edges_to_seed_sum = 0;
		
		forall_items(it, tmp_ec_array)
		  //edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
		edges_to_seed_sum+=G[tmp_ec_array[it]].w;
		
		node s, t;
		forall_items(it, tmp_ec_array) {
		  //G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw
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
	  } else {
		two_tuple<node, double> t(v, 1.0);
		G[v].iw_links.push_back(t);
	  }
	}
  }
  
  /*
  double stable_vertex(node v, TGraphC & G) {
	
	double v_H = 0;
	edge e;
	forall_adj_edges(e, v)
	  v_H = G[e].H;
	return v_H/(double)G.degree(v);
  }
  */
  void mloga_alg_calc_iw_by_rw(TGraphC & G) {

	edge e;
	node v, u;
/*
	double w_sum = 0;
	double rw_sum = 0;
	forall_edges(e, G) {
		w_sum+=G[e].w;
		rw_sum+=G[e].rw;

		u = G.source(e);
		v = G.target(e);
		G[u].w_sum+=G[e].w;
		G[u].rw_sum+=G[e].rw;
		G[v].w_sum+=G[e].w;
		G[v].rw_sum+=G[e].rw;
	}
	int wc = 0; int rwc = 0;
	forall_nodes(v, G) {
		if(G[v].w_sum/w_sum > G[v].rw_sum/rw_sum)
			wc++;
		else
			rwc++;
	}
	cout << "wc = " << wc << ", rwc = " << rwc << endl;
*/
	cerr << "Interpolation weights: sort by algebraic distance" << endl;

	// Recalculation od interpolation order
	int interpolation_order = interpol_order(G);
	int cc=0;
	int N = G.number_of_nodes();

	forall_edges(e, G)
	G[e].iw = 0;

	forall_nodes(v, G) {
		if (G[v].status!=seed) {

			list<edge> tmp_ec_array;
			tmp_ec_array.clear();
			forall_adj_edges(e, v) {
				u = second_adj_for_edge(e, v, G);
				if ((G[u].status==seed))//&&(G[u].initial_id<401))
					tmp_ec_array.push_back(e);
			}


			list_item it;

			TMP_CMP_GRAPHC = &G;
			//if(wc>rwc)
			//	tmp_ec_array.sort(&cmp_edges);
			//else
			tmp_ec_array.sort(&cmp_edges_revrw);
			//tmp_ec_array.sort(&cmp_edges_approx_w1_rw2);
			//tmp_ec_array.sort(&cmp_edges_approx_rw1_w2);
			/*
			forall_items(it, tmp_ec_array)
			{
				cerr << G[tmp_ec_array[it]].rw << " ";
			}
			cerr << endl;
			*/
			double half_size = (double)tmp_ec_array.size()*0.5;


			if((int)half_size>0)
			{
				//list<edge> tmp_ec_array2;
				int c1 = 0; int c2 = 0;

				while(G[tmp_ec_array[tmp_ec_array.first()]].rw> 2 * G[tmp_ec_array[tmp_ec_array.last()]].rw)
				{
					tmp_ec_array.erase(tmp_ec_array.last());
					c2++;
				}
				
				if(c2<1) {

					while((tmp_ec_array.size()>half_size)&&(tmp_ec_array.size()>interpolation_order))
					{
						tmp_ec_array.erase(tmp_ec_array.last());
						//	c1++;
					}
				}
				
				//cerr << "half size " << c1 << "\t i potom eshe " << c2 << endl;
				//cerr << "half weigh " << c2 << "\t i potom eshe " << c1 << endl;
				
				TMP_CMP_GRAPHC = &G;
				tmp_ec_array.sort(&cmp_edges);
			}



			while (tmp_ec_array.size()>interpolation_order/*Params.get_maximum_of_connections_between_fine_and_its_coarses(N)*/)
				tmp_ec_array.erase(tmp_ec_array.last());
			double edges_to_seed_sum = 0;
			node s, t;

			forall_items(it, tmp_ec_array)
			{
				//cout << "taking " << G[tmp_ec_array[it]].w << endl;
				edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
			}

			forall_items(it, tmp_ec_array) {
				//			if(v_H > .1)
				//				G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;
				//			else
				G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;

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
		} else {
			two_tuple<node, double> t(v, 1.0);
			G[v].iw_links.push_back(t);
		}
	}

}

  void alg_calc_iw_by_rw(TGraphC & G, TGraphC & H) {
	
	edge e;
	node v, u;
	
	//double pivot_rw = find_algdist_pivot(G, 0);
	/*
	double AV_CSIZE = 0;
	double IMB_FACTOR = (1+Params.imbalance) * (1.0 + (S+1)*0.1);
	forall_nodes(v, G)
		  AV_CSIZE+=G[v].w;
	AV_CSIZE = AV_CSIZE / (double)G.number_of_nodes();
	*/
	/*
	double w_sum = 0;
	double rw_sum = 0;
	forall_edges(e, G) {
		  w_sum+=G[e].w;
		  rw_sum+=G[e].rw;
		  
		  u = G.source(e);
		  v = G.target(e);
		  G[u].w_sum+=G[e].w;
		  G[u].rw_sum+=G[e].rw;
		  G[v].w_sum+=G[e].w;
		  G[v].rw_sum+=G[e].rw;
  }
  */
	//int wc = 0; int rwc = 0;
	/*
	double tv = 0;
	forall_nodes(v, G) {
		  tv+=G[v].w;
		  //if(G[v].w_sum/w_sum > G[v].rw_sum/rw_sum)
		  //wc++;
		  //else
		  //			rwc++;
  }
  */
	//cout << "wc = " << wc << ", rwc = " << rwc << endl;
	//
	cerr << "alg_calc_iw_by_rw: Interpolation weights: sort by algebraic distance" << endl;
	
	// Recalculation od interpolation order
	int interpolation_order = interpol_order(G);
	int cc=0;
	int N = G.number_of_nodes();
	
	forall_edges(e, G)
	  G[e].iw = 0;
	
	forall_nodes(v, G) {
	  if (G[v].status!=seed) {
		
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if ((G[u].status==seed))//&&(G[u].initial_id<401))
		  tmp_ec_array.push_back(e);
		}
		
		
		list_item it;
		
		TMP_CMP_GRAPHC = &G;
		//if(wc>rwc)
		//tmp_ec_array.sort(&cmp_edges);
		//else
		tmp_ec_array.sort(&cmp_edges_revrw);
		//tmp_ec_array.sort(&cmp_edges_approx_w1_rw2);
		//tmp_ec_array.sort(&cmp_edges_approx_rw1_w2);
		/*
		forall_items(it, tmp_ec_array)
		{
		  cerr << G[tmp_ec_array[it]].rw << " ";
	  }
	  cerr << endl;
	  */
		/*
		double half_size = (double)tmp_ec_array.size()*0.5;
		
		
		if((int)half_size>0)
		{
		  //list<edge> tmp_ec_array2;
		  int c1 = 0; int c2 = 0;
		  
		  while(G[tmp_ec_array[tmp_ec_array.first()]].rw> 2 * G[tmp_ec_array[tmp_ec_array.last()]].rw)
		  {
			tmp_ec_array.erase(tmp_ec_array.last());
			c2++;
	  }
	  
	  if(c2<1) {
		  
		  while((tmp_ec_array.size()>half_size)&&(tmp_ec_array.size()>interpolation_order))
		  {
			tmp_ec_array.erase(tmp_ec_array.last());
			//	c1++;
	  }
	  }
	  
	  //cerr << "half size " << c1 << "\t i potom eshe " << c2 << endl;
	  //cerr << "half weigh " << c2 << "\t i potom eshe " << c1 << endl;
	  
	  TMP_CMP_GRAPHC = &G;
	  tmp_ec_array.sort(&cmp_edges);
	  }
	  */
		
		/*
		forall_items(it, tmp_ec_array) {
		  u = second_adj_for_edge(tmp_ec_array[it], v, G);
		  
		  if((G[u].w > AV_CSIZE*IMB_FACTOR)) {
		  tmp_ec_array.erase(it);
	  }
	  }
	  if(tmp_ec_array.size()==0) {
		  Cnode v_tmp;
		  v_tmp.ArrId = G[v].ArrId;
		  v_tmp.G_prev_ptr= v;
		  v_tmp.status = fine;
		  v_tmp.M = G[v].M;
		  v_tmp.initial_id = G[v].initial_id;
		  node n = H.new_node(v_tmp);
		  G[v].status=seed;
		  G[v].ptr_to_coarse = n;
		  
		  
	  }
	  else {	*/
		//cerr  << "init " << tmp_ec_array.size() << " ";
		//TMP_CMP_GRAPHC = &G;
		//tmp_ec_array.sort(&cmp_edges);
		//tmp_ec_array.sort(&cmp_edges_revrw);
		
		while(tmp_ec_array.size()>interpolation_order) {
		  tmp_ec_array.erase(tmp_ec_array.last());
		  
		}
		
		//TMP_CMP_GRAPHC = &G;
		//tmp_ec_array.sort(&cmp_edges);
		//tmp_ec_array.sort(&cmp_edges_revrw);
		
		//cerr  << " io " << tmp_ec_array.size() << " ";		
		//int xxx = tmp_ec_array.size();
		
		while((tmp_ec_array.size()>1)&&(G[tmp_ec_array[tmp_ec_array.last()]].rw < .5 * G[tmp_ec_array[tmp_ec_array.first()]].rw)) {
		  tmp_ec_array.erase(tmp_ec_array.last());
		}
		
		/*
		while((tmp_ec_array.size()>1)&&
		  (G[second_adj_for_edge(tmp_ec_array[tmp_ec_array.last()], v, G)].w >= 1.5*tv/(double)G.number_of_nodes())) {
		  //(G[tmp_ec_array[tmp_ec_array.last()]].rw < .2 * G[tmp_ec_array[tmp_ec_array.first()]].rw)) {
		  cout << "overloaded "  << G[v].w << " -> " << G[second_adj_for_edge(tmp_ec_array[tmp_ec_array.last()], v, G)].w  << endl;
		tmp_ec_array.erase(tmp_ec_array.last());
		
	  }
	  */
		//cerr  << " strong " << tmp_ec_array.size() << " ";
		//if(xxx!=tmp_ec_array.size() ) cerr << " ----------\n"; 
		//else cerr << "\n";
		//int ss = tmp_ec_array.size();
		//cout << ss << " "  << G[tmp_ec_array[tmp_ec_array.first()]].rw << " " << G[tmp_ec_array[tmp_ec_array.last()]].rw << endl;
		//int ss=0;
		//			while((tmp_ec_array.size()>1)&&(G[tmp_ec_array[tmp_ec_array.last()]].rw < pivot_rw)) {
		  //			tmp_ec_array.erase(tmp_ec_array.last());
		  //		ss++;
		  //}
		  //if(ss > 0)
		  //	cout  << ss << " moves " << tmp_ec_array.size() << endl;
		  /*
		  if((tmp_ec_array.size() >= 7)&&(H.number_of_nodes()>Params.external_part_size)){
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
		  cerr << "added " << endl;
		  two_tuple<node, double> t(v, 1.0);
		  G[v].iw_links.push_back(t);
		  
	} else*/
		  
		  {			
			double edges_to_seed_sum = 0;
			node s, t;
			
			forall_items(it, tmp_ec_array)
			{
			  //cout << "taking " << G[tmp_ec_array[it]].w << endl;
			  //edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
			  edges_to_seed_sum+=G[tmp_ec_array[it]].w;
			}
			
			forall_items(it, tmp_ec_array) {
			  //			if(v_H > .1)
			  //				G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;
			  //			else
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
			  //}
		  }
	}
	
  } else {
	two_tuple<node, double> t(v, 1.0);
	G[v].iw_links.push_back(t);
  }
  }
  
  }
  
  void alg_calc_iw_by_rw_balanced(TGraphC & G, TGraphC & H) {
	
	edge e;
	node v, u;
	
	cerr << "Balanced interpolation weights: sort by algebraic distance" << endl;
	
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
		
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if (G[u].status==seed)
			tmp_ec_array.push_back(e);
		}
		
		
		list_item it;
		
		TMP_CMP_GRAPHC = &G;
		//if(wc>rwc)
		//tmp_ec_array.sort(&cmp_edges);
		//else
		tmp_ec_array.sort(&cmp_edges_revrw);
		
		while(tmp_ec_array.size()>interpolation_order) {
		  tmp_ec_array.erase(tmp_ec_array.last());				
		}
		
		//TMP_CMP_GRAPHC = &G;
		//tmp_ec_array.sort(&cmp_edges);
		//tmp_ec_array.sort(&cmp_edges_revrw);
		//cout << tmp_ec_array.size();
		bool was_attached = false;
		/*
		if((was_attached==false)&&(tmp_ec_array.size()==0)) {
		  cerr << "opsa " << endl;
		  exit(1);
	  }*/
		if(Params.CURRENT_LEVEL<99) {
		  edge best_e1, best_e2;
		  double best_cost = -1;
		  
		  list_item eit1, eit2;
		  forall_items(eit1, tmp_ec_array) {
			forall_items(eit2, tmp_ec_array) {
			  if(eit1<eit2) {
				edge e1 = tmp_ec_array[eit1];
				edge e2 = tmp_ec_array[eit2];
				double cost = G[e1].rw + G[e2].rw;
				node u1 = second_adj_for_edge( e1, v, G);				
				node u2 = second_adj_for_edge( e2, v, G);				
				if((G[u1].wd2 + G[v].w*G[e1].rw/cost <= AVG_AGG_SIZE)&&(G[u2].wd2 + G[v].w*G[e2].rw/cost <= AVG_AGG_SIZE)&&(cost > best_cost)) {
				  best_e1 = e1;
				  best_e2 = e2;
				  best_cost = cost;
				}
			  }
			}
		  }
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
		
		list<edge> tmp_ec_array2 = tmp_ec_array;
		while((was_attached==false)&&(tmp_ec_array.size()>0)) {				
		  u = second_adj_for_edge( tmp_ec_array[tmp_ec_array.first()], v, G);				
		  if(G[u].wd2 + G[v].w <= AVG_AGG_SIZE) {
			while(tmp_ec_array.size()>1)
			  tmp_ec_array.erase(tmp_ec_array.last());
			was_attached = true;
			G[u].wd2 += G[v].w;
		  }
		  else {
			if((tmp_ec_array.size()==1)&&(H.number_of_nodes() > 0.85 * G.number_of_nodes())) {
			  u = second_adj_for_edge( tmp_ec_array2[tmp_ec_array2.first()], v, G);				
			  was_attached = true;
			  G[u].wd2 += G[v].w;				
			  break;
			}
			tmp_ec_array.erase(tmp_ec_array.first());
		  }
		}
		//cerr << tmp_ec_array.size() << endl;
		
		/*			
		while((was_attached==false)&&(tmp_ec_array.size()>0)) {				
		  u = second_adj_for_edge( tmp_ec_array[tmp_ec_array.first()], v, G);				
		  if(G[u].wd2 + G[v].w <= AVG_AGG_SIZE) {
		  while(tmp_ec_array.size()>1)
		  tmp_ec_array.erase(tmp_ec_array.last());
		  was_attached = true;
		  G[u].wd2 += G[v].w;
	  }
	  else {
		if((tmp_ec_array.size()==1)&&(H.number_of_nodes() > 0.7 * G.number_of_nodes())) {
		  was_attached = true;
		  G[u].wd2 += G[v].w;				
		  break;
	  }
	  tmp_ec_array.erase(tmp_ec_array.first());
	  }
	  }
	  */	
		
		
		//cout << " " << tmp_ec_array.size() << endl;
		//			while((tmp_ec_array.size()>1)&&(G[tmp_ec_array[tmp_ec_array.last()]].rw < .5 * G[tmp_ec_array[tmp_ec_array.first()]].rw)) {
		  //				tmp_ec_array.erase(tmp_ec_array.last());
		  //			}
		  
		  //cerr  << " strong " << tmp_ec_array.size() << " ";
		  //if(xxx!=tmp_ec_array.size() ) cerr << " ----------\n"; 
		  //else cerr << "\n";
		  //int ss = tmp_ec_array.size();
		  //cout << ss << " "  << G[tmp_ec_array[tmp_ec_array.first()]].rw << " " << G[tmp_ec_array[tmp_ec_array.last()]].rw << endl;
		  //int ss=0;
		  //			while((tmp_ec_array.size()>1)&&(G[tmp_ec_array[tmp_ec_array.last()]].rw < pivot_rw)) {
		  //			tmp_ec_array.erase(tmp_ec_array.last());
		  //		ss++;
		  //}
		  //if(ss > 0)
		  //	cout  << ss << " moves " << tmp_ec_array.size() << endl;
		  if(was_attached==false) {
			//if((tmp_ec_array.size() >= 7)&&(H.number_of_nodes()>Params.external_part_size)){
			  
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
			double edges_to_seed_sum = 0;
			node s, t;
			
			forall_items(it, tmp_ec_array)
			{
			  //cout << "taking " << G[tmp_ec_array[it]].w << endl;
			  //edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
			  edges_to_seed_sum+=G[tmp_ec_array[it]].rw;
			}
			/*
			if(edges_to_seed_sum==0) {
			  cout << ">>> " << endl;
			  forall_items(it, tmp_ec_array)
			  {
				//cout << "taking " << G[tmp_ec_array[it]].pw << endl;
				cout << "taking " << G[tmp_ec_array[it]].w << endl;
		  }
		  cout << "--------\n";
		  }
		  */
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
	G[v].wd2 = 0;
  }
  cerr << "Interpolation weights: added " << added_fine_nodes << " fine nodes; AVG_AGG_SIZE " << AVG_AGG_SIZE  << endl;
  }
  
  
  
  void clear_f2c_links(TGraphC & G, TGraphC & H, node & v) {
	edge e;
	node w, u;
	
	two_tuple<node, double> t;
	list_item it, it2;
	
	forall_items(it, G[v].iw_links) {
	  w = G[v].iw_links[it].first();
	  //cerr << G[w].iw_links.size() << " ";
	  if(G[w].status!=seed) {
		cerr << "clear_f2c_links: error" << endl;
		exit(-1);
	  }
	  else {
		it2 = G[w].iw_links.first();
		//u = G[w].iw_links[it2].first();
		while((u = G[w].iw_links[it2].first()) != v) {
		  it2 = G[w].iw_links.succ(it2);
		  
		}
		G[w].iw_links.del_item(it2);
	  }
	  //cerr << G[w].iw_links.size() << " ";
	}
	G[v].iw_links.clear();
  }
  
void recalculate_vertex_fvol(TGraphC & G, TGraphC & H) {
	node v, w, u;
	edge e;
	
	forall_nodes(v, G) {
	  G[v].wd = 0;
	}
	double av = 0;
	double c_av = 0;
	forall_nodes(w, G)
	{
	  double edges_to_fine_sum = 0;
	  forall_adj_edges(e, w)
	  {
		u = second_adj_for_edge(e, w, G);
		if(G[u].status!=seed)
		{
		  edges_to_fine_sum+=G[e].iw*G[u].w;
		}
	  }
	  G[w].wd = G[w].w + edges_to_fine_sum;
	  /*
	  if(G[w].status==seed) {
		av+=G[w].wd;
		c_av++;
	}
	*/
	  
	}
	//cerr << "average wdeg " << av / c_av << endl;
}
  
void refine_alg_calc_iw_by_rw(TGraphC & G, TGraphC & H) {
	
	edge e;
	node v, u;
	
	double AV_CSIZE = 0;
	//double IMB_FACTOR = (1+Params.imbalance) * (1.0 + (S+1)*0.1);
	forall_nodes(v, G)
	  AV_CSIZE+=G[v].w;
	AV_CSIZE = 1 * AV_CSIZE / (double)H.number_of_nodes();//Params.min_graph_size;
	cerr << "max csize " << AV_CSIZE << endl;
	
	cerr << "Refine interpolation weights: sort by algebraic distance" << endl;
	
	recalculate_vertex_fvol(G, H);
	
	// Recalculation od interpolation order
	int interpolation_order = interpol_order(G);
	int cc=0;
	int N = G.number_of_nodes();
	
	forall_edges(e, G)
	  G[e].iw = 0;
	
	forall_nodes(v, G) {
	  if (G[v].status!=seed) {
		
		clear_f2c_links(G, H, v);
		
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if ((G[u].status==seed))//&&(G[u].initial_id<401))
			tmp_ec_array.push_back(e);
		}
		
		
		list_item it, it2;
		
		
		TMP_CMP_GRAPHC = &G;
		tmp_ec_array.sort(&cmp_edges_revrw);
		
		
		while(tmp_ec_array.size()>interpolation_order) {
		  tmp_ec_array.erase(tmp_ec_array.last());
		  
		}
		
		//TMP_CMP_GRAPHC = &G;
		//tmp_ec_array.sort(&cmp_edges);
		//tmp_ec_array.sort(&cmp_edges_revrw);
		
		while((tmp_ec_array.size()>1)&&(G[tmp_ec_array[tmp_ec_array.last()]].rw < .5 * G[tmp_ec_array[tmp_ec_array.first()]].rw)) {
		  tmp_ec_array.erase(tmp_ec_array.last());
		}
		
		
		it = tmp_ec_array.last();
		//cerr << tmp_ec_array.size() << " ";
		while(it!=tmp_ec_array.first()) {
		  it2 = it;
		  it = tmp_ec_array.pred(it);
		  node c_neigh = second_adj_for_edge(tmp_ec_array[it2], v, G);
		  //cerr << AV_CSIZE << " " << G[c_neigh].wd << " " << G[c_neigh].w << " " << G.degree(c_neigh) ;
		  if(G[c_neigh].wd > AV_CSIZE) {
			tmp_ec_array.erase(it2);
			//cerr << "-";
		  }				
		  //cerr << endl;
		}
		//cerr << " " << tmp_ec_array.size() << endl;
		
		
		double edges_to_seed_sum = 0;
		node s, t;
		
		forall_items(it, tmp_ec_array)
		{
		  //edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
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
		  
		  if (((G[s].status==seed)&&(G[t].status==seed))||((G[s].status!=seed)&&(G[t].status!=seed))) {
			cerr << "jopa" << endl;
			exit(1);
		  }
		}
		
		//cerr << G[v].iw_links.size() << endl;
		
	  } else {
		two_tuple<node, double> t(v, 1.0);
		G[v].iw_links.push_back(t);
	  }
	}
	
}
  
  
void minla_alg_calc_iw_by_rw(TGraphC & G) {
	edge e;
	node v, u;
	
	cerr << "minla_alg_calc_iw_by_rw: Interpolation weights: sort by algebraic distance" << endl;
	
	// Recalculation od interpolation order
	int interpolation_order = interpol_order(G);
	int cc=0;
	int N = G.number_of_nodes();
	
	forall_edges(e, G)
	  G[e].iw = 0;
	
	forall_nodes(v, G) {
	  if (G[v].status!=seed) {
		
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if ((G[u].status==seed))//&&(G[u].initial_id<401))
			tmp_ec_array.push_back(e);
		}
		
		
		list_item it;
		
		TMP_CMP_GRAPHC = &G;
		//	tmp_ec_array.sort(&cmp_edges);
		tmp_ec_array.sort(&cmp_edges_revrw);
		/*
		double half_size = (double)tmp_ec_array.size()*0.5;
		
		if((int)half_size>0)
		  {
			//list<edge> tmp_ec_array2;
			int c1 = 0; int c2 = 0;
			
			while(G[tmp_ec_array[tmp_ec_array.first()]].rw> 2 * G[tmp_ec_array[tmp_ec_array.last()]].rw)
		  {
			tmp_ec_array.erase(tmp_ec_array.last());
			c2++;
	  }
	  
	  if(c2<1) {
			
			while((tmp_ec_array.size()>half_size)&&(tmp_ec_array.size()>interpolation_order))
		  {
			tmp_ec_array.erase(tmp_ec_array.last());
			//	c1++;
	  }
	  }
	  
	  //cerr << "half size " << c1 << "\t i potom eshe " << c2 << endl;
	  //cerr << "half weigh " << c2 << "\t i potom eshe " << c1 << endl;
	  TMP_CMP_GRAPHC = &G;
	  tmp_ec_array.sort(&cmp_edges);
	  }
	  */
		
		
		while (tmp_ec_array.size()>interpolation_order/*Params.get_maximum_of_connections_between_fine_and_its_coarses(N)*/)
		  tmp_ec_array.erase(tmp_ec_array.last());
		double edges_to_seed_sum = 0;
		node s, t;
		
		forall_items(it, tmp_ec_array)
		{
		  //cout << "taking " << G[tmp_ec_array[it]].w << endl;
		  //edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
		  edges_to_seed_sum+=G[tmp_ec_array[it]].w;
		}
		
		forall_items(it, tmp_ec_array) {
		  //			if(v_H > .1)
		  //				G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;
		  //			else
		  G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].w/edges_to_seed_sum;
		  //G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw/edges_to_seed_sum;
		  
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
	  } else {
		two_tuple<node, double> t(v, 1.0);
		G[v].iw_links.push_back(t);
	  }
	}
	
}
  
  
  
void alg_calc_iw_by_rw_and_w(TGraphC & G) {
	
	cerr << "Interpolation weights: sort by algebraic distance and edge weights" << endl;
	
	// Recalculation od interpolation order
	int interpolation_order = interpol_order(G);
	int cc=0;
	
	edge e;
	node v, u;
	int N = G.number_of_nodes();
	
	forall_edges(e, G)
	  G[e].iw = 0;
	
	forall_nodes(v, G) {
	  if (G[v].status!=seed) {
		
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if ((G[u].status==seed))//&&(G[u].initial_id<401))
			tmp_ec_array.push_back(e);
		}
		
		//			TMP_CMP_GRAPHC = &G;
		//			tmp_ec_array.sort(&cmp_rnd_edges);
		
		list_item it;
		double rw_sum = 0;
		double w_sum = 0;
		forall_items(it, tmp_ec_array) {
		  rw_sum+=G[tmp_ec_array[it]].rw;
		  w_sum+=G[tmp_ec_array[it]].w;//G[tmp_ec_array[it]].nw;
		}
		forall_items(it, tmp_ec_array)
		  G[tmp_ec_array[it]].double_io_score=.5*G[tmp_ec_array[it]].rw/rw_sum+.5*G[tmp_ec_array[it]].w/w_sum;
		
		TMP_CMP_GRAPHC = &G;
		tmp_ec_array.sort(&cmp_rnd_edges);
		tmp_ec_array.sort(&cmp_edges_double_io_score);
		
		/*
		forall_items(it, tmp_ec_array)
		  {
			cout << G[tmp_ec_array[it]].rw << "(" << G[tmp_ec_array[it]].pw << "), " ;
	  }
	  cout << endl;
	  */
		while (tmp_ec_array.size()>interpolation_order/*Params.get_maximum_of_connections_between_fine_and_its_coarses(N)*/)
		  tmp_ec_array.erase(tmp_ec_array.last());
		double edges_to_seed_sum = 0;
		
		forall_items(it, tmp_ec_array) {
		  edges_to_seed_sum+=G[tmp_ec_array[it]].double_io_score;
		}
		
		node s, t;
		forall_items(it, tmp_ec_array) {
		  G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].double_io_score/edges_to_seed_sum;
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
	  } else {
		two_tuple<node, double> t(v, 1.0);
		G[v].iw_links.push_back(t);
	  }
	}
}
  
  
void wag_calc_iw(TGraphC & G) {
	cerr << "Interpolation weights: sort by weighted aggregation (classical)"
	<< endl;
	
	// Recalculation od interpolation order
	int interpolation_order = /*Params.maximum_of_connections_between_fine_and_its_coarses;*/
	interpol_order(G);
	//(int)((double)Params.maximum_of_connections_between_fine_and_its_coarses+log(E0/(double)G.number_of_edges())+0.5);
	//  if(interpolation_order < Params.maximum_of_connections_between_fine_and_its_coarses)
	//    interpolation_order = Params.maximum_of_connections_between_fine_and_its_coarses;
	//  if(interpolation_order > 100)
	//    interpolation_order = 100;
	//  if(S==0)
	//    interpolation_order = Params.maximum_of_connections_between_fine_and_its_coarses;
	//  cerr << "Recalculated interpolation order = " << interpolation_order << endl;
	//////////////////////////////////////////
	
	edge e;
	node v, u;
	int N = G.number_of_nodes();
	
	forall_edges(e, G)
	  G[e].iw = 0;
	
	forall_nodes(v, G) {
	  //G[v].iw_links.clear();
	  if (G[v].status!=seed) {
		list<edge> tmp_ec_array;
		tmp_ec_array.clear();
		forall_adj_edges(e, v) {
		  u = second_adj_for_edge(e, v, G);
		  if (G[u].status==seed)
			tmp_ec_array.push_back(e);
		}
		
		TMP_CMP_GRAPHC = &G;
		tmp_ec_array.sort(&cmp_rnd_edges); // to remove
		tmp_ec_array.sort(&cmp_edges);
		/*
		if(S==0)
		{
		  list_item ppp;
		  forall_items(ppp, tmp_ec_array)
		  {
			e = tmp_ec_array[ppp];
			u = second_adj_for_edge(e, v, G);
			//		  if((G[u].initial_id>=403)&&(tmp_ec_array.size()>1))
			if((G.degree(u)>2000)&&(tmp_ec_array.size()>1))
			tmp_ec_array.erase(ppp);
	  }
	  }
	  */
		while (tmp_ec_array.size()>interpolation_order/*Params.get_maximum_of_connections_between_fine_and_its_coarses(N)*/)
		  tmp_ec_array.erase(tmp_ec_array.last());
		
		//while (G[tmp_ec_array[tmp_ec_array.first()]].pw/G[tmp_ec_array[tmp_ec_array.last()]].pw>100)
		  while (G[tmp_ec_array[tmp_ec_array.first()]].w/G[tmp_ec_array[tmp_ec_array.last()]].w>100)
			// changed from .w - to check!
			tmp_ec_array.erase(tmp_ec_array.last());
		  
		  list_item it;
		  double edges_to_seed_sum = 0;
		  forall_items(it, tmp_ec_array)
			edges_to_seed_sum+=G[tmp_ec_array[it]].w;
		  //edges_to_seed_sum+=G[tmp_ec_array[it]].pw;
		  node s, t;
		  forall_items(it, tmp_ec_array) {
			//G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].pw
			G[tmp_ec_array[it]].iw = G[tmp_ec_array[it]].w/edges_to_seed_sum;
			s = G.source(tmp_ec_array[it]);
			t = G.target(tmp_ec_array[it]);
			
			two_tuple<node, double> t1(s, G[tmp_ec_array[it]].iw);
			two_tuple<node, double> t2(t, G[tmp_ec_array[it]].iw);
			//	      two_tuple<node, double> t1(s,  G[tmp_ec_array[it]].rw);
			//              two_tuple<node, double> t2(t,  G[tmp_ec_array[it]].rw);
			
			G[s].iw_links.push_back(t2);
			G[t].iw_links.push_back(t1);
			/*
			if((G[v].initial_id==311)||(G[v].initial_id==167))
		  cerr << G[s].initial_id << "-" <<  G[t].initial_id << " iw=" << G[tmp_ec_array[it]].iw << " w=" << G[tmp_ec_array[it]].w << " rw=" << G[tmp_ec_array[it]].rw << endl;
		  */
			
			if ((G[s].status==seed)&&(G[t].status==seed)) {
			  cerr << "jopa" << endl;
			  exit(1);
			}
			
			if ((G[s].status!=seed)&&(G[t].status!=seed)) {
			  cerr << "jopa" << endl;
			  exit(1);
			}
			
			//              G[s].iw_links.push_back(tmp_ec_array[it]);
			//              G[t].iw_links.push_back(tmp_ec_array[it]);
		  }
	  } else {
		two_tuple<node, double> t(v, 1.0);
		G[v].iw_links.push_back(t);
	  }
	  //      cerr << G.degree(v) << ", iw links len=" << G[v].iw_links.size() << endl;
	}
}
  
