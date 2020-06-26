#include "externs.h"
#include <LEDA/graph/node_partition.h>
#include <LEDA/core/array.h>
#include <LEDA/system/stream.h>
#include <LEDA/core/tuple.h> 

static int cmp_two_tuple(const two_tuple<edge, node> & n1, const two_tuple<edge, node> & n2)
{
  if((*TMP_CMP_GRAPHC)[n1.second()].interpol_coord>(*TMP_CMP_GRAPHC)[n2.second()].interpol_coord)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1.second()].interpol_coord<(*TMP_CMP_GRAPHC)[n2.second()].interpol_coord)
    return -1;
  else
    return 0;
}

static int cmp(const edge & e1, const edge & e2)
{
  if((*TMP_CMP_GRAPH)[e1]<(*TMP_CMP_GRAPH)[e2])
    return 1;
  else if((*TMP_CMP_GRAPH)[e1]>(*TMP_CMP_GRAPH)[e2])
    return -1;
  else
    return 0;
}

static int cmp_keys(const leda::string & s1, const leda::string & s2)
{
  int n1=0; int n2=0;
  for(int i=0; i<s1.length(); i++)
    {
      if(s1[i]=='1') n1++;
      if(s2[i]=='1') n2++;  
    }
  
  if(n1>n2)
    return 1;
  else if(n1<n2)
    return -1;
  else
    return 0;
}

static int cmp_Pair_LOfInt_Double(const Pair_LOfInt_Double & e1, const Pair_LOfInt_Double & e2)
{
  if(e1.c<e2.c)
    return -1;
  else if(e1.c>e2.c)
    return 1;
  else
    return 0;
}
static int cmp_Pair_LOfInt_Double_by_diff_mark(const Pair_LOfInt_Double & e1, const Pair_LOfInt_Double & e2)
{
  if(e1.diff_mark<e2.diff_mark)
    return -1;
  else if(e1.diff_mark>e2.diff_mark)
    return 1;
  else
    return 0;
}

static int cmp_LA(const LA & la1, const LA & la2)
{
  if(la1.weight<la2.weight)
    return -1;
  else if(la1.weight>la2.weight)
    return 1;
  else
    return 0;
}
/*
static int cmp_tmp_double_one_scope(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].tmp_double_one_scope>(*TMP_CMP_GRAPHC)[n2].tmp_double_one_scope)
    return -1;
  else if((*TMP_CMP_GRAPHC)[n1].tmp_double_one_scope<(*TMP_CMP_GRAPHC)[n2].tmp_double_one_scope)
    return 1;
  else
    return 0;
}
*/


static int cmp_triple(const triple_res & la1, const triple_res & la2)
{
  if(la1.cost<la2.cost)
    return -1;
  else if(la1.cost>la2.cost)
    return 1;
  else
    return 0;
}

static int cmp_ArrId(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].ArrId>(*TMP_CMP_GRAPHC)[n2].ArrId)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].ArrId<(*TMP_CMP_GRAPHC)[n2].ArrId)
    return -1;
  else
    return 0;
}

static int cmp_color(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].col>(*TMP_CMP_GRAPHC)[n2].col)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].col<(*TMP_CMP_GRAPHC)[n2].col)
    return -1;
  else
    return 0;
}
/*
static int cmp_fiedler(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].Fiedler_coord>(*TMP_CMP_GRAPHC)[n2].Fiedler_coord)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].Fiedler_coord<(*TMP_CMP_GRAPHC)[n2].Fiedler_coord)
    return -1;
  else
    return 0;
}
*/
static int cmp_e_w(const edge & n1, const edge & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].w>(*TMP_CMP_GRAPHC)[n2].w)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].w<(*TMP_CMP_GRAPHC)[n2].w)
    return -1;
  else
    return 0;
}

static int cmp_dArrId(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].dArrId>(*TMP_CMP_GRAPHC)[n2].dArrId)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].dArrId<(*TMP_CMP_GRAPHC)[n2].dArrId)
    return -1;
  else
    return 0;
}

static int cmp_ew(const edge & e1, const edge & e2)
{
  if((*TMP_CMP_GRAPHC)[e1].w>(*TMP_CMP_GRAPHC)[e2].w)
    return -1;
  else if((*TMP_CMP_GRAPHC)[e1].w<(*TMP_CMP_GRAPHC)[e2].w)
    return 1;
  else
    return 0;
}

static int cmp_rnd(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].M>(*TMP_CMP_GRAPHC)[n2].M)
    return 1;
  else
    
    return -1;
  /*
  float a = rand(); float b = rand();
  
  if(a/b>0.5)
    return 1;
  else
    return -1;
  */
}
static int cmp_rnd_edges(const edge & n1, const edge & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].M>(*TMP_CMP_GRAPHC)[n2].M)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].M<(*TMP_CMP_GRAPHC)[n2].M)
    
    return -1;
  else
    return 0;
}

static int cmp_lcc_tuples_by_length(const two_tuple<node, node> & t1, const  two_tuple<node, node> & t2)
{
  int l1 = (*TMP_CMP_GRAPHC)[t1.second()].ArrId - (*TMP_CMP_GRAPHC)[t1.first()].ArrId;
  int l2 = (*TMP_CMP_GRAPHC)[t2.second()].ArrId - (*TMP_CMP_GRAPHC)[t2.first()].ArrId;
  if(l1>l2)
    return 1;
  else if(l1<l2)
    return -1;
  else
    return 0;
}

static int cmp_lcc_ArrId(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].lcc_ArrId>(*TMP_CMP_GRAPHC)[n2].lcc_ArrId)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].lcc_ArrId<(*TMP_CMP_GRAPHC)[n2].lcc_ArrId)
    return -1;
  else
    return 0;
}
/*
static int cmp_save_bestLCCArrID(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].save_lcc_ArrId>(*TMP_CMP_GRAPHC)[n2].save_lcc_ArrId)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].save_lcc_ArrId<(*TMP_CMP_GRAPHC)[n2].save_lcc_ArrId)
    return -1;
  else
    return 0;
}
*/
static int cmp_lcc_A(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].lcc_A>(*TMP_CMP_GRAPHC)[n2].lcc_A)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].lcc_A<(*TMP_CMP_GRAPHC)[n2].lcc_A)
    return -1;
  else
    return 0;
}

static int cmp_initial_id(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].initial_id>(*TMP_CMP_GRAPHC)[n2].initial_id)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].initial_id<(*TMP_CMP_GRAPHC)[n2].initial_id)
    return -1;
  else
    return 0;
}


static int cmp_Sval(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].S_value>(*TMP_CMP_GRAPHC)[n2].S_value)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].S_value<(*TMP_CMP_GRAPHC)[n2].S_value)
    return -1;
  else
    return 0;
}

static int cmp_node_w(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].w<(*TMP_CMP_GRAPHC)[n2].w)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].w>(*TMP_CMP_GRAPHC)[n2].w)
    return -1;
  else
    return 0;
}
/*
static int cmp_nodes_rval(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].rval<(*TMP_CMP_GRAPHC)[n2].rval)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].rval>(*TMP_CMP_GRAPHC)[n2].rval)
    return -1;
  else
    return 0;
}
*/
/*
static int cmp_ev2(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].ev2>(*TMP_CMP_GRAPHC)[n2].ev2)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].ev2<(*TMP_CMP_GRAPHC)[n2].ev2)
    return -1;
  else
    return 0;
}
*/
static int cmp_interpol_coord(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].interpol_coord>(*TMP_CMP_GRAPHC)[n2].interpol_coord)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].interpol_coord<(*TMP_CMP_GRAPHC)[n2].interpol_coord)
    return -1;
  else
    return 0;
}

static int cmp_coord_C(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].interpol_coord>(*TMP_CMP_GRAPHC)[n2].interpol_coord)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].interpol_coord<(*TMP_CMP_GRAPHC)[n2].interpol_coord)
    return -1;
  else
    return 0;
}

static int cmp_weighted_degree_C(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].wd>(*TMP_CMP_GRAPHC)[n2].wd)
    return -1;
  else if((*TMP_CMP_GRAPHC)[n1].wd<(*TMP_CMP_GRAPHC)[n2].wd)
    return 1;
  else
    return 0;
}

static int cmp_weighted_degree_rev(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].wd<(*TMP_CMP_GRAPHC)[n2].wd)
    return -1;
  else if((*TMP_CMP_GRAPHC)[n1].wd>(*TMP_CMP_GRAPHC)[n2].wd)
    return 1;
  else
    return 0;
}


static int cmp_imaginary_vertex_move(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].adj_diff_in_coords>(*TMP_CMP_GRAPHC)[n2].adj_diff_in_coords)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].adj_diff_in_coords<(*TMP_CMP_GRAPHC)[n2].adj_diff_in_coords)
    return -1;
  else
    return 0;
}
/*
static int cmp_nCFloat(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].ArrIdFloat>(*TMP_CMP_GRAPHC)[n2].ArrIdFloat)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].ArrIdFloat<(*TMP_CMP_GRAPHC)[n2].ArrIdFloat)
    return -1;
  else
    return 0;
}
*/
static int cmp_nodes_degree(const node & n1, const node & n2)
{
  if((*TMP_CMP_GRAPHC).degree(n1)<(*TMP_CMP_GRAPHC).degree(n2))
    return 1;
  else if((*TMP_CMP_GRAPHC).degree(n1)>(*TMP_CMP_GRAPHC).degree(n2))
    return -1;
  else
    return 0;
}

static int cmp_edges(const tmp_edge_cost & n1, const tmp_edge_cost & n2)
{
  if(n1.w<n2.w)
    return 1;
  else if(n1.w>n2.w)
    return -1;
  else
    return 0;
}

static int cmp_edges_w1_rw2(const edge & n1, const edge & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].w<(*TMP_CMP_GRAPHC)[n2].w)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].w>(*TMP_CMP_GRAPHC)[n2].w)
    return -1;
  else
  {
	  if((*TMP_CMP_GRAPHC)[n1].rw<(*TMP_CMP_GRAPHC)[n2].rw)
    		return 1;
  	  else if((*TMP_CMP_GRAPHC)[n1].rw>(*TMP_CMP_GRAPHC)[n2].rw)
    		return -1;
  	  else
    		return 0;
  }
}

static int cmp_edges_approx_w1_rw2(const edge & n1, const edge & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].w*2<(*TMP_CMP_GRAPHC)[n2].w)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].w>(*TMP_CMP_GRAPHC)[n2].w*2)
    return -1;
  else
  {
	  if((*TMP_CMP_GRAPHC)[n1].rw<(*TMP_CMP_GRAPHC)[n2].rw)
    		return 1;
  	  else if((*TMP_CMP_GRAPHC)[n1].rw>(*TMP_CMP_GRAPHC)[n2].rw)
    		return -1;
  	  else
    		return 0;
  }
}

static int cmp_edges_approx_rw1_w2(const edge & n1, const edge & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].rw*2<(*TMP_CMP_GRAPHC)[n2].rw)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].rw>(*TMP_CMP_GRAPHC)[n2].rw*2)
    return -1;
  else
  {
	  if((*TMP_CMP_GRAPHC)[n1].w<(*TMP_CMP_GRAPHC)[n2].w)
    		return 1;
  	  else if((*TMP_CMP_GRAPHC)[n1].w>(*TMP_CMP_GRAPHC)[n2].w)
    		return -1;
  	  else
    		return 0;
  }
}


static int cmp_erw(const edge & n1, const edge & n2)
{
  if((*TMP_CMP_GRAPHC)[n1].rw<(*TMP_CMP_GRAPHC)[n2].rw)
    return 1;
  else if((*TMP_CMP_GRAPHC)[n1].rw>(*TMP_CMP_GRAPHC)[n2].rw)
    return -1;
  else
    return 0;
}



static int cmp_edges(const edge & e1, const edge & e2)
{
  if((*TMP_CMP_GRAPHC)[e1].w<(*TMP_CMP_GRAPHC)[e2].w)
    return 1;
  else if((*TMP_CMP_GRAPHC)[e1].w>(*TMP_CMP_GRAPHC)[e2].w)
    return -1;
  else
    return 0;
}
static int cmp_edges_rw(const edge & e1, const edge & e2)
{
  if((*TMP_CMP_GRAPHC)[e1].rw<(*TMP_CMP_GRAPHC)[e2].rw)
    return -1;
  else if((*TMP_CMP_GRAPHC)[e1].rw>(*TMP_CMP_GRAPHC)[e2].rw)
    return 1;
  else
    return 0;
}
/*
static int cmp_edges_bamg_iw(const edge & e1, const edge & e2)
{
  if((*TMP_CMP_GRAPHC)[e1].bamg_iw<(*TMP_CMP_GRAPHC)[e2].bamg_iw)
    return 1;
  else if((*TMP_CMP_GRAPHC)[e1].bamg_iw>(*TMP_CMP_GRAPHC)[e2].bamg_iw)
    return -1;
  else
    return 0;
}
static int cmp_edges_rev_bamg_iw(const edge & e1, const edge & e2)
{
  if((*TMP_CMP_GRAPHC)[e1].bamg_iw<(*TMP_CMP_GRAPHC)[e2].bamg_iw)
    return -1;
  else if((*TMP_CMP_GRAPHC)[e1].bamg_iw>(*TMP_CMP_GRAPHC)[e2].bamg_iw)
    return 1;
  else
    return 0;
}
*/
static int cmp_edges_rw2(const edge & e1, const edge & e2)
{
  node u1, v1, u2, v2, u, v;
  u1 = (*TMP_CMP_GRAPHC).source(e1);
  v1 = (*TMP_CMP_GRAPHC).target(e1);
  u2 = (*TMP_CMP_GRAPHC).source(e2);
  v2 = (*TMP_CMP_GRAPHC).target(e2);

  if(u1==u2)
     u = u1; 
  else if(u1==v2)
    { u = u1; v2 = u2; }
  else if(v1==v2)
    { u = v1; v1=u2; v2 = u2; }
  else if(v1==u2)
    { u = v1; v1 = u1; v2 = u1; }

  if((*TMP_CMP_GRAPHC)[v1].future_edges_sum<(*TMP_CMP_GRAPHC)[v2].future_edges_sum)
    return 1;
  else if((*TMP_CMP_GRAPHC)[v1].future_edges_sum>(*TMP_CMP_GRAPHC)[v2].future_edges_sum)
    return -1;
  else
    return 0;
}
static int cmp_edges_revrw(const edge & e1, const edge & e2) // reverse rw
{
  if((*TMP_CMP_GRAPHC)[e1].rw<(*TMP_CMP_GRAPHC)[e2].rw)
    return 1;
  else if((*TMP_CMP_GRAPHC)[e1].rw>(*TMP_CMP_GRAPHC)[e2].rw)
    return -1;
/*
	else if((*TMP_CMP_GRAPHC)[e1].H<(*TMP_CMP_GRAPHC)[e2].H)
    return -1;
	else if((*TMP_CMP_GRAPHC)[e1].H>=(*TMP_CMP_GRAPHC)[e2].H)
    return 1;
*/
  else
    return 0;
  
}
/*
static int cmp_edges_w_rw(const edge & e1, const edge & e2) // reverse rw
{
  if((*TMP_CMP_GRAPHC)[e1].w_rw<(*TMP_CMP_GRAPHC)[e2].w_rw)
    return 1;
  else if((*TMP_CMP_GRAPHC)[e1].w_rw>(*TMP_CMP_GRAPHC)[e2].w_rw)
    return -1;
  else
    return 0;
  
}
*/
static int cmp_edges_w_rw_new(const edge & e1, const edge & e2) // reverse rw
{
	double wrw1 = (*TMP_CMP_GRAPHC)[e1].w * (*TMP_CMP_GRAPHC)[e1].rw;
	double wrw2 = (*TMP_CMP_GRAPHC)[e2].w * (*TMP_CMP_GRAPHC)[e2].rw;
  if(wrw1/wrw2>1)
    return -1;
  else if(wrw1/wrw2<1)
    return 1;
  else
  	return 0;
}

static int cmp_edges_io_score(const edge & e1, const edge & e2) // reverse rw
{
  if((*TMP_CMP_GRAPHC)[e1].io_score<(*TMP_CMP_GRAPHC)[e2].io_score)
    return -1;
  else if((*TMP_CMP_GRAPHC)[e1].io_score>(*TMP_CMP_GRAPHC)[e2].io_score)
    return 1;
  else
    return 0;
  
}
static int cmp_edges_double_io_score(const edge & e1, const edge & e2) // reverse rw
{
  if((*TMP_CMP_GRAPHC)[e1].double_io_score<(*TMP_CMP_GRAPHC)[e2].double_io_score)
    return 1;
  else if((*TMP_CMP_GRAPHC)[e1].double_io_score>(*TMP_CMP_GRAPHC)[e2].double_io_score)
    return -1;
  else
    return 0;
  
}


static int cmp_elen(const edge & e1, const edge & e2)
{
  node x1 = (*TMP_CMP_GRAPHC).source(e1); node y1 = (*TMP_CMP_GRAPHC).target(e1);
  node x2 = (*TMP_CMP_GRAPHC).source(e2); node y2 = (*TMP_CMP_GRAPHC).target(e2);
  double l1 = (*TMP_CMP_GRAPHC)[e1].w * fabs((*TMP_CMP_GRAPHC)[x1].S_value - (*TMP_CMP_GRAPHC)[y1].S_value);
  double l2 = (*TMP_CMP_GRAPHC)[e2].w * fabs((*TMP_CMP_GRAPHC)[x2].S_value - (*TMP_CMP_GRAPHC)[y2].S_value);
    
  if(l1<l2)
    return 1;
  else if(l1>l2)
    return -1;
  else
    return 0;
}
