#include <cmath>
#include <LEDA/graph/ugraph.h>
#include <LEDA/system/stream.h>
#include <LEDA/graph/node_array.h>
#include <LEDA/graph/edge_array.h>
#include <LEDA/core/array.h>
#include <LEDA/core/tuple.h>

using namespace leda;
using namespace std;

#define xdmin(X, Y)  ((X) < (Y) ? (double)(X) : (double)(Y))
#define xdmax(X, Y)  ((X) > (Y) ? (double)(X) : (double)(Y))

//#define quad(X) ((X)*(X))

#define DOUBLE_DIFF 0.00000000000000001
#define STR_LEN_1K 1024
#define MAX_BLOCK_SIZE 4
//#define MIN_GRAPH_SIZE 8
#define ONLY_FINE true
#define ALL false
#define SIZE 256

#define START_MAX_TEMP_SA 10

#define all_pairs_improve_bound 10

typedef UGRAPH<int, double> TGraph;

enum TStatus {seed = 0, fine = 2}; // green = choosen for a cluster
enum fvol_sort_type {nosort = 0, alg = 1, wd = 2, mix_alg_wag = 3}; 
enum alg_dist_type {jacobian = 0, gauss_seidel = 1, noalgdist = 2, wag_and_gs = 3, maximum_matching = 4, asym_gauss_seidel = 5}; 
enum sym_type {symmetric = 0, nonsymmetric = 1};
enum T_iw_calc {iwbamg = 0, iw_alg = 1, iw_wag = 2, iw_alg_and_wag = 3}; 

class Bandwidth
{
 public:
  double c;

  list<edge> l; 

  friend ostream& operator<<(ostream& o, const Bandwidth& s)
    { 
      return o;
    }

  friend istream& operator>>(istream& i, Bandwidth& s)
    {
      return i;
    }


  Bandwidth()
    {
            c = 0;
            l.clear();
    }
  /*
  friend bool& operator<(const Bandwidth& a, const Bandwidth& b)
    { 
      return a.c<b.c;
    }
  
  friend bool& operator>(const Bandwidth& a, const Bandwidth& b)
    { 
      return a.c>b.c;
    }
  
  friend bool& operator<=(const Bandwidth& a, const Bandwidth& b)
    { 
      return a.c<=b.c;
    }
  
  friend bool& operator>=(const Bandwidth& a, const Bandwidth& b)
    { 
      return a.c>=b.c;
    }
  friend bool& operator==(const Bandwidth& a, const Bandwidth& b)
    { 
      return a.c==b.c;
    }
  */
 
};

class edge_to_node {
 public:
  node v;
  double w;
  
  friend ostream& operator<<(ostream& o, const  edge_to_node& s)
    {
      return o;
    }
  
  friend istream& operator>>(istream& i,  const edge_to_node& s)
    {
      return i;
    }
};

class Cnode
{
 public:
   int sol;
int lemon_id;
   double T_max;
   double contr_edges;
 double wd2;
 int metis_initial_id;
 int part_sol;
 double boundary;
 bool visited;
bool is_he;
 double hedge_rw;

  double bval;
 leda::string name;
 double w_sum;
 double rw_sum;
 //double part_rw;
 double min_v;
 //list<node> wfn;
 // int wfnc;
int longest_redge;

 bool matched;
 int col;
 //int ordNum;
	//double twd;
	//double min_rw;
 double maxrw;

  int number_of_Cneighbors;
  bool bamg_insignificant;

  double all_edges_sum;
  double w_all_edges_sum;

  double allrw;
 // double jac_rv;
//double jac_rv2;
  //bool tflag; // to remove;

  //list<node> wf;
  //array<double> rv;//, prev_rv;
  //array<double> rv2;
  //double Fiedler_coord;
  //double rval;

  //list<edge_to_node> deleted_edges;
  //double total_deleted;
  
  //bool tmp_in_window;
  //node ptr_to_induced_wind;

  
  //double tmp_double_one_scope;
  //bool coarse_zero_degree;
  //bool bfs_visited;
  
  //int window1;
  //int window2;

  //double before_correction;
  //double scaled_before_correction;
  //double correction;
  
  //double ev2;
  //  array<int> ArrIds; // from different basis solutions
  list<two_tuple<node, double> > iw_links;

  //  double lh_pts;
  //  double rh_pts;
  
  double double_dummy;
  double future_edges_sum;
  double w_future_edges_sum;
  double future_edges_toC_sum;
  
  double diff_in_coords;
  double adj_diff_in_coords;
  
  int touched;
  int lcc_depth;

  int seg_min_save_ArrId;
  int seg_min_mainsave_ArrId;
  node was_flipped_with;
  
  int tmp_order;
  //  bool has_big_wd; // more than mid_wd+sd_wd
  double wd;
  //  double weighted_deg; // Sort by weighted degree before choosing coarse nodes
  double w;

  //  list<node> my_fine_nodes;                           // for coarse node
  //  list<int> positions_of_my_fine_nodes;               // for coarse node
  //  list<node> my_coarse_nodes;                         // for fine node
  //  list<int> positions_in_clusters_of_coarse_nodes;    // for fine node
  
  double S_value;
  double save_S_value;
  double interpol_coord;
  // double save_interpol_coord_on_local_window_gs;
  //  double prev_interpol_coord;
  bool captured;
  //  int weight;
  //  double total_edges_weight;
  int save_ArrId;
  int after_interpol_ArrId;
  int ArrId;
  double dArrId;
  //  float ArrIdFloat;
  TStatus status;
  
  //  int ID;
  node G_prev_ptr; // in matching - ptr to previous
  node ptr_to_coarse; // if this node is green then it will be a pointer to its coarse node
  int initial_id;
  int last_adj_num;
  //  list<node> orange_adj; // only for green nodes. stores its orange neighb.

  // LCC variables
  double M;
  //list<node> contr_nodes;
  void save(file_ostream o)
    {
      o << initial_id << endl << ArrId << endl << w << endl;
      if(status==seed)
        o << 1 << endl;
      else
        o << 0 << endl;
    }
  void load(file_istream i)
    {
      int col;
      i >> initial_id;
      i >> ArrId;
      i >> w;
      i >> col;
      if(col==1)
        status=seed;
      else
        status=fine;
    }
    
  double SM;  
  bool       lcc_flag;
  node       lcc_left_adj;
  node       lcc_right_adj;
  int        lcc_ArrId;
  double     lcc_SM;
  double     lcc_A;
  double     lcc_S_value;
  //  double     adj_edges_w;
  
  node left_real_neighbor;
  node right_real_neighbor;

  //  double bestLCC_SM;
  //  double bestLCC_LCC_M;
  //  int bestLCC_ArrId;


  //  double LCC_A;

  //  int save_lcc_ArrId;
  
  Cnode()
    {
      contr_edges = 0;
	visited = false;
is_he = false;
hedge_rw = 0;

    w_sum = 0;
    rw_sum = 0;
    matched = false;
      bamg_insignificant = false;
      number_of_Cneighbors = 0;
      longest_redge = 0;
      M = pow(((double)random())/((double)RAND_MAX), 5);
      bval = M*100;      
      status = fine;
      last_adj_num = 0;
      lcc_flag = false;
      left_real_neighbor = nil;
      right_real_neighbor = nil;
      lcc_depth = 0;
      iw_links.clear();
      future_edges_sum = 0;

      double_dummy = 0;
      wd = 0;
    }

  void reset()
    {
	visited = false;
is_he = false;
hedge_rw = 0;

    w_sum = 0;
    rw_sum = 0;
    matched = false;
      bamg_insignificant = false;
      number_of_Cneighbors = 0;
      longest_redge = 0;
      M = pow(((double)random())/((double)RAND_MAX), 5);
      bval = M*100;      
      status = fine;
      last_adj_num = 0;
      lcc_flag = false;
      left_real_neighbor = nil;
      right_real_neighbor = nil;
      lcc_depth = 0;
      iw_links.clear();
      future_edges_sum = 0;

      double_dummy = 0;
      wd = 0;
    }
	  
  friend ostream& operator<<(ostream& o, const Cnode& s)
    { 
      return o;
    }

  friend istream& operator>>(istream& i, Cnode& s)
    {
      return i;
    }
};

class CEdge {
 public:
           double ps_from_t;
        double pt_from_s;
        bool in_mst;

	//list<double> stdev_vals;
	//array<int> rw_distr;
	 //double nw;
  //bool w_change;
  bool imaginary;
  int io_score;
  double double_io_score;
 // double old_W;
  //double old_rw;
  double rw;
  //double rw2;  
  //double w_rw;
  //double H;
  //double SDV;
  //double before_normalization;
  //double prev_w;
  double M;
  double w;
  double real_w;
  //double pw;
  double iw,pw;
  //double bamg_iw;
  //double power_coeff;
  int visited;
  
  bool is_new;
  bool be_deleted;
  

  //  double w3;

  bool InFlipFlag;

  CEdge()
    {
	  
	  imaginary = false;
      M = rand();
      rw =0 ;
      //old_rw = 0;
      w = 0;
     // prev_w = 0;
      //      w3 = 0;
      InFlipFlag = false;
      be_deleted = false;
      is_new = false;
      visited = -1;
    }
  void reset()
    {
      M = rand();
      rw =0 ;
      //old_rw = 0;
      
     // prev_w = 0;
      //      w3 = 0;
      InFlipFlag = false;
      be_deleted = false;
      is_new = false;
      visited = -1;
    }
  void save(file_ostream o)
    {
      //o << w << endl << pw << endl << iw << endl;
    }
  
  friend ostream& operator<<(ostream& o, const CEdge& s)
    { 
      return o;
    }
  
  friend istream& operator>>(istream& i, CEdge& s)
    {
      return i;
    }
  
};

typedef UGRAPH<Cnode, CEdge> TGraphC;

class tmp_edge_cost {
 public:
  edge e;
  double w;
  
  friend ostream& operator<<(ostream& o, const  tmp_edge_cost& s)
    {
      return o;
    }
  
  friend istream& operator>>(istream& i,  const tmp_edge_cost& s)
    {
      return i;
    }
};


class triple_res {
public:
  int u_pos;
  int v_pos;
  int w_pos;

  double cost;
  
  friend ostream& operator<<(ostream& o, const triple_res& s)
  { 
    return o;
  }
  
  friend istream& operator>>(istream& i, const triple_res& s)
  {
    return i;
  }
};

class LCCNode {
 public:
  int pos;
  double M;
  double SM;
};

class LCC {
 public:
  node_array<LCCNode> order;
};
class CS_Window {
public :
  double start;
  double finish;
  list<node> l;

  friend ostream& operator<<(ostream& o, const CS_Window& s)
    { 
      return o;
    }
  
  friend istream& operator>>(istream& i, CS_Window& s)
    {
      return i;
    }

  CS_Window() { start = 0; finish = 0;}

  void reset()
  {start = 0; finish = 0;l.clear();
  }
  void insert(node v, TGraphC & G)
  {
    if(G[v].interpol_coord - G[v].w/2.0<start)
      start = G[v].interpol_coord - G[v].w/2.0;
    if(G[v].interpol_coord + G[v].w/2.0>finish)
      finish = G[v].interpol_coord + G[v].w/2.0;
    l.push_back(v);
  }   
};
