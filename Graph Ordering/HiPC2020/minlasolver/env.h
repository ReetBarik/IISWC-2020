#include "common.h"

using namespace std;

extern double gen_prob();

enum TSparsCoeff {abs_coeff = 0, rel_coeff = 1};

enum TImprove {all_pairs = 0, rnd_all = 1, sa = 2, no = 3, saR = 4};
enum TMethod {METHOD_DCLUST_MAX = 0,
              METHOD_DCLUST_AVG = 1,
              METHOD_ECLUST = 2,
              METHOD_SRG = 3,
              METHOD_S2 = 4,
              METHOD_S3 = 5};

enum RAISED_POWER_WIND_TYPE {adaptive = 0, non_adaptive = 1};
enum ADD_HEAVY_TYPE {local_heavy = 0, global_heavy = 1, no_heavy = 2};
enum TProblem {mloga = 0, kpart = 1, netgen = 2, nothing = 3, netresp = 4, minla = 5};

enum TAlgDistVolumes {future_wdeg = 0, wdeg = 1, novol = 2};
enum UNSAT_ENFORCE_TYPE {all_levels = 0, finest_level = 1};
enum FRAMEWORK_TYPE {multiscale = 0, relaxation = 1};

enum GRAPH_INPUT_TYPE {graph = 0, rmf = 1};

class StatInfo {
 public:
 TSparsCoeff sparsification_coeff_type;
 double sparsification_coeff;
   bool balanced_coarsening;
int io2uptolevel;
	FRAMEWORK_TYPE framework;
	UNSAT_ENFORCE_TYPE unsat_enforce;

   int E0, V0;
   int CURRENT_LEVEL;
 int min_part_size;
 int zero_degree_nodes;
 int external_part_size;
 double max_node_vol;
 bool fmg;
 double imbalance;
 double imb_incr;
  //bool add_heavy;
 //bool logsum;
 TProblem problem;
 
 double threshold_heavy_nodes;
bool output_ordering;
  alg_dist_type alg_dist;
  int alg_distance_iter;
  int alg_distance_rvec;
  TAlgDistVolumes alg_dist_volumes;
  bool alg_dist_normalize;
  ADD_HEAVY_TYPE add_heavy;
  T_iw_calc iw_calc_type;
  sym_type sym_nonsym;

  RAISED_POWER_WIND_TYPE raised_power_window_type;
  int raised_power_window_size;
  int raised_power_power_step;
  int original_correction_scheme_window_size; // window size start



  int correction_scheme_window_size; // window size start
  int correction_scheme_iterations;
  //  int correction_scheme_iterations;

  leda::string gen_graph;
  int segment_minimization;
  int solution_number;
  int seed;
  int prog_id;
  int hot_cold_steps;
  double threshold_edge_weight;
  bool improved_la_function;
  fvol_sort_type future_volume_sort;
  bool individual_minimization;
  bool many_cycles_in_min;
  int gs_relaxation_sweeps;
  int nlp_relaxation_sweeps;
  int comp_relaxation_sweeps;
  bool random_reorder;
  bool const_number_of_hc;
  TImprove improve;
  int relax_sweeps;
  leda::string statInfo;
  //  int improve_bound;
  int min_graph_size;
  int basis_solutions_number;
  int basis_solutions_type; // 0 - first X b.s., 1 - C/D best X b.s., 2 - first X b.s. with step
  int pair_dist;
  int correction_scheme_delta_dist;
  int insert_minimization_dist;
  int insert_sa_dist;
  bool gml_out;
  //  bool do_sa;
  bool do_sa_1, do_sa_2, do_sa_3;
  bool use_coord_in_sa;
  bool use_edges_abs_filter;
  bool use_lcc;
  int  max_lcc_length;
  bool reverse_sa;
  bool use_sa_in_relax;
  bool use_sa_only_in_last_two_levels;
  bool check_best_relax;
  bool use_compatible_relaxation;
  bool use_gauss_seidel_relaxation;
  int number_of_sweeps_in_one_hc;
  int strict_minimization; // 0 - neighbors, 1 - triples, 2 - wmm
  double sa_coor;
  TMethod method;
  
  double sa_alpha; // Simulated annealing decreasing factor
  int start_sa1perc; // Start perm. percent for sa, dist=1
  int finish_sa1perc;
  int start_sa2perc; // Start perm. percent for sa, dist=2
  int finish_sa2perc;
  int sa_permition_percent;
  //  int start_sa_triples_perc;
  //  int finish_sa_triples_perc;
  bool const_number_of_sa1pr;
  
  double threshold_coarse_edges;
  double threshold_coarse_nodes_start;
  //  float threshold_coarse_nodes_finish;
  //  float threshold_coarse_nodes_step;
  double threshold_coarse_nodes;
  double threshold_coarse_nodes_current;
  double amount_of_coarse_nodes; // not more than ... of the number of input nodes
  int power_of_connection_between_fine_and_its_coarses; // In %
  int maximum_of_connections_between_fine_and_its_coarses;
  int V_Cycle_Iter_Num; // The number of V-cycles

  /*
  int get_maximum_of_connections_between_fine_and_its_coarses(long N)
    {
      
      if(N<500)
        return maximum_of_connections_between_fine_and_its_coarses + 10;
      else if((N>=500)&&(N<1000))
        return maximum_of_connections_between_fine_and_its_coarses + 7;
      else if((N>=1000)&&(N<5000))
        return maximum_of_connections_between_fine_and_its_coarses + 5;
      else if((N>=5000)&&(N<10000))
        return maximum_of_connections_between_fine_and_its_coarses + 3;
      else
      
        return maximum_of_connections_between_fine_and_its_coarses;
      
    }
  
  double rel_filter(long N)
    {
      
      if(N<500)
        return threshold_edge_weight*0.01;
      else if((N>=500)&&(N<1000))
        return threshold_edge_weight*0.05;
      else if((N>=1000)&&(N<5000))
        return threshold_edge_weight*0.1;
      else if((N>=5000)&&(N<10000))
        return threshold_edge_weight*0.5;
        else
      
        return threshold_edge_weight;
 
    }
  */
};

class LA {
public:
  node_array<int> V;
  double weight;

  friend ostream& operator<<(ostream& o, const  LA& s)
    {
      return o;
    }
  
  friend istream& operator>>(istream& i,  const LA& s)
    {
      return i;
    }
};

class Pair_LOfInt_Double {
public:
  array<double> coords;
  list<int> L;
  double c;
  double diff_mark;
  

  friend ostream& operator<<(ostream& o, const  Pair_LOfInt_Double& s)
    {
      return o;
    }
  
  friend istream& operator>>(istream& i,  const Pair_LOfInt_Double& s)
    {
      return i;
    }
};

