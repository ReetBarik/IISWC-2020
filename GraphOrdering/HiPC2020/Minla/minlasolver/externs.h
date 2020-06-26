#include "env.h"
#include <time.h>

#define dblmax(X, Y)  ((X) > (Y) ? (double)(X) : (double)(Y))
#define dblmin(X, Y)  ((X) < (Y) ? (double)(X) : (double)(Y))
#define mymax(X, Y)  ((X) > (Y) ? (X) : (Y))
#define mymin(X, Y)  ((X) < (Y) ? (X) : (Y))
extern double SolveRelax(TGraphC & G);
extern double Solve_byHighestDegreeRemoval(TGraphC & G);

extern void print_G_lemon(TGraphC& G);
extern void dot_write(TGraphC& G, leda::string f);

extern void squize_nodes(TGraphC& G);
extern void sparsify_graph(TGraphC& G);

extern void balanced_coarsening_operator_by_algdist(TGraphC & G, TGraphC & H);
extern void balanced_coarsening_operator_by_edgeweights(TGraphC & G, TGraphC & H);

extern void Jacobi_get_best_rw_value(int rit, TGraphC& G, node_array<array<double> > & rv);
extern void mloga_Jacobi_get_best_rw_value(int rit, TGraphC& G, node_array<array<double> > & rv);
extern void mloga_Jacobi_correct_rw2(TGraphC & G);
extern void Jacobi_correct_rw2(TGraphC & G);
extern void calc_wdegree_futurevolume(TGraphC & G);

extern void chris_output(TGraphC& G, TGraphC& H);
extern void chris_coarsest_output(TGraphC& G);

extern array<list<Pair_LOfInt_Double> > basis_solutions; // minla.cpp
extern array<node> dummyNodesArr;
extern bool char_search(leda::string s, char c); // utils.cpp
extern bool move_with_min_prob(double q, double p);
extern clock_t prev_clock;
extern double apply_assignment_by_initId(TGraphC & G, leda::string s);
extern double calc_all_max_weights(node u, node v, TGraphC & G);
extern double calc_all_weights(node u, node v, TGraphC & G);
extern double calc_cut(TGraphC & G);
extern double calcd_bamg_iw(TGraphC & G, list<edge> & le, node fnode); //linmat.cpp
extern double calc_energy_by_interpol_coords(TGraphC & G);
extern double calc_laC(TGraphC & G); // calc_la.cpp
extern double calc_la_w(UGRAPH<int, double>& U, node_array<int>& Vu);
extern double calc_loglaC(TGraphC & G); // calc_la.cpp
extern double calc_logsum(TGraphC & G); // calc_la.cpp
extern double calc_print_cut(TGraphC & G, bool & ret_violate);
extern double calc_real_laC(TGraphC & G); // calc_la.cpp
extern double check_seeds_by_compatible_relaxation(TGraphC & G);
extern double coarse(TGraphC & G, bool first_time_down, bool in_fmg); // s3_clust.cpp
extern double coarse(TGraphC & G); // s3_clust.cpp
extern double coords_by_linear_order(TGraphC & G); //relaxation.cpp
extern double cpu_time;
extern double E0;

extern double fiedler_relaxation(TGraphC & G, int dist);
extern double find_min_laC(TGraphC & U);
extern double fix_edge_relative_weight(TGraphC & G, edge e);
extern double flip_two_lcc_nodesC(node v, node u, TGraphC & G, double current_la); // flips.cpp
extern double flip_two_nodesC(node v, node u, TGraphC & G, double current_la); // flips.cpp
extern double gen_prob(); //utils.cpp
extern double insert_minimization(TGraphC & G, bool  for_fines_only); // relaxation.cpp
extern double insert_sa(TGraphC & G, double lin_arr);
extern double lcc_triples_minimization(TGraphC & G); //relaxation.cpp
extern double light_minimization(TGraphC & G); //relaxation.cpp
extern double linear_order_by_coord(TGraphC & G); // relaxation.cpp
extern double method_s3(TGraphC & G, bool first_time_down); // s3_clust.cpp

extern double minimization(TGraphC & G); //relaxation.cpp
extern double minimize_dist(TGraphC& G, int dist); // perm_relax.cpp
extern double new_sa_on_distance(TGraphC & G, int dist, bool on_fines_only); //improve.cpp
extern double one_pass_matching_minimization(TGraphC & G, double lin_arr);
extern double only_neighbors_relaxation(TGraphC & G); //relaxation.cpp
extern double print_graph_info(TGraphC & G);
extern double quad(double x);
extern double restore_lcc_arrangement(TGraphC & G, array<node> & gap, double lin_arr);
extern double s2d(leda::string s); // utils.cpp
//extern double sa_improveC(TGraphC & G);
extern double sa_on_distance(TGraphC & G, int dist); //improve.cpp
extern double sa_on_distance(TGraphC & G, int dist, int HC); //improve.cpp
//extern double sa_on_fines_distance(TGraphC & G, int dist); //improve.cpp
extern double sa_only_neighbors(TGraphC & G); //improve.cpp
extern double sa_on_triples(TGraphC & G); //sa_triples.cpp
//extern double saRb_improveC(TGraphC & G, list<node> & green_nodes);
//extern double saR_improveC(TGraphC & G);
extern double segment_minimization(TGraphC & G); // relaxation.cpp
//extern double segment_relaxation(TGraphC & G); //relaxation.cpp
extern double solve_exact_minla_C(TGraphC & G);
extern double sort_fine_nodes_by_wdegree(TGraphC & G);
extern double sort_fine_nodes_by_weighted_degree_s3(TGraphC & G);
extern double start_threshold_on_nodes;
extern double triples_minimization(TGraphC & G); //relaxation.cpp
extern double uncoarse(TGraphC & G, TGraphC & H);
extern double V_Cycle_Iterations(TGraphC & G); // s3_clust.cpp
extern double wmm_minimization(TGraphC & G); //relaxation.cpp
extern edge edge_test(node v, node w, TGraphC & G);
extern edge is_edge(node v, node w,  TGraphC & G);
extern int algdist_Jacobi_underrelax(TGraphC & G);
extern int algdist_minla_relax(TGraphC & G);
extern int algdist_nonsymJacobi_underrelax(TGraphC & G);
extern int algdist_symJacobi_underrelax_stage2(TGraphC & G);
extern int algdist_symJacobi_underrelax(TGraphC & G);

extern int basis_solutions_passed; // s3_clust.cpp
extern int calc_HC_steps(TGraphC & G); // sa.cpp
extern int cossim_algdist_symJacobi_underrelax(TGraphC & G);
extern int interpol_order(TGraphC & G);
extern int MAX_LEVEL;
extern int nodes_algdist_compatible_relaxation2d(TGraphC & G);
extern int nodes_algdist_relaxation2d(TGraphC & G);
extern int nodes_algdist_relaxation(TGraphC & G);

extern int s2i(leda::string s); // utils.cpp
extern int sourceG_NofNodes;
extern int V_Cycle_Number;
extern LA best_previous_V_Cycle_LA;
extern LA find_min_la(UGRAPH<int, double>& U, node_array<int>& Vu);
extern leda::string d2s(double); // utils.cpp
extern leda::string getTime (); // utils.cpp
extern leda::string i2s(int); // utils.cpp
extern leda::string SI;
extern list<TGraphC> AllLevelsGraphs; // s3_clust.cpp
extern node second_adj_for_edge(edge e, node v, TGraphC & G);
extern StatInfo Params;
extern TGraphC * TMP_CMP_GRAPHC;
extern TGraph * TMP_CMP_GRAPH;
//extern void add_heavy_nodes(TGraphC & G, TGraphC & H, double & sum_fines_wd);
extern void add_heavy_nodes(TGraphC & G, TGraphC & H, double & sum_fines_wd);
extern void add_locally_heavy_nodes(TGraphC & G, TGraphC & H);
extern void alg_calc_iw_by_rw_and_w(TGraphC & G);
extern void alg_calc_iw_by_rw_balanced_new(TGraphC & G, TGraphC & H);
extern void alg_calc_iw_by_rw_balanced(TGraphC & G, TGraphC & H);
extern void alg_calc_iw_by_rw(TGraphC & G);
extern void alg_calc_iw_by_rw(TGraphC & G, TGraphC & H);
extern void alg_calc_iw(TGraphC & G);
extern void  algdist_matching(TGraphC & G, TGraphC & H/*, list<node> & ret*/);
extern void bamg_calc_iw(TGraphC & G);
extern void calc_bamg_iw(TGraphC & G, list<edge> & le, node fnode); //linmat.cpp
extern void calc_preliminary_bamg_iw_without_prefiltering(TGraphC & G);
extern void calc_preliminary_bamg_iw_with_prefiltering(TGraphC & G);
extern void calc_pw(TGraphC & G);
extern void clear_iw_links(TGraphC & G);
extern void clust_step(TGraph & G);
extern void coarse_nodes_alg_and_wag_sum(TGraphC & G, TGraphC & H/*, list<node> & ret*/);
extern void coarse_nodes_alg_max(TGraphC & G, TGraphC & H/*, list<node> & ret*/);
extern void coarse_nodes_alg_sum(TGraphC & G, TGraphC & H/*, list<node> & ret*/);
extern void color_and_reorder_vertices(TGraphC & G);
extern void construct_new_graph_dclust(TGraphC & H, TGraphC & G);
extern void construct_new_graph_max_dclust(TGraphC & H, TGraphC & G);
extern void cut_relax(bool compatible_only, TGraphC & G, bool fine_tuning_only);
extern void cut_relax_kl(TGraphC & G);
extern void cut_relax_resolve_violations(TGraphC & G);
extern void cut_relax_rw(bool compatible_only, TGraphC & G, bool fine_tuning_only);
extern void decrease_lin_arr_by_gap(TGraphC & G, array<node> & gap, double & lin_arr); // utils.cpp
extern void decrease_lin_arr_by_nlist(TGraphC & G, list<node> & gap, double & lin_arr); // utils.cpp
extern void define_S_values_after_swap(node u, node v, TGraphC & G); //calc_la.cpp
extern void define_S_values(TGraphC & G); //calc_la.cpp
extern void dummy_flip_two_nodesC(node v, node u, TGraphC & G); // flips.cpp
extern void get_imb_factor(array<double> & csize, double & AV_CSIZE,  double & ret_avg_imb, double & ret_max_imb);
extern void gml_ordered_write(leda::string f, TGraphC & G);
extern void gml_weights_write(leda::string f, TGraphC & G,edge_array<double> & new_costs);
extern void gml_write(leda::string f, TGraphC & G);
extern void gml_write_with_iw(leda::string f, TGraphC & G);
extern void gml_write_with_pw(leda::string f, TGraphC & G);
extern void gml_write_with_w(leda::string f, TGraphC & G);
extern void graph_format_print( TGraphC & G, leda::string title);
extern void elist_format_print( TGraphC & G, leda::string title);
extern void graph_print_4netresp( TGraphC & G, leda::string title);
extern void graph_print(TGraphC & G, leda::string title);
extern void hmetis_format_print( TGraphC & G, leda::string title);
extern void increase_lin_arr_by_gap(TGraphC & G, array<node> & gap, double & lin_arr); // utils.cpp
extern void increase_lin_arr_by_nlist(TGraphC & G, list<node> & gap, double & lin_arr); // utils.cpp
//extern void interpolate_C(edge_array<double> & Etmp, list<node> & green_nodes, TGraphC & G, TGraphC & H);
extern void interpolate_s2(TGraphC & G, TGraphC & H);
extern void interpolate_s3(edge_array<double> & new_edges_costs, TGraphC & G, TGraphC & H);
extern void lcc_init(TGraphC & G); // lcc.cpp
extern void lcc_update(TGraphC & G, double & lin_arr, double & lcc_cost); // lcc.cpp
extern void load_metis_solution( TGraphC & G, leda::string title);
extern void local_save_arrangement(TGraphC & G, node u, node v); // utils.cpp
extern void matrix_print( TGraphC & G, leda::string title);
extern void maxmatching_coarsening(TGraphC & G, TGraphC & H/*, list<node> & ret*/);
extern void metis_format_print( TGraphC & G, leda::string title);
extern void minla_alg_calc_iw_by_rw(TGraphC & G);
extern void mloga_alg_calc_iw_by_rw(TGraphC & G);
extern void MY_error(char *fmt, ...);
extern void new_bamg_alg_calc_iw_by_rw(TGraphC & G);
extern void nodes_gs(TGraphC & G);
extern void nodes_jac(TGraphC & G);
extern void /*old_good_*/minimize_vertex_coordinate(node v, TGraphC & G);
extern void Perm (int n, int NN, array<int> & p, array<int> & pi, array<int> & dir, list<list<int> > & L);
void output_visMatlab(TGraphC& G, leda::string fname);
extern void print_diff_save_ArrId_and_ArrId(TGraphC & G); // utils.cpp
extern void print_Laplacian(TGraphC& G, leda::string f);
extern void print_part_solution( TGraphC & G, leda::string title);
extern void put_arrangement_into_G(TGraphC & G, list<int> & A);
extern void read_arrangement(TGraphC & G, leda::string filename);
extern void read_graph(TGraphC & G, leda::string filename);
extern void read_lemon(TGraphC & G, leda::string filename);
extern void read_mtx(TGraphC & G, leda::string filename);
extern void read_rmf2(TGraphC & G, leda::string filename);
extern void read_rmfext(TGraphC & G, leda::string filename);
extern void read_rmf(TGraphC & G, leda::string filename);
extern void read_ine(TGraphC & G, leda::string filename);
extern void refine_alg_calc_iw_by_rw(TGraphC & G, TGraphC & H);
extern void relaxation_gauss_seidel(TGraphC & G);
extern void relaxation_of_fines(TGraphC & G); // relaxation.cpp
extern void relaxation_of_gs(TGraphC & G); // relaxation.cpp
extern void restore_arrangement(TGraphC & G/*, array<int> & MinArr*/); // utils.cpp
extern void restore_best_lcc_arrangement(TGraphC & G);
extern void row_graph_print( TGraphC & G, leda::string title);
extern void save_arrangement(TGraphC & G/*, array<int> & MinArr*/); // utils.cpp
extern void save_best_lcc_arrangement(TGraphC & G);
extern void save_lcc_arrangement(TGraphC & G, array<node> & gap);
extern void store_graph_arrangement(TGraphC & G, leda::string filename);
extern void update_time();
extern void wag_calc_iw(TGraphC & G);
extern void wag_coarse_nodes(TGraphC & G, TGraphC & H/*, list<node> & ret*/);
void check_same_coordinates(TGraphC & G); // utils.cpp
