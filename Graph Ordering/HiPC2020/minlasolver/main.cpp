//#include "externs.h"
#include "cmpfuncs.h"
//#include "minla.h"
#include <time.h>
#include <stdio.h>

leda::string SI;

leda::string read_params(char * infile, StatInfo& p, char *);
void init();
extern double sa_only_neighbors_old(TGraphC & G);

int sourceG_NofNodes;
double start_threshold_on_nodes;

double calc_cost_by_initial_id(TGraphC & G)
{
  node v;
  forall_nodes(v, G)
    G[v].ArrId = G[v].initial_id;

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

  cerr << "Calculated cost (according to initial id) = " << calc_laC(G) << endl;
  exit(1);
}

void calc_qp(TGraphC & G)
{
  double s=0;edge e;
  forall_edges(e, G)
    {
      s+=pow(fabs(G[G.source(e)].ArrId -G[G.target(e)].ArrId),2);

    }
  s = sqrt(s);
  cerr << "QP : " << s << endl;
}
leda::string instance_file;

double apply_assignment_by_initId(TGraphC & G, leda::string s1)
{
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_initial_id);

  node v = G.first_node();
  leda::string int_str = "";
  //cerr << "strlen = " << s1.length() << endl;
  for(int i=0; i<s1.length()-1; i++)
    {
      if((s1[i]=='0')||(s1[i]=='1')||(s1[i]=='2')||(s1[i]=='3')||(s1[i]=='4')||(s1[i]=='5')||(s1[i]=='6')||(s1[i]=='7')||(s1[i]=='8')||(s1[i]=='9'))
        int_str+=s1[i];
      else
        {
          //cout << G[v].initial_id << " gets " << s2i(int_str) << endl;
          int node_pos = s2i(int_str);
          G[v].ArrId = node_pos;

          int_str = "";
          if(v!=G.last_node())
          v = G.succ_node(v);//cerr << G[v].ArrId << endl;
        }
    }

  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_ArrId);
  define_S_values(G);

  return calc_laC(G);
}

int main(int argc, char **argv)
{
  instance_file = argv[1];
  if((leda::string)argv[1]=="-lcc")
    {
      int L=1000000;
      TGraphC G;
      if(argv[2][strlen(argv[2])-1]=='f')
        read_rmf(G, argv[2]);
      else if(argv[2][strlen(argv[2])-1]=='r')
          read_rmf2(G, argv[2]);

      else
        read_graph(G, argv[2]);

      leda::string arr_file = (leda::string)argv[3];
      FILE * inFile;
      inFile = fopen(argv[3], "r");
      //      cerr << "arrangements file opened" << endl;
      char s[L];
      fgets(s, L, inFile);
      double la_cost = apply_assignment_by_initId(G, s);
      //      cerr << "first arrangement applied " << la_cost << endl;
      lcc_init(G);
      double lcc_cost = la_cost;
      cerr << "Initial cost : " << la_cost << endl;

       while(!feof(inFile))
         {
           fgets(s, L, inFile);
           la_cost = apply_assignment_by_initId(G, s);
           cerr << "arrangement readed " << la_cost << endl;
           lcc_update(G, la_cost, lcc_cost);
           if(la_cost<lcc_cost)
             {
               lcc_init(G);
               lcc_cost = la_cost;
             }
           cerr << "LCC NEW COST " << lcc_cost << endl;
         }
       fclose(inFile);

      return 0;
    }

  clock_t start = clock();
  long t = time(0);
  init();
  leda::string all_params = read_params(argv[2], Params, argv[3]);


      Params.threshold_coarse_nodes_current=Params.threshold_coarse_nodes;

      // SI = all_params;
      SI += getTime();
      SI += "InFile: " + (leda::string)argv[1] + "\n";

      cerr << argv[3] << endl;
      TGraphC G;
      if(argv[1][strlen(argv[1])-1]=='f')
        read_rmf(G, argv[1]);
      else if(argv[1][strlen(argv[1])-1]=='t')
        read_rmfext(G, argv[1]);
      else if(argv[1][strlen(argv[1])-1]=='x')
        read_mtx(G, argv[1]);
      else if(argv[1][strlen(argv[1])-1]=='r')
        {
          cerr <<"gr" << endl;
          read_rmf2(G, argv[1]);
        }
      else
        read_graph(G, argv[1]);

      sourceG_NofNodes = G.number_of_nodes();

      //      calc_cost_by_initial_id(G);

      /*
      S=0;
      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_nC);
      define_S_values(G);
      triples_minimization(G);
      double min = calc_laC(G);
      cerr << "do sa     " <<  min << endl;

      min = sa_on_distance(G, 1);
      min = sa_on_distance(G, 2);
      cerr << "posle sa  " <<  min << endl;
      exit(1);
      */


      /*
      TMP_CMP_GRAPHC = &G;
      G.sort_nodes(&cmp_nC);
      define_S_values(G);
      minimize_dist(G,3);
      //      cerr << calc_laC(G) << endl;
      //triples_minimization(G);
      cerr << calc_laC(G) << endl;
      //      cerr << sa_on_distance(G, 1) << endl;
      exit(1);
      
      
      node x;
      forall_nodes(x, G) {
	if(G.degree(x)==0)
	  cerr << G[x].initial_id << "(" << G[x].ArrId << "), ";
      }
      cerr << endl;
*/
     /* 
      // adding communication graph
      forall_nodes(x, G) {
	if((G[x].initial_id > 1024)&&(G[x].initial_id < 1041)) {
	  for(int k=0; k<8; k++) {
	    Cnode n;
	    n.ArrId = G.number_of_nodes()+1;
	    n.initial_id = G.number_of_nodes()+1;
	    n.name = i2s(G.number_of_nodes()+1);
	    n.w = 1;
	    G.new_node(n);
	    
	    CEdge ne; ne.w = 2000;
	    ne.real_w = 2000;
	    G.new_edge(x, n, ne); 
	    
	    for(int l=0; l<7; l++) {
	      Cnode n1;
	      n1.ArrId = G.number_of_nodes()+1;
	      n1.initial_id = G.number_of_nodes()+1;
	      n1.name = i2s(G.number_of_nodes()+1);
	      n1.w = 1;
	      G.new_node(n1);
	    
	      CEdge ne; ne.w = 256;
	      ne.real_w = 256;
	      G.new_edge(n1, n, ne); 	      
	    }

	  }
	}
      }
      */
      double res;
      edge x;
     /* 
      forall_edges(x, G) {
	node p = G.source(x);
	node q = G.target(x);
	
	G[x].w = 1;
	if((G[p].initial_id > 1040) || (G[q].initial_id > 1040))
	  G[x].w = 1000;//G[x].w * 10000;
	
	if(((G[p].initial_id > 1024) && (G[p].initial_id < 1041) && (G[q].initial_id > 1040) && (G[q].initial_id > 1040)) || 
	   ((G[q].initial_id > 1024) && (G[q].initial_id < 1041) && (G[p].initial_id > 1040) && (G[p].initial_id > 1040)))
	  G[x].w = 1000;//G[x].w * 10000;

	if(((G[p].initial_id > 1024) && (G[p].initial_id < 1041) && (G[q].initial_id < 1025) && (G[q].initial_id < 1025)) || 
	   ((G[q].initial_id > 1024) && (G[q].initial_id < 1041) && (G[p].initial_id < 1025) && (G[p].initial_id < 1025)))
	  G[x].w = 1000;//G[x].w * 10000;

	//if((G[p].initial_id < 1025) && (G[q].initial_id < 1025))
	 // G[x].w = 1;//G[x].w / 100.0;
	
	
	G[x].real_w = G[x].w;
      }
      */
      res = V_Cycle_Iterations(G);
      TMP_CMP_GRAPHC = &G;
      
      
      //G.sort_nodes(&cmp_ArrId);

      //cout << "p ghct " << G.number_of_nodes() << " " << G.number_of_edges() << endl;
      node v;
      /*
      forall_nodes(v, G)
	cout << G[v].ArrId << " ";
      cout << endl;
      */
      //print_Laplacian(G, "lap.m");
      //exit(1);

      //      graph_print(G, "Final : ");

      char s[20];
      sprintf(s, "%f", res);

      cerr << "FILE=" << (leda::string)argv[1] <<",COST=" << s << endl;
      //      SI += "Threshold on coarse nodes=" + d2s(Params.threshold_coarse_nodes_current) + "\n";
      //SI = SI + "MinLA=" + s + "\n";
      //      sprintf(s, "%f", G.number_of_nodes()*G.number_of_edges()/3.0);
      //      SI = SI + "Expected randomly generated order cost = " + s + "\n";
      //      sprintf(s, "%d", S);
      //      SI += "Levels=" + (leda::string)s + "\n";
      //cerr << SI;

      //      cerr << "Total computation CPU time = " << ((double)(clock() - start)) / (double)CLOCKS_PER_SEC << " s.\n" << endl;
      //      cerr << "Total computation CPU time = " << ((double)(clock() - start)) / (double)CLOCKS_PER_SEC/60.0 << " m.\n" << endl;

      cerr << "REAL TOTAL CPU TIME = " << cpu_time << " M." << endl;

      //      cerr << "TIME (H) " << (double)(time(0)-t)/3600.0 << endl;
      //      cerr << "TIME (M) " << (double)(time(0)-t)/60.0 << endl;
      //      cerr << "TIME (S) " << time(0)-t << endl;

      //calc_qp(G);
      if(Params.output_ordering==true) {
        G.sort_nodes(&cmp_initial_id);
        output_visMatlab(G, "visgraph.m");

		graph_print(G, "nodes_order.dat");
	  }

      //      graph_print(G, "final : ");

      //  store_graph_arrangement(G, "ga.dat");
}

void init()
{
  basis_solutions_passed = 0;
  srand(Params.seed);
  //  srandom(Params.seed);
  Params.threshold_coarse_edges = 0;
  Params.do_sa_1 = false;
  Params.do_sa_2 = false;
  Params.do_sa_3 = false;

  cpu_time = 0;
  prev_clock = clock();
  update_time();

  Params.prog_id = rand();
}

leda::string read_params(char * infile, StatInfo& p, char * possible_seed)
{
  p.V_Cycle_Iter_Num = 1;

  leda::string all_params = "";
  ifstream fl;
  char s[512];
  fl.open(infile, ios::in);
  while(fl.getline(s, 512, '\n'))
    {
      leda::string str = s;
      all_params += str;
      all_params +="\n";

      cerr << s << endl;

      if(str(0,5) == "[SEED]")
        {
          leda::string ss = str(6, str.length());
          if(ss=="TIME")
            {
              p.seed = time(0);
            }
          else
            {
              int res = s2i(ss);
              p.seed = res;
            }
        }

      if(str(0,5) == "[SEED]")
        {
          leda::string ss = str(6, str.length());
          if(ss=="TIME")
            {
              p.seed = time(0);
            }
          else if(ss=="CL")
            {
              p.seed = s2i(possible_seed);
            }
          else
            {
              int res = s2i(ss);
              p.seed = res;
            }
          cerr << "GOT SEED " << p.seed << endl;
        }

      if(str(0,5) == "[IMPR]")
        {
          leda::string ss = str(6, str.length());
          int res = s2i(ss);
          if(res==0)
            p.improve = all_pairs;
          else if(res==2)
            p.improve = sa;
          else if(res==4)
            p.improve = saR;
          else if(res==3)
            p.improve = no;
          //          SI += "Improve=" + i2s(p.improve) + "\n";
        }
      /*
      if(str(0,7) == "[METHOD]")
        {
          leda::string ss = str(8, str.length());
          int res = s2i(ss);
          if(res==0)
            p.method = METHOD_DCLUST_MAX;
          else if(res==1)
            p.method = METHOD_DCLUST_AVG;
          else if(res==2)
            p.method = METHOD_ECLUST;
          else if(res==3)
            p.method = METHOD_SRG;
          else if(res==4)
            p.method = METHOD_S2;
          else if(res==5)
            p.method = METHOD_S3;
          //          SI += "Method=" + i2s(p.method) + "\n";
        }
      */
      /*
      if(str(0,6) == "[IMPRB]")
        {
          leda::string ss = str(7, str.length());
          p.improve_bound = s2i(ss);
          SI += "Improve bound=" + i2s(p.improve_bound) + "\n";
        }
      */
      if(str(0,13) == "[RELAX_SWEEPS]")
        {
          leda::string ss = str(14, str.length());
          p.relax_sweeps = s2i(ss);
          //          SI += "Num. of Relaxation sweeps=" + i2s(p.relax_sweeps) + "\n";
        }
      if(str(0,14) == "[START_SA1PERC]")
        {
          leda::string ss = str(15, str.length());
          p.start_sa1perc = s2i(ss);
        }
      if(str(0,21) == "[SA_PERMITION_PERCENT]")
        {
          leda::string ss = str(22, str.length());
          p.sa_permition_percent = s2i(ss);
        }
      if(str(0,8) == "[SA_TYPE]")
        {
          leda::string ss = str(9, str.length());
          p.use_coord_in_sa = char_search(ss, 'C'); // use coordinates
        }
      if(str(0,13) == "[EDGES_FILTER]")
        {
          leda::string ss = str(14, str.length());
          p.use_edges_abs_filter = char_search(ss, 'A'); // use coordinates
        }
      if(str(0,15) == "[FINISH_SA1PERC]")
        {
          leda::string ss = str(16, str.length());
          p.finish_sa1perc = s2i(ss);
        }
      /*
      if(str(0,17) == "[FINISH_SATRIPLES]")
        {
          leda::string ss = str(18, str.length());
          p.finish_sa_triples_perc = s2i(ss);
        }
      */
      if(str(0,14) == "[START_SA2PERC]")
        {
          leda::string ss = str(15, str.length());
          p.start_sa2perc = s2i(ss);
        }
      if(str(0,15) == "[FINISH_SA2PERC]")
        {
          leda::string ss = str(16, str.length());
          p.finish_sa2perc = s2i(ss);
        }

      if(str(0,13) == "[V_CYCLES_NUM]")
        {
          leda::string ss = str(14, str.length());
          p.V_Cycle_Iter_Num = s2i(ss);
          //          SI += "V_Cycle_Iter_Num = " + i2s(p.V_Cycle_Iter_Num) + "\n";
        }
      if(str(0,23) == "[BASIS_SOLUTIONS_NUMBER]")
        {
          leda::string ss = str(24, str.length());
          p.basis_solutions_number = s2i(ss);
          //          SI += "Number of basis solutions =" + i2s(p.basis_solutions_number) + "\n";
        }
      if(str(0,21) == "[BASIS_SOLUTIONS_TYPE]")
        {
          leda::string ss = str(22, str.length());
          p.basis_solutions_type = s2i(ss);
        }
      if(str(0,9) == "[SA_ALPHA]")
        {
          leda::string ss = str(10, str.length());
          p.sa_alpha = s2d(ss);
          //          SI += "SA Alpha=" + d2s(p.sa_alpha) + "\n";
        }

      if(str(0,23) == "[THRESHOLD_COARSE_EDGES]")
        {
          leda::string ss = str(24, str.length());
          p.threshold_edge_weight = s2d(ss);
          //          SI += "Threshold on coarse edges=" + d2s(p.threshold_coarse_edges) + "\n";
        }
      if(str(0,23) == "[AMOUNT_OF_COARSE_NODES]")
        {
          leda::string ss = str(24, str.length());
          p.amount_of_coarse_nodes = s2d(ss);
          start_threshold_on_nodes = p.amount_of_coarse_nodes;

          //          SI += "Amount of coarse nodes=" + d2s(p.amount_of_coarse_nodes) + "\n";
        }
      if(str(0, 49) == "[POWER_OF_CONNECTION_BETWEEN_FINE_AND_ITS_COARSES]")
        {
          leda::string ss = str(50, str.length());
          p.power_of_connection_between_fine_and_its_coarses = s2i(ss);
          //          SI += "POWER_OF_CONNECTION_BETWEEN_FINE_AND_ITS_COARSES = " + i2s(p.power_of_connection_between_fine_and_its_coarses) + "\n";
        }
      if(str(0, 52) == "[MAXIMUM_OF_CONNECTIONS_BETWEEN_FINE_AND_ITS_COARSES]")
        {
          leda::string ss = str(53, str.length());
          p.maximum_of_connections_between_fine_and_its_coarses = s2i(ss);
          //          SI += "MAXIMUM_OF_CONNECTIONS_BETWEEN_FINE_AND_ITS_COARSES = " + i2s(p.maximum_of_connections_between_fine_and_its_coarses) + "\n";
        }
      if(str(0,13) == "[HOT_COLD_STP]")
        {
          leda::string ss = str(14, str.length());
          p.hot_cold_steps = s2i(ss);
          //          SI += "Num. of hot-cold sweeps=" + i2s(p.hot_cold_steps) + "\n";
        }
      if(str(0,15) == "[MIN_GRAPH_SIZE]")
        {
          leda::string ss = str(16, str.length());
          p.min_graph_size = s2i(ss);
        }
      if(str(0,15) == "[MAX_LCC_LENGTH]")
        {
          leda::string ss = str(16, str.length());
          p.max_lcc_length = s2i(ss);
        }
      if(str(0,23) == "[THRESHOLD_COARSE_NODES]")
        {
          leda::string ss = str(24, str.length());
          p.threshold_coarse_nodes = s2d(ss);

          /*
          leda::string ss = str(24, str.length());
          int separator_1_pos = ss.pos(":",0);
          int separator_2_pos = ss.pos(":",separator_1_pos+1);
          //          cerr << separator_1_pos << "\t" <<  separator_2_pos<<endl;
          p.threshold_coarse_nodes_start = s2d(ss(0,separator_1_pos-1));
          //          cerr << ss(separator_1_pos+1,separator_2_pos-separator_1_pos+1);
          p.threshold_coarse_nodes_finish = s2d(ss(separator_1_pos+1,separator_2_pos-1));
          p.threshold_coarse_nodes_step = s2d(ss(separator_2_pos+1,ss.length()));
          */
        }
      if(str(0,9) == "[PAIRDIST]")
        {
          leda::string ss = str(10, str.length());
          p.pair_dist = s2i(ss);
          //          SI += "Pair distance=" + i2s(p.pair_dist) + "\n";
        }

      if(str(0,8) == "[SA_COOR]")
        {
          leda::string ss = str(9, str.length());
          //          cerr << ss << endl;
          //          cerr << ss << endl;
          p.sa_coor = s2d(ss);
          //          cerr << p.path_length<<endl;
          //          SI += "SA_COOR=" + d2s(p.sa_coor) + "\n";
          //          cerr << p.sa_coor << endl;
          //          exit(1);
        }
      else if(str(0,6) == "[GMLWR]")
        {
          leda::string ss = str(7, str.length());
          if(ss=="YES")
              p.gml_out = true;
          else
            p.gml_out = false;
        }
     else if(str(0,21) == "[IMPROVED_LA_FUNCTION]")
        {
          leda::string ss = str(22, str.length());
          if(ss=="YES")
              p.improved_la_function = true;
          else
              p.improved_la_function = false;
        }
      else if(str(0,19) == "[CONST_NUMBER_OF_HC]")
        {
          leda::string ss = str(20, str.length());
          if(ss=="YES")
              p.const_number_of_hc = true;
          else
              p.const_number_of_hc = false;
        }
      else if(str(0,19) == "[CONST_NUMBER_SA1PR]")
        {
          leda::string ss = str(20, str.length());
          if(ss=="YES")
              p.const_number_of_sa1pr = true;
          else
              p.const_number_of_sa1pr = false;
        }
      else if(str(0,18) == "[ALG_DISTANCE_RVEC]")
        {
          leda::string ss = str(19, str.length());
	  p.alg_distance_rvec = s2i(ss);
        }
      else if(str(0,18) == "[ALG_DISTANCE_ITER]")
        {
          leda::string ss = str(19, str.length());
	  p.alg_distance_iter = s2i(ss);
        }

      else if(str(0,20) == "[STRICT_MINIMIZATION]")
        {
          leda::string ss = str(21, str.length());
          p.strict_minimization = s2i(ss);
        }
      else if(str(0,25) == "[INSERT_MINIMIZATION_DIST]")
        {
          leda::string ss = str(26, str.length());
          p.insert_minimization_dist = s2i(ss);
        }
      else if(str(0,29) == "[CORRECTION_SCHEME_DELTA_DIST]")
        {
          leda::string ss = str(30, str.length());
          p.correction_scheme_delta_dist = s2i(ss);
        }
      else if(str(0,29) == "[CORRECTION_SCHEME_ITERATIONS]")
        {
          leda::string ss = str(30, str.length());
          p.correction_scheme_iterations = s2i(ss);
        }
      else if(str(0,29) == "[CORRECTION_SCHEME_WIND_START]")
        {
          leda::string ss = str(30, str.length());
          p.correction_scheme_window_size = s2i(ss);
        }
      else if(str(0,15) == "[INSERT_SA_DIST]")
        {
          leda::string ss = str(16, str.length());
          p.insert_sa_dist = s2i(ss);
        }
      else if(str(0,15) == "[GS_RELAXATIONS]")
        {
          leda::string ss = str(16, str.length());
          p.gs_relaxation_sweeps = s2i(ss);
        }
     else if(str(0,15) == "[SM_RELAXATIONS]")
        {
          leda::string ss = str(16, str.length());
          p.segment_minimization = s2i(ss);
        }
      else if(str(0,17) == "[COMP_RELAXATIONS]")
        {
          leda::string ss = str(18, str.length());
          p.comp_relaxation_sweeps = s2i(ss);
        }
      else if(str(0,16) == "[SOLUTION_NUMBER]")
        {
          leda::string ss = str(17, str.length());
          p.solution_number = s2i(ss);
        }

      else if(str(0,27) == "[NUMBER_OF_SWEEPS_IN_ONE_HC]")
        {
          leda::string ss = str(28, str.length());
          p.number_of_sweeps_in_one_hc = s2i(ss);
        }
      else if(str(0,16) == "[RELAXATION_TYPE]")
        {
          leda::string ss = str(17, str.length());
          p.use_compatible_relaxation = char_search(ss, 'C');
          p.use_gauss_seidel_relaxation = char_search(ss, 'G');
        }
      else if(str(0,9) == "[GENGRAPH]")
        {
          leda::string ss = str(10, str.length());
          p.gen_graph = ss;
        }
      else if(str(0,18) == "[USE_LCC_ALGORITHM]")
        {
          leda::string ss = str(19, str.length());
          if(ss=="YES")
              p.use_lcc = true;
          else
              p.use_lcc = false;
        }
      else if(str(0,16) == "[OUTPUT_ORDERING]")
        {
          leda::string ss = str(17, str.length());
          if(ss=="YES")
              p.output_ordering = true;
          else
              p.output_ordering = false;
        }
        
      else if(str(0,11) == "[SYM_NONSYM]")
        {
          leda::string ss = str(12, str.length());
	  if(ss=="SYM")
	    p.sym_nonsym = symmetric;
	  else if(ss=="NONSYM")
	    p.sym_nonsym = nonsymmetric;
	  else {
	    cerr << "wrong parameter SYM_NONSYM" << endl;
	    exit(1);
	  }
	}
      
      else  if(str(0,22) == "[THRESHOLD_HEAVY_NODES]")
        {
          leda::string ss = str(23, str.length());
          p.threshold_heavy_nodes = s2d(ss);
        }
      else if(str(0,8) == "[PROBLEM]")
        {

          leda::string ss = str(9, str.length());
          if(ss=="MLOGA")
              p.problem = mloga;
          else if(ss=="KPART")
              p.problem = kpart;
          else if(ss=="NETGEN")
              p.problem = netgen;
	  else if(ss=="MINLA")
	    p.problem = minla;
		  else {
			cerr << "bad problem type" << endl;
			exit(-1);
		  }
        }
 
      else if(str(0,18) == "[ALG_DISTANCE_TYPE]")
        {
          leda::string ss = str(19, str.length());
          if(ss=="JACOBIAN")
              p.alg_dist = jacobian;
	  else if(ss=="GAUSS-SEIDEL")
              p.alg_dist = gauss_seidel;
	else if(ss=="MIX_ALG_WAG")
		p.alg_dist = wag_and_gs;
	else if(ss=="MM")
		p.alg_dist = maximum_matching;
          else
              p.alg_dist = noalgdist;
        }
      else if(str(0,19) == "[FUTURE_VOLUME_SORT]")
        {
          leda::string ss = str(20, str.length());
          if(ss=="WD")
              p.future_volume_sort = wd;
	else if(ss=="MIX_ALG_WAG")
              p.future_volume_sort = mix_alg_wag;
	  else if(ss=="ALG")
              p.future_volume_sort = alg;
          else
              p.future_volume_sort = nosort;
        }
        else if(str(0,13) == "[IW_CALC_TYPE]")
        {
          leda::string ss = str(14, str.length());
          if(ss=="BAMG")
              p.iw_calc_type = iwbamg;
	  else if(ss=="ALG")
              p.iw_calc_type = iw_alg;
	  else if(ss=="MIX_ALG_WAG")
              p.iw_calc_type = iw_alg_and_wag;
          else
              p.iw_calc_type = iw_wag;
        }


      else if(str(0,16) == "[USE_SA_IN_RELAX]")
        {
          leda::string ss = str(17, str.length());
          if(ss=="YES")
              p.use_sa_in_relax = true;
          else
              p.use_sa_in_relax = false;
        }
      else if(str(0,16) == "[USE_SA_IN_LAST2]")
        {
          leda::string ss = str(17, str.length());
          if(ss=="YES")
              p.use_sa_only_in_last_two_levels = true;
          else
              p.use_sa_only_in_last_two_levels = false;
        }

      else if(str(0,17) == "[CHECK_BEST_RELAX]")
        {

          leda::string ss = str(18, str.length());
          if(ss=="YES")
              p.check_best_relax = true;
          else
              p.check_best_relax = false;
        }
      else if(str(0,11) == "[REVERSE_SA]")
        {

          leda::string ss = str(12, str.length());
          if(ss=="YES")
              p.reverse_sa = true;
          else
              p.reverse_sa = false;
        }
      else if(str(0,15) == "[RANDOM_REORDER]")
        {

          leda::string ss = str(16, str.length());
          if(ss=="YES")
              p.random_reorder = true;
          else
              p.random_reorder = false;
        }
      else if(str(0,15) == "[INDIVIDUAL_MIN]")
        {

          leda::string ss = str(16, str.length());
          if(ss=="1C") // 1 cycle individual
            { p.individual_minimization = true; p.many_cycles_in_min = false; }
          else if(ss=="MC") // many cycles individual
            { p.individual_minimization = true; p.many_cycles_in_min = true; }
          else // regular minimization
            { p.individual_minimization = false;p.many_cycles_in_min = true; }
        }

     else if(str(0,3) == "[SA]")
        {
          leda::string ss = str(4, str.length());
          p.do_sa_1 = char_search(ss, '1');
          p.do_sa_2 = char_search(ss, '2');
          p.do_sa_3 = char_search(ss, '3');
        }

   }

  fl.close();
  srand(Params.seed);
  return all_params;
}
