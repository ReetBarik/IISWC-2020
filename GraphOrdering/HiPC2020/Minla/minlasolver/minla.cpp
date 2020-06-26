#include "cmpfuncs.h"
#include <time.h>

void Perm (int n, int NN,
           array<int> & p,
           array<int> & pi,
           array<int> & dir,
           list<list<int> > & L)
{ 
  int i, z, j, jd;

  if (n > NN)
    {
      list<int> l_int;      
      for (i=1; i <= NN; ++i)
        l_int.push_back(p[i]);
      L.push_back(l_int);
    }
  else {
       Perm (n+1, NN, p, pi, dir, L);
       for (i=1; i <= n-1; ++i)
         {
           j = (int)pi[n];
           jd = (int)j+dir[n];
           z = (int)p[jd];
           p[j] = z;
           p[jd] = n;
           pi[z] = j;
           pi[n] = jd;
           Perm (n+1, NN, p, pi, dir, L); 
           }
       dir[n] = -dir[n];
  }
}

void AllPerm(list<list<int> > & L, list<int> LR, list<int> C)
{
  int i, j;
  list_item it;

  if(LR.size()==0)
    {
      //      list_item k;
      //      forall_items(k, C)
      //        cerr << C[k] << ", ";
      //      cerr << endl;
      L.push_back(C);
    }
  else
    {
      forall_items(it, LR)
        {
          //          cerr << LR[it] << endl;
          C.push_back(LR[it]);
          list<int> NewLR = LR;
          NewLR.remove(LR[it]);
          AllPerm(L, NewLR, C);
          C.del(C.last());
        }
    }
}

array<list<Pair_LOfInt_Double> >  basis_solutions;

class res_pair {
public:
  double lin_arr;
  list_item ptr;

  friend ostream& operator<<(ostream& o, const res_pair& s)
  { 
    return o;
  }
  
  friend istream& operator>>(istream& i, res_pair& s)
  {
    return i;
  }
};

static int cmp_pair_res(const res_pair& a, const res_pair & b)
{
  if(a.lin_arr < b.lin_arr)
    return -1;
  else
    return 1;
}

double getdiff(Pair_LOfInt_Double basis, Pair_LOfInt_Double other, double total_segment_len)
{
  int i,j,k;
  double diff=0;
  double reverse_diff=0;

  //  double total_segment_len = ;
  
  
  list_item itb = basis.L.first_item();
  for(i=0; i<basis.L.size(); i++)
    {      
      list_item ito = other.L.first_item();
      for(j=0; j<other.L.size(); j++)
        {
          if(basis.L[itb]==other.L[ito])
            {
              diff+=fabs(basis.coords[i]-other.coords[j]);
              reverse_diff+=fabs(total_segment_len-basis.coords[i]+1-other.coords[j]);
            }
          ito = other.L.succ(ito);
        }
      itb = basis.L.succ(itb);
    }


  diff = fabs(basis.c - other.c)/mymin(diff,reverse_diff);
  
  return diff;
}

double find_min_laC(TGraphC & G)
{
  node Q;
  double total_segment_len = 0; // was Params.min_graph_size
  forall_nodes(Q, G)
    total_segment_len += G[Q].w;
    
  //  cerr << "Solving exact problem" << endl;
  int i;
  list<list<int> > L;
  L.clear();
  int n = G.number_of_nodes();

  //  cerr << "\nPermutations of " << n << endl;
  list<int> C; C.clear();
  list<int> LR;
  for(i=1; i<n+1; i++)
    LR.push_back(i);
  AllPerm(L, LR, C);
  
  //  cout << "generated" << endl;
  //  G.write();
  //   cerr << "List generated size=" << L.size() << endl;
  //   exit(1);

  
  LA MinLA;
  MinLA.weight = 99999999;
  MinLA.V.init(G);

  //  array<res_pair> results(L.size());
  //  i=0;
  list_item it, it2;
  forall_items(it, L)
    {
      
      it2 = L[it].first();
      
      node v;
      //      cerr << "Perm: ";
      TMP_CMP_GRAPHC = &G;  
      G.sort_nodes(&cmp_initial_id);
      
      forall_nodes(v, G)
        {
          //          cerr << L[it][it2] << " ";
          G[v].ArrId = L[it][it2];
          
          it2 = L[it].succ(it2);
         
        }

      TMP_CMP_GRAPHC = &G;  
      G.sort_nodes(&cmp_ArrId);
      define_S_values(G);
      /*
      cerr << "BASIS : " ;
      forall_nodes(v, G)
        cerr << "(" << G[v].S_value << ", " << G[v].w << ")" << endl;
      cerr << endl;
      */
      double w = calc_laC(G);
      //      cerr << "----------------- " << w << endl;
      //      results[i].lin_arr = w;
      //      results[i].ptr = it;
      Pair_LOfInt_Double P;
      P.coords.resize(Params.min_graph_size);
      P.c = w;
      P.L = L[it];
      int pos = 0;
      forall_nodes(v, G)
        {
          P.coords[pos] = G[v].S_value;
          pos++;
        }
      
        
      basis_solutions[Params.CURRENT_LEVEL].push_back(P);
      
      if(MinLA.weight > w)
        {
          forall_nodes(v, G)
            MinLA.V[v] = G[v].ArrId;
          MinLA.weight = w;
        }
    }
  /*
  cerr << "-----------------------" << endl;
  node v;
  forall_nodes(v, G)
      G[v].ArrId = MinLA.V[v];
  
  TMP_CMP_GRAPHC = &G;
  G.sort_nodes(&cmp_nC);
  cerr << "<<<<" << endl;
  forall_nodes(v, G)
    cerr << G[v].initial_id << "; ";
  cerr << ">>>>" << endl;
  define_S_values(G);
  cerr << calc_laC(G) << endl;
  cerr << "-----------------------" << endl;
  list<int> a; a.clear();
  a.push_back(3);
  a.push_back(2);
  a.push_back(1);
  a.push_back(4);
  
  put_arrangement_into_G(G, a);
  define_S_values(G);
  cerr << calc_laC(G) << endl;
  */
  basis_solutions[Params.CURRENT_LEVEL].sort(&cmp_Pair_LOfInt_Double);

  //  cerr << basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].first()].c << endl;
  //  cerr << basis_solutions[Params.CURRENT_LEVEL][basis_solutions[Params.CURRENT_LEVEL].succ(basis_solutions[Params.CURRENT_LEVEL].first())].c << endl;

  //  cerr << basis_solutions[Params.CURRENT_LEVEL].size() << "\t" << Params.basis_solutions_number<<endl;
  if(Params.basis_solutions_number > basis_solutions[Params.CURRENT_LEVEL].size())
    Params.basis_solutions_number = basis_solutions[Params.CURRENT_LEVEL].size();
  else
    {
      list<Pair_LOfInt_Double> tmp_list;
      
      if(Params.basis_solutions_type==2)
        {
          // Take 10 first with step 100      
          
          list<Pair_LOfInt_Double> tmp_list;
          list_item it = basis_solutions[Params.CURRENT_LEVEL].first_item();
          
          for(int k=0; k<Params.basis_solutions_number; k++)
            {          
              tmp_list.push_back(basis_solutions[Params.CURRENT_LEVEL][it]);
              
              for(int l=0; l<100; l++)
                it = basis_solutions[Params.CURRENT_LEVEL].succ(it);            
            }
          basis_solutions[Params.CURRENT_LEVEL] = tmp_list;
          
        }
      else if(Params.basis_solutions_type==1)
        {
          // Take 10 most different from 1st solution
          
          
          // it = 3rd item
          list_item it = basis_solutions[Params.CURRENT_LEVEL].first_item();
          list_item first_it = it;
          it = basis_solutions[Params.CURRENT_LEVEL].succ(it);   it = basis_solutions[Params.CURRENT_LEVEL].succ(it);

          //          cerr << "ogogo : " << basis_solutions[Params.CURRENT_LEVEL].size() << endl;
          
          for(int k=0; k<basis_solutions[Params.CURRENT_LEVEL].size()/2-10; k++)
            {
              double diff = getdiff(basis_solutions[Params.CURRENT_LEVEL][first_it],basis_solutions[Params.CURRENT_LEVEL][it], total_segment_len);
              Pair_LOfInt_Double basis;
              basis.L = basis_solutions[Params.CURRENT_LEVEL][it].L;
              basis.c = diff;
              tmp_list.push_back(basis);
              
              it = basis_solutions[Params.CURRENT_LEVEL].succ(it); // jump over reverse of the order
              it = basis_solutions[Params.CURRENT_LEVEL].succ(it);   
            }
          tmp_list.sort(&cmp_Pair_LOfInt_Double);
          basis_solutions[Params.CURRENT_LEVEL] = tmp_list;
          if(Params.basis_solutions_number == 1)
            {
              while(basis_solutions[Params.CURRENT_LEVEL].size()>Params.solution_number)
                basis_solutions[Params.CURRENT_LEVEL].erase(basis_solutions[Params.CURRENT_LEVEL].last());
              while(basis_solutions[Params.CURRENT_LEVEL].size()>1)
                basis_solutions[Params.CURRENT_LEVEL].erase(basis_solutions[Params.CURRENT_LEVEL].first());
            }
          else
            {
              while(basis_solutions[Params.CURRENT_LEVEL].size()>Params.basis_solutions_number)
                basis_solutions[Params.CURRENT_LEVEL].erase(basis_solutions[Params.CURRENT_LEVEL].last());
            }
        }
      else if(Params.basis_solutions_type==0)
        {
          // Take 10 first solutions
          while(basis_solutions[Params.CURRENT_LEVEL].size()>Params.basis_solutions_number)
            basis_solutions[Params.CURRENT_LEVEL].erase(basis_solutions[Params.CURRENT_LEVEL].last());
        }
    }

  cerr << "Remained bs = " << basis_solutions[Params.CURRENT_LEVEL].size() << endl;
  return MinLA.weight;
}

double solve_exact_minla_C(TGraphC & G)
{
  //  cerr << Params.basis_solutions_number << endl;
  //  exit(1);
  node v;
  
  int i = 1;
  forall_nodes(v, G)
    {
      G[v].ArrId=i;
      i++;
    }

  long t = time(0);
  double min_la = find_min_laC(G);
  cerr << "Basis solving time (s) = " << (double)(time(0)-t)/1.0 << endl;
  cerr << "Solving exactly at S=" << Params.CURRENT_LEVEL << endl;
}
