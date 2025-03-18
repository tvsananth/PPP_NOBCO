#include "selfish_routing_leblanc.h"
#include "func_leblanc.h"
using namespace std;

void do_sr_leblanc(SelfishRoutingSolver &solver1, NetworkDesignProblem &problem1, int &pnedges_undir, vector<int> &pcur_fr_bool, double &pcur_time, double &pcur_time_evac, vector<double> &pcur_sol_sr_double)
     {
      int pnedges_dir;
      pnedges_dir = (pnedges_undir)/2;
      
      std::set<std::pair<int,int> > cur_FrEdges;
      for (int i1=0; i1< pnedges_dir; i1++)
         {
          if (((pcur_fr_bool)[i1] == 1) || ((pcur_fr_bool)[i1 + pnedges_dir] == 1))
             {
              cur_FrEdges.insert(problem1.edgeIndices[i1]);
             }
         }
      solver1.solve(cur_FrEdges, FunctionType::USER_OPTIMAL, pcur_time, pcur_time_evac, pcur_sol_sr_double);
     }


void do_sr_leblanc_bi(SelfishRoutingSolver &solver1, NetworkDesignProblem &problem1, int &pnedges_undir, vector<int> &pcur_fr_bool, double &pcur_time, double &pcur_time_evac, vector<double> &pcur_sol_sr_double)
     {
      int pnedges_dir;
      pnedges_dir = (pnedges_undir)/2;

      std::set<std::pair<int,int> > cur_FrEdges;
      for (int i1=0; i1< pnedges_dir; i1++)
         {
          if (((pcur_fr_bool)[2*i1] == 1) || ((pcur_fr_bool)[2*i1 + 1] == 1))
             {
              cur_FrEdges.insert(problem1.edgeIndices[i1]);
             }
         }
      solver1.solve(cur_FrEdges, FunctionType::USER_OPTIMAL, pcur_time, pcur_time_evac, pcur_sol_sr_double);
     }

