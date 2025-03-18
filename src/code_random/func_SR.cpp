#include <iostream>
#include <cstdlib>
#include <cmath>
#include "func_cost.h"
#include "func_SR.h"
using namespace std;

void do_sr(int &pnedges_undir, int &pngraver, vector<int> &pnedges_graver_ij, vector<vector<int>> &pgraver_ij_nonzero_ids, vector<vector<int>> &pgraver_ij, vector<int> &pfwalk_sol, vector<int> &pcur_sol, double &pcur_time, double &pcur_time_evac, vector<double> &pcur_time_edges, vector<double> &pcur_time_evac_edges, vector<int> &pcur_fr_bool, vector<double> &pedge_sol_lookup, vector<double> &pedge_sol_lookup_evac, vector<double> &pcur_cij_undir, vector<double> &ptij_undir, vector<double> &pfwalk_time_edges, vector<double> &pfwalk_time_evac_edges, double &palpha, double &pbeta, double &pt1, int &pn_counter, int &pSR_aug_counter)
     {
      double pcur_time_prev_1, pcur_time_prev_2, tol_12;
      pcur_time_prev_1 = 0.0;
      pcur_time_prev_2 = 0.0;
      tol_12 = 1.0;

      int pchk_bool, pcc, pnc, pnc_bool, pe_label;
      pchk_bool = 1;
      pcc = 0;
      double pfwalk_time, pfwalk_time_evac;
      int pj, pg_id;
      pfwalk_time = 0;
      pfwalk_time_evac = 0;
      pj = 0;
      pg_id = 0;
      while (pchk_bool >= 1)
         {
          pnc = 0;
          for (int pi=0; pi < pngraver; pi++)
              {
               pnc_bool = 1;
               pj = 0;
               pg_id = 0;
               pfwalk_time = 0;
               pfwalk_time_evac = 0;
               pfwalk_time += pcur_time;
               pfwalk_time_evac += pcur_time_evac;
              
               while ((pnc_bool >= 1) && (pg_id < (pnedges_graver_ij)[pi]))
                   {
                    pj = ((pgraver_ij_nonzero_ids)[pi])[pg_id];
                    (pfwalk_sol)[pj] = (pcur_sol)[pj] + ((pgraver_ij)[pi])[pj];
                    if ((pfwalk_sol)[pj] < 0)
                       {
                        pnc_bool = 0;
                       }
                    else
                       {
                        pfwalk_time -= (pcur_time_edges)[pj];
                        pfwalk_time_evac -= (pcur_time_evac_edges)[pj];
                        pe_label = ((pfwalk_sol)[pj])*(pnedges_undir)*2 + (pcur_fr_bool)[pj]*(pnedges_undir) + pj;
                        
                        if ((pedge_sol_lookup)[pe_label] == -1)
                           {
                             t_func_int_j(pfwalk_sol, pcur_cij_undir, ptij_undir, pj, pfwalk_time, pfwalk_time_edges, pfwalk_time_evac, pfwalk_time_evac_edges, palpha, pbeta, pt1, pn_counter);
                           
                            (pedge_sol_lookup)[pe_label] = (pfwalk_time_edges)[pj];
                            (pedge_sol_lookup_evac)[pe_label] = (pfwalk_time_evac_edges)[pj];
                           }
                        else
                           {
                            (pfwalk_time_edges)[pj] = (pedge_sol_lookup)[pe_label];
                            pfwalk_time += (pfwalk_time_edges)[pj];
                            (pfwalk_time_evac_edges)[pj] = (pedge_sol_lookup_evac)[pe_label];
                            pfwalk_time_evac += (pfwalk_time_evac_edges)[pj];
                           } 
                       } 
                    pg_id += 1;
                   }

               pj = 0;
               if (pnc_bool >= 1)
                  {
                   if (pfwalk_time < pcur_time)
                      {
                       for (pg_id=0; pg_id < (pnedges_graver_ij)[pi]; pg_id ++)
                           {
                            pj = ((pgraver_ij_nonzero_ids)[pi])[pg_id];
                            (pcur_sol)[pj] = (pfwalk_sol)[pj];
                            (pcur_time_edges)[pj] = (pfwalk_time_edges)[pj];
                            (pcur_time_evac_edges)[pj] = (pfwalk_time_evac_edges)[pj];
                           }
                        pcur_time = pfwalk_time;
                        pcur_time_evac = pfwalk_time_evac;
			pSR_aug_counter += 1;
                      }
                   else
                      {
                       pnc += 1;
                      }
                  }
               else
                  {
                   pnc += 1;
                  }
              }

          tol_12 = (pcur_time_prev_1/(pcur_time)) - 1.0;
          pcur_time_prev_2 = pcur_time_prev_1;
          pcur_time_prev_1 = pcur_time;

          if (pnc >= pngraver)
             {
              pchk_bool = 0;
             }
          pcc += 1;
         }
    }
