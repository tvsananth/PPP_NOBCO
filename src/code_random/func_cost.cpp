#include <cmath>
#include "func_cost.h"

void t_func_int(vector<int> &x_sol, vector<double> &cij_undir, vector<double> &tij_undir, int &nedges_undir, double &tsol, vector<double> &tsol_edges, double &tsol_evac, vector<double> &tsol_evac_edges)
     {
      tsol = 0;
      tsol_evac = 0;
      double alpha=0.15;
      double beta=4.0;
      double t1 = 0;
      for (int i=0; i< nedges_undir; i++)
          {
           t1 = alpha*(pow((x_sol)[i]/(cij_undir)[i], beta));
           (tsol_edges)[i] = (tij_undir)[i]*((x_sol)[i])*(1 + t1/(beta + 1.0));
           tsol += (tsol_edges)[i];
           (tsol_evac_edges)[i] = (tij_undir)[i]*((x_sol)[i])*(1 + t1);
           tsol_evac += (tsol_evac_edges)[i];
          }
     }


void t_func_int_j(vector<int> &x_sol, vector<double> &cij_undir, vector<double> &tij_undir, int &j, double &tsol, vector<double> &tsol_edges, double &tsol_evac, vector<double> &tsol_evac_edges, double &alpha, double &beta, double &t1, int &counter)
     {
      t1 = (alpha)*(pow((x_sol)[j]/(cij_undir)[j], (beta)));
      (tsol_edges)[j] = (tij_undir)[j]*((x_sol)[j])*(1 + t1/((beta) + 1.0));
      tsol += (tsol_edges)[j];
      (tsol_evac_edges)[j] = (tij_undir)[j]*((x_sol)[j])*(1 + t1);
      tsol_evac += (tsol_evac_edges)[j];
      t1 = 0;
      counter += 1;
     }

//
void t_func_double(vector<double> &x_sol, vector<double> &cij_undir, vector<double> &tij_undir, int &nedges_undir, double &tsol, vector<double> &tsol_edges, double &tsol_evac, vector<double> &tsol_evac_edges)
     {
      tsol = 0;
      tsol_evac = 0;
      double alpha=0.15;
      double beta=4.0;
      double t1 = 0;
      for (int i=0; i< nedges_undir; i++)
          {
           t1 = alpha*(pow((x_sol)[i]/(cij_undir)[i], beta));
           (tsol_edges)[i] = (tij_undir)[i]*((x_sol)[i])*(1 + t1/(beta + 1.0));
           tsol += (tsol_edges)[i];
           (tsol_evac_edges)[i] = (tij_undir)[i]*((x_sol)[i])*(1 + t1);
           tsol_evac += (tsol_evac_edges)[i];
          }
     }

//
void t_func_double_j(vector<double> &x_sol, vector<double> &cij_undir, vector<double> &tij_undir, int &j, double &tsol, vector<double> &tsol_edges, double &tsol_evac, vector<double> &tsol_evac_edges, double &alpha, double &beta, double &t1, int &counter)
     {
      t1 = (alpha)*(pow((x_sol)[j]/(cij_undir)[j], (beta)));
      (tsol_edges)[j] = (tij_undir)[j]*((x_sol)[j])*(1 + t1/((beta) + 1.0));
      tsol += (tsol_edges)[j];
      (tsol_evac_edges)[j] = (tij_undir)[j]*((x_sol)[j])*(1 + t1);
      tsol_evac += (tsol_evac_edges)[j];
      t1 = 0;
      counter += 1;
     }

