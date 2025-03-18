#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>

#include <vector>

#include "func_cost.h"
#include "func_SR.h"
#include "selfish_routing_leblanc.h"
#include "func_leblanc.h"
using namespace std;


void get_fsize(ifstream &, streampos &);
void get_paths(std::string &, std::string &, std::string &, std::string &, std::string &, std::string &, std::string &, std::string &);
void get_paths_ini_sol_fr(std::string &, std::string &, int &, std::string &);
void get_nedges(std::string &, int &, int &);
void get_ngraver(std::string &, int &, int &);
void load_data_txt(std::string &, vector<double> &, int &);
void load_data_txt_int(std::string &, vector<int> &, int &);
void load_data_bin_int(std::string &, vector<int> &, int &, int &);


void get_paths(std::string &fname_graver, std::string &fname_graver_sr, std::string &fname_cij, std::string &fname_tij, std::string &fname_ini_sol, std::string &fname_ini_fr, std::string &fname_lanes, std::string &path_dir_name){
     fname_graver = path_dir_name + "/fsol_graver_r.npy";
     fname_graver_sr = path_dir_name + "/fsol_graver_r.npy";
     fname_cij = path_dir_name + "/cij_edges_undir.txt";
     fname_tij = path_dir_name + "/tij_edges_undir.txt";
     fname_ini_sol = path_dir_name + "/ini_sol_rand_testfr.npy";
     fname_ini_fr = path_dir_name + "/ini_fr_rand_testfr.npy";
     fname_lanes = path_dir_name + "/lanes_edges_undir.txt";

     cout << "filename : " << fname_graver << " ... graver_fr ..."<< "\n";
     cout << "filename : " << fname_graver_sr << " ... graver_evacuee ..."<< "\n";
     }


void get_paths_ini_sol_fr(string &fname_ini_sol, string &fname_ini_fr, int &id_sol, std::string &path_dir_name){
     char buffer1[200];
     sprintf(buffer1, "/rand_ini_sols_dir/ini_sol_rand_testfr%d.npy", id_sol);
     fname_ini_sol = path_dir_name + buffer1;
     cout << "filename : " << fname_ini_sol << "\n";
     char buffer2[200];
     sprintf(buffer2, "/rand_ini_fr_dir/ini_fr_rand_testfr%d.npy", id_sol);
     fname_ini_fr = path_dir_name + buffer2;
     cout << "filename : " << fname_ini_fr << "\n";
     cout << "done0" << '\n';
     }

 
void get_nedges(std::string &fname_cij, int &nedges_dir, int &nedges_undir){
     std::string line_cij;
     ifstream file_cij(fname_cij);
     if (file_cij.is_open())
       {
        while(getline(file_cij, line_cij))
           {
             nedges_undir += 1;
           }
        file_cij.close(); 
       } 
     nedges_dir = (nedges_undir)/2;
     }


void get_fsize(ifstream &file, streampos &size){
     streampos file_ini, file_fin;
     file_ini = (file).tellg();
     (file).seekg(0, ios::end);
     file_fin = (file).tellg();
     (file).close();
     size = file_fin - file_ini;
     }


void load_data_txt(std::string &fname_cij, vector<double> &data, int &ndata){
     int count = 0;
     std::string line_cij;
     ifstream file_cij(fname_cij);
     if (file_cij.is_open())
        {
         while(getline(file_cij, line_cij))
            {
             (data)[count] = std::stod(line_cij);
             count += 1;
            }
         file_cij.close();              
        }
     }


void load_data_txt_int(std::string &fname_cij, vector<int> &data, int &ndata){
     int count = 0;
     std::string line_cij;
     ifstream file_cij(fname_cij);
     if (file_cij.is_open())
        {
         while(getline(file_cij, line_cij))
            {
             (data)[count] = std::stoi(line_cij);
             count += 1;
            }
         file_cij.close();
        }
     }

void load_data_bin_int(std::string &fname_cij, vector<int> &data, int &ndata, int &rand_index){
     if (rand_index == 0)
       {
        data = vector<int>(ndata, 0);
       }
     ifstream file_cij(fname_cij, ios::binary);
     int d_i = 0;
     for (int i=0; i < ndata; i++)
        {
         file_cij.read((char*) &d_i, sizeof(int));
         (data)[i] = d_i;
        }
     }


void load_data_graver(std::string &fname_graver, vector<vector<int>> &data_g, int &ngraver, int &nedges_undir){

     data_g = vector<vector<int>>(ngraver);

     for (int i=0; i < ngraver; i++)
        {
         (data_g)[i] = vector<int>(nedges_undir, 0);
        }
     ifstream file_graver(fname_graver, ios::binary);
     int g_ij = 0;
     for (int i=0; i < ngraver; i++)
        {
         for (int j=0; j< nedges_undir; j++)
            {
             file_graver.read((char*) &g_ij, sizeof(int));
             ((data_g)[i])[j] = g_ij;
            }
        }
     }


void get_ngraver(std::string &fname_graver, int &ngraver, int &nedges_undir){
     streampos file_size;
     ifstream file_graver(fname_graver, ios::binary);
     get_fsize(file_graver, file_size);
     ngraver = file_size/sizeof(int)/(nedges_undir);
     cout << file_size << " " << sizeof(int) << " " << ngraver << '\n';
     }



//Function to return vector with id's of correspoding opposite edge
void opposite_edges(vector<int> &id_e_opp, int &nedges_dir){
    for (int i=0; i<nedges_dir; i++)
        {
         id_e_opp[i] = i + nedges_dir;
         id_e_opp[i + nedges_dir] = i;
        }
}

//
void initialize_graver_nedges(int &nedges_undir, int &ngraver, vector<vector<int>> &graver_ij, vector<int> &nedges_graver_ij){
    for (int i=0; i<ngraver; i++)
       {
        for (int j=0; j<nedges_undir; j++)
           {
            if (graver_ij[i][j] != 0)
              {
               nedges_graver_ij[i] += 1;
              }
           }
       }

    }

//
void initialize_graver_nonzero_ids(int &nedges_undir, int &ngraver, vector<vector<int>> &graver_ij, vector<int> &nedges_graver_ij, vector<vector<int>> &graver_ij_nonzero_ids, std::chrono::duration<double> &elapsed){

    for (int i=0; i<ngraver; i++)
       {
        graver_ij_nonzero_ids[i] = vector<int>(nedges_graver_ij[i], 0);
       }
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    int g_count;
    for (int i=0; i<ngraver; i++)
       {
        g_count = 0;
        for (int j=0; j<nedges_undir; j++)
           {
            if (graver_ij[i][j] != 0)
              {
               graver_ij_nonzero_ids[i][g_count] = j;
               g_count += 1;
              }
           }
       }
}

//
void initialize_FR_SR_data(int &nedges_dir, int &nedges_undir, vector<int> &id_e_opp, double &cur_time, double &cur_time_evac, double &cur_time_sr, double &cur_time_evac_sr, double &fwalk_time, double &fwalk_time_evac, vector<double> &cij_undir, vector<int> &lanes_undir, vector<double> &cij_undir_fr, vector<int> &cur_fr, vector<int> &cur_fr_bool, vector<double> &cur_cij_undir, vector<int> &ini_fr, vector<int> &ini_sol){
    cur_time = 0.0;
    fwalk_time = 0.0;

    cur_time_evac = 0.0;
    fwalk_time_evac = 0.0;

    cur_time_sr = 0.0;
    cur_time_evac_sr = 0.0;

    for (int i=0; i<nedges_undir; i++){cij_undir_fr[i] = cij_undir[i]*(lanes_undir[i] - 1.0)/lanes_undir[i];}
    for (int i=0; i<nedges_undir; i++)
            {
             cur_fr[i] = ini_fr[i];
            }
    for (int i=0; i<nedges_dir; i++)
            {
             if ((cur_fr[i] >= 1) || (cur_fr[id_e_opp[i]] >= 1))
                {
                 cur_fr_bool[i] = 1;
                 cur_fr_bool[id_e_opp[i]] = 1;
                }
             else
               {
                cur_fr_bool[i] = 0;
                cur_fr_bool[id_e_opp[i]] = 0;
               }
            }

    for (int i=0; i<nedges_dir; i++)
            {
             if ((cur_fr[i] >= 1) || (cur_fr[id_e_opp[i]] >= 1))
                {
                 cur_cij_undir[i] = cij_undir_fr[i];
                 cur_cij_undir[id_e_opp[i]] = cij_undir_fr[id_e_opp[i]];
                }
             else
                {
                 cur_cij_undir[i] = cij_undir[i];
                 cur_cij_undir[id_e_opp[i]] = cij_undir[id_e_opp[i]];
                }
            }
}


//Function checking FR convergence conditions
void is_FR_graver_step_feasible(int &nc_bool_fr, int &g_id_fr, int &nedges_undir, vector<int> &nedges_graver_ij, vector<vector<int>> &graver_ij_nonzero_ids, vector<int> &fwalk_fr, vector<int> &cur_fr, vector<int> &cur_sol, vector<vector<int>> &graver_ij, vector<double> &cij_undir_fwalk, vector<double> &cij_undir_fr, vector<double> &cij_undir, vector<double> &tij_undir, vector<int> &id_e_opp, vector<double> &cur_time_edges, vector<double> &cur_time_evac_edges, double &fwalk_time, vector<double> &fwalk_time_edges, double &fwalk_time_evac, vector<double> &fwalk_time_evac_edges, int &e_label, vector<double> &edge_sol_lookup, vector<double> &edge_sol_lookup_evac, int &fwalk_fr_bool, int &i_fr, int &j, int &j1, double &alpha, double &beta, double &t1, int &n_counter){

          while ((nc_bool_fr >= 1) && (g_id_fr < nedges_graver_ij[i_fr]))
                  {
                   j = graver_ij_nonzero_ids[i_fr][g_id_fr];
                   fwalk_fr[j] = cur_fr[j] + graver_ij[i_fr][j];
                   fwalk_fr[id_e_opp[j]] = cur_fr[id_e_opp[j]] + graver_ij[i_fr][id_e_opp[j]];
                   if ((fwalk_fr[j] >= 1) || (fwalk_fr[id_e_opp[j]] >= 1))
                     {
                      cij_undir_fwalk[j] = cij_undir_fr[j];
                      cij_undir_fwalk[id_e_opp[j]] = cij_undir_fr[id_e_opp[j]];
                     }
                   else
                     {
                      cij_undir_fwalk[j] = cij_undir[j];
                      cij_undir_fwalk[id_e_opp[j]] = cij_undir[id_e_opp[j]];
                     }
                   if (fwalk_fr[j] < 0)
                      {
                       nc_bool_fr = 0;
                      }
                   g_id_fr += 1;
                   }

          g_id_fr = 0;
          fwalk_fr_bool = 0;
          while ((nc_bool_fr >= 1) && (g_id_fr < nedges_graver_ij[i_fr]))
                  {
                   j = graver_ij_nonzero_ids[i_fr][g_id_fr];
                   if ((fwalk_fr[j] >= 1) || (fwalk_fr[id_e_opp[j]] >= 1))
                     {fwalk_fr_bool = 1;}
                   else
                     {fwalk_fr_bool = 0;}

                   if (fwalk_fr[j] < 0)
                      {
                       nc_bool_fr = 0;
                      }
                   else
                      {
                       fwalk_time -= cur_time_edges[j];
                       fwalk_time_evac -= cur_time_evac_edges[j];
                       e_label = cur_sol[j]*nedges_undir*2 + fwalk_fr_bool*nedges_undir + j;
                       if (edge_sol_lookup[e_label] == -1)
                          {
                           t_func_int_j(cur_sol, cij_undir_fwalk, tij_undir, j, fwalk_time, fwalk_time_edges, fwalk_time_evac, fwalk_time_evac_edges, alpha, beta, t1, n_counter);
                           edge_sol_lookup[e_label] = fwalk_time_edges[j];
                           edge_sol_lookup_evac[e_label] = fwalk_time_evac_edges[j];
                          }
                       else
                          {
                           fwalk_time_edges[j] = edge_sol_lookup[e_label];
                           fwalk_time += fwalk_time_edges[j];

                           fwalk_time_evac_edges[j] = edge_sol_lookup_evac[e_label];
                           fwalk_time_evac += fwalk_time_evac_edges[j];
                           }

                       j1 = id_e_opp[j];

                       fwalk_time -= cur_time_edges[j1];
                       fwalk_time_evac -= cur_time_evac_edges[j1];
                       e_label = cur_sol[j1]*nedges_undir*2 + fwalk_fr_bool*nedges_undir + j1;
                       if (edge_sol_lookup[e_label] == -1)
                          {
                           t_func_int_j(cur_sol, cij_undir_fwalk, tij_undir, j1, fwalk_time, fwalk_time_edges, fwalk_time_evac, fwalk_time_evac_edges, alpha, beta, t1, n_counter);
                           edge_sol_lookup[e_label] = fwalk_time_edges[j1];
                           edge_sol_lookup_evac[e_label] = fwalk_time_evac_edges[j1];
                          }
                       else
                          {
                           fwalk_time_edges[j1] = edge_sol_lookup[e_label];
                           fwalk_time += fwalk_time_edges[j1];

                           fwalk_time_evac_edges[j1] = edge_sol_lookup_evac[e_label];
                           fwalk_time_evac += fwalk_time_evac_edges[j1];
                           }

                      }                      
                   g_id_fr += 1;
                  }	      
}

//Function to initialize SR flow values and updated capacities in a graver step
void update_SR_initial_input(int &i_fr, int &nedges_dir, int &nedges_undir, vector<int> &id_e_opp, vector<vector<int>> &graver_ij, vector<int> &cur_sol, vector<int> &cur_fr, vector<double> &cur_time_edges, vector<double> &cur_time_evac_edges, vector<double> &fwalk_time_edges, vector<double> &fwalk_time_evac_edges, vector<int> &cur_sol_sr, vector<int> &cur_fr_sr, vector<double> &cur_time_edges_sr, vector<double> &cur_time_evac_edges_sr, vector<double> &cij_undir, vector<double> &cij_undir_fr, vector<double> &cij_undir_fwalk_sr, vector<int> &cur_fr_bool_sr, double &fwalk_time, double &fwalk_time_evac, double &cur_time_sr, double &cur_time_evac_sr){

             for (int j_sr=0; j_sr<nedges_undir; j_sr++)
                {
                 cur_sol_sr[j_sr] = cur_sol[j_sr];
                 cur_fr_sr[j_sr] = cur_fr[j_sr] + graver_ij[i_fr][j_sr];
                 cur_time_edges_sr[j_sr] = cur_time_edges[j_sr];
                 cur_time_evac_edges_sr[j_sr] = cur_time_evac_edges[j_sr];
                 }

             for (int j_sr=0; j_sr<nedges_undir; j_sr++)
                {
                 if (graver_ij[i_fr][j_sr] != 0)
                   {
                    cur_time_edges_sr[j_sr] = fwalk_time_edges[j_sr];
                    cur_time_edges_sr[id_e_opp[j_sr]] = fwalk_time_edges[id_e_opp[j_sr]];
                    cur_time_evac_edges_sr[j_sr] = fwalk_time_evac_edges[j_sr];
                    cur_time_evac_edges_sr[id_e_opp[j_sr]] = fwalk_time_evac_edges[id_e_opp[j_sr]];
                   }
                }

             for (int j_sr=0; j_sr<nedges_dir; j_sr++)
                {
                 if ((cur_fr_sr[j_sr] >= 1) || (cur_fr_sr[id_e_opp[j_sr]] >= 1))
                   {
                    cij_undir_fwalk_sr[j_sr] = cij_undir_fr[j_sr];
                    cij_undir_fwalk_sr[id_e_opp[j_sr]] = cij_undir_fr[id_e_opp[j_sr]];
                    cur_fr_bool_sr[j_sr] = 1;
                    cur_fr_bool_sr[id_e_opp[j_sr]] = 1;
                   }
                 else
                   {
                    cij_undir_fwalk_sr[j_sr] = cij_undir[j_sr];
                    cij_undir_fwalk_sr[id_e_opp[j_sr]] = cij_undir[id_e_opp[j_sr]];
                    cur_fr_bool_sr[j_sr] = 0;
                    cur_fr_bool_sr[id_e_opp[j_sr]] = 0;
                   }
            }
             cur_time_sr = fwalk_time;
             cur_time_evac_sr = fwalk_time_evac;
}

//Function to update current best FR path, SR flow values and evacuation times
void update_current_best_SR_soution(int &nedges_dir, int &nedges_undir, vector<int> &id_e_opp, vector<int> &cur_sol, vector<int> &cur_fr, vector<int> &cur_fr_bool, vector<double> &cur_time_edges, vector<double> &cur_time_evac_edges, vector<int> &cur_sol_sr, vector<int> &cur_fr_sr, vector<double> &cur_time_edges_sr, vector<double> &cur_time_evac_edges_sr, double &cur_time, double &cur_time_evac, double &cur_time_sr, double &cur_time_evac_sr){
             for (int j_sr=0; j_sr<nedges_undir; j_sr++)
                {
                 cur_sol[j_sr] = cur_sol_sr[j_sr];
                 cur_fr[j_sr] = cur_fr_sr[j_sr];
                 cur_time_edges[j_sr] = cur_time_edges_sr[j_sr];
                 cur_time_evac_edges[j_sr] = cur_time_evac_edges_sr[j_sr];
                }

             for (int j_sr=0; j_sr<nedges_dir; j_sr++)
                {
                 if ((cur_fr[j_sr] >= 1) || (cur_fr[id_e_opp[j_sr]] >= 1))
                    {
                     cur_fr_bool[j_sr] = 1;
                     cur_fr_bool[id_e_opp[j_sr]] = 1;
                    }
                 else
                    {
                     cur_fr_bool[j_sr] = 0;
                     cur_fr_bool[id_e_opp[j_sr]] = 0;
                    }
                }


             cur_time = cur_time_sr;
             cur_time_evac = cur_time_evac_sr;
}



void update_output_data(int &rand_ini_count, int &nedges_undir, vector<double> &sr_sol, double &cur_time_evac, vector<double> &time_sol, double &time_sol_gaga_rand_ini, vector<double> &sr_sol_ini, vector<int> &cur_fr, vector<int> &cur_sol, vector<double> &cij_undir_fwalk_sr, vector<double> &sr_sol_leblanc, double &cur_time_evac_leblanc, vector<double> &time_sol_leblanc, double &time_sol_leblanc_rand_ini, vector<double> &cur_sol_leblanc, vector<int> &ncount_aug_fr, vector<int> &ncount_aug_sr, int &FR_counter_aug, int &SR_counter_aug, std::string &path_dir_name){

    sr_sol[rand_ini_count] = cur_time_evac;
    time_sol[rand_ini_count] = time_sol_gaga_rand_ini;
    sr_sol_leblanc[rand_ini_count] = cur_time_evac_leblanc;
    time_sol_leblanc[rand_ini_count] = time_sol_leblanc_rand_ini;
    ncount_aug_fr[rand_ini_count] = FR_counter_aug;
    ncount_aug_sr[rand_ini_count] = SR_counter_aug;
    rand_ini_count += 1;


    ofstream file_write1(path_dir_name + "/Graver_walk/e21time_sol_srallFR_op1000.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write1.write((char *) &time_sol[i1], sizeof(double));
        }
    ofstream file_write2(path_dir_name + "/Graver_walk/e21sr_sol_srallFR_op1000.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write2.write((char *) &sr_sol[i1], sizeof(double));
        }
    ofstream file_write3(path_dir_name + "/Graver_walk/e21sr_sol_ini_srallFR_op1000.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write3.write((char *) &sr_sol_ini[i1], sizeof(double));
        }

    ofstream file_write1l(path_dir_name + "/Graver_walk/e21sr_time_nogama_op1000.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write1l.write((char *) &time_sol_leblanc[i1], sizeof(double));
        }
    ofstream file_write2l(path_dir_name + "/Graver_walk/e21sr_sol_nogama_op1000.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write2l.write((char *) &sr_sol_leblanc[i1], sizeof(double));
        }
    
    ofstream file_write1_fr_aug(path_dir_name + "/Graver_walk/ncount_aug_FR.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write1_fr_aug.write((char *) &ncount_aug_fr[i1], sizeof(int));
        }

    ofstream file_write1_sr_aug(path_dir_name + "/Graver_walk/ncount_aug_SR.npy", ios::out | ios::binary);
    for (int i1=0; i1<1000; i1++)
        {
         file_write1_sr_aug.write((char *) &ncount_aug_sr[i1], sizeof(int));
        }

    std::string path_fname_fr;
    std::string path_fname_sr_sol, path_fname_cij_sr_sol;
    std::string rand_ini_count_str;
    rand_ini_count_str = std::to_string(rand_ini_count);
    path_fname_fr = path_dir_name + "/Graver_walk/FR_" + rand_ini_count_str + ".npy";
    path_fname_sr_sol = path_dir_name + "/Graver_walk/SR_" + rand_ini_count_str + ".npy";
    path_fname_cij_sr_sol = path_dir_name + "/Graver_walk/cijSR_" + rand_ini_count_str + ".npy";
    ofstream file_write4(path_fname_fr, ios::out | ios::binary);
    for (int j_sr=0; j_sr<nedges_undir; j_sr++)
        {
         file_write4.write((char *) &cur_fr[j_sr], sizeof(int));
        }
    ofstream file_write5(path_fname_sr_sol, ios::out | ios::binary);
    for (int j_sr=0; j_sr<nedges_undir; j_sr++)
        {
         file_write5.write((char *) &cur_sol[j_sr], sizeof(int));
        }
    ofstream file_write6(path_fname_cij_sr_sol, ios::out | ios::binary);
    for (int j_sr=0; j_sr<nedges_undir; j_sr++)
        {
         file_write6.write((char *) &cij_undir_fwalk_sr[j_sr], sizeof(double));
        }

    std::string path_fname_sr_leblanc_sol;
    path_fname_sr_leblanc_sol = path_dir_name + "/Graver_walk/SR_sol_double_nogama_" + rand_ini_count_str + ".npy";
    ofstream file_write5l(path_fname_sr_leblanc_sol, ios::out | ios::binary);
    for (int j_sr=0; j_sr<nedges_undir; j_sr++)
        {
         file_write5.write((char *) &cur_sol_leblanc[j_sr], sizeof(double));
        }
}


void check_nan_capacities(vector<double> &data_cij, int &ndata_cij){
	for (int i=0; i < ndata_cij; i++){
		if (data_cij[i] == 0){
			data_cij[i] = 1;
		}
	}
}


void do_main_FR(std::string path_dir_name, std::string nodes_fname_csv, std::string edges_fname_csv){

    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

  
    vector<vector<int>> graver_ij;
    double cur_time, fwalk_time, cur_time_sr;
    double cur_time_evac, fwalk_time_evac, cur_time_evac_sr;

    double cur_time_leblanc, cur_time_evac_leblanc; 

    std::string fname_graver, fname_graver_sr, fname_cij, fname_tij, fname_ini_sol, fname_ini_fr, fname_lanes;
    
    int nedges_dir, nedges_undir;
    nedges_dir = 0;
    nedges_undir = 0;

    int ngraver; 
    ngraver = 0;

    double alpha = 0.15;
    double beta = 4.0;
    double t1=0.0;

    int n_counter = 0;
    int FR_counter_aug = 0;
    int SR_counter_aug = 0;

    get_paths(fname_graver, fname_graver_sr, fname_cij, fname_tij, fname_ini_sol, fname_ini_fr, fname_lanes, path_dir_name);
    get_nedges(fname_cij, nedges_dir, nedges_undir);
    get_ngraver(fname_graver, ngraver, nedges_undir);
    cout << " nedges: " << nedges_dir << " nedges_undir: " << nedges_undir << " ngraver: " << ngraver << '\n';

    vector<int> id_e_opp(nedges_undir, 0);
    opposite_edges(id_e_opp, nedges_dir);

    vector<double> cij_undir(nedges_undir, 0);
    load_data_txt(fname_cij, cij_undir, nedges_undir);
    check_nan_capacities(cij_undir, nedges_undir);

    vector<double> tij_undir(nedges_undir, 0);
    load_data_txt(fname_tij, tij_undir, nedges_undir);

    int rand_index = 0;
    get_paths_ini_sol_fr(fname_ini_sol, fname_ini_fr, rand_index, path_dir_name);

    vector<int> ini_sol(nedges_undir, 0);
    load_data_bin_int(fname_ini_sol, ini_sol, nedges_undir, rand_index);

    vector<int> ini_fr(nedges_undir, 0);
    load_data_bin_int(fname_ini_fr, ini_fr, nedges_undir, rand_index);

    vector<int> lanes_undir(nedges_undir, 0);
    load_data_txt_int(fname_lanes, lanes_undir, nedges_undir);
    
    load_data_graver(fname_graver, graver_ij, ngraver, nedges_undir);

    int ngraver_sr;
    ngraver_sr = 0;
    get_ngraver(fname_graver_sr, ngraver_sr, nedges_undir);
    cout << " nedges: " << nedges_dir << " nedges_undir: " << nedges_undir << " ngraver_SR: " << ngraver_sr << '\n';
    
    vector<vector<int>> graver_sr_ij;
    load_data_graver(fname_graver_sr, graver_sr_ij, ngraver_sr, nedges_undir);

    vector<double> edge_sol_lookup(100000*nedges_undir*2, -1);
    vector<double> edge_sol_lookup_evac(100000*nedges_undir*2, -1);
    vector<int> cur_sol(nedges_undir, 0);
    vector<int> fwalk_sol(nedges_undir, 0);

    vector<double> cur_time_edges(nedges_undir, 0);
    vector<double> cur_time_evac_edges(nedges_undir, 0);
    vector<double> fwalk_time_edges(nedges_undir, 0);
    vector<double> fwalk_time_evac_edges(nedges_undir, 0);
    
    vector<int> cur_fr(nedges_undir, 0);
    vector<int> cur_fr_bool(nedges_undir, 0);
    vector<int> fwalk_fr(nedges_undir, 0);

    vector<double> cij_undir_fr(nedges_undir, 0);
    vector<double> cur_cij_undir(nedges_undir, 0);
    vector<double> cij_undir_fwalk(nedges_undir, 0);

    vector<int> cur_sol_sr(nedges_undir, 0);

    vector<double> cur_time_edges_sr(nedges_undir, 0);
    vector<double> cur_time_evac_edges_sr(nedges_undir, 0);

    vector<int> cur_fr_sr(nedges_undir, 0);
    vector<int> cur_fr_bool_sr(nedges_undir, 0);

    vector<double> cij_undir_fwalk_sr(nedges_undir, 0);

    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";


    vector<int> nedges_graver_ij(ngraver, 0);
    initialize_graver_nedges(nedges_undir, ngraver, graver_ij, nedges_graver_ij);
    
    vector<vector<int>> graver_ij_nonzero_ids(ngraver);
    initialize_graver_nonzero_ids(nedges_undir, ngraver, graver_ij, nedges_graver_ij, graver_ij_nonzero_ids, elapsed);
       
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    vector<int> nedges_graver_sr_ij(ngraver_sr, 0);
    initialize_graver_nedges(nedges_undir, ngraver_sr, graver_sr_ij, nedges_graver_sr_ij);

    vector<vector<int>> graver_sr_ij_nonzero_ids(ngraver_sr);
    initialize_graver_nonzero_ids(nedges_undir, ngraver_sr, graver_sr_ij, nedges_graver_sr_ij, graver_sr_ij_nonzero_ids, elapsed);   
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    int rand_ini_count = 0;
    vector<double> sr_sol_ini(1000, 0);
    vector<double> time_sol(1000, 0);
    vector<double> sr_sol(1000, 0);
    vector<double> time_sol_leblanc(1000, 0);
    vector<double> sr_sol_leblanc(1000, 0);
    vector<int> ncount_aug_fr(1000, 0);
    vector<int> ncount_aug_sr(1000, 0);

    double time_sol_0 = 0;
    
    
    while (rand_ini_count < 100)
    {

    n_counter = 0;
    FR_counter_aug = 0;
    SR_counter_aug = 0;
    
    for (int i=0; i<nedges_undir; i++){cur_sol[i] = ini_sol[i];}

    rand_index = rand_ini_count;
    get_paths_ini_sol_fr(fname_ini_sol, fname_ini_fr, rand_ini_count, path_dir_name);
    load_data_bin_int(fname_ini_sol, ini_sol, nedges_undir, rand_ini_count);
    load_data_bin_int(fname_ini_fr, ini_fr, nedges_undir, rand_ini_count);

    initialize_FR_SR_data(nedges_dir, nedges_undir, id_e_opp, cur_time, cur_time_evac, cur_time_sr, cur_time_evac_sr, fwalk_time, fwalk_time_evac, cij_undir, lanes_undir, cij_undir_fr, cur_fr, cur_fr_bool, cur_cij_undir, ini_fr, ini_sol);

    for (int i=0; i<nedges_undir; i++){cur_sol[i] = ini_sol[i];}

    t_func_int(cur_sol, cur_cij_undir, tij_undir, nedges_undir, cur_time, cur_time_edges, cur_time_evac, cur_time_evac_edges);
    std::cout.precision(16);

    cout << "ini " << cur_time << " " << cur_time_evac << '\n';
    do_sr(nedges_undir, ngraver_sr, nedges_graver_sr_ij, graver_sr_ij_nonzero_ids, graver_sr_ij, fwalk_sol, cur_sol, cur_time, cur_time_evac, cur_time_edges, cur_time_evac_edges, cur_fr_bool, edge_sol_lookup, edge_sol_lookup_evac, cur_cij_undir, tij_undir, fwalk_time_edges, fwalk_time_evac_edges, alpha, beta, t1, n_counter, SR_counter_aug);
    cout << "SR ini " << cur_time << " " << cur_time_evac << '\n';  

    int chk_bool = 1;
    int cc = 0;
    int nc = 0;
    int nc_bool = 1;
    int e_label = 0;

    int nc_fr = 0;
    int nc_bool_fr = 0;
     
    sr_sol_ini[rand_ini_count] = cur_time_evac;

    int chk_bool_fr = 1;
    
    while (chk_bool_fr >= 1)
    {
     nc_fr = 0;
     int g_id_fr, fwalk_fr_bool;
     
     for (int i_fr=0; i_fr<ngraver; i_fr++)
         {
          nc_bool_fr = 1;
          int j = 0;
          int j1 = 0;
          g_id_fr = 0;
          fwalk_time = 0;
          fwalk_time += cur_time;
          fwalk_time_evac = 0;
          fwalk_time_evac += cur_time_evac;

	  is_FR_graver_step_feasible(nc_bool_fr, g_id_fr, nedges_undir, nedges_graver_ij, graver_ij_nonzero_ids, fwalk_fr, cur_fr, cur_sol, graver_ij, cij_undir_fwalk, cij_undir_fr, cij_undir, tij_undir, id_e_opp, cur_time_edges, cur_time_evac_edges, fwalk_time, fwalk_time_edges, fwalk_time_evac, fwalk_time_evac_edges, e_label, edge_sol_lookup, edge_sol_lookup_evac, fwalk_fr_bool, i_fr, j, j1, alpha, beta, t1, n_counter);

         if (nc_bool_fr >= 1)
            {
             update_SR_initial_input(i_fr, nedges_dir, nedges_undir, id_e_opp, graver_ij, cur_sol, cur_fr, cur_time_edges, cur_time_evac_edges, fwalk_time_edges, fwalk_time_evac_edges, cur_sol_sr, cur_fr_sr, cur_time_edges_sr, cur_time_evac_edges_sr, cij_undir, cij_undir_fr, cij_undir_fwalk_sr, cur_fr_bool_sr, fwalk_time, fwalk_time_evac, cur_time_sr, cur_time_evac_sr);

	     
             do_sr(nedges_undir, ngraver_sr, nedges_graver_sr_ij, graver_sr_ij_nonzero_ids, graver_sr_ij, fwalk_sol, cur_sol_sr, cur_time_sr, cur_time_evac_sr, cur_time_edges_sr, cur_time_evac_edges_sr, cur_fr_bool_sr, edge_sol_lookup, edge_sol_lookup_evac, cij_undir_fwalk_sr, tij_undir, fwalk_time_edges, fwalk_time_evac_edges, alpha, beta, t1, n_counter, SR_counter_aug);


             if (cur_time_evac_sr < cur_time_evac)
               {
		       update_current_best_SR_soution(nedges_dir, nedges_undir, id_e_opp, cur_sol, cur_fr, cur_fr_bool, cur_time_edges, cur_time_evac_edges, cur_sol_sr, cur_fr_sr, cur_time_edges_sr, cur_time_evac_edges_sr, cur_time, cur_time_evac, cur_time_sr, cur_time_evac_sr);
		       FR_counter_aug += 1;
               }
             else
             {
             nc_fr += 1;
             }

            }
          else
            {
            nc_fr += 1;
            }
         }
      chk_bool_fr += 1;
      
     if (nc_fr >= 1*ngraver)
       {
        chk_bool_fr = 0;
       }
     cc += 1;
    }


    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "number of function calls: " << n_counter << '\n';
    std::cout << "Elapsed time: " << elapsed.count() << " s\n"; 

    for (int i=0; i<nedges_undir; i++){cout << cur_fr[i] << " ";}
    cout << '\n';
    for (int i=0; i<nedges_undir; i++){cout << cur_sol[i] << " ";}
    cout << '\n';

    double time_sol_gaga_rand_ini = 0;
    time_sol_gaga_rand_ini = elapsed.count() - time_sol_0;
    time_sol_0 = elapsed.count();

    cout << "time_before_leblanc: " << time_sol_0 << '\n';
    cur_time_leblanc = 0;
    cur_time_evac_leblanc = 0;
    vector<double> cur_sol_leblanc(nedges_undir, 0);
    NetworkDesignProblem problem_t(nodes_fname_csv, edges_fname_csv);
    SelfishRoutingSolver solver_t(problem_t, true);
    do_sr_leblanc(solver_t, problem_t, nedges_undir, cur_fr_bool, cur_time_leblanc, cur_time_evac_leblanc, cur_sol_leblanc);

    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    
    double time_sol_leblanc_rand_ini = 0;
    time_sol_leblanc_rand_ini = elapsed.count() - time_sol_0;
    time_sol_0 = elapsed.count();
    cout << "time_after_leblanc: " << time_sol_0 << '\n';

    update_output_data(rand_ini_count, nedges_undir, sr_sol, cur_time_evac, time_sol, time_sol_gaga_rand_ini, sr_sol_ini, cur_fr, cur_sol, cij_undir_fwalk_sr, sr_sol_leblanc, cur_time_evac_leblanc, time_sol_leblanc, time_sol_leblanc_rand_ini, cur_sol_leblanc, ncount_aug_fr, ncount_aug_sr, FR_counter_aug, SR_counter_aug, path_dir_name);
}
}

int main(int argc, char **argv){
    std::string path_dir_name;
    std::string nodes_fname_csv;
    std::string edges_fname_csv;
    path_dir_name = argv[1];
    nodes_fname_csv = argv[2];
    edges_fname_csv = argv[3];
    do_main_FR(path_dir_name, nodes_fname_csv, edges_fname_csv);
}
