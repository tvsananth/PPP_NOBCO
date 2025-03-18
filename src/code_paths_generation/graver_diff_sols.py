import numpy as np

def compute_kernel_feas_sols(feas_sols):

    ker_sols = []
    for i in range(0, feas_sols.shape[0]):
        for j in range(i+1, feas_sols.shape[0]):
            ker_sols.append(feas_sols[i] - feas_sols[j])
            ker_sols.append(feas_sols[j] - feas_sols[i])
            if (np.sum(np.abs(feas_sols[i] - feas_sols[j])) == 1):
               print(i, j)
    ker_sols = np.array(ker_sols)
    return ker_sols


def normal_form(s_i, g_basis, nc):
    g_basis_list = g_basis.copy()
    r_i = s_i
    while (len(g_basis_list) > 0):
          g_i_chk = g_basis_list[0]
          f_chk1 = ((g_i_chk * r_i) >= 0).all()
          f_chk2 = (np.abs(g_i_chk) <= np.abs(r_i)).all()
          if (f_chk1 & f_chk2):
             r_i = r_i - g_i_chk
          else:
             g_basis_list = np.delete(g_basis_list, 0, axis=0)
    return r_i


def compute_gbasis_par(ker_par):

    ker_el_ones_copy = ker_par.copy()
    tst1 = True
    cc = 0
    while (tst1):
       tst1 = False
       del_id_rows = []
       for i in range(ker_el_ones_copy.shape[0]):
           g1_i = ker_el_ones_copy[i]
           g1_i_del = np.delete(ker_el_ones_copy, i, axis=0)
           g1_r_i = normal_form(g1_i, g1_i_del, len(g1_i_del))
           if ((g1_r_i != g1_i).any()):
              if ((g1_r_i != 0).any()):
                 ker_el_ones_copy[i] = g1_r_i
                 tst1 = True
                 print(i, cc)
              else:
                 del_id_rows.append(i)
       if (len(del_id_rows) > 0):
          ker_el_ones_copy = np.delete(ker_el_ones_copy, np.array(del_id_rows), axis=0)
       cc += 1
    return ker_el_ones_copy
