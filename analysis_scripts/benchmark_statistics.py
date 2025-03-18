import json
import argparse
import numpy as np
import random
import scipy
from scipy.stats import lognorm, CensoredData
import matplotlib
import matplotlib.pyplot

def load_json_file(fname):

    with open(fname, 'r') as f_input:
        data = json.load(f_input)
    n_ins = len(data)

    return data, n_ins


def bstrap_Tq(t1, t1_bool, T95):
    N = t1.shape[0]
    M = N
    T95_M = np.zeros(M)
    for i1 in range(M):
        t1_i1 = np.asarray(random.choices(t1, k=N))
        t1_i1_argsort = np.argsort(t1_i1)
        t1_i1_sort = t1_i1[t1_i1_argsort]
        T95_M[i1] = t1_i1_sort[94]
    T95_M_err = np.sqrt(((T95_M - T95)**2.0).sum()/(M))

    return T95_M_err


def get_stats_Tq(t1, t1_bool):

    t1_argsort = np.argsort(t1)
    t1_sort = t1[t1_argsort]
    t1_bool_sort = t1_bool[t1_argsort]

    T95 = t1_sort[94]
    T95_err = 2.0*bstrap_Tq(t1, t1_bool, T95)
    T95_interval_01 = [T95 - 0.5*T95_err, T95 + 0.5*T95_err]
    T95_interval0 = round(T95_interval_01[0], 2)
    T95_interval1 = round(T95_interval_01[1], 2)
    T95_interval = [T95_interval0, T95_interval1]

    return T95, T95_err, T95_interval


def bstrap_samples_Tq(t1, t1_bool):
    N = t1.shape[0]
    M = N
    T95_M = np.zeros(M)
    for i1 in range(M):
        t1_i1 = np.asarray(random.choices(t1, k=N))
        t1_i1_argsort = np.argsort(t1_i1)
        t1_i1_sort = t1_i1[t1_i1_argsort]
        T95_M[i1] = t1_i1_sort[94]

    return T95_M


def samples_Tq(t1, t1_bool):

    t1_argsort = np.argsort(t1)
    t1_sort = t1[t1_argsort]
    t1_bool_sort = t1_bool[t1_argsort]

    T95 = t1_sort[94]
    T95_M = bstrap_samples_Tq(t1, t1_bool)

    return T95, T95_M


def get_stats_Tr(t1, t1_bool):

    n_total = t1.shape[0]
    n_succ  = t1_bool.sum()
    n_fail  = (1 - t1_bool).sum()

    p_succ = n_succ*1.0/n_total

    if (p_succ < 1):
        T95 = t1.max()*np.log(1 - 0.95)/np.log(1 - p_succ)
        p_succ_l = scipy.stats.beta.ppf(1.0 - 0.05/2.0, 0.5 + n_succ, 0.5 + n_fail)
        p_succ_u = scipy.stats.beta.ppf(0.05/2.0, 0.5 + n_succ, 0.5 + n_fail)
        T95_l = t1.max()*np.log(1 - 0.95)/np.log(1 - p_succ_l)
        T95_u = t1.max()*np.log(1 - 0.95)/np.log(1 - p_succ_u)
        T95_err = (T95_u - T95_l)
        T95_interval = [round(T95_l, 2), round(T95_u, 2)]
    else:
        T95 = 0
        T95_err = 0
        T95_interval = [0, 0]
    return T95, T95_err, T95_interval


def fit_stats(t1, t1_bool):

    try:
        n_total = t1.shape[0]
        n_succ  = t1_bool.sum()
        n_fail  = (1 - t1_bool).sum()
        p_succ = n_succ*1.0/n_total
        ft1 = t1.copy()
        if (n_succ == n_total):
            fit_pars = lognorm.fit(ft1)
        else:
            Tr99 = t1.max()*np.log(1 - 0.99)/np.log(1 - p_succ)
            ft1[t1_bool == 0] = Tr99
            ft1_d = CensoredData.right_censored(ft1, ft1 > t1.max())
            fit_pars = lognorm.fit(ft1_d)
        T95 = lognorm.ppf(0.95, fit_pars[0], fit_pars[1], fit_pars[2])
    except:
        T95 = 0

    return T95


def get_stats_Tmle(t1, t1_bool):

    n_total = t1.shape[0]
    n_succ  = t1_bool.sum()
    n_fail  = (1 - t1_bool).sum()
    p_succ = n_succ*1.0/n_total
    T95 = fit_stats(t1, t1_bool)
    T95_bstrap = np.zeros(n_total)
    t1_ids = np.arange(n_total)
    for i1 in range(n_total):
        ids_rand = np.asarray(random.choices(t1_ids, k=n_total))
        T95_i1 = fit_stats(t1[ids_rand], t1_bool[ids_rand])
        T95_bstrap[i1] = T95_i1
    T95_lu = np.quantile(T95_bstrap, 0.025) 
    T95_ru = np.quantile(T95_bstrap, 0.975)
    T95_err = 1.0*(T95_ru - T95_lu)
    T95_interval = [T95_lu, T95_ru]

    if (p_succ < 1):
        T95 = 0
        T95_err = 0
        T95_interval = 0

    return T95, T95_err, T95_interval


def do_stats(i, data_i, T_bnb, metric_type):

    if (metric_type == 'TTS'):
        t1 = np.asarray(data_i['runtime_seconds'])

    if (metric_type == 'ETS'):
        t1 = np.asarray(data_i['energy'])

    t1_bool = 1 - np.asarray(data_i['n_unsat_clauses'])
    
    Tq95, Tq95_err, Tq95_interval = get_stats_Tq(t1, t1_bool)

    Tr95, Tr95_err, Tr95_interval = get_stats_Tr(t1, t1_bool)

    Tmle95, Tmle95_err, Tmle95_interval = get_stats_Tmle(t1, t1_bool)

    stats = {}
    stats['Min. Obj. (% fractional difference from cSOA)'] = round(np.array(data_i['configurations']).min(), 2)
    stats['Max. Obj. (% fractional difference from cSOA)'] = round(np.array(data_i['configurations']).max(), 2)
    stats['Tq95'] = {}
    stats['Tr95'] = {}
    stats['Tmle95'] = {}
    stats['Tq95'][metric_type] = round(Tq95, 2)  
    stats['Tq95']['Error'] = round(Tq95_err, 2)
    stats['Tq95']['Confidence Interval'] = Tq95_interval
    stats['Tq95']['Digital Efficiency'] = round(T_bnb/Tq95, 2)
    stats['Tr95'][metric_type] = round(Tr95, 2)
    stats['Tr95']['Error'] = round(Tr95_err, 2)
    stats['Tr95']['Confidence Interval'] = Tr95_interval
    if (Tr95 != 0):
        stats['Tr95']['Digital Efficiency'] = round(T_bnb/Tr95, 2)
    else:
        stats['Tr95']['Digital Efficiency'] = '-'
    stats['Tmle95'][metric_type] = round(Tmle95, 2)
    stats['Tmle95']['Error'] = round(Tmle95_err, 2)
    stats['Tmle95']['Confidence Interval'] = Tmle95_interval
    if (Tmle95 != 0):
        stats['Tmle95']['Digital Efficiency'] = round(T_bnb/Tmle95, 2)
    else:
        stats['Tmle95']['Digital Efficiency'] = '-'
    print('Instance:', i, 'Success fraction (%): ', t1_bool.sum())

    print(i, stats['Min. Obj. (% fractional difference from cSOA)'], stats['Max. Obj. (% fractional difference from cSOA)'], t1_bool.sum(), stats['Tq95'][metric_type], stats['Tq95']['Confidence Interval'], stats['Tq95']['Digital Efficiency'], stats['Tr95'][metric_type], stats['Tr95']['Confidence Interval'], stats['Tr95']['Digital Efficiency'])
    return stats
    

def do_stats_quantile(tts, tts_M, p):

    n = tts.shape[0]
    M = tts_M.shape[1]
    q_all = np.quantile(tts, p)
    q_M = np.zeros(M)
    for i in range(M):
        q_M[i] = np.quantile(tts_M[:, i], p)
    q_M_err = np.sqrt(((q_M - q_all)**2.0).sum()/(M))
    q_M_mean = np.mean(q_M)

    return q_all, q_M_err, q_M_mean


def batch_metrics(n_ins, data, metric_type):

    tts = []
    tts_M = []

    for i in range(n_ins):
        
        if (metric_type == 'TTS'):
            t1 = np.asarray(data[i]['runtime_seconds'])

        if (metric_type == 'ETS'):
            t1 = np.asarray(data[i]['energy'])

        t1_bool = 1 - np.asarray(data[i]['n_unsat_clauses'])

        T95, T95_M = samples_Tq(t1, t1_bool)

        tts.append(T95)
        tts_M.append(T95_M)

    tts = np.array(tts)
    tts_M = np.array(tts_M)

    q50, q50_err, q50_mean = do_stats_quantile(tts, tts_M, 0.5)
    q75, q75_err, q75_mean = do_stats_quantile(tts, tts_M, 0.75)
    q90, q90_err, q90_mean = do_stats_quantile(tts, tts_M, 0.9)

    cdf = np.linspace(0, 1, 101)
    qcdf_all = np.zeros(cdf.shape[0])
    qcdf_err = np.zeros(cdf.shape[0])
    qcdf_mean = np.zeros(cdf.shape[0])
    for i in range(cdf.shape[0]):
        qp, qp_err, qp_mean = do_stats_quantile(tts, tts_M, cdf[i])
        qcdf_all[i] = qp
        qcdf_err[i] = qp_err
        qcdf_mean[i] = qp_mean

    l_tts = np.log10(tts)
    print('Batch Metric:', l_tts.mean(), l_tts.std())
    print(q50, q50_err, q50_mean)
    print(q75, q75_err, q75_mean)
    print(q90, q90_err, q90_mean)

    return q50, q50_err, q75, q75_err, q90, q90_err, cdf, qcdf_all, qcdf_err


def main_compute():

    parser = argparse.ArgumentParser()

    parser.add_argument('-path_results_dir', '--path_results_dir_cinput', type=str, help='path to the results directory')

    parser.add_argument('-json_filename', '--json_filename_cinput', type=str, help='json file name')

    parser.add_argument('-graph_size', '--graph_size_cinput', type=int, default=0, help='graph node size')

    parser.add_argument('-graph_den', '--graph_den_cinput', type=float, default=0, help='graph edge density')

    parser.add_argument('-decomp_str', '--decomp_str_cinput', type=str, default='', help='Type of decomposition used')

    parser.add_argument('-metric_type', '--metric_type_cinput', type=str, help='Metric type: TTS or ETS')

    parser.add_argument('-T_bnb', '--T_bnb_cinput', type=float, help='Run time (in seconds) of Branch-and-Bound')

    args = parser.parse_args()

    path_results_dir = args.path_results_dir_cinput

    path_plots_dir = path_results_dir + '/plots/'

    fname = path_results_dir + '/' + args.json_filename_cinput

    graph_size = args.graph_size_cinput

    graph_den = args.graph_den_cinput

    decomp_str = args.decomp_str_cinput

    metric_type = args.metric_type_cinput

    T_bnb = args.T_bnb_cinput


    if (graph_size == 20):
        decomp_str = '_' + decomp_str

    data, n_ins = load_json_file(fname)

    print('n_instances:', n_ins, 'keys:', data[0].keys())

    fstats = []
    for i in range(n_ins):
        stats = do_stats(i, data[i], T_bnb, metric_type)
        fstats.append(stats)
        #print(stats)

    with open(path_results_dir + '/metrics_' + metric_type + '_' + str(graph_size) + '_' + str(graph_den) + decomp_str + '.json', 'w') as f_out:
        json.dump(fstats, f_out)
    with open(path_results_dir + '/metrics_' + metric_type + '_' + str(graph_size) + '_' + str(graph_den) + decomp_str + '.json', 'r') as f_in:
        fstats_jfile = json.load(f_in)

    q50, q50_err, q75, q75_err, q90, q90_err, cdf, qcdf_all, qcdf_err = batch_metrics(n_ins, data, metric_type)

    matplotlib.pyplot.figure(1)
    matplotlib.pyplot.plot(qcdf_all, cdf)
    matplotlib.pyplot.xlabel(metric_type + ' (seconds)', fontsize=20)
    matplotlib.pyplot.ylabel(r'$\mathrm{CDF}$', fontsize=20)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(path_plots_dir + '/CDF_' + metric_type + '_' + str(graph_size) + '_' + str(graph_den) + decomp_str + '.png')
    matplotlib.pyplot.show()

    matplotlib.pyplot.figure(2)
    matplotlib.pyplot.errorbar(cdf*100, qcdf_all, yerr=qcdf_err)
    matplotlib.pyplot.xlabel(r'Quantile', fontsize=20)
    matplotlib.pyplot.ylabel(metric_type + ' (seconds)', fontsize=20)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(path_plots_dir + '/Quantiles_' + metric_type + '_' + str(graph_size) + '_' + str(graph_den) + decomp_str + '.png')
    matplotlib.pyplot.show()


if __name__ == '__main__':
     main_compute()

