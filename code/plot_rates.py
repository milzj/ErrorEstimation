#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import ast
import sys
import lsqs_label

def open_file(file_name):
    with open(file_name, 'r') as f:
        data = f.read()
    dictionary = ast.literal_eval(data)
    return dictionary

def find_nref(file_name):
    start_position = file_name.find('nref_')+5
    end_position = file_name.find('_',start_position,-1)
    return int(file_name[start_position:end_position])

def convert_nref_base2(nref):
    x = int(np.log2(nref))
    assert 2**x == nref
    return x

def plot_estimate_error(xpoints,canonical_map,normal_map,gap,rgap,nref,plot_number=None):
    log_xpoints = np.log2(xpoints)
    log_can_map = np.log2(canonical_map)
    log_nor_map = np.log2(normal_map)
    log_gap = np.log2(gap)
    log_rgap = np.log2(rgap)

    X = np.vstack([np.ones(log_xpoints.size), log_xpoints]).T
    b1 = np.linalg.lstsq(X, log_can_map, rcond=None)[0]
    b2 = np.linalg.lstsq(X, log_nor_map, rcond=None)[0]
    b3 = np.linalg.lstsq(X, log_gap, rcond=None)[0]
    b4 = np.linalg.lstsq(X, log_rgap, rcond=None)[0]
    log_can_fit = X @ b1
    #log_nor_fit = X @ b2
    log_gap_fit = X @ b3
    #log_rgap_fit = X @ b4
    b = [b1,b2,b3,b4]
    c = np.zeros(4)
    r = np.zeros(4)
    lsqs_base = "n"
    for i in range(4):
        c[i] = np.exp(b[i][0])
        r[i] = b[i][1]

    global_color = '#000000'
    plt.plot([], [], ' ', label=lsqs_label.lsqs_label_base(base=2, rate=convert_nref_base2(nref)))
    plt.plot(xpoints,canonical_map, 'o', markerfacecolor='none', markeredgecolor=global_color,markersize=12,label = r"$\chi_{\mathrm{can}}(u; \tau)$")
    plt.plot(xpoints, np.exp2(log_can_fit), '--', label=lsqs_label.lsqs_label(rate=r[0], constant=c[0], base=lsqs_base), color=global_color)
   
    plt.plot(xpoints, normal_map, 'd', label=r"$\chi_{\mathrm{nor}}(v; \tau)$")  
   #plt.plot(xpoints, np.exp2(log_nor_fit), ':', label=lsqs_label.lsqs_label(rate=r[1], constant=c[1], base=lsqs_base),color=global_color)
    
    plt.plot(xpoints, gap, '^', markerfacecolor='none', markeredgecolor=global_color,markersize=12, label=r"$\chi_{\mathrm{gap}}(u)$")  
    plt.plot(xpoints, np.exp2(log_gap_fit), '-.', label=lsqs_label.lsqs_label(rate=r[2], constant=c[2], base=lsqs_base),color=global_color)
    
    #plt.plot(xpoints, rgap, 'D', markerfacecolor='none', markeredgecolor=global_color,markersize=7, label=r"$\chi_{\mathrm{rgap}}(u;\tau)$")  
    #plt.plot(xpoints, np.exp2(log_rgap_fit), '', label=lsqs_label.lsqs_label(rate=r[3], constant=c[3], base=lsqs_base),color=global_color)
    
    plt.xlabel(r"${}$".format(lsqs_base))
    plt.yscale('log',base=10)
    plt.xscale('log',base=2)
    plt.legend(loc='lower left')
    
    #now = sys.argv[1] 
    
    """
    if plot_number is None:
        #timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        base_filename = "Linear_Problem_{}".format(now)
    elif plot_number == 1:
        base_filename = "Bilinear_Problem_{}".format(now)
    else:
        base_filename = "Semilinear_Problem_{}".format(now)
    for f in ["png", "pdf"]:
        plt.savefig(f"{base_filename}.{f}")
    """
    if plot_number is None:
        Problem = 'LinearProblem'
    elif plot_number == 1:
        Problem = 'BilinearProblem'
    else:
        Problem = 'SemilinearProblem'
    
    #name = Problem().__str__()
    outdir = "output/{}/{}/".format(now,Problem)

    for f in ["png", "pdf"]:
        #plt.savefig(f"{base_filename}.{f}")
        plt.savefig(outdir+"{}_criticality_measures_{}.{}".format(Problem,now,f), bbox_inches='tight')
    plt.close()

def criticality_measure_rate_plot(file_name,plot_number=None):
    dict_data = open_file(file_name)
    nref = find_nref(file_name)
    xpoints = dict_data['n']
    canonical_map = dict_data['canonical_map']
    normal_map = dict_data['normal_map']
    gap = dict_data['gap']
    rgap = dict_data['rgap']
    plot_estimate_error(xpoints,canonical_map,normal_map,gap,rgap,nref,plot_number)

criticality_measure_rate_plot('criticality_measures_nref_524288_18-Feb-2024-14-08-31.txt',plot_number=None)
criticality_measure_rate_plot('criticality_measures_nref_2048_18-Feb-2024-14-08-31.txt',plot_number=1)
criticality_measure_rate_plot('criticality_measures_nref_2048_18-Feb-2024-14-08-31 (1).txt',plot_number=2)