#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import ast
import sys
import lsqs_label

from figure_style import *

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

def plot_estimate_error(xpoints,canonical_map,normal_map,gap,rgap,nref,timestamp,problem):
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
    log_nor_fit = X @ b2
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
   
    # canonical criticality measure
    plt.plot(xpoints,canonical_map, 'o', markerfacecolor='none', markeredgecolor=global_color,markersize=12,label = r"$\chi_{\mathrm{can}}(u_h; 1)$")
    #plt.plot(xpoints, np.exp2(log_can_fit), '--', label=lsqs_label.lsqs_label(rate=r[0], constant=c[0], base=lsqs_base), color=global_color)
   
    # normal map
    plt.plot(xpoints, normal_map, 'd', label=r"$\chi_{\mathrm{nor}}(v_h; 1)$", color=global_color,  markeredgecolor=global_color)
    plt.plot(xpoints, np.exp2(log_nor_fit), '--', label=lsqs_label.lsqs_label(rate=r[1], constant=c[1], base=lsqs_base), color=global_color)

    # normal gap function
    plt.plot(xpoints, gap, '^', markerfacecolor='none', markeredgecolor=global_color,markersize=12, label=r"$\chi_{\mathrm{gap}}(u_h)$")
    plt.plot(xpoints, np.exp2(log_gap_fit), '-.', label=lsqs_label.lsqs_label(rate=r[2], constant=c[2], base=lsqs_base),color=global_color)
    
    plt.xlabel(r"${}$".format(lsqs_base))
    plt.yscale('log',base=10)
    plt.xscale('log',base=2)
    plt.legend(loc='best')
    
    outdir = "output/{}/{}/".format(timestamp,problem)

    for f in ["png", "pdf"]:
        plt.savefig(outdir+"{}_criticality_measures_{}.{}".format(problem,timestamp,f), bbox_inches='tight')
    plt.close()

def criticality_measure_rate_plot(file_name,timestamp, problem):
    dict_data = open_file(file_name + ".txt")
    nref = find_nref(file_name)
    xpoints = dict_data['n']
    canonical_map = dict_data['canonical_map']
    normal_map = dict_data['normal_map']
    gap = dict_data['gap']
    rgap = dict_data['rgap']
    plot_estimate_error(xpoints,canonical_map,normal_map,gap,rgap,nref,timestamp, problem)


if __name__ == "__main__":

    timestamp = sys.argv[1]
    problem = sys.argv[2]
    filename = np.loadtxt("output/{}/{}/criticality_measures_filename.txt".format(timestamp, problem), dtype="str")
    criticality_measure_rate_plot(str(filename), timestamp, problem)
