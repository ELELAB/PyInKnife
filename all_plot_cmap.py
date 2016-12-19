#!/usr/bin/env python

__author__ = 'Juan Salamanca Viloria'

#    PyinKnife script to create dataframes and plots for the connected components.
#    Copyright (C) 2016 Juan Salamanca Viloria <juan.salamanca.viloria@gmail.com>, Elena Papaleo <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

##
import re
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(description=
"""
Plot the total number of hubs for the whole trajectory and each resampling. 
""")
parser.add_argument("files", nargs= "*")
parser.add_argument("-x", metavar="xvg", help= "xvg file")
parser.add_argument("-hb" , type=str, choices=['all','mc-mc','mc-sc','sc-sc'],nargs= "*")
args = parser.parse_args()


last_line = open(args.x).readlines()[-1]
time_length = last_line.split()[0]
s_time = int(round(float(time_length), -1))

def calculate_contact(file):
    component = []
    for line in file:
        line = line.strip()
        if line == "Connected component 1:":
            break
    for line in file:
        line = line.strip()
        if "Connected component" not in line:
            number = re.search(r'([0-9]+)',line)
            #print "number", number.group()
            #print "LINEA", line
            component.append(int(number.group(0)))

    return component

def list_to_matrix(lista):
    
    matrix_hubs= []
    for x in lista:
        matrix_hubs.append((x,lista.count(x)))
    matrix_hubs = sorted(set(matrix_hubs))
    matrix_hubs = np.array(matrix_hubs)
    return matrix_hubs


n= 0
list_files= []
for filename in args.files[0:]:
    n+=1
 
    list_files.append(open(filename, "r"))

traj_length =  round(s_time/float(1000),-1)

steps= traj_length/float(10)

contact_map = [calculate_contact(x) for x in list_files]


number_ccont=[len(c) for c in contact_map]
max_ccont = [max(d) for d in contact_map]
    
top5_ccont = [list(enumerate(sorted(e, reverse=True),1))[:5] for e in contact_map]

#print top5_ccont

final = np.hstack(top5_ccont)

#print final
final=  np.delete(final, [2,4,6,8,10,12,14,16,18,20], axis = 1)



frames_position= range(0,int(traj_length + steps),int(steps))



if any(i in ["all", "mc-mc", "sc-sc", "mc-sc"] for i in args.hb):
    plt.plot(frames_position, number_ccont,marker='o', linestyle='--', color='r')
    plt.title("Number of connected components")
    plt.xlim(-1,)
    plt.ylim(min(number_ccont)-5,max(number_ccont) +5)
    plt.xlabel("$ns$")

    for i in args.hb: 
        if i == "all":
            plt.savefig(os.getcwd() + "/total_connectedcomp.pdf", format='pdf')
            
        if i == "sc-sc":
            plt.savefig(os.getcwd() + "/total_connectedcomp_sc.pdf", format='pdf')
        if i == "mc-mc":
            plt.savefig(os.getcwd() + "/total_connectedcomp_mc.pdf", format='pdf')
        if i == "mc-sc":
            plt.savefig(os.getcwd() + "/total_connectedcomp_sc-mc.pdf", format='pdf')
    plt.close()
    
    plt.plot(frames_position, max_ccont,marker='o', linestyle='--', color='r')
    plt.xlim(-1,)
    plt.ylim(min(max_ccont)-5,max(max_ccont)+5)
    plt.xlabel("$ns$")
    plt.ylabel("Maximum number of connected components")
    for i in args.hb: 
        if i == "all":
            plt.savefig(os.getcwd() + "/total_connectedcomp_max.pdf", format='pdf')
            np.savetxt("final_cmap.txt",final)
        if i == "sc-sc":
            plt.savefig(os.getcwd() + "/total_connectedcomp_max_sc.pdf", format='pdf')
            np.savetxt("final_cmap_sc.txt",final)
        if i == "mc-mc":
            plt.savefig(os.getcwd() + "/total_connectedcomp_max_mc.pdf", format='pdf')
            np.savetxt("final_cmap_mc.txt",final)
        if i == "mc-sc":
            plt.savefig(os.getcwd() + "/total_connectedcomp_max_mc-sc.pdf", format='pdf')
            np.savetxt("final_cmap_mc-sc.txt",final)
    plt.close()
    
#else:
#    plt.plot(frames_position, number_ccont,marker='o', linestyle='--', color='r')
#    plt.title("Number of connected components")
#    plt.xlim(-1,)
#    plt.ylim(min(number_ccont)-5,max(number_ccont) +5)
#    plt.xlabel("$ns$")
#    plt.savefig(os.getcwd() + "/total_connectedcomp.pdf", format='pdf')
#    plt.close()
#
#
#    plt.plot(frames_position, max_ccont,marker='o', linestyle='--', color='r')
#    plt.xlim(-1,)
#    plt.ylim(min(max_ccont)-5,max(max_ccont)+5)
#    plt.xlabel("$ns$")
#    plt.ylabel("Maximum number of connected components")
#    plt.savefig(os.getcwd() + "/total_connectedcomp_max.pdf", format='pdf')
#    plt.close()
#    np.savetxt("final_cmap.txt",final)
