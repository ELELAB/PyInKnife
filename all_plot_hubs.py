#!/usr/bin/env python

__author__ = 'Juan Salamanca Viloria'



###
##
##

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
#print s_time
def calculate_hubs(file):
    hubs = []
    i = 0
    for line in file:
        line = line.strip()
        if line == "Hubs:":
            break
        elif "No hubs" in line:
            print "None Hubs Found. Please Check your results. Replaced by 0"
            hubs.append(0)
    for line in file:
        line = line.strip()
        i+= 1
        if not line.startswith("Node"):
            #print line.split("\t")
            #hubs += line.split("\t")[1]
            hubs.append(int(line.split("\t")[1]))
    return hubs

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
        if x == 0 :
            matrix_hubs.append((0,0))
        else:
            matrix_hubs.append((x,lista.count(x)))
    matrix_hubs = sorted(set(matrix_hubs))
    matrix_hubs = np.array(matrix_hubs)
    return matrix_hubs


n= 0
list_files= []
#print args.files
for filename in args.files[0:]:
    #print filename
    n+=1
    locals()["hub_list"+ str(n)] = open(filename, "r")
    list_files.append(open(filename, "r"))
    
traj_length = round(s_time/float(1000),-1)
steps= traj_length/float(10)


hubs = []
for x in list_files:
    hubs.append(calculate_hubs(x))

matrix_hubs= []
for m in hubs:
    matrix_hubs.append(list_to_matrix(m))


def transform_hubs(array1, array2):

    if array1.shape != array2.shape:
        if set(array2[:,0]) == None : 
            print "None hubs found"
        #print "array1", set(array1[:,0])
        #print "array2", set(array2[:,0])
        dif_number= set(array1[:,0]) - set(array2[:,0])
        #print "diff number", dif_number
        if len(dif_number) ==  1:
            array3= np.vstack([array2, [list(dif_number)[0], 0]])
            #print "posible resultado", array3
            array3= array3[np.argsort(array3[:,0])]
            return array3
        elif len(dif_number) > 1:
            array3= np.vstack([array2, [list(dif_number)[0], 0]])
            for i in range(1, len(dif_number)):
                array3= np.vstack([array3, [list(dif_number)[i], 0]])
            array3 = array3[np.argsort(array3[:,0])]
            return array3
        
        elif len(dif_number) < 1:
            print "None hubs found"
            array3 = np.vstack([array2, [list(dif_number)[0], 0]])
            array3= array3[np.argsort(array3[:,0])]
            return array3
    else:
        dif_number= set(array1[:,0]) - set(array2[:,0])
        #print "DIff number", dif_number
        if len(dif_number) ==  1:
            array3= np.vstack([array2, [list(dif_number)[0], 0]])
            #print "posible resultado", array3
            array3= array3[np.argsort(array3[:,0])]
            return array3
        elif len(dif_number) > 1:
            array3= np.vstack([array2, [list(dif_number)[0], 0]])
            for i in range(1, len(dif_number)):
                array3= np.vstack([array3, [list(dif_number)[i], 0]])
            array3 = array3[np.argsort(array3[:,0])]
            return array3

    return array2


def shape_list(list):
    shapes= []
    for x in list:
        shapes.append(x.shape[0])    
    return shapes

dim_list= []
for x in matrix_hubs:
    dim_list.append(x.shape[0])

reference=  np.argmax(dim_list)
#print "input", matrix_hubs
#print "reference", reference
new_array_list= []
for i in range(0, len(matrix_hubs)):
    if i != reference:
        #print "i", i, matrix_hubs[i]
        #print "transformation", transform_hubs(matrix_hubs[reference], matrix_hubs[i])
        new_array_list.append((transform_hubs(matrix_hubs[reference], matrix_hubs[i])))

new_array_list.insert(reference, matrix_hubs[reference])
#print "output", new_array_list
#
#final = np.hstack(new_array_list)
#
#print final
#final=  np.delete(final, [2,4,6,8,10,12,14,16,18,20], axis = 1)
#np.savetxt("distribution_hubs_sc.txt",final)

dim_list2= shape_list(new_array_list)
while len(set(dim_list2)) > 1:
### Check it again in case one element was missed
    #print "PRUEBA SET", set(dim_list2)

    reference2=  np.argmax(dim_list2)
#print "REFERENCE ARRAY", new_array_list[reference2]
#print "reference", reference2
    new_array_list2= []
    for i in range(0, len(new_array_list)):
        if i != reference2:
        #print "i", i, new_array_list[i]
        #print "transformation", transform_hubs(new_array_list[reference2], new_array_list[i])
            new_array_list2.append((transform_hubs(new_array_list[reference2], new_array_list[i])))

    new_array_list2.insert(reference2,new_array_list[reference2])
    dim_list2 = shape_list(new_array_list2)
    new_array_list= new_array_list2


final = np.hstack(new_array_list)
#print "FINAL", final
final=  np.delete(final, [2,4,6,8,10,12,14,16,18,20], axis = 1)


#print len(hubs)
#print len(hubs2)

##next time give as argument the lenght of the simulation and the step
frames_position= range(0,int(traj_length + steps),int(steps))

hub_lengths = [len(a) for a in hubs]

    
if any(i in ["all", "mc-mc", "sc-sc", "mc-sc"] for i in args.hb):
    plt.plot(frames_position,hub_lengths, marker='o', linestyle='--', color='r')
    plt.title("Total number of hubs per frame")
    plt.xlim(-1,)

    for i in args.hb: 
        if i == "all":
            plt.savefig(os.getcwd() + "/total_hubs.pdf",format='pdf')
            np.savetxt("distribution_hubs.txt",final[~(final==0).all(1)])
            
        if i == "sc-sc":
            plt.savefig(os.getcwd() + "/total_hubs_sc.pdf",format='pdf')
            np.savetxt("distribution_hubs_sc.txt",final[~(final==0).all(1)])
        if i == "mc-mc":
            plt.savefig(os.getcwd() + "/total_hubs_mc.pdf",format='pdf')
            np.savetxt("distribution_hubs_mc.txt",final[~(final==0).all(1)])
        if i == "mc-sc":
            plt.savefig(os.getcwd() + "/total_hubs_mc-sc.pdf",format='pdf')
            np.savetxt("distribution_hubs_mc-sc.txt",final[~(final==0).all(1)])
    plt.close()
#    
#else:
#    plt.plot(frames_position,hub_lengths, marker='o', linestyle='--', color='r')
#    plt.title("Total number of hubs per frame")
#    plt.xlim(-1,)
#    plt.savefig(os.getcwd() + "/total_hubs.pdf",format='pdf')
#    plt.close()
#    np.savetxt("distribution_hubs.txt",final[~(final==0).all(1)])