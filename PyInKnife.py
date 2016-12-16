#!/usr/bin/env python

#    Pipeline for graph analysis, a software suite to analyze interactions and interaction network in structural 
#    ensembles.
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

import argparse
import sys
import os
import time

start = time.time()


parser = argparse.ArgumentParser(description=
"""
Pipeline for graph analysis with PYINTERAPH. It creates two folders, contact and h-bond
Required -f trajectory no jumps (filtered) -s tpr file   -cut_off 1st cut off for the contact map -cut_off2 2nd cut off
""")


parser.add_argument("-f", metavar="xtc", required=True, help="xtc trajectory no jumps")
parser.add_argument("-s", metavar="tpr", required=True, help="tpr file")


parser.add_argument("-cut_off",type=str,help="1st cut off for the contact map", default=5)
parser.add_argument("-cut_off2",type=str,help="cut off for the contact map")
parser.add_argument("-filter", type =int, help="rate of persistance", default = 20)
parser.add_argument("-ff", type=str, help ="Force field to be used (for masses calculation only)")
parser.add_argument("-k", type=int, help="Minimum number of connections for hubs")
parser.add_argument("-hb" , type=str, choices=['all','mc-mc','mc-sc','sc-sc', 'no'],nargs= "*")
parser.add_argument("-kbp", action ="store_true", help = "Output file for the side-chains interactions energy in weighted adjacency matrix format")
parser.add_argument("-x", metavar="xvg", help= "xvg file")
parser.add_argument("-jackknife", action ="store_true", help ="To perform the convergence method (jackknife)")
parser.add_argument("-pbc_nojump", action ="store_true", help="Remove PBC artifacts")
#parser.add_argument("-steps", type= int, help = "Number that will divide the whole trajectory in frames of those length")
args = parser.parse_args()

## Store the path
path = os.getcwd()
## Create new directory call pyinteraph and move to it
os.mkdir(path + "/pyinteraph")
os.chdir(path + "/pyinteraph")
##Create the index for pyinteraph (keep 1 and q)
os.system("make_ndx -f ../%s -o prot.ndx" % (args.s))
pbc = False
if args.pbc_nojump:

    pbc = True
    os.system("trjconv -f ../%(1)s -s ../%(2)s -pbc mol -ur compact -o traj_mol_ur.xtc  -n prot.ndx" % {"1": args.f, "2": args.s})
    os.system("trjconv -f traj_mol_ur.xtc -s ../%s -pbc nojump -o traj_mol_ur_nojump.xtc" % args.s)

## Create nohup file
os.system("touch nohup.out")
## Create a child process to run all in the background
pid = os.fork()
if pid == 0:
    # Child
    os.setsid()  # This creates a new session
    print "In background:", os.getpid()

    if pbc:

        os.system("trjconv -f traj_mol_ur_nojump.xtc -s ../%s  -fit rot+trans -sep -o reference.fixed.pdb -b 0 -e 0 -n prot.ndx" % args.s)

    ## Extract the 1st frame to get the pdb structure and convert the trajectory file into a pdb
    else:
        os.system("trjconv -f ../%(1)s -s ../%(2)s  -fit rot+trans -sep -o reference.fixed.pdb -b 0 -e 0 -n prot.ndx" % {"1": args.f, "2": args.s})

## Extract the topologies, need to format gromacs to pdb
    os.system("editconf -f reference.fixed0.pdb -o reference0.gro")
    ## Store the new path
    path1 = os.getcwd()
    ## Create a new directory and move to it
    os.mkdir(path1 + "/contact")
    os.chdir(path1 + "/contact/")
    ## Create a directory for the cut off and move to it
    os.mkdir(str(args.cut_off))
    os.chdir(path1 + "/contact/%s/" % args.cut_off)
    ## Create a nohup file
    os.system("touch nohup.out")

    print("Pyinteraph analysis 1st cut off")
    ## Command line used to work on the pyinteraph environment
    execfile("/usr/local/envs/pyinteraph/bin/activate_this.py", dict(__file__="/usr/local/envs/pyinteraph/bin/activate_this.py"))

    ## Contact map analysis (using all the aminoacids except GLY)
    os.system("nohup pyinteraph -v -s ../../reference0.gro -t ../../../%(1)s  -r ../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(2)s --ff-masses %(3)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": args.f, "2": args.cut_off, "3": args.ff})
    if pbc:
        os.system("nohup pyinteraph -v -s ../../reference0.gro -t ../../traj_mol_ur_nojump.xtc  -r ../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(1)s --ff-masses %(2)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": args.cut_off, "2": args.ff})

    ## Estimate significance threshold for the graph and use it to filter it. ( 20 by default)
    os.system("filter_graph -d contact.dat -o filter_contact.dat -t %s" % args.filter)


    ##Calculate the component contact and list of hubs.  Make a copy of the pdb because replace columns
    os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact.dat -c -cb reference_cb.pdb > component_contact.out 2>&1")
    os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact.dat -u -ub reference_hub.pdb -k %s > hubs_list.out 2>&1" % args.k)
    ## In case we want to use a 2nd cut off it will make the same calculations but in a new folder.
    if args.cut_off2 > 0:

      os.chdir(path1 + "/contact/")
      os.mkdir(str(args.cut_off2))
      os.chdir(path1 + "/contact/%s/" % args.cut_off2)

      os.system("touch nohup.out")
      print("Pyinteraph analysis 2nd cut off")
      os.system("nohup pyinteraph -v -s ../../reference0.gro -t ../../../%(1)s  -r ../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(2)s --ff-masses %(3)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": args.f, "2": args.cut_off2, "3": args.ff})
      if pbc:
          os.system("nohup pyinteraph -v -s ../../reference0.gro -t ../../traj_mol_ur_nojump.xtc  -r ../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(1)s --ff-masses %(2)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": args.cut_off2, "2": args.ff})

      os.system("filter_graph -d contact.dat -o filter_contact.dat -t %s" % args.filter)

      os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact.dat -c -cb reference_cb.pdb > component_contact.out 2>&1")
      os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact.dat -u -ub reference_hub.pdb -k %s > hubs_list.out 2>&1" % args.k)

    ## Create new folder for the h-bond calculations

    #print any(i in ["all", "mc-mc", "sc-sc", "mc-sc"] for i in args.hb)
    if any(i in ["all", "mc-mc", "sc-sc", "mc-sc"] for i in args.hb):
        os.chdir(path1)
        os.mkdir("h-bond")
        os.chdir(path1 + "/h-bond/")


        for i in args.hb:
            if i == "all":
                print("Pyinteraph analysis h-bonds all")


                os.system("touch nohup.out")
                print("Pyinteraph analysis h-bonds")
        ##as a default it will use all hbonds
                os.system("nohup pyinteraph -s ../reference0.gro -t ../../%s -r ../reference.fixed0.pdb -y --hb-graph hb-graph.dat  --ff-masses %s --hb-dat hb-dat.test > nohup.out 2>&1" %  (args.f, args.ff))
                if pbc:

                    os.system("nohup pyinteraph -s ../reference0.gro -t ../traj_mol_ur_nojump.xtc -r ../reference.fixed0.pdb -y --hb-graph hb-graph.dat  --ff-masses %s --hb-dat hb-dat.test > nohup.out 2>&1" %  args.ff)


                os.system("filter_graph -d hb-graph.dat -o hb_filter_contact.dat -t %s" % args.filter)
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact.dat -c -cb reference_cb.pdb > component_contact.out 2>&1")
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact.dat -u -ub reference_hub.pdb -k %s > hubs_list.out 2>&1" % args.k)

        ##if you want to do also sc-sc
        ## Last change: remove reference.gro due to problems in the analysis, now whit the reference.pdb seems to be ok

            if i == "sc-sc":
                print("Pyinteraph analysis h-bonds sc-sc")
                os.system(" nohup pyinteraph -s ../reference.fixed0.pdb  -t ../../%s -r ../reference.fixed0.pdb -y --hb-graph hb-graph_sc.dat  --ff-masses %s --hb-dat hb-dat_sc.test --hb-class sc-sc > nohup_sc.out 2>&1" % (args.f, args.ff))
                if pbc:

                    os.system(" nohup pyinteraph -s ../reference.fixed0.pdb  -t ../traj_mol_ur_nojump.xtc -r ../reference.fixed0.pdb -y --hb-graph hb-graph_sc.dat  --ff-masses %s --hb-dat hb-dat_sc.test --hb-class sc-sc > nohup_sc.out 2>&1" % args.ff)

                os.system("filter_graph -d hb-graph_sc.dat -o hb_filter_contact_sc.dat -t %s" % args.filter)
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact_sc.dat -c -cb reference_cb_sc.pdb > component_contact_sc.out 2>&1")
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact_sc.dat -u -ub reference_hub_sc.pdb -k %s > hubs_list_sc.out 2>&1" % args.k)


            if i == "mc-mc":
                print("Pyinteraph analysis h-bonds mc-mc")

                os.system(" nohup pyinteraph -s ../reference.fixed0.pdb  -t ../../%s -r ../reference.fixed0.pdb -y --hb-graph hb-graph_mc.dat  --ff-masses %s --hb-dat hb-dat_mc.test --hb-class mc-mc > nohup_mc.out 2>&1" % (args.f, args.ff))
                if pbc:
                    os.system(" nohup pyinteraph -s ../reference.fixed0.pdb  -t ../traj_mol_ur_nojump.xtc -r ../reference.fixed0.pdb -y --hb-graph hb-graph_mc.dat  --ff-masses %s --hb-dat hb-dat_mc.test --hb-class mc-mc > nohup_mc.out 2>&1" % args.ff)


                os.system("filter_graph -d hb-graph_mc.dat -o hb_filter_contact_mc.dat -t %s" % args.filter)
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact_mc.dat -c -cb reference_cb_mc.pdb > component_contact_mc.out 2>&1")
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact_mc.dat -u -ub reference_hub_mc.pdb -k %s > hubs_list_mc.out 2>&1" % args.k)


            if  i == "mc-sc":
                print("Pyinteraph analysis h-bonds mc-sc")

                os.system(" nohup pyinteraph -s ../reference.fixed0.pdb  -t ../../%s -r ../reference.fixed0.pdb -y --hb-graph hb-graph_sc.dat  --ff-masses %s --hb-dat hb-dat_mc-sc.test --hb-class mc-sc > nohup_mc-sc.out 2>&1" % (args.f, args.ff))
                if pbc:
                    os.system(" nohup pyinteraph -s ../reference.fixed0.pdb  -t ../traj_mol_ur_nojump.xtc -r ../reference.fixed0.pdb -y --hb-graph hb-graph_sc.dat  --ff-masses %s --hb-dat hb-dat_mc-sc.test --hb-class mc-sc > nohup_mc-sc.out 2>&1" % args.ff)


                os.system("filter_graph -d hb-graph_mc-sc.dat -o hb_filter_contact_mc-sc.dat -t %s" % args.filter)
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact_mc-sc.dat -c -cb reference_cb_mc-sc.pdb > component_contact_mc-sc.out 2>&1")
                os.system("nohup graph_analysis -r ../reference.fixed0.pdb -a hb_filter_contact_mc-sc.dat -u -ub reference_hub_mc-sc.pdb -k %s > hubs_list_mc-sc.out 2>&1" % args.k)

    else:
        print "No h-bonds calculated"
    ##  Knowledge-based potential calculation
    if args.kbp:

        os.chdir(path1)
        os.mkdir("kbp")
        os.chdir(path1 + "/kbp/")
        print("Knowledge-based potential calculation")
        os.system(" nohup pyinteraph -s ../reference0.gro  -t ../../%s -r ../reference.fixed0.pdb -p --kbp-graph kbp-graph.dat > nohup.out 2>&1" % args.f)
        if pbc:
            os.system(" nohup pyinteraph -s ../reference0.gro  -t ../traj_mol_ur_nojump.xtc -r ../reference.fixed0.pdb -p --kbp-graph kbp-graph.dat > nohup.out 2>&1")


      ### SPlitting the trajectory
    if args.jackknife and args.x:
        os.chdir(path)
        last_line = open(args.x).readlines()[-1]
        time_length = last_line.split()[0]
        s_time = int(round(float(time_length), -1))
        step= s_time/10
        n= 10
        for trajectories in range(step, s_time + step ,step):
            n-= 1
            if trajectories == step:
                print "resampling%s" % (n)
                os.mkdir(path1 + "/resampling%s" % (n))
                os.chdir(path1 + "/resampling%s" % (n))
                if pbc:
                    os.system("trjcat -f ../traj_mol_ur_nojump.xtc -b %s -e %s -o trajectories_frames%s-%s.xtc" % (0, s_time-step, 0, s_time-step))
                else:
                    os.system("trjcat -f ../../%s -b %s -e %s -o trajectories_frames%s-%s.xtc" % (args.f,0, s_time-step, 0, s_time-step))

            elif trajectories == s_time:
                print "trajectories_frames%s" % (n)
                os.mkdir(path1 + "/resampling%s" % (n))
                os.chdir(path1 + "/resampling%s" % (n))
                if pbc:
                    os.system("trjcat -f ../traj_mol_ur_nojump.xtc -b %s -e %s -o trajectories_frames%s-%s.xtc" % (step, s_time, step, s_time))
                else: 
                    os.system("trjcat -f ../../%s -b %s -e %s -o trajectories_frames%s-%s.xtc" % (args.f, step, s_time, step, s_time))

            else:
                os.mkdir(path1 + "/resampling%s" % (n))
                os.chdir(path1 + "/resampling%s" % (n))
                if pbc:
                    os.system("trjcat -f ../traj_mol_ur_nojump.xtc -b %s -e %s -o trajectories_frames%s-%s.xtc" % (0, s_time-trajectories, 0, s_time-trajectories))
                    os.system("trjcat -f ../traj_mol_ur_nojump.xtc -b %s -e %s -o trajectories_frames%s-%s.xtc" % (s_time-trajectories+ step,  s_time, s_time-trajectories+ step, s_time))
                else:
                    os.system("trjcat -f ../../%s -b %s -e %s -o trajectories_frames%s-%s.xtc" % (args.f, 0, s_time-trajectories, 0, s_time-trajectories))
                    os.system("trjcat -f ../../%s -b %s -e %s -o trajectories_frames%s-%s.xtc" % (args.f, s_time-trajectories+ step,  s_time, s_time-trajectories+ step, s_time))

                os.system("trjcat -f trajectories_frames%s-%s.xtc trajectories_frames%s-%s.xtc -o trajectories_frames%s-%s.xtc" % (0, s_time-trajectories, s_time-trajectories+ step, s_time,0, s_time-trajectories ))

            path2 = os.getcwd()
            os.mkdir(path2 + "/contact")
            os.chdir(path2 + "/contact/")
            os.mkdir(args.cut_off)
            os.chdir(path2 + "/contact/%s/" % args.cut_off)
            os.system("touch nohup.out")

            print("Creating the pdb with all the trajectories")
            execfile("/usr/local/envs/pyinteraph/bin/activate_this.py", dict(__file__="/usr/local/envs/pyinteraph/bin/activate_this.py"))


  #            os.system("nohup pyinteraph -v -s reference0.gro -t traj_fixed.xtc  -r reference.fixed0.pdb -f --hc-graph contact.dat --hc-co 5 --ff-masses charmm27 --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS > nohup.out", shell = True)
            if trajectories == s_time:
                os.system("nohup pyinteraph -v -s ../../../reference0.gro -t ../../trajectories_frames%(1)d-%(2)s.xtc -r ../../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(3)s --ff-masses %(4)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": step, "2": s_time, "3":args.cut_off,"4":args.ff})
            else:
                os.system("nohup pyinteraph -v -s ../../../reference0.gro -t ../../trajectories_frames%(1)d-%(2)s.xtc -r ../../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(3)s --ff-masses %(4)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": 0, "2": s_time-trajectories, "3":args.cut_off,"4":args.ff})

            os.system("filter_graph -d contact.dat -o filter_contact.dat -t %s" % args.filter)
            

          ###make a copy of the pdb becasue replace columns
              #os.system("nohup graph_analysis -r ../../../reference.fixed0.pdb -a filter_contact.dat -c -cb reference_hub.pdb > nohup.out")
            os.system("nohup graph_analysis -r ../../../reference.fixed0.pdb -a filter_contact.dat -c -cb reference_cb.pdb > component_contact.out 2>&1")
            os.system("nohup graph_analysis -r ../../../reference.fixed0.pdb -a filter_contact.dat -u -ub reference_hub.pdb -k %s > hubs_list.out 2>&1" % args.k)
          #
            if args.cut_off2 > 0:
                os.chdir(path2 + "/contact/")
                os.mkdir(str(args.cut_off2))
                os.chdir(path2 + "/contact/%s/" % args.cut_off2)
          #
                os.system("touch nohup.out")

                if trajectories == s_time:
                    os.system("nohup pyinteraph -v -s ../../../reference0.gro -t ../../trajectories_frames%(1)d-%(2)s.xtc  -r ../../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(3)s --ff-masses %(4)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": step, "2": s_time, "3":args.cut_off2,"4":args.ff})
                else:
                    os.system("nohup pyinteraph -v -s ../../../reference0.gro -t ../../trajectories_frames%(1)d-%(2)s.xtc  -r ../../../reference.fixed0.pdb -f --hc-graph contact.dat --hc-co %(3)s --ff-masses %(4)s --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN > nohup.out 2>&1" % {"1": 0, "2": s_time-trajectories, "3":args.cut_off2,"4":args.ff})
          #
                os.system("filter_graph -d contact.dat -o filter_contact.dat -t %s" % args.filter)

                os.system("nohup graph_analysis -r ../../../reference.fixed0.pdb -a filter_contact.dat -c -cb reference_cb.pdb > component_contact.out 2>&1")
                os.system("nohup graph_analysis -r ../../../reference.fixed0.pdb -a filter_contact.dat -u -ub reference_hub.pdb -k %s > hubs_list.out 2>&1" % args.k)

            if any(i in ["all", "mc-mc", "sc-sc", "mc-sc"] for i in args.hb):
                os.chdir(path2)
                os.mkdir("h-bond")
                os.chdir(path2 + "/h-bond/")
         #
                os.system("touch nohup.out")
         ##as a default it will use all hbonds
                for i in args.hb:
                    if i == "all":
                        print("Pyinteraph analysis h-bonds all")
                        if trajectories == s_time:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph.dat  --ff-masses %(3)s --hb-dat hb-dat.test > nohup.out 2>&1" % {"1":step, "2":s_time,"3":args.ff})
                        else:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph.dat  --ff-masses %(3)s --hb-dat hb-dat.test > nohup.out 2>&1" % {"1":0, "2":s_time-trajectories,"3":args.ff})

                        os.system("filter_graph -d hb-graph.dat -o filter_contact.dat -t %s" % args.filter)
         #
             ###make a copy of the pdb becasue replace columns
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact.dat -c -cb reference_cb.pdb > component_contact.out 2>&1")
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact.dat -u -ub reference_hub.pdb -k %s > hubs_list.out 2>&1" % args.k)

                    if i == "sc-sc":
                        print("Pyinteraph analysis h-bonds sc-sc")
                        if trajectories == s_time:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph_sc.dat  --ff-masses %(3)s --hb-dat hb-dat_sc.test  --hb-class sc-sc > nohup_sc.out 2>&1" % {"1":step, "2":s_time,"3":args.ff})
                        else:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph_sc.dat  --ff-masses %(3)s --hb-dat hb-dat_sc.test --hb-class sc-sc > nohup_sc.out 2>&1" % {"1":0, "2":s_time-trajectories,"3":args.ff})

                        os.system("filter_graph -d hb-graph_sc.dat -o filter_contact_sc.dat -t %s" % args.filter)
         #
             ###make a copy of the pdb becasue replace columns
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact_sc.dat -c -cb reference_cb_sc.pdb > component_contact_sc.out 2>&1")
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact_sc.dat -u -ub reference_hub_sc.pdb -k %s > hubs_list_sc.out 2>&1" % args.k)

                    if i == "mc-mc":
                        print("Pyinteraph analysis h-bonds mc-mc")
                        if trajectories == s_time:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph_mc.dat  --ff-masses %(3)s --hb-dat hb-dat_mc.test  --hb-class mc-mc > nohup_mc.out 2>&1" % {"1":step, "2":s_time,"3":args.ff})
                        else:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph_mc.dat  --ff-masses %(3)s --hb-dat hb-dat_mc.test  --hb-class mc-mc > nohup_mc.out 2>&1" % {"1":0, "2":s_time-trajectories,"3":args.ff})

                        os.system("filter_graph -d hb-graph_mc.dat -o filter_contact_mc.dat -t %s" % args.filter)
         #
             ###make a copy of the pdb becasue replace columns
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact_mc.dat -c -cb reference_cb_mc.pdb > component_contact_mc.out 2>&1")
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact_mc.dat -u -ub reference_hub_mc.pdb -k %s > hubs_list_mc.out 2>&1" % args.k)

                    if i == "mc-sc":
                        print("Pyinteraph analysis h-bonds mc-sc")
                        if trajectories == s_time:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph_mc-sc.dat  --ff-masses %(3)s --hb-dat hb-dat_mc-sc.test --hb-class mc-sc > nohup_mc-sc.out 2>&1" % {"1":step, "2":s_time,"3":args.ff})
                        else:
                            os.system("nohup pyinteraph -s ../../reference0.gro -t ../trajectories_frames%(1)d-%(2)s.xtc -r ../../reference.fixed0.pdb -y --hb-graph hb-graph_mc-sc.dat  --ff-masses %(3)s --hb-dat hb-dat_mc-sc.test --hb-class mc-sc > nohup_mc-sc.out 2>&1" % {"1":0, "2":s_time-trajectories,"3":args.ff})

                        os.system("filter_graph -d hb-graph_mc-sc.dat -o filter_contact_mc-sc.dat -t %s" % args.filter)
         #
             ###make a copy of the pdb becasue replace columns
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact_mc-sc.dat -c -cb reference_cb_mc-sc.pdb > component_contact_mc-sc.out 2>&1")
                        os.system("nohup graph_analysis -r ../../reference.fixed0.pdb -a filter_contact_mc-sc.dat -u -ub reference_hub_mc-sc.pdb -k %s > hubs_list_mc-sc.out 2>&1" % args.k)
            os.chdir(path1 + "/resampling%s" % (n))
            if trajectories == step:
                os.remove("trajectories_frames%s-%s.xtc" % (0, s_time-step))
            elif trajectories == s_time:
                os.remove("trajectories_frames%s-%s.xtc" % (step, s_time))
            else:
                os.remove("trajectories_frames%s-%s.xtc" % (0, s_time-trajectories))
                os.remove("trajectories_frames%s-%s.xtc" % (s_time-trajectories+ step,  s_time))


elapsed = (time.time() - start)
        #print "Seconds %s" % (elapsed)
print "%d:%02dmn" % (elapsed / 60, elapsed % 60)

print "Process finished"
sys.exit()


