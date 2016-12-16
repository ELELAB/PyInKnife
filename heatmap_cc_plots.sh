#!/bin/bash

## Script created to plot the heatmap figures from the PyInKnife
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

awk 'FNR==1{f++}{a[f,FNR]=$10}END{for(x=1;x<=FNR;x++){for(y=1;y<ARGC;y++)printf("%s ",a[y,x]);print ""}}' resampling*/contact/$1/reference_cb.pdb > test_all_cc_$2.txt

awk '{print $4, $5, $10}' contact/$1/reference_cb.pdb > test_cc_traj_$2.txt

paste -d ' ' test_cc_traj_$2.txt test_all_cc_$2.txt > cc_merged_$2.txt

Rscript heatmap_cc_traj.R cc_merged_$2.txt heatmap_cc_$2_$3_$4.pdf $3 $4

Rscript heatmap_cc_traj.R cc_merged_$2.txt heatmap_cc_$2_$5_$6.pdf $5 $6

#Rscript /data/user/juan/pyinteraph_elena/final_version/heatmap_cc_traj.R cc_merged_$2.txt heatmap_cc_$2_$7_$8.pdf $7 $8