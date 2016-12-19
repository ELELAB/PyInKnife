#!/bin/bash

## Script created to fix the error for heatmap figures from the PyInKnife
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

awk '{s= substr($9, 1, 4)}{g=substr($9,5, length($9))}{print $1,$2,$3,$4,$5,$6,$7,$8,s,g,$10, $11}' contact/4.5/reference_cb.pdb > contact/4.5/reference_cb.txt 
for i in 0 1 2 3 4 5 6 7 8 9 
do 
  awk '{s= substr($9, 1, 4)}{g=substr($9,5, length($9))}{print $1,$2,$3,$4,$5,$6,$7,$8,s,g,$10, $11}' resampling"$i"/contact/4.5/reference_cb.pdb > resampling"$i"/contact/4.5/reference_cb.txt
done

