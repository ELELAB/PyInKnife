#!/bin/bash

## Script created to plot the distribution figures from the PyInKnife
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

cd contact/$1/

all_plot_hubs.py hubs_list.out ../../resampling*/contact/$1/hubs_list.out -x ../../$3 -hb $4


all_plot_cmap.py component_contact.out ../../resampling*/contact/$1/component_contact.out -x ../../$3 -hb $4

all_hubs_distribution.R distribution_hubs.txt

all_cmap_distribution.R final_cmap.txt 

cd ../$2/

all_plot_hubs.py hubs_list.out ../../resampling*/contact/$2/hubs_list.out -x ../../$3 -hb $4


all_plot_cmap.py component_contact.out ../../resampling*/contact/$2/component_contact.out -x ../../$3 -hb $4

all_hubs_distribution.R distribution_hubs.txt

all_cmap_distribution.R final_cmap.txt 

if [ -d "../../h-bond/" ]; then

	cd ../../h-bond
	all_plot_hubs.py hubs_list.out ../resampling*/h-bond/hubs_list.out -x ../$3 -hb $4


	all_plot_cmap.py component_contact.out ../resampling*/h-bond/component_contact.out -x ../$3 -hb $4

	all_hubs_distribution.R distribution_hubs.txt

	all_cmap_distribution.R final_cmap.txt 

	if [ $5 ]; then

		all_plot_hubs.py hubs_list_sc.out ../resampling*/h-bond/hubs_list_sc.out -x ../$3 -hb $5

		all_plot_cmap.py component_contact_sc.out ../resampling*/h-bond/component_contact_sc.out -x ../$3 -hb $5

		all_hubs_distribution.R distribution_hubs_sc.txt $5

		all_cmap_distribution.R final_cmap_sc.txt $5
	fi
fi
