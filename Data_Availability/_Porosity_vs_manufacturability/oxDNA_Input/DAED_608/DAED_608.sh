#!/bin/sh

current_dir=$(pwd)
oxDNA /home/mmbl/Documents/avetturi/run_oxDNA/08-17_23-15/DAED_608/DAED_608_first_min
mv "$current_dir/last_conf.dat" /home/mmbl/Documents/avetturi/run_oxDNA/08-17_23-15/DAED_608/first_min_DAED_608.conf
oxDNA /home/mmbl/Documents/avetturi/run_oxDNA/08-17_23-15/DAED_608/DAED_608_first_relax
mv "$current_dir/last_conf.dat" /home/mmbl/Documents/avetturi/run_oxDNA/08-17_23-15/DAED_608/first_relax_DAED_608.conf
oxDNA /home/mmbl/Documents/avetturi/run_oxDNA/08-17_23-15/DAED_608/DAED_608_simulation
mv "$current_dir/last_conf.dat" /home/mmbl/Documents/avetturi/run_oxDNA/08-17_23-15/DAED_608/simulated_DAED_608.conf
