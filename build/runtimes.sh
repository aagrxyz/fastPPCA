#!/bin/bash
cd /home/aman/Desktop/Intern/fastPPCA/build
cmake ..
make
echo "Running eff"
time (./em_eff -p par_eff.txt > eff_log.txt) 2>> eff_log.txt

echo "Running not eff"
time (./fastppca -p par.txt > log.txt) 2>> log.txt

echo "Done"

