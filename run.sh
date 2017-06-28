#!/bin/bash
cd bin
echo "Running em_mailman on small"
time (./em_mailman -g /home/aman/Desktop/Intern/data_em/timing_comparison/small/genotype.txt -k 2 -m 20 -v -o /home/aman/Desktop/Intern/data_em/timing_comparison/small/ > logfile_mailman_small.txt) 2>> logfile_mailman_small.txt
echo "Running em_naive on small"
time (./em_naive -g /home/aman/Desktop/Intern/data_em/timing_comparison/small/genotype.txt -k 2 -m 20 -v -a -o /home/aman/Desktop/Intern/data_em/timing_comparison/small/ > logfile_naive_small.txt) 2>> logfile_naive_small.txt

echo "Running em_mailman on medium"
time (./em_mailman -g /home/aman/Desktop/Intern/data_em/timing_comparison/medium/genotype.txt -k 2 -m 20 -v -o /home/aman/Desktop/Intern/data_em/timing_comparison/medium/ > logfile_mailman_medium.txt) 2>> logfile_mailman_medium.txt
echo "Running em_naive on medium"
time (./em_naive -g /home/aman/Desktop/Intern/data_em/timing_comparison/medium/genotype.txt -k 2 -m 20 -v -o /home/aman/Desktop/Intern/data_em/timing_comparison/medium/ > logfile_naive_medium.txt) 2>> logfile_naive_medium.txt

# echo "Running em_mailman on big"
# time (./em_mailman -g /home/aman/Desktop/Intern/data_em/timing_comparison/big/genotype.txt -k 2 -m 20 -v -o /home/aman/Desktop/Intern/data_em/timing_comparison/big/ > logfile_mailman_big.txt) 2>> logfile_mailman_big.txt
# echo "Running em_naive on big"
# time (./em_naive -g /home/aman/Desktop/Intern/data_em/timing_comparison/big/genotype.txt -k 2 -m 5 -v -o /home/aman/Desktop/Intern/data_em/timing_comparison/big/ > logfile_naive.txt) 2>> logfile_naive.txt

# echo "Done"