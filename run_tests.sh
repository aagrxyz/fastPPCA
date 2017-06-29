cd /home/aman/Desktop/Intern/fast_em_pca/build

rm -rf *
cmake ..
make

echo "Running on Small Data"

echo "Running for k=5"

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s1.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s1_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s1_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s2.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s2_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s2_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s3.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s3_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s3_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s4.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s4_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s4_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s5.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s5_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s5_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s6.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s6_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s6_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s7.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s7_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s7_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s8.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s8_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s8_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s9.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s9_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s9_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s10.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s10_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s10_logfile.txt


echo "Running for k=10"

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s1.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s1_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s1_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s2.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s2_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s2_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s3.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s3_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s3_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s4.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s4_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s4_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s5.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s5_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s5_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s6.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s6_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s6_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s7.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s7_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s7_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s8.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s8_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_10/s8_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s9.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s9_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s9_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/s10.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s10_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5/s10_logfile.txt



