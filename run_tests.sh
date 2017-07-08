
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




echo "Running on Big Data"

echo "Running for k=5"

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b1.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b1_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b1_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b2.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b2_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b2_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b3.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b3_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b3_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b4.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b4_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b4_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b5.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b5_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b5_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b6.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b6_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b6_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b7.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b7_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b7_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b8.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b8_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b8_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b9.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b9_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b9_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b10.geno -k 5 -l 5 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b10_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b10_logfile.txt


echo "Running for k=10"

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b1.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b1_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b1_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b2.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b2_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b2_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b3.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b3_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b3_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b4.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b4_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b4_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b5.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b5_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b5_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b6.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b6_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b6_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b7.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b7_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b7_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b8.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b8_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_10/b8_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b9.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b9_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b9_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/b10.geno -k 10 -l 10 -m 1000 -a -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b10_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5/b10_logfile.txt

echo "Done"
