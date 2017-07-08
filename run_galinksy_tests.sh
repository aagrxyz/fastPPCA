
cd /home/aman/Desktop/Intern/fast_em_pca/build

rm -rf *
cmake ..
make

# echo "Running on Small Data"

# echo "Running for k=5 , m =20"

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s1.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s1_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s1_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s2.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s2_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s2_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s3.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s3_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s3_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s4.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s4_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s4_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s5.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s5_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s5_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s6.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s6_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s6_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s7.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s7_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s7_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s8.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s8_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s8_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s9.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s9_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s9_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s10.geno -k 5 -l 15 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s10_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_20/s10_logfile.txt




echo "Running on Big Data"

echo "Running for k=5 , m=20"

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b1.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b1_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b1_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b2.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b2_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b2_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b3.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b3_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b3_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b4.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b4_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b4_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b5.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b5_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b5_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b6.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b6_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b6_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b7.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b7_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b7_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b8.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b8_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b8_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b9.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b9_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b9_logfile.txt

time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b10.geno -k 5 -l 5 -m 20 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b10_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_20/b10_logfile.txt





# echo "Running on Small Data"

# echo "Running for k=5 , m =40"

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s1.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s1_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s1_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s2.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s2_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s2_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s3.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s3_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s3_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s4.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s4_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s4_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s5.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s5_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s5_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s6.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s6_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s6_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s7.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s7_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s7_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s8.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s8_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s8_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s9.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s9_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s9_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s10.geno -k 5 -l 15 -m 40 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s10_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_15_m_40/s10_logfile.txt




# echo "Running on Big Data"

# echo "Running for k=5 , m=50"

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b1.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b1_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b1_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b2.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b2_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b2_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b3.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b3_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b3_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b4.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b4_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b4_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b5.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b5_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b5_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b6.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b6_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b6_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b7.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b7_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b7_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b8.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b8_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b8_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b9.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b9_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b9_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/big-data/geno/b10.geno -k 5 -l 5 -m 50 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b10_ > /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/big-data/em/k_5_m_50/b10_logfile.txt




# echo "Running on Small Data"

# echo "Running for k=5 , m =100"

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s1.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s1_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s1_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s1_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s2.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s2_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s2_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s2_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s3.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s3_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s3_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s3_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s4.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s4_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s4_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s4_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s5.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s5_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s5_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s5_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s6.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s6_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s6_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s6_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s7.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s7_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s7_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s7_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s8.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s8_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s8_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s8_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s9.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s9_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s9_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s9_logfile.txt

# time (./em_mailman -g /home/aman/Desktop/Intern/simulator/small-data/geno/s10.geno -k 5 -l 15 -m 100 -cl 0.001 -o /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s10_ > /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s10_logfile.txt ) 2>> /home/aman/Desktop/Intern/simulator/small-data/em/k_5_m_100/s10_logfile.txt
