./compare_ppm cd.txt >results/log_ppm.log
./compare_agk >results/log_agk.log
./compare_cpd >results/log_cpd.log
./compare_sd >results/log_sd.log

cd results
cat table_pos_ppm.csv > table_pos.csv
cat table_pos_agk.csv | grep -v "decl" >> table_pos.csv
cat table_pos_cpd.csv | grep -v "decl" >> table_pos.csv
cat table_pos_sd.csv | grep -v "decl" >> table_pos.csv

cat table_mag_ppm.csv > table_mag.csv
cat table_mag_agk.csv | grep -v "decl" >> table_mag.csv
cat table_mag_sd.csv | grep -v "decl" >> table_mag.csv
