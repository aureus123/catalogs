
# Makefile - Made in 2024 by Daniel E. Severin.

CC = g++

CCOPT = -Wall
CCLNFLAGS = -L wcstools-3.9.7/libwcs/ -lwcs

all: compare_ppm compare_agk compare_cpd compare_ppm_bd compare_gc compare_sd compare_cd find_coord mag_cd mag_bd

mag_cd: mag_cd.o read_cd.o read_ppm.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

mag_cd.o: mag_cd.cpp
	$(CC) $(CCFLAGS) -c $<

compare_agk: compare_agk.o read_cd.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_agk.o: compare_agk.cpp
	$(CC) $(CCFLAGS) -c $<

compare_gc: compare_gc.o read_cd.o read_old.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_gc.o: compare_gc.cpp
	$(CC) $(CCFLAGS) -c $<

compare_cpd: compare_cpd.o read_cd.o read_cpd.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_cpd.o: compare_cpd.cpp
	$(CC) $(CCFLAGS) -c $<

compare_ppm: compare_ppm.o read_cd.o read_ppm.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_ppm.o: compare_ppm.cpp
	$(CC) $(CCFLAGS) -c $<

compare_sd: compare_sd.o read_cd.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_sd.o: compare_sd.cpp
	$(CC) $(CCFLAGS) -c $<

compare_cd: compare_cd.o read_cd.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_cd.o: compare_cd.cpp
	$(CC) $(CCFLAGS) -c $<

compare_ppm_bd: compare_ppm_bd.o read_bd.o read_ppm.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

compare_ppm_bd.o: compare_ppm_bd.cpp
	$(CC) $(CCFLAGS) -c $<

find_coord: find_coord.o read_cd.o read_old.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

find_coord.o: find_coord.cpp
	$(CC) $(CCFLAGS) -c $<

read_cd.o: read_cd.cpp
	$(CC) $(CCFLAGS) -c $<

read_old.o: read_old.cpp
	$(CC) $(CCFLAGS) -c $<

mag_bd: mag_bd.o read_bd.o read_ppm.o trig.o misc.o
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)

mag_bd.o: mag_bd.cpp
	$(CC) $(CCFLAGS) -c $<

read_bd.o: read_bd.cpp
	$(CC) $(CCFLAGS) -c $<

read_ppm.o: read_ppm.cpp
	$(CC) $(CCFLAGS) -c $<

read_cpd.o: read_cpd.cpp
	$(CC) $(CCFLAGS) -c $<

trig.o: trig.cpp
	$(CC) $(CCFLAGS) -c $<

misc.o: misc.cpp
	$(CC) $(CCFLAGS) -c $<

.PHONY: clean

clean:
	rm -f *.o
	rm -f compare_ppm compare_agk compare_cpd compare_ppm_bd compare_gc compare_sd compare_cd find_coord mag_cd mag_bd
