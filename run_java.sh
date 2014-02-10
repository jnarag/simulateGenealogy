#!/bin/tcsh
#
#$ -S /bin/tcsh -cwd
#$ -o run1.out -j y
#$ -l mem_free=10G
java -Xmx10000m -jar simulateGenealogy.jar
