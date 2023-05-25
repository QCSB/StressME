#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-laurence
#SBATCH --mem=6G

#python foldme_test.py
python simulation.py $1 $2 $3
