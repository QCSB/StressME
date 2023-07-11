# StressME using Linux (clusters)
## Installation
### To set up virtual environment for python 3.6

virtualenv -p python MePython363 
source MePython363/bin/activate
### install dependencies
Accessories		Version			Installation
(a)	cython			0.28.2				pip
(b)	sympy			1.1.1				pip
(c)	numpy			1.14.3				pip
(d)	scipy			1.1.0				pip
(e)	pytest			3.5.1				pip
(f)	pandas			0.22.0				pip
(g)	cycler			0.11.0				pip
(h)	matplotlib		2.2.2				pip
(i)	biopython		1.76				pip
(j)	qMINOS* 		5.6 				https://github.com/SBRG/solvemepy
(k)	cobrame		StressME 1.1			https://github.com/QCSB/StressME	
[l] ecolime			StressME 1.1 			https://github.com/QCSB/StressME
[m] oxidizeme		StressME 1.1			https://github.com/QCSB/StressME
[n] acidifyme		StressME 1.1			https://github.com/QCSB/StressME
[o] meuser 		StressME 1.1	 		https://github.com/QCSB/StressME

*See https://github.com/SBRG/solvemepy for qMINOS installation

## Simulations on Linux clusters by slurm

salloc --time=2:0:0 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --account=<your_account> python StressME_wildtype.py 42 5.0 10
or: 
sbatch --mem=8G --account=<your_account> --time=2:00:00 --output StressME_wildtype StressME_wildtype.sh
where StressME_wildtype.sh is coded as: 
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=<your_account>
python StressME_wildtype.py 42 5.0 10

Here “42 5.0 10” refers to the triple stress conditions at temperature 42 ℃, pH 5.0 and ROS 10X of the basal level. 
