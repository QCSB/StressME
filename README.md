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

salloc --time=2:0:0 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --account=<your_account> python StressME_wildtype.py 42 5.0 10, 
or: 
sbatch --mem=8G --account=<your_account> --time=2:00:00 --output StressME_wildtype StressME_wildtype.sh, 
where StressME_wildtype.sh is coded as: 
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=<your_account>
python StressME_wildtype.py 42 5.0 10

Here “42 5.0 10” refers to the triple stress conditions at temperature 42 ℃, pH 5.0 and ROS 10X of the basal level. 

# StressME using Docker

Docker allows ME-model users to run StressME locally without going through the complicated processes of installing solvers and dependencies that may be incompatible with each other due to their version update.  
The Docker image of StressME can be found on Docker Hub (queensysbio/stressme:v1.1). This build was developed from the modified version of COBRAme and EcoliME kernel to integrate FoldME and AcidifyME with OxidizeME. This build also includes qMINOS solver that users can use to solve StressME using solvemepy.
The installed COBRAme (StressME version), ECOLIme (StressME version), AcidifyME, OxidizeME and solvemepy packages can be found in /source/ after a new container has been created from the image of StressME. The working directory is /home/meuser, where ME-model users can run simulations and export output from the StressME container to the host. 

## Installation on Windows Subsystem for Linux (WSL2)

Steps to run a Docker container from the image of StressME with everything required to run the StressME using qMINOS:

(a)	Install Ubuntu on WSL2 on Windows 10/11. 
https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview
(b)	Download Docker Desktop for Windows https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe
(c)	Start Docker Desktop from the Windows Start menu. From the Docker menu, select Settings and then General. Select the Use WSL 2 based engine check box. Select Apply & Restart
(d)	Pull the latest StressME image from DockerHub with:
docker pull queensysbio/stressme:v1.1.

### Simulations on Windows Subsystem for Linux (WSL2)

(a)	Start Docker Desktop from the Windows Start menu
(b)	Start another PowerShell window from the Windows Start menu. In the command line, type “ubuntu”
(c)	Type “cd /mnt/c/Users/<your user name in this local computer>/Desktop” and then type “mkdir StressME_wild_simulations” and “mkdir StressME_heat_simulations” to create host directories to receive output from the Docker StressME container
(d)	Type “cd StressME_wild_simulations” or “cd StressME_heat_simulations” and then type “docker run -p 8888:8888 --rm -i -v $(pwd):/mount_point/ -t queensysbio/stressme:v1.1 bash” . This will start the Docker StressME container (virtual machine) in the working directory: /home/meuser. 
(e)	Type “ls” to check files available in the working directory.
(f)	Type “python StressME_wildtype.py 42 5.0 10” or “python StressME_heatevolved.py 42 5.0 10” to run simulations for the wild type or the heat evolved strains. Here “42 5.0 10” refers to the triple stress conditions at temperature 42 ℃, pH 5.0 and ROS 10X of the basal level. 
(g)	After simulations are done, three csv files (TripleStressME_proteome.csv, TripleStressME_phenotypes.csv, and TripleStressME_fluxes.csv) are found in the working directory /home/meuser, representing the protein mass fractions, metabolic fluxes and phenotypes under the test stress conditions. Type “cp *.csv /mount_point/ StressME_wild_simulations” or “cp *.csv /mount_point/ StressME_heat_simulations” to export results to the corresponding host directories. 
(h)	Use “CTRL+D” to exit Docker StressME container to go back to the host directory. 
(i)	Type “explorer.exe .” to open File Explorer to view (e.g. EXCEL) and process (e.g., R for Rstudio or Matlab) the csv data. 

### Simulations on Jupiter notebook

(a)	Start Docker Desktop from the Windows Start menu
(b)	In WSL2 or powershell, run ‘jupyter notebook --ip=0.0.0.0 --port=8888 --allow-root’
(c)	Copy the random token that is generated on the screen.
(d)	Paste the token on http://localhost:8888/ to launch jupyter notebook
(e)	Run StressME_heatevolved.ipynb and StressME_wildtype.ipynb for simulations.
