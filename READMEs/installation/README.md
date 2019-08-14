# Installing snakemake

### Well... first you need to install conda...

Step 1. Get conda
```
curl -O https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh
bash Anaconda3-5.3.1-Linux-x86_64.sh
```
Then log out and in.

Step 2. Update conda
```
conda update conda
conda update anaconda
```

Now you're set!!!

### Installing snakemake with conda

In short:
```
conda create -n snakemake
conda activate snakemake
conda install -c bioconda -c conda-forge snakemake
```
When finished
```
conda deactivate
```
* note: the first time activating a conda environment, you may need to specify the language that is being used. If this is the case, select `bash`
* note: this means that you will need to activate the snakemake conda environement every time you use snakemake with `conda activate snakemake`

### Installing snakemake with pip

If you have root access:
```
pip3 install snakemake
```
if not...
```
pip3 install snakemake --user
```

### If these options don't work...
Read the official documentation... 

The official documentation can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

# Installing singularity
The official documentation can be found [here](https://singularity.lbl.gov/install-linux) and [here](https://sylabs.io/guides/3.0/user-guide/installation.html).

On Centos
```
yum install singularity
```
On Ubuntu
```
apt-get install -y singularity-container
```
If you don't have root or sudo privileges, ask your administrator to install singularity.
