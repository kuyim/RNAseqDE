# RNAseq project
Opened VSCode terminal through Ubuntu by running:
```
mkdir -p ~/rnaseqde
cd ~/rnaseqde
code .
```

## Setup
1. Downloaded Nextflow using:
```
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```
Nextflow is a JVM tool and requires Java 17+. \
2. Installed, created, and activated venv:
```
sudo apt install python3-venv -y
python3 -m venv .venv
activate .venv/bin/activate
```
Downloading the nf-core CLI within virtual environment is preferred as global package installation can potentially break system tools. Not necessarily a huge issue in this case since most of the packages used here are not Python. \
3. Install nf-core:
```
pip install nf-core
```

## Dataset
# Wrong Dataset! Don't use this section
Realized that the GSE86337 dataset doesn't have treated vs controls, it is just 5 different cell lines, so DE analysis is not useful
### Dataset Identity: GSE86337 from NIH SRA Run Selector
Load the ```matrix/ minimal/ soft/ ``` and ``` suppl/``` data from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86337/ using Filezilla. These contain the metadata and counts matrix, not the raw FASTQ files.
### Matrix
```GSExxx_series_matrix.txt.gz or GSExxx-GPLxxx_series_matrix.txt.gz``` \
gzipped Series-matrix files

Series_matrix files are summary text files that include a tab-delimited value-matrix table generated from the 'VALUE' column of each Sample record, headed by Sample and Series metadata. These files include SOFT attribute labels. Data generated from multiple Platforms are contained in separate files. It is recommended to view Series_matrix files in a spreadsheet application like Excel.

CAUTION: value data are extracted directly from the original records with no consideration as to whether the values are directly comparable.
### Minimal
```GSExxx_family.xml.tgz``` \
tarred gzipped MINiML files by Series (GSE)
GSExxx_family files contain MINiML-formatted data for all Platforms (GPL) and Samples (GSM) associated with one Series (GSE).
### Soft
```GSExxx_family.soft.gz``` \
gzipped SOFT files by Series (GSE)
GSExxx_family files contain SOFT-formatted data for all Platforms (GPL)
and Samples (GSM) associated with one Series (GSE).
### Suppl
``` GSExxx_RAW.tar ``` \
tarred files for all Sample supplementary files corresponding to a Series, as well as any additional files the submitter wants make available.

All submitters have been asked to provide supplementary data (for example, Affymetrix
.CEL files or cDNA array .GPR files) to accompany their GEO records.  If you are unable
to locate supplementary data for your experiment of interest, we suggest that you contact
the submitter directly to encourage that they supply raw data files to GEO so that we may
make them available to the scientific community.

If you are interested in locating all instances of a particular file type, we
suggest that you use Entrez GEO DataSets at
http://www.ncbi.nlm.nih.gov/gds/.  For example, to locate all .cel files
corresponding to Affymetrix HG-U133A array that has GEO accession GPL96, search with:
GPL96 AND "cel"[Supplementary Files]

### Unzip files
```
gunzip data/GSE86337_series_matrix.txt.gz
gunzip data/GSE86337_family.soft.gz
gunzip data/GSE86337_reverse.stranded.unfiltered.count.matrix.txt.gz
tar -xzvf  data/GSE86337_family.xml.tgz -C data/
```

### Located raw FASTQ data
Used the following link with the SRA study accession number on the NCBI website for this study
https://www.ebi.ac.uk/ena/browser/text-search?query=SRP083954

Downloaded the fastq.gz files using the website, but they are also accessible via wget. Alternatively, nf-core supports feeding in the SRA study accession number and it will fetch the fastq files.
# Good after this line
# ***************************************************

### Created accession list for test human airway dataset
```
mkdir -p assets
touch assets/airway_ids.csv
```
This dataset GSE52778 with SRA study number SRP033351 contains 4 primary cell lines of the human airway muscle transcriptome at baseline and under treatment with dexamethasone 1 uM (corticosteroid) for 18 hours. This is a subset of the full study, which uses 16 libraries and contains albuterol arms.

### Created config file
```
mkdir -p conf
touch conf/run.config
```

## Resolved some issues with Docker
Saved Docker config to a backup
```
cp ~/.docker/config.json ~/.docker/config.json.bak
```
Docker config was showing:
```
{
  "credsStore": "desktop.exe"
}
```
Replaced with:
```
{"auths": {}}
```

## Nextflow Version Issues
Within ```run_pipeline.sh```, Nextflow version ```export NXF_VER=24.10.5``` worked for fetchngs but not for rnaseq. The fix was to update the line containing Nextflow version to ```export NXF_VER=25.04.6```.

## Increased Memory Allocation for WSL to Allow STAR Alignment
Created ```.wslconfig``` file in ```%UserProfile%```:
```
[ws12]
memory=25GB
processors=8
```
Then reset wsl:
```
wsl --shutdown
```
Confirm if Linux sees the higher limit:
```
free -g
nproc
```
Updated ```resourceLimits``` to 18GB memory in ```run.config``` to reflect the higher RAM. Initially was set to 26GB but was having issues with crashing WSL, so lowered it to 18.

## Even with Increased Memory Allocation, STAR not working, Added the Following Lines to run.config
```
params {
    pseudo_aligner = 'salmon'
    skip_alignment = true
}
```
Kept getting crashes and infinite reconnecting to WSL:Ubuntu loops.
Reset to 6 CPUs and 14GB RAM in config file, 24GB, 8 processors, 12GB swap in .wslconfig.
Main fix was moving Ubuntu to D drive and freeing up 190GB of space from ext4.vhdx file that was not clearing properly.

After migrating Ubuntu to D drive, was able to run the pipeline on all 8 samples but somehow overnight, the ext4 file was corrupted and forced the file system to read only. Copied important files back to Windows to complete R differential expression analysis.