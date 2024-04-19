# RNAcompete Analysis Pipeline
A new and efficient Python pipeline to analyze RNAcompete data.

## Instructions

### Step 1. Create the Conda environment
**Conda** is required to install  the dependencies of the pipeline.
```
conda env create -f scripts/environment.yml
conda activate rnacompete
```

### Step 2. Prepare the input files
RNAcompete experiments are **normalized in batches**. The input files for a sample batch `HybID00025_00035` is provided as an example.

#### Step 2.1 Create the batch folder
```
mkdir HybID00025_00035
```
#### Step 2.2 Create the metadata file
Populate the metadata file `rnacompete_metadata.csv` in `scripts/metadata`.
- Each row represents an RNAcompete experiment.
- Please consult the sample metadata file for details.

#### Step 2.3 Create the probe intensity file
Prepare the probe intensity file `probe_intensity.csv` and place it under the batch folder.
- Each row represents an RNA probe.
- Each experiment is represented by two columns. For example, the column `HybID00025` contains float values representing median probe intensities, whereas the column `HybID00025_flag` contains either `0` or `1`, where `1` represents probes that are flagged as artifacts. 


### Step 3. Run the analysis pipeline

The pipeline could either be run **locally or remotely** through the LSF scheduler:
```
cd scripts

# LOCAL
python3 analyze_data.py HybID00025_00035

# REMOTE
mkdir logs
bsub < analyze_data.sh
```

After running the pipeline, the following files should be present in the batch folder:
```
├── logo
│   ├── HybID*_*.png
├── pwm
│   ├── HybID*_*.txt
├── scatter
│   ├── HybID*.txt
├── feature.csv
├── kmer_zscore.csv
├── metadata.csv
├── probe_intensity.csv
├── probe_zscore.csv
└── summary.csv
```

### Step 4. Generate the HTML reports
To generate **HTML reports** and an **output summary** of **all experiments**:
```
python3 scripts/generate_html.py
```
To generate an **output summary** of **all experiments** only:
```
python3 scripts/generate_html.py --nohtml
```
To perform the above for **selected experiments**, provide the filename of a file containing the names of selected experiments as an argument:
```
printf 'HybID00025\nHybID00027\n' > selection.txt
python3 generate_html.py selection.txt --nohtml
```


## Differences from [the original pipeline](https://github.com/morrislab/RNAcompete)

|  | Original pipeline | New pipeline |
| ------------- | ------------- | ------------- | 
| **Implementation**  | A mix of MATLAB, Perl, R, Python scripts, and third-party software | Python only |
| **Processing**  | Contains a small bug in k-mer normalization and alignment  | Bugs fixed |
| **Output** | | Includes classifier output

