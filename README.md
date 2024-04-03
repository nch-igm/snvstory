# SNVstory
Rapid and accurate ancestry inference using single nucleotide variants.

## Citation

If you use SNVstory in your research, please consider citing our paper:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05703-y


## Dependencies

SNVstory requires Docker to run. Download and install [Docker Desktop](https://docs.docker.com/get-docker/).

You will also need the AWS Command Line Interface to download the resources needed to run the models. Please follow the installation instructions for [AWS CLI Version 2](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

## Installation

Build SNVstory.
```bash
cd snvstory
docker-compose build
```

Copy the resources into a location the container can find. This step will take some time, but only needs to be copied once.
```bash
aws s3 sync s3://igm-public-dropbox/snvstory/resource_dir/ dev/data/resource_dir/ --no-sign-request
```

## Execution

SNVstory is executed with Docker by running the following on the terminal. Make sure Docker Destop is running or this command will not work.
```bash
docker-compose run ancestry <arguments>
```


Run to see all possible arguments:
```bash
docker-compose run ancestry --help
```

Example run with VCF on s3.
```bash
docker-compose run ancestry \
    --path s3://path-to-input-file \
    --resource "/data/resource_dir" \
    --output-dir s3://path-to-output-directory \
    --genome-ver 38 \
    --mode WES
```

Example run with VCF on local computer. Make sure both input and output files are location in the ```/data/``` directory so that the Docker container can find them.
```bash
docker-compose run ancestry \
    --path "/data/path-to-input-file" \
    --resource "/data/resource_dir" \
    --output-dir "/data/path-to-output-directory" \
    --genome-ver 38 \
    --mode WES
```

## Output

### Ancestry Report
SNVstory returns a .csv report which includes the probabilities of each label. A .pdf is also returned, which summarizes these model probabilities in dot plots. The subcontinental model probabilities are weighted by the corresponding continental model result. 

E.g., in the following example case the gnomAD continental probability for the 'eas' label is 0, so the gnomAD East Asian subcontinental model probabilities are multiplied by 0 in the dot plot.


![Example Report](assets/ExampleAncestryReport.svg)


### UMAP
SNVstory also outputs a UMAP transformation of the user input sample (in black) on each set of training samples (color labeled by continent). The interactive plots are saved to .html files (see ./assets). A hover tool is used to display the country and population of nearby training samples.

![Example Report](assets/Example_1kGP_umap.png)



## Feature Importance

The feature importance anaysis is executed separately, and requires micromamba/conda for installation.

### Installation

```bash
cd Feature_Importance
micromamba create -f config/snvstoryfeats.yml
micromamba activate snvstoryfeats
```

### Execution

Feature importance requires your single/multi-sample VCF as input and will run on all samples in the VCF.

```bash
(snvstoryfeats)$ python Feature_importance.py --help
usage: Feature_importance.py [-h] [-b] [-s] [-i] vcf output resource_dir {hg19,hg38}

Calculate feature importance aggregated to genes and cytolocations. Returns two .npz files with shap values. Optionally create summary plots.

positional arguments:
  vcf                   Path to input VCF file.
  output                Output folder to write results. Will exit if folder exists.
  resource_dir          Path to resource directory.
  {hg19,hg38}           Genome version (hg19 or hg38) of input vcf.

optional arguments:
  -h, --help            show this help message and exit
  -b, --bar-plot        Flag to indicate whether to create mean(|SHAP val|) bar plot. All samples in the VCF are aggregated together.
  -s, --stacked-bar-plot
                        Flag to indicate whether to create stacked bar plot. All samples in the VCF are aggregated together.
  -i, --ideogram-plot   Flag to indicate whether to create mean(|SHAP val|) bar plot. This will create a separate plot for each sample in the VCF.

```

Run the following example with the single sample VCF in the resources directory.

```bash
python Feature_importance.py ../dev/data/resource_dir/feature_importance/HG00096.example.vcf.gz test_output/ ../dev/data/resource_dir/ hg38 -b -s -i
```

The following samples should be created in the output directory:
```
test_output/
├── HG00096_ideogram.png
├── Top_20_genes.png
├── Top_20_genes_stacked_label.png
├── shap_cyto.npz
└── shap_gene.npz
```

See ```assets/``` for the .png results.






