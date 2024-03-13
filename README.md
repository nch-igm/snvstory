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

## Using Docker Image

To use the `heritagehelix/snvstory:3.0.2` Docker image directly, follow these steps to run SNVstory. Ensure Docker Desktop is running before executing these commands.

### Prerequisites
- Docker must be installed and running on your system.
- Your data and resources should be organized in known directories for volume mounting.

### Command

```bash
docker run \
  -v /path/to/your/data:/data \
  -v /path/to/your/output:/output \
  -t "heritagehelix/snvstory:3.0.2" \
  --path /data/your_input_file.vcf.gz \
  --resource /data/resource_dir \
  --output-dir /output \
  --genome-ver 38 \
  --mode WES
  ```

### Parameters Explained
- ```-v /path/to/your/data:/data```: This mounts the directory containing your input VCF file and resource_dir into the container. Replace /path/to/your/data with the local directory containing your VCF files and resource_dir.
- ```-v /path/to/your/output:/output```: Specifies the directory where output files will be saved. Replace /path/to/your/output with the local directory where you want to save the output.
- ```-t "heritagehelix/snvstory:3.0.2"```: Specifies the Docker image to use.
- ```--path /data/your_input_file.vcf.gz```: The path to your input VCF file within the Docker container. Replace your_input_file.vcf.gz with the name of your VCF file.
- ```--resource /data/resource_dir```: Points to the directory containing necessary resources for running SNVstory, as seen by the Docker container.
- ```--output-dir /output```: Designates the output directory within the Docker container.
- ```--genome-ver 38```: Specifies the genome version.
- ```--mode WES```: Indicates the mode of operation, which can be Whole Exome Sequencing (WES) in this case.

## Using Docker Compose 

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

## Using Singularity


This guide walks you through running SNVstory using Singularity. Singularity is a container platform that allows you to run applications in isolated environments.

Prerequisites
- Singularity must be installed on your system. You can find installation instructions for Singularity on the official Singularity website.
- Ensure you have your input VCF file and the required resource_dir ready in known directories on your local machine.

### Step-by-Step Command
Execute the following command in your terminal, making sure to replace the placeholders with your actual data and output paths:

```bash
singularity exec -B $PWD:$PWD,/path/to/resource_dir/:$PWD/resource_dir \
   --env PYTHONPATH=/opt/Ancestry/:$PYTHONPATH --env TMP_DIR=$PWD/tmp \
    /path/to/snvstory.sif python3 -m igm_churchill_ancestry  \
   --path your_input_file.vcf.gz --resource $PWD/resource_dir --genome-ver 38  \
   --mode WES --output-dir /output_dir
```


### Parameters Explained
- ```-B $PWD:$PWD,/path/to/resource_dir/:$PWD/resource_dir```: This binds directories from the host system to the container. Replace ```/path/to/resource_dir/``` with the path to your resource_dir on the host system.
- ```--env PYTHONPATH=/opt/Ancestry/:$PYTHONPATH```: Sets the ```PYTHONPATH``` environment variable inside the container to include the ```/opt/Ancestry/``` directory.
- ```--env TMP_DIR=$PWD/tmp```: Sets the ```TMP_DIR``` environment variable to the tmp directory in the current working directory.
- ```/path/to/snvstory.sif```: The path to the Singularity image file for SNVstory. Replace this with the actual path to your ```snvstory.sif``` file.
- ```python3 -m igm_churchill_ancestry```: Runs the igm_churchill_ancestry module using Python 3 inside the container.
- ```--path your_input_file.vcf.gz```: The path to your input VCF file. Replace ```your_input_file.vcf.gz``` with the name of your VCF file.
- ```--resource $PWD/resource_dir```: Points to the directory containing necessary resources for running SNVstory, as seen by the container.
- ```--genome-ver 38```: Specifies the genome version.
- ```--mode WES```: Indicates the mode of operation, which can be Whole Exome Sequencing (WES) in this case.
- ```--output-dir /output_dir```: Designates the output directory where the results will be saved. Replace ```/output_dir``` with the desired output directory path.

Make sure to replace the placeholders (```/path/to/resource_dir/```, ```/path/to/snvstory.sif```, ```your_input_file.vcf.gz```, ```/output_dir```) with the actual paths and file names relevant to your setup.

## Output

### Ancestry Report
SNVstory returns a .csv report which includes the probabilities of each label. A .pdf is also returned, which summarizes these model probabilities in dot plots. The subcontinental model probabilities are weighted by the corresponding continental model result. 

E.g., in the following example case the gnomAD continental probability for the 'eas' label is 0, so the gnomAD East Asian subcontinental model probabilities are multiplied by 0 in the dot plot.


![Example Report](assets/ExampleAncestryReport.svg)


### UMAP
SNVstory also outputs a UMAP transformation of the user input sample (in black) on each set of training samples (color labeled by continent). The interactive plots are saved to .html files (see ./assets). A hover tool is used to display the country and population of nearby training samples.

![Example Report](assets/Example_1kGP_umap.png)








