# SNVstory
Rapid and accurate ancestry inference using SNVs.


## Installation

SNVstory requires Docker to run. Download and install [Docker Desktop](https://docs.docker.com/get-docker/).

Build SNVstory.
```bash
cd snvstory
docker-compose build
```

Copy the resources into a location the container can find. This step will take some time, but only needs to be copied once.
```bash
aws s3 sync s3://path-to-directory dev/data/resource_dir/
```

## Execution

```bash
docker-compose run ancestry <arguments>
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

Example run with VCF on local computer.
```bash
docker-compose run ancestry \
    --path "/data/path-to-input-file" \
    --resource "/data/resource_dir" \
    --output-dir "/data/path-to-output-directory" \
    --genome-ver 38 \
    --mode WES
```





