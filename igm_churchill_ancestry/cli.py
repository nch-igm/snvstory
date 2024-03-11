import os
import uuid
import logging
import tempfile
import argparse
import concurrent.futures

from igm_churchill_ancestry.pipelines.variables import variables
from igm_churchill_ancestry.utilities.flex import flex_input, flex_output
from igm_churchill_ancestry.pipelines.ancestry_prediction import run_ancestry_pipeline
from igm_churchill_ancestry.utilities.utilities import get_extension, filter_extension, is_vcf_multisample, check_resources





def setup_workspace():
    """Create the workspace environment for this job."""
    tmp_dir = os.environ.get("TMP_DIR")
    if not tmp_dir:
        tmp_dir = tempfile.gettempdir()
    job_id = os.environ.get("AWS_BATCH_JOB_ID")
    if not job_id:
        job_id = uuid.uuid4().hex
    wrk_dir = f"{tmp_dir}/{job_id}"
    try:
        os.makedirs(wrk_dir)
    except OSError:
        raise RuntimeError(f"Could not create working directory: {wrk_dir}")
    return wrk_dir


def run_ancestry():
    """Parse cli args, download from s3, run the normal pipeline, upload to s3."""
    parser = argparse.ArgumentParser(description='Ancestry Prediction v1.0')
    parser.add_argument('--path', dest="path", required=True, type=str, help="<REQUIRED> specify the vcf/multi-sample-vcf s3 or local file path or the path to directory containing the file(s)")
    parser.add_argument('--resource', dest="resource", required=True, type=str, help="<REQUIRED> specify the location of the resource folder")
    parser.add_argument('--sample_pos', dest="sp", required=False, nargs='+', default='all', help="if the input is in a multi-sample vcf format specify which sample to select designated by position in the MultiSample vcf e.g. 1,2,3 ect or all.\nAlternatively you can specify the name of the sample e.g. mother, father, proband")
    parser.add_argument('--logging', dest='logging', type=str, help="<OPTIONAL> provide the output path for the logging file")
    parser.add_argument('-j', '--vcf_json', required=False, nargs='+', default=None, help="<REQUIRED> provide a json format of the data e.g. 2_238272966_TC:1")
    parser.add_argument('--output-dir', dest='output_dir', type=str, required=True, help="<REQUIRED> provide the dir path")
    parser.add_argument('--genome-ver', dest='genome_ver', type=str, required=True, choices=['37', '38'], default='38', help="<REQUIRED> select a human genome version")
    parser.add_argument('--mode', dest='mode', type=str, required=True, nargs='+', default='WES', help="<REQUIRED> Mode that sequence allocation analyses were run in. Provide a value for each VCF if multiple VCFs are being submitted.")
    parser.add_argument('--output_filename', dest='output_filename', type=str, default=None, help="<REQUIRED> File name used for the prediction csv")
    args = parser.parse_args()
    """
    Argument Paths:
    Example command:
    in src dir
    pip install .
    docker build -t latest .


    docker-compose -f docker-compose.yml run ancestry \
    --path /dev/data/test2.vcf \
    --resource /dev/data/resource_dir  \
    --output-dir /dev/data/output/ \
    --genome-ver 37 \
    --mode WES

    path -> is dir or file?
         dir -> is s3 or local
             s3 -> download
             local -> reassign variable
         file -> is s3 or local
              s3 -> download
              local -> reassign variable

    resources -> is s3 or local?
              s3 -> download
              local -> reassign variable

    sample_comp -> is dir or file?
              dir -> LOOP -> file
                  file -> is single or multi samples
                       single -> execute ancestry_single_sample
                       multi -> execute ancestry_multi_sample
              file -> is single or multi samples
                  single -> execute ancestry_single_sample
                  multi -> execute ancestry_multi_sample


                        run
                       //  \\
                      //    \\
                 single      multi
                //   \\     //   \\
               gvcf   vcf  gvcf   vcf
              //|\\       //|\\
             0 - all     0 - all

    """

    # setup directories
    DATA_DIR = setup_workspace()
    INPUT_DIR = f"{DATA_DIR}/input/"
    LOGGING_FILE = f"{DATA_DIR}/logging.txt"
    logging.basicConfig(filename=LOGGING_FILE, level=logging.DEBUG, format="%(asctime)s:%(levelname)s:%(message)s")
    RSRC_DIR = f"{DATA_DIR}/resources/"
    os.makedirs(RSRC_DIR)
    OUT_DIR = f"/{DATA_DIR}/output/"
    os.makedirs(OUT_DIR)

    logging.info("Executing Ancestry Pipeline")

    # determine extension
    file_extz = get_extension(args.path)

    # placeholders
    local_vcf_file = ''
    local_vcf_dir = ''

    if file_extz != 'dir':
        # Attempt to download single input from s3
        try:
            local_vcf_file = flex_input(args.path)
            logging.debug(f"Input VCF success. File size: {os.path.getsize(local_vcf_file)}")
            print(f"Input VCF success. File size: {os.path.getsize(local_vcf_file)}")
        except Exception:
            logging.debug(f"Input file: {args.path} does not appear to be local or s3 input.")
            print(f"Input file: {args.path} does not appear to be local or s3 input.")
    else:
        # Attempt to download directory from s3
        try:
            local_vcf_dir = flex_input(args.path, INPUT_DIR, directory=True)
            logging.debug(f"Input VCF directory success. Directory size: {os.path.getsize(INPUT_DIR)}") # this doesn't give the size of the files in directory rn
            print(f"Input VCF directory success. Directory size: {os.path.getsize(INPUT_DIR)}")
        except Exception:
            logging.debug(f"Directory: {args.path} does not appear to be local or an s3 input")
            print(f"Directory: {args.path} does not appear to be local or an s3 input")
            raise RuntimeError

    # Download resource folder e.g. models
    try:
        RSRC_DIR = flex_input(args.resource, RSRC_DIR, directory=True)
    except Exception:
        logging.debug(f"Resource folder: {args.resource} does not appear to be local or an s3 input")
        print(f"Resource folder: {args.resource} does not appear to be local or an s3 input")
        raise RuntimeError

    var = variables(RSRC_DIR)
    check_resources(var)

    # Determine if the file is multi-sampled or single
    if local_vcf_dir:
        print(f'Expected input: {args.path} is a directory')
        local_vcf_dir = [os.path.join(INPUT_DIR, x) for x in os.listdir(local_vcf_dir)]
        local_vcf_dir = filter_extension(local_vcf_dir)
        if len(args.sp) != len(local_vcf_dir):
            if args.sp != 'all':
                print(f"Length of 'sample_pos' is not equal to the number of relevant input VCFs. Analyzing all samples.")
            sp = ['all'] * len(local_vcf_dir)
        else:
            sp = args.sp
        for i, f in enumerate(local_vcf_dir):
            multi_sample_status, sample_name = is_vcf_multisample(f, True)
            if args.output_filename is None:
                ofn = f"{os.path.splitext(f)[0].split('/')[-1]}_{sample_name[0]}.csv"
            else:
                ofn = args.output_filename
            print(f'Multisample: {multi_sample_status} and sample name: {sample_name}')

            run_ancestry_pipeline(vcf_path=f, multi_sample_status=multi_sample_status,
                                  sample=sample_name, sample_position=sp, var=var,
                                  outdir=OUT_DIR, genome_ver=args.genome_ver, mode=args.mode[i],
                                  ofn=ofn)
        flex_output(OUT_DIR, args.output_dir)

    elif local_vcf_file:
        print(f'Expected input: {local_vcf_file} is a file')
        multi_sample_status, sample_name = is_vcf_multisample(local_vcf_file, True)
        print(f'Multisample: {multi_sample_status} and sample name: {sample_name[0]}')
        if args.output_filename is None:
            ofn = f"{os.path.splitext(local_vcf_file)[0].split('/')[-1]}_{sample_name[0]}.csv"
        else:
            ofn = args.output_filename

        run_ancestry_pipeline(vcf_path=local_vcf_file, multi_sample_status=multi_sample_status,
                              sample=sample_name, sample_position=args.sp, var=var,
                              outdir=OUT_DIR, genome_ver=args.genome_ver, mode=args.mode[0],
                              ofn=ofn)
        flex_output(OUT_DIR, args.output_dir)

    logging.debug(f"Completed")
    try:
        if args.logging:
            flex_output(LOGGING_FILE, args.logging)
    except Exception:
        pass


