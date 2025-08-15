import os
import uuid
import logging
import tempfile
import argparse
from typing import List, Tuple, Optional

from igm_churchill_ancestry.pipelines.variables import variables
from igm_churchill_ancestry.utilities.flex import flex_input, flex_output, get_dir_size
from igm_churchill_ancestry.pipelines.ancestry_prediction import run_ancestry_pipeline
from igm_churchill_ancestry.utilities.utilities import (
    get_extension,
    filter_extension,
    is_vcf_multisample,
    check_resources,
)





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


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments for the ancestry pipeline."""

    parser = argparse.ArgumentParser(description="Ancestry Prediction v1.0")
    parser.add_argument(
        "--path",
        dest="path",
        required=True,
        type=str,
        help="<REQUIRED> specify the vcf/multi-sample-vcf s3 or local file path or the path to directory containing the file(s)",
    )
    parser.add_argument(
        "--resource",
        dest="resource",
        required=True,
        type=str,
        help="<REQUIRED> specify the location of the resource folder",
    )
    parser.add_argument(
        "--sample_pos",
        dest="sp",
        required=False,
        nargs="+",
        default=["all"],
        help=(
            "if the input is in a multi-sample vcf format specify which sample to select "
            "designated by position in the MultiSample vcf e.g. 1,2,3 ect or all.\nAlternatively you can specify the name of the sample e.g. mother, father, proband"
        ),
    )
    parser.add_argument(
        "--logging",
        dest="logging",
        type=str,
        help="<OPTIONAL> provide the output path for the logging file",
    )
    parser.add_argument(
        "-j",
        "--vcf_json",
        required=False,
        nargs="+",
        default=None,
        help="<REQUIRED> provide a json format of the data e.g. 2_238272966_TC:1",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        type=str,
        required=True,
        help="<REQUIRED> provide the dir path",
    )
    parser.add_argument(
        "--genome-ver",
        dest="genome_ver",
        type=str,
        required=True,
        choices=["37", "38"],
        default="38",
        help="<REQUIRED> select a human genome version",
    )
    parser.add_argument(
        "--mode",
        dest="mode",
        type=str,
        required=True,
        nargs="+",
        default="WES",
        help=(
            "<REQUIRED> Mode that sequence allocation analyses were run in. "
            "Provide a value for each VCF if multiple VCFs are being submitted."
        ),
    )
    parser.add_argument(
        "--output_filename",
        dest="output_filename",
        type=str,
        default=None,
        help="<OPTIONAL> File name used for the prediction csv",
    )
    return parser.parse_args()


def resolve_inputs(path: str, sample_positions: List[str], input_dir: str) -> Tuple[List[str], List[str]]:
    """Resolve input path into a list of VCF files and sample positions."""

    file_ext = get_extension(path)
    if file_ext != "dir":
        try:
            local_vcf = flex_input(path)
            logging.info("Input VCF success. File size: %s", os.path.getsize(local_vcf))
            vcf_paths = [local_vcf]
        except ValueError as err:
            logging.error("Input file not found or invalid: %s", path)
            raise FileNotFoundError(f"Input file not found or invalid: {path}") from err
    else:
        try:
            local_dir = flex_input(path, input_dir, directory=True)
            logging.info(
                "Input VCF directory success. Directory size: %s", get_dir_size(local_dir)
            )
            vcf_paths = [os.path.join(local_dir, x) for x in os.listdir(local_dir)]
            vcf_paths = filter_extension(vcf_paths)
        except ValueError as err:
            logging.error("Directory not found or invalid: %s", path)
            raise FileNotFoundError(f"Directory not found or invalid: {path}") from err

    if len(sample_positions) != len(vcf_paths):
        logging.warning(
            "Length of 'sample_pos' is not equal to the number of relevant input VCFs. Analyzing all samples."
        )
        sample_positions = ["all"] * len(vcf_paths)

    return vcf_paths, sample_positions


def download_resources(resource_path: str, resource_dir: str) -> variables:
    """Download and validate required resources for the pipeline."""

    try:
        local_rsrc_dir = flex_input(resource_path, resource_dir, directory=True)
    except ValueError as err:
        logging.error("Resource folder not found or invalid: %s", resource_path)
        raise FileNotFoundError(
            f"Resource folder not found or invalid: {resource_path}"
        ) from err
    var = variables(local_rsrc_dir)
    check_resources(var)
    return var


def execute_pipeline(
    vcf_paths: List[str],
    sample_positions: List[str],
    var: variables,
    out_dir: str,
    genome_ver: str,
    modes: List[str],
    output_filename: Optional[str],
) -> None:
    """Run the ancestry prediction pipeline for each VCF file."""

    if len(modes) == 1 and len(vcf_paths) > 1:
        modes = modes * len(vcf_paths)
    if len(modes) != len(vcf_paths):
        raise ValueError(
            f"Number of --mode values ({len(modes)}) must match number of input VCFs ({len(vcf_paths)})"
        )

    for i, vcf in enumerate(vcf_paths):
        logging.info("Processing file: %s", vcf)
        multi_sample_status, sample_name = is_vcf_multisample(vcf, True)
        logging.info(
            "Multisample: %s and sample name(s): %s", multi_sample_status, sample_name
        )
        if output_filename is None:
            ofn = f"{os.path.splitext(vcf)[0].split('/')[-1]}_{sample_name[0]}.csv"
        else:
            ofn = output_filename

        run_ancestry_pipeline(
            vcf_path=vcf,
            multi_sample_status=multi_sample_status,
            sample=sample_name,
            sample_position=sample_positions[i],
            var=var,
            outdir=out_dir,
            genome_ver=genome_ver,
            mode=modes[i],
            ofn=ofn,
        )


def run_ancestry() -> None:
    """Parse CLI args, manage resources, and run the ancestry pipeline."""

    args = parse_arguments()

    data_dir = setup_workspace()
    input_dir = f"{data_dir}/input/"
    logging_file = f"{data_dir}/logging.txt"
    logging.basicConfig(
        filename=logging_file,
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s:%(message)s",
    )
    rsrc_dir = f"{data_dir}/resources/"
    os.makedirs(rsrc_dir)
    out_dir = f"/{data_dir}/output/"
    os.makedirs(out_dir)

    logging.info("Executing Ancestry Pipeline")

    vcf_paths, sample_positions = resolve_inputs(args.path, args.sp, input_dir)
    var = download_resources(args.resource, rsrc_dir)
    execute_pipeline(
        vcf_paths,
        sample_positions,
        var,
        out_dir,
        args.genome_ver,
        args.mode,
        args.output_filename,
    )
    flex_output(out_dir, args.output_dir)

    logging.info("Completed")
    try:
        if args.logging:
            flex_output(logging_file, args.logging)
    except (ValueError, OSError) as err:
        logging.error("Failed to copy log file to %s: %s", args.logging, err)



