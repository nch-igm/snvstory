import numpy as np
from igm_churchill_ancestry.pipelines.variables import variables


def parse_vcf(o, gz_file):
    f"""
    A function to read a VCF and store lines with actionable
    information. Acceptable formats are: {variables.EXTENSIONS}. Assumes
    only single sample.

    args
    ----
    o - vcf file read into memory and broken into lines
    gz_file - boolean represents whether file is gzipped

    returns
    -------
    parsed_vcf - Returns non-comment lines as a list.

    """
    parsed_vcf = []
    for num, line in enumerate(o):
        if gz_file:
            line = line.decode('utf-8').strip()
        if not line.startswith('#'):
            line = line.strip()
            parsed_vcf.append(line)
    return parsed_vcf


def parse_multisample_vcf_sample(o, gz_file, sample):
    """
    Extract the sample position from a vcf or gvcf

    args
    ----
    o - vcf file read into memory and broken into lines
    gz_file - boolean represents whether file is gzipped
    sample - sample position or name (str)

    returns
    -------
    ext - returns index position(s) of samples to be analyzed

    """
    numeric = False
    try:
        pos = int(sample)
        numeric = True
    except Exception:
        pass
    for num, line in enumerate(o):
        if gz_file:
            line = line.decode('utf-8').strip()
        if line.startswith('#CHROM'):
            line = line.strip()
            break
    aline = np.asarray(line.split('\t'))
    if numeric is True:
        # we assume samples always come after the FORMAT field
        index = np.where(aline == 'FORMAT')[0][0]
        sample_names = aline[index + pos]
        return index + pos, sample_names
    elif sample == 'all':
        index = np.where(aline == 'FORMAT')[0][0]
        indices = np.arange(index + 1, len(aline)).tolist()
        sample_names = aline[indices]
        return indices, sample_names
    else:
        try:
            index = np.where(aline == sample)[0][0]
            if index.size != 0:
                sample_names = aline[index]
                return index, sample_names
            else:
                print(f'Unable to find index using sample: {sample}')
        except Exception as e:
            print(f'Unable to find index using sample: {sample} error: {e}')
            return


def parse_multisample_vcf(o, gz_file, sample):
    f"""
    A function to read a gVCF and store lines with actionable information.
    Acceptable formats are {variables.EXTENSIONS}


    args
    ----
    o - vcf file read into memory and broken into lines
    gz_file - boolean represents whether file is gzipped
    sample - sample position or name (str)

    returns
    -------
    parsed_vcf - Returns non-comment lines as a list.

    """
    parsed_vcf = []
    for num, line in enumerate(o):
        if gz_file:
            line = line.decode('utf-8').strip()
        if not line.startswith('#'):
            line = line.strip()
            # assume the sample has been converted to the numeric column
            try:
                sline = line.split('\t')
                sample_gt = sline[sample]
                nline = '\t'.join(sline[:9] + [sample_gt])  # I assume that the samples start after the FORMAT field which is the 8th field by zero index.
                parsed_vcf.append(nline)
            except Exception as e:
                print(f'Failed to parse sample {sample} from line {line}. {e}')
                return
    return parsed_vcf
