from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List

# File extension support
_ext = ['.g.vcf', '.gvcf', 'vcf']
EXTENSIONS = tuple(_ext + [f'{x}.gz' for x in _ext])

# Display information for labels
LABEL_DISPLAY = {
    'afr': {'abbr': 'afr', 'color': '#B847A3'},
    'amr': {'abbr': 'amr', 'color': '#FBDF6C'},
    'asj': {'abbr': 'asj', 'color': '#2EDB7E'},
    'eas': {'abbr': 'eas', 'color': '#ED592A'},
    'eur': {'abbr': 'eur', 'color': '#2FA4DC'},
    'sas': {'abbr': 'sas', 'color': '#DC2E31'},
    'fin': {'abbr': 'fin', 'color': '#2FA4DC'},
    'nfe_bgr': {'abbr': 'nfe_bgr', 'color': '#2FA4DC'},
    'nfe_est': {'abbr': 'nfe_est', 'color': '#2FA4DC'},
    'nfe_nwe': {'abbr': 'nfe_nwe', 'color': '#2FA4DC'},
    'nfe_onf': {'abbr': 'nfe_onf', 'color': '#2FA4DC'},
    'nfe_seu': {'abbr': 'nfe_seu', 'color': '#2FA4DC'},
    'nfe_swe': {'abbr': 'nfe_swe', 'color': '#2FA4DC'},
    'eas_jpn': {'abbr': 'eas_jpn', 'color': '#ED592A'},
    'eas_kor': {'abbr': 'eas_kor', 'color': '#ED592A'},
    'eas_oea': {'abbr': 'eas_oea', 'color': '#ED592A'},
    'Puerto Rican in Puerto Rico': {'abbr': 'amr_pri', 'color': '#FBDF6C'},
    'Colombian in Medellin, Colombia': {'abbr': 'amr_col', 'color': '#FBDF6C'},
    'Peruvian in Lima, Peru': {'abbr': 'amr_per', 'color': '#FBDF6C'},
    'Mexican Ancestry in Los Angeles, California': {'abbr': 'amr_mex', 'color': '#FBDF6C'},
    'African Caribbean in Barbados': {'abbr': 'afr_brb', 'color': '#B847A3'},
    'Gambian in Western Division, The Gambia': {'abbr': 'afr_gmb', 'color': '#B847A3'},
    'Esan in Nigeria': {'abbr': 'afr_nga_esan', 'color': '#B847A3'},
    'Mende in Sierra Leone': {'abbr': 'afr_sle', 'color': '#B847A3'},
    'Yoruba in Ibadan, Nigeria': {'abbr': 'afr_nga_yoru', 'color': '#B847A3'},
    'Luhya in Webuye, Kenya': {'abbr': 'afr_ken', 'color': '#B847A3'},
    'African Ancestry in Southwest US': {'abbr': 'afr_sw_usa', 'color': '#B847A3'},
    'Southern Han Chinese, China': {'abbr': 'eas_chn_s_han', 'color': '#ED592A'},
    'Chinese Dai in Xishuangbanna, China': {'abbr': 'eas_chn_dai', 'color': '#ED592A'},
    'Kinh in Ho Chi Minh City, Vietnam': {'abbr': 'eas_vnm', 'color': '#ED592A'},
    'Han Chinese in Bejing, China': {'abbr': 'eas_chn_n_han', 'color': '#ED592A'},
    'Japanese in Tokyo, Japan': {'abbr': 'eas_jpn', 'color': '#ED592A'},
    'British in England and Scotland': {'abbr': 'eur_gbr', 'color': '#2FA4DC'},
    'Finnish in Finland': {'abbr': 'eur_fin', 'color': '#2FA4DC'},
    'Iberian populations in Spain': {'abbr': 'eur_esp', 'color': '#2FA4DC'},
    'Toscani in Italy': {'abbr': 'eur_ita', 'color': '#2FA4DC'},
    'Punjabi in Lahore,Pakistan': {'abbr': 'sas_pak', 'color': '#DC2E31'},
    'Bengali in Bangladesh': {'abbr': 'sas_bgd', 'color': '#DC2E31'},
    'Sri Lankan Tamil in the UK': {'abbr': 'sas_lka', 'color': '#DC2E31'},
    'Indian Telugu in the UK': {'abbr': 'sas_se_ind', 'color': '#DC2E31'},
    'Gujarati Indian in Houston,TX': {'abbr': 'sas_nw_ind', 'color': '#DC2E31'},
    'WestEurasia': {'abbr': 'weur', 'color': '#2FA4DC'},
    'Oceania': {'abbr': 'oce', 'color': '#646464'},
    'America': {'abbr': 'amr', 'color': '#FBDF6C'},
    'Africa': {'abbr': 'afr', 'color': '#B847A3'},
    'EastAsia': {'abbr': 'eas', 'color': '#ED592A'},
    'SouthAsia': {'abbr': 'sas', 'color': '#DC2E31'},
    'CentralAsiaSiberia': {'abbr': 'cas', 'color': '#FFC1C1'},
}

# Titles for plots keyed by model identifier
PLOT_TITLES = {
    'gnomAD_continental': 'gnomAD_continental',
    'gnomAD_eur': 'gnomAD_eur',
    'gnomAD_eas': 'gnomAD_eas',
    '1kGP_amr': '1KGP_amr',
    '1kGP_afr': '1KGP_afr',
    '1kGP_eas': '1KGP_eas',
    '1kGP_eur': '1KGP_eur',
    '1kGP_sas': '1KGP_sas',
    '1kGP_continental': '1KGP_continental',
    'SGDP_continental': 'SGDP_continental',
}


@dataclass
class LabelMappings:
    """Mappings between numeric class indices and label information."""

    continental: Dict[int, Tuple[str, int]]
    subcontinental_eur: Dict[int, Tuple[str, int]]
    subcontinental_eas: Dict[int, Tuple[str, int]]
    genomes_1000_amr: Dict[int, Tuple[str, int]]
    genomes_1000_afr: Dict[int, Tuple[str, int]]
    genomes_1000_eas: Dict[int, Tuple[str, int]]
    genomes_1000_eur: Dict[int, Tuple[str, int]]
    genomes_1000_sas: Dict[int, Tuple[str, int]]
    genomes_1000_continental: Dict[int, Tuple[str, int]]
    sgdp: Dict[int, Tuple[str, int]]


@dataclass
class ResourcePaths:
    """Filesystem paths to resource directories and converters."""

    continental_dir: Path
    subcontinental_dir_eur: Path
    subcontinental_dir_eas: Path
    genomes_1000_continental_dir: Path
    genomes_1000_amr_dir: Path
    genomes_1000_afr_dir: Path
    genomes_1000_eur_dir: Path
    genomes_1000_sas_dir: Path
    genomes_1000_eas_dir: Path
    sgdp_continental_dir: Path
    hg38_json_converter: Path
    wes_b37_json_converter: Path
    wgs_b37_json_converter: Path
    sgdp_b37_json_converter: Path
    json_converts: Dict[str, Dict[str, Dict[str, Path]]]
    matrix_att_dirs: List[Path]
    model_dirs: List[Path]


@dataclass
class ModelConfig:
    """Aggregate configuration for ancestry prediction models."""

    labels: LabelMappings
    paths: ResourcePaths
    n_classes_continental: int
    n_classes_continental_nygc: int
    n_classes_subcontinental_eur: int
    n_classes_subcontinental_eas: int
    n_classes_1000_genomes_afr: int
    n_classes_1000_genomes_amr: int
    n_classes_1000_genomes_sas: int
    n_classes_1000_genomes_eur: int
    n_classes_1000_genomes_eas: int
    n_classes_sgdp_continental: int
    labs_converter: Dict[str, Dict[int, Tuple[str, int]]]
    r_dirs: List[Tuple[str, str, int, str]]
    axis_loc: Dict[str, int]
    titles: Dict[str, str]


def load_variables(resource_root: Path) -> ModelConfig:
    """Build all variables required for the ancestry pipeline.

    Parameters
    ----------
    resource_root : Path
        Root directory containing model resources.

    Returns
    -------
    ModelConfig
        Populated configuration dataclass.
    """

    labels = LabelMappings(
        continental={0: ('afr', 5), 1: ('amr', 1), 2: ('asj', 0), 3: ('eas', 2), 4: ('eur', 3), 5: ('sas', 4)},
        subcontinental_eur={0: ('fin', 0), 1: ('nfe_bgr', 1), 2: ('nfe_est', 2), 3: ('nfe_nwe', 3),
                            4: ('nfe_onf', 4), 5: ('nfe_seu', 5), 6: ('nfe_swe', 6)},
        subcontinental_eas={0: ('eas_jpn', 0), 1: ('eas_kor', 1), 2: ('eas_oea', 2)},
        genomes_1000_amr={0: ('Puerto Rican in Puerto Rico', 0), 1: ('Colombian in Medellin, Colombia', 1),
                          2: ('Peruvian in Lima, Peru', 2), 3: ('Mexican Ancestry in Los Angeles, California', 3)},
        genomes_1000_afr={0: ('African Caribbean in Barbados', 0), 1: ('Gambian in Western Division, The Gambia', 1),
                          2: ('Esan in Nigeria', 2), 3: ('Mende in Sierra Leone', 3),
                          4: ('Yoruba in Ibadan, Nigeria', 4), 5: ('Luhya in Webuye, Kenya', 5),
                          6: ('African Ancestry in Southwest US', 6)},
        genomes_1000_eas={0: ('Southern Han Chinese, China', 0), 1: ('Chinese Dai in Xishuangbanna, China', 1),
                          2: ('Kinh in Ho Chi Minh City, Vietnam', 2), 3: ('Han Chinese in Bejing, China', 3),
                          4: ('Japanese in Tokyo, Japan', 4)},
        genomes_1000_eur={0: ('British in England and Scotland', 0), 1: ('Finnish in Finland', 1),
                          2: ('Iberian populations in Spain', 2), 3: ('Toscani in Italy', 3)},
        genomes_1000_sas={0: ('Punjabi in Lahore,Pakistan', 0), 1: ('Bengali in Bangladesh', 1),
                          2: ('Sri Lankan Tamil in the UK', 2), 3: ('Indian Telugu in the UK', 3),
                          4: ('Gujarati Indian in Houston,TX', 4)},
        genomes_1000_continental={0: ('eur', 3), 1: ('eas', 2), 2: ('amr', 1), 3: ('sas', 4), 4: ('afr', 5)},
        sgdp={0: ('Africa', 0), 1: ('America', 1), 2: ('EastAsia', 2), 3: ('Oceania', 3),
              4: ('WestEurasia', 4), 5: ('CentralAsiaSiberia', 5), 6: ('SouthAsia', 6)},
    )

    paths = ResourcePaths(
        continental_dir=resource_root / 'continental',
        subcontinental_dir_eur=resource_root / 'subcontinental' / 'eur',
        subcontinental_dir_eas=resource_root / 'subcontinental' / 'eas',
        genomes_1000_continental_dir=resource_root / '1k_genomes' / 'wes_v2' / 'continental',
        genomes_1000_amr_dir=resource_root / '1k_genomes' / 'wes_v2' / 'amr',
        genomes_1000_afr_dir=resource_root / '1k_genomes' / 'wes_v2' / 'afr',
        genomes_1000_eur_dir=resource_root / '1k_genomes' / 'wes_v2' / 'eur',
        genomes_1000_sas_dir=resource_root / '1k_genomes' / 'wes_v2' / 'sas',
        genomes_1000_eas_dir=resource_root / '1k_genomes' / 'wes_v2' / 'eas',
        sgdp_continental_dir=resource_root / 'sgdp' / 'continental',
        hg38_json_converter=resource_root / 'genome_ver_converters' / 'hg38tob37' / 'hg38_liftoverAIMs.json',
        wes_b37_json_converter=resource_root / 'genome_ver_converters' / 'b37tohg38' / 'b37.liftover_to_hg38.1kGP.nygc.json',
        wgs_b37_json_converter=resource_root / 'genome_ver_converters' / 'b37tohg38' / 'b37.liftover_to_hg38.1kGP.nygc.json',
        sgdp_b37_json_converter=resource_root / 'genome_ver_converters' / 'b37tohg38' / 'b37_liftover_to_hg38.SGDP.json',
        json_converts={},
        matrix_att_dirs=[],
        model_dirs=[],
    )

    paths.json_converts = {
        '37': {
            'WES': {'1kGP': paths.wes_b37_json_converter, 'gnomAD': None, 'SGDP': paths.sgdp_b37_json_converter},
            'WGS': {'1kGP': paths.wgs_b37_json_converter, 'gnomAD': None, 'SGDP': paths.sgdp_b37_json_converter},
        },
        '38': {
            'WES': {'1kGP': None, 'gnomAD': paths.hg38_json_converter, 'SGDP': None},
            'WGS': {'1kGP': None, 'gnomAD': paths.hg38_json_converter, 'SGDP': None},
        },
    }

    matrix_att = 'matrix_attributes'
    ml_models = 'machine_learning_models'

    sources = [
        paths.continental_dir,
        paths.subcontinental_dir_eur,
        paths.subcontinental_dir_eas,
        paths.genomes_1000_amr_dir,
        paths.genomes_1000_afr_dir,
        paths.genomes_1000_eas_dir,
        paths.genomes_1000_eur_dir,
        paths.genomes_1000_sas_dir,
        paths.genomes_1000_continental_dir,
        paths.sgdp_continental_dir,
    ]
    paths.matrix_att_dirs = [p / matrix_att for p in sources]
    paths.model_dirs = [p / ml_models for p in sources]

    n_classes_continental = 6
    n_classes_continental_nygc = 5
    n_classes_subcontinental_eur = 7
    n_classes_subcontinental_eas = 3
    n_classes_1000_genomes_afr = len(labels.genomes_1000_afr)
    n_classes_1000_genomes_amr = len(labels.genomes_1000_amr)
    n_classes_1000_genomes_sas = len(labels.genomes_1000_sas)
    n_classes_1000_genomes_eur = len(labels.genomes_1000_eur)
    n_classes_1000_genomes_eas = len(labels.genomes_1000_eas)
    n_classes_sgdp_continental = 7

    mode = [
        'gnomAD_continental', 'gnomAD_eur', 'gnomAD_eas', '1kGP_amr', '1kGP_afr',
        '1kGP_eas', '1kGP_eur', '1kGP_sas', '1kGP_continental', 'SGDP_continental'
    ]
    label_order = [
        labels.continental, labels.subcontinental_eur, labels.subcontinental_eas,
        labels.genomes_1000_amr, labels.genomes_1000_afr, labels.genomes_1000_eas,
        labels.genomes_1000_eur, labels.genomes_1000_sas,
        labels.genomes_1000_continental, labels.sgdp
    ]
    labs_converter = dict(zip(mode, label_order))
    n_classes = [
        n_classes_continental, n_classes_subcontinental_eur, n_classes_subcontinental_eas,
        n_classes_1000_genomes_amr, n_classes_1000_genomes_afr, n_classes_1000_genomes_eas,
        n_classes_1000_genomes_eur, n_classes_1000_genomes_sas,
        n_classes_continental_nygc, n_classes_sgdp_continental
    ]
    r_dirs = list(zip(
        [str(p) for p in paths.matrix_att_dirs],
        [str(p) for p in paths.model_dirs],
        n_classes,
        mode
    ))

    axis_order = [0, 2, 4, 6, 7, 5, 3, 8, 1, 9]
    axis_loc = dict(zip(mode, axis_order))
    titles = {m: PLOT_TITLES[m] for m in mode}

    return ModelConfig(
        labels=labels,
        paths=paths,
        n_classes_continental=n_classes_continental,
        n_classes_continental_nygc=n_classes_continental_nygc,
        n_classes_subcontinental_eur=n_classes_subcontinental_eur,
        n_classes_subcontinental_eas=n_classes_subcontinental_eas,
        n_classes_1000_genomes_afr=n_classes_1000_genomes_afr,
        n_classes_1000_genomes_amr=n_classes_1000_genomes_amr,
        n_classes_1000_genomes_sas=n_classes_1000_genomes_sas,
        n_classes_1000_genomes_eur=n_classes_1000_genomes_eur,
        n_classes_1000_genomes_eas=n_classes_1000_genomes_eas,
        n_classes_sgdp_continental=n_classes_sgdp_continental,
        labs_converter=labs_converter,
        r_dirs=r_dirs,
        axis_loc=axis_loc,
        titles=titles,
    )

