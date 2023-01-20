def join_paths(path1, path2):
    path = path1 + path2
    return [path]


class variables:
    # static variables
    ext = ['.g.vcf', '.gvcf', 'vcf']
    ext_gz = [f'{x}.gz' for x in ext]
    EXTENSIONS = tuple(ext + ext_gz)

    def __init__(self, rsrc_root):

        # model labels
        self.LABS_CONT = {0: ('afr', 5), 1: ('amr', 1), 2: ('asj', 0), 3: ('eas', 2), 4: ('eur', 3), 5: ('sas', 4)}
        self.LABS_SUBCONT_EUR = {0: ('fin',0), 1: ('nfe_bgr',1), 2: ('nfe_est',2), 3: ('nfe_nwe',3), 4: ('nfe_onf',4), 5: ('nfe_seu',5), 6: ('nfe_swe',6)}
        self.LABS_SUBCONT_EAS = {0: ('eas_jpn',0), 1: ('eas_kor',1), 2: ('eas_oea', 2)}
        self.LABS_GENOMES_1000_AMR = {0: ('Puerto Rican in Puerto Rico', 0), 1: ('Colombian in Medellin, Colombia', 1), 2: ('Peruvian in Lima, Peru',2), 3: ('Mexican Ancestry in Los Angeles, California', 3)}
        self.LABS_GENOMES_1000_AFR = {0: ('African Caribbean in Barbados', 0), 1: ('Gambian in Western Division, The Gambia', 1), 2: ('Esan in Nigeria', 2), 3: ('Mende in Sierra Leone', 3), 4: ('Yoruba in Ibadan, Nigeria', 4), 5: ('Luhya in Webuye, Kenya', 5), 6: ('African Ancestry in Southwest US', 6)}
        # self.LABS_GENOMES_1000_EAS = {0: ('Southern Han Chinese, China', 4), 1: ('Chinese Dai in Xishuangbanna, China', 1), 2: ('Kinh in Ho Chi Minh City, Vietnam', 2), 3: ('Han Chinese in Bejing, China', 3), 4: ('Japanese in Tokyo, Japan', 0)}
        self.LABS_GENOMES_1000_EAS = {0: ('Southern Han Chinese, China', 0), 1: ('Chinese Dai in Xishuangbanna, China', 1), 2: ('Kinh in Ho Chi Minh City, Vietnam', 2), 3: ('Han Chinese in Bejing, China', 3), 4: ('Japanese in Tokyo, Japan', 4)}
        # self.LABS_GENOMES_1000_EUR = {0: ('British in England and Scotland', 1), 1: ('Finnish in Finland', 0), 2: ('Iberian populations in Spain', 2), 3: ('Toscani in Italy', 3)}
        self.LABS_GENOMES_1000_EUR = {0: ('British in England and Scotland', 0), 1: ('Finnish in Finland', 1), 2: ('Iberian populations in Spain', 2), 3: ('Toscani in Italy', 3)}
        self.LABS_GENOMES_1000_SAS = {0: ('Punjabi in Lahore,Pakistan', 0), 1: ('Bengali in Bangladesh', 1), 2: ('Sri Lankan Tamil in the UK', 2), 3: ('Indian Telugu in the UK', 3), 4: ('Gujarati Indian in Houston,TX', 4)}
        self.LABS_GENOMES_1000_CONTINENTAL = {0: ('eur', 3), 1: ('eas', 2), 2: ('amr',1), 3: ('sas', 4), 4: ('afr', 5)}
        self.LABS_SGDP = {0: ('Africa', 0), 1: ('America', 1), 2: ('EastAsia', 2), 3: ('Oceania', 3), 4: ('WestEurasia', 4), 5: ('CentralAsiaSiberia', 5), 6: ('SouthAsia', 6)}
        

        self.ABBR = {'afr' : ('afr', '#B847A3'), 'amr' : ('amr', '#FBDF6C'), 'asj' : ('asj', '#2EDB7E'), 'eas' : ('eas', '#ED592A'), 
            'eur' : ('eur', '#2FA4DC'), 'sas' : ('sas', '#DC2E31'), 'fin' : ('fin', '#2FA4DC'), 'nfe_bgr' : ('nfe_bgr', '#2FA4DC'), 
            'nfe_est' : ('nfe_est', '#2FA4DC'), 'nfe_nwe' : ('nfe_nwe', '#2FA4DC'), 'nfe_onf' : ('nfe_onf', '#2FA4DC'), 
            'nfe_seu' : ('nfe_seu', '#2FA4DC'), 'nfe_swe' : ('nfe_swe', '#2FA4DC'), 'eas_jpn' : ('eas_jpn', '#ED592A'), 
            'eas_kor' : ('eas_kor', '#ED592A'), 'eas_oea' : ('eas_oea', '#ED592A'), 'Puerto Rican in Puerto Rico' : ('amr_pri', '#FBDF6C'), 
            'Colombian in Medellin, Colombia' : ('amr_col', '#FBDF6C'), 'Peruvian in Lima, Peru' : ('amr_per', '#FBDF6C'), 
            'Mexican Ancestry in Los Angeles, California' : ('amr_mex', '#FBDF6C'), 'African Caribbean in Barbados' : ('afr_brb', '#B847A3'), 
            'Gambian in Western Division, The Gambia' : ('afr_gmb' , '#B847A3'), 'Esan in Nigeria' : ('afr_nga_esan', '#B847A3'), 
            'Mende in Sierra Leone' : ('afr_sle', '#B847A3'), 'Yoruba in Ibadan, Nigeria' : ('afr_nga_yoru', '#B847A3'), 
            'Luhya in Webuye, Kenya' : ('afr_ken', '#B847A3'), 'African Ancestry in Southwest US' : ('afr_sw_usa', '#B847A3'), 
            'Southern Han Chinese, China' : ('eas_chn_s_han', '#ED592A'), 'Chinese Dai in Xishuangbanna, China' : ('eas_chn_dai', '#ED592A'), 
            'Kinh in Ho Chi Minh City, Vietnam' : ('eas_vnm', '#ED592A'), 'Han Chinese in Bejing, China' : ('eas_chn_n_han', '#ED592A'), 
            'Japanese in Tokyo, Japan' : ('eas_jpn', '#ED592A'), 'British in England and Scotland' : ('eur_gbr', '#2FA4DC'), 
            'Finnish in Finland' : ('eur_fin', '#2FA4DC'), 'Iberian populations in Spain' : ('eur_esp', '#2FA4DC'), 
            'Toscani in Italy' : ('eur_ita', '#2FA4DC'), 'Punjabi in Lahore,Pakistan' : ('sas_pak', '#DC2E31'), 
            'Bengali in Bangladesh' : ('sas_bgd', '#DC2E31'), 'Sri Lankan Tamil in the UK' : ('sas_lka', '#DC2E31'), 
            'Indian Telugu in the UK' : ('sas_se_ind', '#DC2E31'), 'Gujarati Indian in Houston,TX' : ('sas_nw_ind', '#DC2E31'), 
            'eur' : ('eur', '#2FA4DC'), 'eas' : ('eas', '#ED592A'), 'amr' : ('amr', '#FBDF6C'), 'sas' : ('sas', '#DC2E31'), 
            'afr' : ('afr', '#B847A3'),
            'WestEurasia' : ('weur', '#2FA4DC'), 'Oceania' : ('oce', '#646464'), 'America' : ('amr', '#FBDF6C'), 'Africa' : ('afr', '#B847A3'),
            'EastAsia' : ('eas', '#ED592A'), 'SouthAsia' : ('sas', '#DC2E31'), 'CentralAsiaSiberia' : ('cas', '#FFC1C1')}


        # static paths
        self.CONTINENTAL_DIR = f'{rsrc_root}/continental/'
        self.SUBCONTINENTAL_DIR_EUR = f'{rsrc_root}/subcontinental/eur/'
        self.SUBCONTINENTAL_DIR_EAS = f'{rsrc_root}/subcontinental/eas/'
        self.GENOMES_1000_CONTINENTAL_DIR = f'{rsrc_root}/1k_genomes/wes_v2/continental/'
        self.GENOMES_1000_AMR_DIR = f'{rsrc_root}/1k_genomes/wes_v2/amr/'
        self.GENOMES_1000_AFR_DIR = f'{rsrc_root}/1k_genomes/wes_v2/afr/'
        self.GENOMES_1000_EUR_DIR = f'{rsrc_root}/1k_genomes/wes_v2/eur/'
        self.GENOMES_1000_SAS_DIR = f'{rsrc_root}/1k_genomes/wes_v2/sas/'
        self.GENOMES_1000_EAS_DIR = f'{rsrc_root}/1k_genomes/wes_v2/eas/'
        self.SGDP_CONTINENTAL_DIR = f'{rsrc_root}/sgdp/continental/'
        self.HG38_JSON_CONVERTER = f'{rsrc_root}/genome_ver_converters/hg38tob37/hg38_liftoverAIMs.json'
        self.WES_b37_JSON_CONVERTER = f'{rsrc_root}/genome_ver_converters/b37tohg38/b37.liftover_to_hg38.1kGP.nygc.json'
        self.WGS_b37_JSON_CONVERTER = f'{rsrc_root}/genome_ver_converters/b37tohg38/b37.liftover_to_hg38.1kGP.nygc.json'
        self.SGDP_b37_JSON_CONVERTER = f'{rsrc_root}/genome_ver_converters/b37tohg38/b37_liftover_to_hg38.SGDP.json'

        self.JSON_CONVERTS = {'37': {'WES': {'1kGP': self.WES_b37_JSON_CONVERTER, 'gnomAD': None, 'SGDP': self.SGDP_b37_JSON_CONVERTER}, 'WGS': {'1kGP': self.WGS_b37_JSON_CONVERTER, 'gnomAD': None, 'SGDP': self.SGDP_b37_JSON_CONVERTER}}, '38': {'WES': {'1kGP': None, 'gnomAD': self.HG38_JSON_CONVERTER, 'SGDP': None}, 'WGS': {'1kGP': None, 'gnomAD': self.HG38_JSON_CONVERTER, 'SGDP': None}}}

        self.MATRIX_ATT = '/matrix_attributes/'
        self.ML_MODELS = '/machine_learning_models/'

        # folder list
        self.MATRIX_ATT_DIRS = join_paths(self.CONTINENTAL_DIR, self.MATRIX_ATT) + join_paths(self.SUBCONTINENTAL_DIR_EUR, self.MATRIX_ATT) + join_paths(self.SUBCONTINENTAL_DIR_EAS, self.MATRIX_ATT) + join_paths(self.GENOMES_1000_AMR_DIR, self.MATRIX_ATT) + join_paths(self.GENOMES_1000_AFR_DIR, self.MATRIX_ATT) + join_paths(self.GENOMES_1000_EAS_DIR, self.MATRIX_ATT) + join_paths(self.GENOMES_1000_EUR_DIR, self.MATRIX_ATT) + join_paths(self.GENOMES_1000_SAS_DIR, self.MATRIX_ATT) + join_paths(self.GENOMES_1000_CONTINENTAL_DIR, self.MATRIX_ATT) + join_paths(self.SGDP_CONTINENTAL_DIR, self.MATRIX_ATT)

        self.MODEL_DIRS = join_paths(self.CONTINENTAL_DIR, self.ML_MODELS) + join_paths(self.SUBCONTINENTAL_DIR_EUR, self.ML_MODELS) + join_paths(self.SUBCONTINENTAL_DIR_EAS, self.ML_MODELS) + join_paths(self.GENOMES_1000_AMR_DIR, self.ML_MODELS) + join_paths(self.GENOMES_1000_AFR_DIR, self.ML_MODELS) + join_paths(self.GENOMES_1000_EAS_DIR, self.ML_MODELS) + join_paths(self.GENOMES_1000_EUR_DIR, self.ML_MODELS) + join_paths(self.GENOMES_1000_SAS_DIR, self.ML_MODELS) + join_paths(self.GENOMES_1000_CONTINENTAL_DIR, self.ML_MODELS) + join_paths(self.SGDP_CONTINENTAL_DIR, self.ML_MODELS)

        # n_class
        self.N_CLASSES_CONTINENTAL = 6
        self.N_CLASSES_CONTINENTAL_NYGC = 5
        self.N_CLASSES_SUBCONTINENTAL_EUR = 7
        self.N_CLASSES_SUBCONTINENTAL_EAS = 3
        self.N_CLASSES_1000_GENOMES_AFR = len(self.LABS_GENOMES_1000_AFR)
        self.N_CLASSES_1000_GENOMES_AMR = len(self.LABS_GENOMES_1000_AMR)
        self.N_CLASSES_1000_GENOMES_SAS = len(self.LABS_GENOMES_1000_SAS)
        self.N_CLASSES_1000_GENOMES_EUR = len(self.LABS_GENOMES_1000_EUR)
        self.N_CLASSES_1000_GENOMES_EAS = len(self.LABS_GENOMES_1000_EAS)
        self.N_CLASSES_SGDP_CONTINENTAL = 7
        #self.N_CLASSES_SIM_REGION = 5

        # indicators
        mode = ['gnomAD_continental', 'gnomAD_eur', 'gnomAD_eas', '1kGP_amr', '1kGP_afr', '1kGP_eas', '1kGP_eur', '1kGP_sas', '1kGP_continental', 'SGDP_continental']
        label_order = [self.LABS_CONT, self.LABS_SUBCONT_EUR, self.LABS_SUBCONT_EAS, self.LABS_GENOMES_1000_AMR, self.LABS_GENOMES_1000_AFR, self.LABS_GENOMES_1000_EAS, self.LABS_GENOMES_1000_EUR, self.LABS_GENOMES_1000_SAS, self.LABS_GENOMES_1000_CONTINENTAL, self.LABS_SGDP]
        self.LABS_CONVERTER = dict(zip(mode, label_order))
        n_classes = [self.N_CLASSES_CONTINENTAL, self.N_CLASSES_SUBCONTINENTAL_EUR, self.N_CLASSES_SUBCONTINENTAL_EAS, self.N_CLASSES_1000_GENOMES_AMR, self.N_CLASSES_1000_GENOMES_AFR, self.N_CLASSES_1000_GENOMES_EAS, self.N_CLASSES_1000_GENOMES_EUR, self.N_CLASSES_1000_GENOMES_SAS, self.N_CLASSES_CONTINENTAL_NYGC, self.N_CLASSES_SGDP_CONTINENTAL]
        self.R_DIRS = list(zip(self.MATRIX_ATT_DIRS, self.MODEL_DIRS, n_classes, mode))

        # plotting axis order the plot is 2,5
        axis_order = [0, 2, 4, 6, 7, 5, 3, 8, 1, 9]
        self.axis_loc = dict(zip(mode, axis_order))

        # change names to something more descriptive
        subplt_titles = ['gnomAD_continental', 'gnomAD_eur', 'gnomAD_eas', '1KGP_amr', '1KGP_afr', '1KGP_eas', '1KGP_eur', '1KGP_sas', '1KGP_continental', 'SGDP_continental']
        self.TITLES = dict(zip(mode, subplt_titles))