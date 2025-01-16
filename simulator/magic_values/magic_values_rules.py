# Country names
AUSTRIA = 'Austria'
BELGIUM = 'Belgium'
CROATIA = 'Croatia'
GERMANY = 'Germany'
HUNGARY = 'Hungary'
LUXEMBOURG = 'Luxembourg'
NETHERLANDS = 'Netherlands'
SLOVENIA = 'Slovenia'

# Center codes
CENTER_VIENNA = 'AWGTP'
CENTER_UPPERAUSTRIA = 'AOETP'

# Immunization statuses
T = 'T'
HU = 'HU'
NT = 'NT'
HI = 'HI'
IMM = 'I'

# Other statuses
FU = 'FU'
R = 'R'
D = 'D'

# Blood group
BG_A = 'A'
BG_B = 'B'
BG_AB = 'AB'
BG_O = 'O'

# Age groups for balances
LT18 = '0-17'
YOUNGADULT = '18-49'
OLDADULT = '50-64'
ELDERLY = '≥65'

# HLA types
HLA_A = 'hla_a'
HLA_B = 'hla_b'
HLA_C = 'hla_c'
HLA_DR = 'hla_dr'
HLA_DR345 = 'hla_drb345'
HLA_DQB = 'hla_dqb'
HLA_DQA = 'hla_dqa'
HLA_DPB = 'hla_dpb'
HLA_DPA = 'hla_dpa'
ALL_HLA_LOCI = (
    HLA_A, HLA_B, HLA_C, HLA_DR, HLA_DR345,
    HLA_DQB, HLA_DQA, HLA_DPB, HLA_DPA
)


PUBLICS = tuple('Cw' + str(i) for i in range(20)) + ('Bw6', 'Bw4')
PUBLICS = set(tuple(x.upper() for x in PUBLICS))

# Unsplittables
DR3 = 'DR3'
DR17 = 'DR17'
DR18 = 'DR18'

IMPORT = 'import'
EXPORT = 'export'
NATIONAL = 'national'
REGIONAL = 'regional'

MMB_HLA_A = 'mmb_hla_a'
MMB_HLA_B = 'mmb_hla_b'
MMB_HLA_DR = 'mmb_hla_dr'

MMS_HLA_A = 'mms_hla_a'
MMS_HLA_B = 'mms_hla_b'
MMS_HLA_DR = 'mms_hla_dr'

NHS_MISMATCH_LEVEL = 'nhs_mm_level'
ET_MISMATCH_LEVEL = 'et_mm_level'
YPB_MISMATCH_LEVEL = 'ypb_mm_level'

FH = 'fh'
DR_plus = 'dr+'
DR0 = 'dr0',
DR1 = 'dr1',
DR2 = 'dr2'

PRF = 'PRF'
DIAG = 'DIAG'
HLA = 'HLA'
UNACC = 'UNACC'
URG = 'URG'
AM = 'AM'
DIAL = 'DIAL'
PRA = 'PRA'
PRG = 'PRG'

ESP = 'ESP'
ETKAS = 'ETKAS'

ESP_TIER_A = 'A'
ESP_TIER_B = 'B'
ESP_TIER_C = 'C'
ESP_TIER_D = 'D'
ESP_TIER_E = 'E'
ESP_TIER_F = 'F'
ESP_TIER_G = 'G'
ESP_TIER_H = 'H'
ESP_TIER_I = 'I'
ESP_TIERS = (
    ESP_TIER_A, ESP_TIER_B, ESP_TIER_C, ESP_TIER_D, ESP_TIER_E,
    ESP_TIER_F, ESP_TIER_G, ESP_TIER_H, ESP_TIER_I
)
ESP_TIER_REVERSAL_DICT = dict(
    zip(ESP_TIERS, ESP_TIERS[::-1])
)
FMR_TIERS_ESP = {
    ESP_TIER_A, ESP_TIER_B,
    ESP_TIER_C, ESP_TIER_D
}


COUNTRIES_SUBREGIONS = {GERMANY, HUNGARY}
