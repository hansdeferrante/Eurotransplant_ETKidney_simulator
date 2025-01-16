
# The ETKidney Simulator

The ETKidney simulator is a discrete event simulator designed for the analysis of ETKAS and ESP, which are the kidney-only allocation programs of Eurotransplant. To this end, the simulator implements the actual rules of ETKAS and ESP used from 2021 onwards. The simulator is written in Python 3.12. In this Git repository, we release the code of the simulator.

Unfortunately, Eurotransplant is not allowed to publicly release the real patient and donor datasets which drive ETKidney simulations. Instead, synthetic (fake) minimal datasets are included in this publicly available repository (fake_data.zip), which all have the prefix fake_. These synthetic datasets do not contain information from real patients or donors from Eurotransplant, although distributions of baseline characteristics should be similar. Additionally, all included HLA data was randomly generated from German haplotype distributions on the ABDR loci only, and the number of patients and donors per center were randomized.

The fake datasets thus cannot be used for research. The sole purpose of the fake datasets is to allow external researchers to work with the ETKidney simulator simulator if they do not have real data from Eurotransplant.

Researchers in simulation of ETKAS / ESP with realistic data may instead send a study proposal to the Eurotransplant Kidney Advisory Committee (ETKAC).

# Installation

Clone the repository and install dependencies using:

```bash
git clone https://github.com/hansdeferrante/etkidney_simulator.git
conda create --name <env_name> --file requirements.txt
conda activate <env_name>
```

# Usage

To run a simulation, use the `simulate_etkas_and_esp.py` file (or `multiprocess_simulations.py` for parallelization). All simulation settings are specified using YAML-files, which contain paths to input data streams (fake data in this repository), the simulation start and end dates, prioritization rules, and parameters for allocation scoring. The format of these YAML-files is discussed at the end of this document. Example YAML-files can be found in the `simulator/sim_yamls` subdirectory.

To simulate allocation according to your settings, update the pointers in Python scripts:

```python
sim_set = read_sim_settings(
    os.path.join(
        es.DIR_SIM_SETTINGS,
        '2025-01-16',
        'CurrentETKAS_validation_2_2.yml'
    )
)
```


# Anticipated workflow for the ETKidney Simulator

The ETKidney simulator can be used to analyze the impact of changes to ETKAS and Eurotransplant allocation rules. Below are examples of how the simulator can be adapted for different policy evaluations.

1. Simulating simple changes to the ETKAS point system

> The ETKAS point system has prioritized HLA matching, the geographical location of the candidate w.r.t. the donor, pediatric candidates, candidates with the High Urgency (HU) status. The points awarded to these components is saved in yml files. For instance, the current ETKAS point system is parametrized by:


```yaml
POINTS_ETKAS:
  YEARS_ON_DIAL: 33.333333
  PREVIOUS_WT_YEARS: 33.333333
  INTERCEPT: 400
  MMB_HLA_A: -66.66666
  MMB_HLA_B: -66.66666
  MMS_HLA_DR: -66.66667
  INTERCEPT:R_PEDIATRIC: 400
  MMB_HLA_A:R_PEDIATRIC: -66.66666
  MMB_HLA_B:R_PEDIATRIC: -66.66666
  MMS_HLA_DR:R_PEDIATRIC: -66.66667
  ET_MMP: 1
  BAL_NAT: 30
  BAL_REG: 1
  MATCH_L: 300
  MATCH_N: 100
  R_PEDIATRIC: 100
  PATIENT_IS_HU: 500
```

> end users of the ETKidney simulator may adapt these points to simulate ETKAS under alternative scoring systems. For instance, for the first case study included in the accompanying publication, we changed the relative weighting placed on ABDR-matching:

```yaml
  INTERCEPT: 400
  MMB_HLA_A: 0
  MMB_HLA_B: -100
  MMS_HLA_DR: -100
  INTERCEPT:R_PEDIATRIC: 400
  MMB_HLA_A:R_PEDIATRIC: 0
  MMB_HLA_B:R_PEDIATRIC: -100
  MMS_HLA_DR:R_PEDIATRIC: -100
```


2. Adding other (transformed) variables to the ETKAS point system

> In the United States, immunized candidates are prioritizes using a sliding scale (an exponential function). This sliding scale is parametrized by a weight (the maximum number of points for a vPRA of 100%) and a base, which controls the steepness of the sliding scale. To simulate such a sliding scale with the ETKidney simulator, one can add the vPRA as a component to the ETKAS point system, and transform this vPRA using a to-be-defined function. For instance, to simulate kidney allocation with a sliding scale with base 2 and weight 100, one has to specify:

```yaml
POINTS_ETKAS:
  VPRA: 100
TRAFOS_ETKAS:
  VPRA: exponential_base2
```

> The transformation `exponential_base2` has to be defined in the file `etkidney_simulator_settings.py`


3. Changes to ESP Allocation Rules

> ESP allocation tiers and age-based criteria are specified in the input yml file and in `magic_values/ESP_match_tiers.yml`. To simulate changes to ESP prioritization, these files can be modified.


# Data

The ETKidney simulator is a data-driven simulator. Users of the ETKidney simulator thus have to provide paths to input streams for the simulator to work. Four input files are of critical importance, which are:
1. the patients file
2. the donors file
3. the status update file
4. the allocation profile file

We describe the requirements for these files here. ETKidney simulations are preferably based on realistic (historic) data from Eurotransplant.
Because Eurotransplant is not allowed to make such data publicly available, synthetic (fake) data is released with the ETKidney simulator.

## Donors File (by default `fake_donors_for_etkidney_simulator.csv`):

This defines static information relating to the donors. This includes:

- The date the donor was reported to Eurotransplant (`d_date`)
- Unique identifier assigned to each donor (`id_donor`)
- The country where the donor is reported from (`d_country`)
- The specific region within the donor’s country (`d_region`)
- The ESP subregion from which the donor was reported (`d_center_subregion_esp`)
- The center responsible for the procurement of the donor’s organ (`d_center`)
- Detailed donor profile information, such as:
  - Age of the donor (`donor_age`)
  - Sex of the donor (`d_sex`)
  - Donor’s blood group (`d_bloodgroup`)
  - Whether the donor tested positive for specific markers like HBsAg (`graft_hbsag`), HCV antibodies (`graft_hcvab`), and HBcAb (`graft_hbcab`)
  - Whether the donor had specific conditions like sepsis (`graft_sepsis`) or meningitis (`graft_meningitis`)
  - Additional information about the donor, such as history of malignancy (`donor_malignancy`), drug abuse (`donor_drug_abuse`), marginality (`donor_marginal_free_text`), and history of tumors (`donor_tumor_history`)
  - Donor HLA typing in ETRL match determinants (`donor_hla`)
  - Cause of death group (`death_cause_group`)
  - Whether the donor was a DCD donor (`graft_dcd`)

## Patients File (by default `fake_patients_for_etkidney_simulator.csv`):

This defines the static information relating to listings of transplantation candidates. Listings are uniquely defined by a registration ID (`id_registration`) and are coupled to a patient by a recipient ID (`id_recipient`). Information required in this file is:

- The date of waiting list activation (`inc_date_time`)
- Time to deregistration (`time_to_dereg`), measured in days from the registration date
- The deregistration reason (`removal_reason`)
- The time since the candidate's previous transplantation at listing (`time_since_prev_txp`), measured in days
- Whether the previous transplantation was a living transplantation (`prev_txp_living`)
- Whether the previous transplantation was pediatric (`prev_txp_ped`)
- How many previous kidney transplantations the candidate has had (`n_previous_transplants`)
- In which center the candidate is registered (`recipient_center`), in which region (`recipient_region`), and in which country (`recipient_country`)
- Candidate information such as the candidate’s blood group (`r_bloodgroup`), date of birth (`r_dob`), and sex (`patient_sex`)
- The candidate’s age at listing (`age_at_listing`)
- The candidate’s HLA typing at listing in ETRL match determinants, if known (`r_hla`)
- A choice for ETKAS or ESP, which are mutually exclusive in Germany (`wlki_programme`)
- The first date of dialysis, which is necessary to calculate accrued dialysis time (`date_first_dial`)
- Returned dialysis time measured in days (`previous_t`). Candidates may receive such returned dialysis time in case they were re-listed for a repeat kidney transplantation shortly after the primary transplantation.


## Status Updates File (by default `fake_patstat1_for_etkidney_simulator.csv`):

This file includes all changes to a candidate's status. Required columns for this file are:

- Registration ID (`id_registration`), used to couple status updates to patient listings (see the patients file)
- The listing date (`inc_date_time`)
- `tstart`, which indicates when an update was reported relative to the listing date (in days)
- Update type (`type_update`), which indicates the type of status update (see next section)
- Additional details in `variable_value`, `variable_detail`, and `variable_detail2`, depending on update type

### Update Types for the ETKidney Simulator:

- **PRA**: Updates to panel reactive antibody levels.
	+ variable_value: valid (1) or invalid (0) antibody screening
	+ variable_detail: PRA (%)
	+ variable_detail2: type of screening
- **HLA**: Updates to HLA typing
	- variable_value: full HLA typing, in HLA match determinants
	- variable_detail: 1-ABDR HLA mismatch frequency)
- **DIAG**: Updates to the candidate’s primary disease group
	- variable_value: diagnosis
	- variable_detail: diagnosis group)
- **DIAL**: Changes to candidate's reported date of dialysis
	- variable_value: dialysis date
- **URG**: Changes to the candidate's urgency status.
	+ variable_value: urgency code (transplantable `T`, non-transplantable `NT`, immunized `I`, highly immunized `HI`, removed `R`, death `D`, or transplanted `FU`)
- **AM**: Changes to AM status:
	- variable_value: `A` for active, or `N` for ineligible
- **PRF**: update to the candidate's allocation profile
	+ variable_value:
- **UNACC**: Changes to a candidate's unacceptable antigens.
	+ variable_value: space-separated string of unacceptables (ETRL match determinants)
	+ variable_detail: vPRA (%)

## The allocation profile and HLA mismatch criteria file (by default `fake_profiles_for_etkidney_simulator.csv`):

This file contains minimal HLA mismatch criteria and allocation profile specified for patients. Minimal HLA mismatch criteria are used by centers to indicate that their candidate does not want to be considered for transplantation in case of a particular HLA match quality (for instance, `222` indicates 2 mismatches on all HLA-ABDR loci). The allocation profile is used to exclude offers from donors with certain characteristics (for instance, age and virology).

## Optional Input Files:

### Donor pool file (by default `./fake_donorpool_for_etkidney_simulator.csv`):

This file should contain 10,000 donors with their HLA-typings. The HLA system module can calculate the vPRA and mismatch probability against this pool, in case such information is not available from input streams.

Ideally this dataset is similar to the dataset used by the ETRL to calculate vPRAs.
Eurotransplant cannot release such data publicly, so a fake data was created based on German ABDR haplotype distributions.

### Historic international transplantations (by default `./fake_donors_for_balances_for_etkidney_simulator.csv`):

This file contains a list of international transplantations which occurred in the past. This file can optionally be supplied to initialize the Eurotransplant balance systems, and schedule balance update events.

Synthetic (fake) example datasets are provided in the repository for demonstration purposes.

# YAML Files

## Configuration Requirements for the ETKidney Simulator

The YAML file used in the ETKidney simulator must specify paths to the required datasets:

- **`PATH_DONORS`**: Path to the donor input dataset.
- **`PATH_PATIENTS`**: Path to static information on patients.
- **`PATH_STATUS_UPDATES`**: Path to the patient status update input stream used for simulation.
- **`PATH_PROFILES`**: Path to allocation profiles and HLA mismatch criteria.

- **`PATH_BALANCES`**: Path to donor balance information for the simulation.
- **`PATH_PROGRAM_UPDATES`**: Path to program choices for German candidates.
- **`PATH_DONOR_POOL`**: Path to the donor pool, against which mismatch probabilities and the vPRA can be calculated.

The YAML file also contains paths to general simulation settings, which can be modified by end users of the ETKidney simulator.

- **`PATH_MATCH_TABLE`**: Path to the ETRL match determinants. Glyphs defined in this file are recognized as HLAs by the simulator. If a typing is only partially provided (for instance, only the allele), the simulator translates this typing into the full typing (allele-split-broad) based on ETRL match tables.
- **`PATH_ETKAS_TIERS`**: Path to ETKAS match tiers file, used to rank candidates in ETKAS
- **`PATH_ESP_TIERS`**: Path to ESP match tiers file, used to rank candidates in ESP.
- **`PATH_ALLELE_FREQUENCIES_BROAD`**: Path to broad allele frequency data, used to calculate mismatch probabilities
- **`PATH_ALLELE_FREQUENCIES_SPLIT`**: Path to split allele frequency data, used to calculate mismatch probabilities
- **`PATH_BG_FREQUENCIES`**: Path to blood group frequency data, used to calculate mismatch probabilities
- **`PATH_GEOGRAPHIC_DEFS`**: Geographic clustering of transplantation centers within Eurotransplant, used to priorities in ESP and non-standard allocation.
- **`PATH_MATCH_POTENTIALS`**: Path to a file with pre-calculated HLA 1-ABDR mismatch frequencies, which can speed up simulations.


Settings which control how the simulation modules are initialized are:

- **`SEED`**: Sets the random seed for simulation.
- **`SIM_START_DATE`** and **`SIM_END_DATE`**: Define the simulation time window.

> Graft offering module:

- **`SIMULATE_RESCUE`**: Boolean flag, whether to simulate extended/rescue allocation.
- **`CENTER_ACC_POLICY`** and **`PATIENT_ACC_POLICY`**: Define the acceptance policies for centers and patients. Options include:
  - `LR`: Simulate with logistic regression-based acceptance.
  - `always`: Simulate always accepting offers.
- **`ALLOW_DISCARDS`**: Boolean flag to allow/disallow discards of unaccepted kidneys during the simulation.
- **`SIMULATE_RANDOM_EFFECTS`**: Whether to simulate random effects for acceptance predictions.
- **`VARCOMPS_RANDOM_EFFECTS`**: Variance components for random effects in graft offer acceptance.
- **`MAX_OFFERS_PER_CENTER`**: Limits on the number of offers per center for different allocation systems (separately for ESP and ETKAS).

> Post-transplant module:

- **`LOAD_RETXS_FROM`** and **`LOAD_RETXS_TO`**: Controls which repeat kidney transplantation candidates are loaded from input streams. These repeat kidney transplantations are used by the post-transplant module to create synthetic re-listings.

> Balance system:
- **`WINDOW_FOR_BALANCE`**: Number of days after which a to-be-balanced transplantation expires. This was 365 days until April 1st 2019, after which to-be-balanced transplantations no longer expire.
- **`STARTTIME_FOR_BALANCE`**: Date back to which balances are calculated (2019-04-01 in reality).
- **`REMOVE_BALANCES_BEFORE_STARTTIME`**: Boolean flag, whether to ignore balances before the starttime.
- **`BALANCE_GROUP_VARS`**: Tuple of variables on which the balances are stratified. Since April 1st 2019, donor age group. Since March 2023, DBD/DCD donation as well.

ETKAS and ESP rules which are directly specifiable in the yml file are:

- **`POINTS_ETKAS`**: Settings for ETKAS scoring.
- **`POINTS_ESP`**: Settings for ESP scoring.
- **`PEDIATRIC_CANDIDATE_AGE`** and **`PEDIATRIC_DONOR_AGE`**: Age thresholds for pediatric classification of candidates and donors. Patient age was 16 until March 2021. Donor 16 until March 2023. Both are currently 18.
- **`PEDIATRIC_REQUIRES_DIALYSIS`**: Boolean, indicating whether a candidate is only pediatric if they receive dialysis before turning 18. This was the case until March 2023.
- **`DONOR_AGE_ESP_ELIGIBLE`**: Specifies ESP donor age eligibility by country.
- **`CAND_AGE_ESP_ELIGIBLE`**: Specifies ESP candidate age eligibility by country.

- **`LOCI_ZERO_MISMATCH`**: HLA-loci which must have 0 mismatches for a candidate to be zero-mismatched.
- **`LOCI_MMP_SPLIT`** and **`LOCI_MMP_BROAD`**: Loci for calculating the mismatch probabilities on the split and broad level, resp.

The YAML file also specifies settings and filenames for simulation outputs:
- **`RESULTS_FOLDER`**: Specifies the output folder for simulation results.
- **`SAVE_MATCH_LISTS`**: Boolean flag indicating whether to save match lists to output files.
- **`STORE_SCORE_COMPONENTS`**: Boolean flag, indicating whether to also save the number of points awarded per component in match lists.
- **`PATH_TRANSPLANTATIONS`**: File to which transplantation data is written.
- **`PATH_EXITS`**: File to record patient exits from the waitlist.
- **`PATH_FINAL_PATIENT_STATUS`**: File for final patient statuses.
- **`PATH_DISCARDS`**: File for discarded grafts (if applicable).
- **`PATH_MATCH_LISTS`**: File to store match lists (if `SAVE_MATCH_LISTS` is `True`).

Example YAML templates for the case studies conducted in the accompanying paper can be found in the `simulator/sim_yamls` directory.


# YAML-templates

We anticipate that end-users of the ETKidney simulator will want to simulate a scenario multiple times to assess variability in outcomes. To keep simulation settings the same across yaml-files, it is helpful to use yaml-templates as is illustrated in the `sim_yamls/templates/` subfolder. Attribute values in these templates are parametrized with double brackets. An example of such a template is:


	SEED: {{seed}}
	PATH_TRANSPLANTATIONS: transplantations_k{{kth_stat}}_s{{sim_seed}}.csv

For a seed of 1 and the first status update file, this template will lead to the following yaml-file:

	SEED: 1
	PATH_TRANSPLANTATIONS: transplantations_k1_s1.csv

There is an R-markdown file in the `sim_yamls/` directory, which construct YAML-files necessary for the validation of the ETKidney simulator. Also included are three templates for scenarios evaluated as part of the case studies.

# Output

Simulation results include:

- All simulated transplantations
- Waiting list outcomes for all candidates
- List of discarded kidneys
- List of all exits (deaths, removals, transplantations)


# License

This project is licensed under the GNU General Public License.

# Contact

For questions and support, please open an issue in the repository, or contact Eurotransplant.
