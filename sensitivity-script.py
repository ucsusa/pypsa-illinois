import os
import re
import itertools
import shutil

# Path to the base config file
base_config_file = "config.yml"
backup_config_file = "config_backup.yml"

# Backup the original config file to restore it later
shutil.copyfile(base_config_file, backup_config_file)

# Define the variables and their possible values
variable_ranges = {
    'fuel_cost_year': [2023],
    'load_growth': [0, 0.01, 0.02, 0.025],
    'total_demand': [140e6, 185e6],
    'atb_scenario': ['Moderate']  # 'Conservative', 'Advanced'
}

# Create a string of all the rules that must be forced to run
# NOTE: add retrieve_costs when also modifying atb scenarios
forcerun = "retrieve_fuel_costs add_electricity"

# Generate all possible combinations of the variables
variable_combinations = list(itertools.product(
    variable_ranges['fuel_cost_year'],
    variable_ranges['load_growth'],
    variable_ranges['total_demand'],
    variable_ranges['atb_scenario']
))

# ATB scenarios must be handled differently because sub-parameters can't be
# passed to snakemake in the command line
atb_scenario_targets = [
    "scenario: \'Moderate\'",
    "scenario: \'Conservative\'",
    "scenario: \'Advanced\'"
]


def replace_strings_in_config(file_path, target_strings, new_value):
    """
    This function looks for a list of target strings and replaces them with a
    new value (from ChatGPT).
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    pattern = '|'.join(map(re.escape, target_strings))
    updated_content = re.sub(pattern, new_value, content)

    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(updated_content)


run_no = 0
for idx, (fuel_cost_year,
          load_growth,
          total_demand,
          atb_scenario) in enumerate(variable_combinations):

    # Name the scenario
    scenario = (
        f"cost-{fuel_cost_year}_"
        f"growth-{load_growth}_"
        f"demand-{total_demand:.2E}_"
        f"atb-{atb_scenario}"
    )

    # Combine all the top-line config variables
    config = (
        f"time_res=1 "
        f"fuel_cost_year={fuel_cost_year} "
        f"load_growth={load_growth} "
        f"total_demand={total_demand} "
        f"scenario={scenario}"
    )

    # Modify the config file to update the ATB scenario.
    # This is disabled because the model won't run --
    # alternate scenarios don't exist for nuclear
    # atb_scenario_text = f"scenario: \'{atb_scenario}\'"
    # replace_strings_in_config(
    #    base_config_file,
    #    atb_scenario_targets,
    #    atb_scenario_text
    # )

    # Run Snakemake, passing the new variables as configs
    print(f"Running Snakemake with config {run_no} and the following command:")
    print(f"snakemake --cores=4 --forcerun {forcerun} --config {config}")
    os.system(f"snakemake --cores=4 --forcerun {forcerun} --config {config}")
    run_no += 1

    # for debugging:
    # if run_no == 1: break

# Restore the original base config after running all scenarios
shutil.copyfile(backup_config_file, base_config_file)
os.remove(backup_config_file)
print("Restored the original config file.")
