# Variables to change for sensitivities

- **Coal and methane prices.** Fuel prices are modeled as 12 values from a historical year, derived from this PUDL database: https://data.catalyst.coop/pudl/core_eia923__monthly_fuel_receipts_costs.
- **ATB cases.** The model uses the 2022 ATB, which has three scenarios: "Moderate", "Conservative", and "Advanced". For sensitivity analysis we can ignore "Moderate" and just use "Conservative" and "Advanced".
- **Demand growth.** This will be modeled as either 0.00 or 1.00.
- **Exports.** There's no toggle or setting for this -- rather the starting load determines whether or not the model is matching the in-state load or the in-state generation. For no export the starting value should be 136e6; for export it should be 185e6.

| Parameter   | Variable in `config.yml`| Scripts                                       | Values                           |
|-------------|-------------------------|-----------------------------------------------|----------------------------------|
|Fuel prices  |`fuel_cost_year`         |`retrieve_fuel_prices.py`, `add_electricity.py`|2018, 2019, 2020, 2021, 2022, 2023|
|ATB cases    |`scenario`               |`retrieve_costs.py`                            |"Conservative", "Advanced"        |
|Demand growth|`load_growth`            |`add_electricity.py`                           |0.00, 1.00                        |
|Exports      |`total_demand`           |`add_electricity.py`                           |136e6, 185e6                      |

Total possibilities = 6x2x2x2 = 48

The order in the snakemake workflow of the rules we need to modify:

| Rule | Script                  |
|------|-------------------------|
|3     |`retrieve_costs.py`      |
|4     |`retrieve_fuel_prices.py`|
|11    |`add_electricity.py`     |