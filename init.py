# load system packages
import os
from pathlib import Path

# load custom packages (requirements)
from dotenv import load_dotenv
import pandas as pd
import matplotlib.pyplot as plt

# load local packages
from downloader import get_era5_monthly, get_EIA_plant_information, get_EIA_plant_generation

# load environment variables
load_dotenv()
EIA_API_KEY = os.environ["EIA_API_KEY"]

# get EIA plant information for wind assets
plant_info_file=Path("data/plant_list.csv").resolve()
if not plant_info_file.is_file():
    df_plant = get_EIA_plant_information(EIA_API_KEY)
    df_plant.to_csv(plant_info_file)
else:
    df_plant = pd.read_csv(plant_info_file)

print(df_plant)

# select an EIA plant for analysis
selected_plant = df_plant.iloc[1]
print(selected_plant)

# get the EIA generation for the selected plant
plant_data_file=Path(f"data/plant_generation_{selected_plant['plantid']}.csv").resolve()
if not plant_data_file.is_file():
    df_generation = get_EIA_plant_generation(EIA_API_KEY,plantCode=selected_plant["plantid"])
    df_generation.to_csv(plant_data_file)
else:
    df_generation = pd.read_csv(plant_data_file)

# set index to converted datetime
df_generation = df_generation.set_index("period")
df_generation.index = pd.to_datetime(df_generation.index)
print(df_generation)

# get era5 data for the selected plant
df_era5 = get_era5_monthly(
    lat=selected_plant["latitude"],
    lon=selected_plant["longitude"]%360,
    save_pathname="data",
    save_filename=selected_plant["plantid"],
    data_pathname="data/era5_monthly",
    data_filename="era5_monthly",
    start_date="2000-01",
)

print(df_era5)

# combine wind speed and generation into single dataframe
df = pd.concat([df_era5[["windspeed_ms"]],df_generation[["generation","gross-generation"]]],axis=1)

print(df.tail(20))

# add days in month column
df["days-in-month"] = df.index.to_series().dt.days_in_month

# add availability column
df["availability-reported"] = df["generation"]/df["gross-generation"]

# add normalised generation column
df["generation-normalised"] = df["gross-generation"]*30/df["days-in-month"]

print(df.tail(20))

# fig, ax = plt.subplots(figsize=(20,10)) 
# df.plot(y = "generation", ax = ax) 
# df.plot(y = "gross-generation", ax = ax) 
# df.plot(y = "windspeed_ms", ax = ax, secondary_y = True) 
# plt.show()

fig, ax = plt.subplots(figsize=(20,10)) 
df.plot.scatter(x="windspeed_ms",y="generation-normalised",ax=ax)
plt.show()

# Next steps:
# 1) DONE - normalise generation (days in month, gross generation, degradation)
# 2) plot correlation
# 3) remove outliers
# 4) linear fit and synthesise generation
# 5) LT corrections
# 6) P50 output with availability shown

# Improvements:
# Done - ERA5 download for whole of world to speed up results
# EIA download for all plants to speed up results
# MERRA2 included as well
# Cover wind and solar
# Degradation analysis
# Simple front-end
# Examples
# Documentation
# Extend to other open generation datasets 