# load system packages
import os

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
df_plant = get_EIA_plant_information(EIA_API_KEY)
print(df_plant)

# select an EIA plant for analysis
selected_plant = df_plant.iloc[0]
print(selected_plant)

# get the EIA generation for the selected plant
df_generation = get_EIA_plant_generation(EIA_API_KEY,plantCode=selected_plant["plantid"])

# set index to converted datetime
df_generation = df_generation.set_index("period")
df_generation.index = pd.to_datetime(df_generation.index)
print(df_generation)

# get era5 data for the selected plant
df_era5 = get_era5_monthly(
    lat=selected_plant["latitude"],
    lon=selected_plant["longitude"],
    save_pathname="data/",
    save_filename=selected_plant["plantid"],
    start_date="2000-01",
)

print(df_era5)

df = pd.concat([df_era5[["windspeed_ms"]],df_generation[["generation","gross-generation"]]],axis=1)

print(df)

fig, ax = plt.subplots(figsize=(20,10)) 

df.plot(y = 'generation', ax = ax) 
df.plot(y = 'gross-generation', ax = ax) 
df.plot(y = 'windspeed_ms', ax = ax, secondary_y = True) 

plt.show()