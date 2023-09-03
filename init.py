import requests
import os
from dotenv import load_dotenv

load_dotenv()

EIA_API_KEY = os.environ["EIA_API_KEY"]

base_url = "https://api.eia.gov/v2/electricity/facility-fuel/data/?"

params = {"frequency":"monthly",
          "data[0]":"generation",
          "facets[fuel2002][]":"WND",
          "facets[plantCode][]":"56810",
          "facets[primeMover][]":"ALL",
          "start":"2023-01",
          "sort[0][column]":"period",
          "sort[0][direction]":"desc",
          "offset":"0",
          "length":"5000",
          "api_key":EIA_API_KEY}
          
r = requests.get(base_url, params=params)

print(r.status_code)
print(r.text)