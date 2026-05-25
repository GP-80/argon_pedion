"""
Fetch flood-season precipitation and reference evapotranspiration (ET₀)
for Argon Pedion from the Open-Meteo archive API (ERA5-Land reanalysis).

WHAT THIS SCRIPT DOES
─────────────────────
For each flood season (defined as October of year Y-1 through April of year Y)
it queries the Open-Meteo archive for daily precipitation and ET₀ at the
centroid of the study area, sums both over the season, and computes the
net water balance (rainfall - ET₀).

WHY NET BALANCE AND NOT JUST RAINFALL?
───────────────────────────────────────
Raw rainfall alone does not determine how much water accumulates in the plain.
In warm, sunny spells between rain events, evapotranspiration removes moisture
from the soil and the water surface, reducing effective flood duration.
Net balance (rainfall - ET₀) is a first-order proxy for the water available
to pond on the plain and is more informative than rainfall alone when comparing
years with different temperature and sunshine patterns (e.g. 2026, where rapid
drying between storms moderated flood extent despite high rainfall totals).

ET₀ (FAO Penman-Monteith reference evapotranspiration) represents potential
evaporation from a short green grass reference surface — it captures the
atmospheric demand for water under the prevailing temperature, humidity,
wind, and solar radiation conditions.

OUTPUT
──────
Outputs/CSV/argon_pedion_rainfall.csv
Columns: year | rain_mm | et0_mm | balance_mm | season
"""

import pathlib   # standard library: file path handling
import requests  # third-party: HTTP requests (pip install requests)
import pandas as pd  # third-party: tabular data (pip install pandas)


# ── Paths ────────────────────────────────────────────────────────
# __file__ is the path of this script file.
# .parent goes one level up (the Python/ folder).
# .parent again goes one more level up (the project root).
# This makes the script work regardless of where you run it from.
ROOT = pathlib.Path(__file__).parent.parent
OUT  = ROOT / 'Outputs' / 'CSV' / 'argon_pedion_rainfall.csv'


# ── Study area & years ───────────────────────────────────────────
# Centroid of the Argon Pedion AOI in WGS84 decimal degrees.
# Open-Meteo returns the nearest ERA5-Land grid point (~9 km resolution),
# so sub-kilometre precision does not matter here.
LAT   = 37.636
LON   = 22.447

# Match the years available in the GEE comparison script.
# Each year Y represents the flood season Oct(Y-1)–Apr(Y).
YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026]


# ── API endpoint ─────────────────────────────────────────────────
# Open-Meteo's historical archive endpoint.
# No API key or account needed — it is free and anonymous.
# Documentation: https://open-meteo.com/en/docs/historical-weather-api
URL = 'https://archive-api.open-meteo.com/v1/archive'


# ── Fetch data for each flood season ────────────────────────────
rows = []

for yr in YEARS:
    # Build ISO 8601 date strings for the flood season window.
    # e.g. for yr=2019: start='2018-10-01', end='2019-04-30'
    start = f'{yr - 1}-10-01'
    end   = f'{yr}-04-30'

    # Make an HTTP GET request to the Open-Meteo API.
    # 'params' is a dictionary of query-string parameters:
    #   latitude / longitude  — study area centroid
    #   start_date / end_date — date range in YYYY-MM-DD format
    #   daily                 — comma-separated list of variables to fetch
    #   timezone              — local timezone for daily aggregation
    #                           (important: a 'day' runs midnight–midnight
    #                            in the specified timezone, not UTC)
    #   timeout=30            — raise an error if the server takes > 30 s
    resp = requests.get(URL, params={
        'latitude':   LAT,
        'longitude':  LON,
        'start_date': start,
        'end_date':   end,
        'daily':      'precipitation_sum,et0_fao_evapotranspiration',
        'timezone':   'Europe/Athens',
    }, timeout=30)

    # Raise a Python exception if the server returned an error code
    # (e.g. 404 Not Found, 500 Server Error). Without this, a failed
    # request would silently produce wrong results.
    resp.raise_for_status()

    # Parse the JSON response. The structure is:
    # {
    #   "daily": {
    #     "time": ["2018-10-01", "2018-10-02", ...],
    #     "precipitation_sum": [3.2, 0.0, 11.4, ...],
    #     "et0_fao_evapotranspiration": [1.1, 1.3, ...]
    #   }
    # }
    daily = resp.json()['daily']

    # Sum each variable over the season, skipping any missing days (None).
    # The 'v for v in list if v is not None' pattern is a Python generator
    # expression — it loops over the list and passes only non-None values
    # to sum(). Missing values occur when ERA5-Land data is not yet available
    # for recent dates (e.g. the tail end of the current year).
    rain    = round(sum(v for v in daily['precipitation_sum']          if v is not None), 1)
    et0     = round(sum(v for v in daily['et0_fao_evapotranspiration'] if v is not None), 1)
    balance = round(rain - et0, 1)

    rows.append({
        'year':       yr,
        'rain_mm':    rain,
        'et0_mm':     et0,
        'balance_mm': balance,
        'season':     f'{start} / {end}'   # human-readable date range
    })

    # Print a one-line summary per year so you can see progress and
    # spot obvious errors (e.g. implausibly low values) while running.
    # The '+' in '{balance:+.0f}' forces a sign prefix (+ or −).
    print(f'{yr}  rain={rain:.0f} mm  ET0={et0:.0f} mm  balance={balance:+.0f} mm')


# ── Save to CSV ──────────────────────────────────────────────────
# pd.DataFrame converts the list of row-dictionaries into a table.
# .to_csv(..., index=False) writes it without a row-number column.
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f'\nSaved: {OUT}')
