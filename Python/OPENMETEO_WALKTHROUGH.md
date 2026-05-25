# Fetching Climate Data with Open-Meteo
## Argon Pedion · ERA5-Land Precipitation & ET₀

This guide explains how to pull historical climate data for any location using
the Open-Meteo archive API, what the data represents, and how the script works.
No API key or account needed.

---

## What you will produce

`Outputs/CSV/argon_pedion_rainfall.csv` — a table with one row per flood season:

| year | rain_mm | et0_mm | balance_mm | season |
|---|---|---|---|---|
| 2019 | 621.3 | 303.1 | +318.2 | 2018-10-01 / 2019-04-30 |
| 2024 | 389.4 | 533.2 | −143.8 | 2023-10-01 / 2024-04-30 |

---

## Open-Meteo vs Copernicus CDS — which to use?

Both services deliver ERA5-Land reanalysis data (the same underlying model),
but they have very different interfaces:

| | **Open-Meteo** | **Copernicus CDS** |
|---|---|---|
| Account / API key | None — anonymous | Registration + key required |
| Setup time | Zero | 10–15 minutes |
| Response format | JSON (works with `requests`) | NetCDF or GRIB (needs `xarray`, `cfgrib`, `cdsapi`) |
| Variables available | ~30 common weather variables | ERA5, ERA5-Land, CMIP6, satellite products, hundreds of datasets |
| Spatial queries | Point only (nearest grid cell) | Point or spatial bounding box |
| Bulk downloads | Not designed for it | Yes — download regional grids |
| Best for | Quick point queries, prototyping | Large-area analysis, non-ERA5 datasets |

**Rule of thumb:** If you need a time series at one point and ERA5-Land covers
your variable, use Open-Meteo. If you need spatial coverage, hourly data over
a large region, or data from other models (e.g. CMIP6 climate projections), use CDS.

---

## What is ERA5-Land?

ERA5-Land is a global reanalysis dataset produced by ECMWF (European Centre
for Medium-Range Weather Forecasts). "Reanalysis" means historical weather
observations (station data, radiosondes, satellite retrievals) are combined
with a numerical weather model to produce a consistent, gap-free gridded
record going back to 1950.

Key properties:
- **Resolution:** ~9 km grid (0.1° × 0.1°)
- **Temporal coverage:** 1950 – present (~5-day lag for recent data)
- **Timestep:** hourly (Open-Meteo aggregates to daily for us)
- **Temporal resolution caveat:** ERA5-Land captures seasonal and multi-day
  patterns well but can miss the intensity of individual convective storms

---

## Prerequisites

Install the required Python packages if you haven't already:

```
pip install requests pandas
```

Both are standard data-science libraries:
- `requests` — makes HTTP calls to web APIs
- `pandas` — reads/writes CSV and provides the DataFrame table structure

---

## Step 1 — Understand the API URL

Open-Meteo's historical archive endpoint is a simple HTTP GET request:

```
https://archive-api.open-meteo.com/v1/archive?latitude=37.636&longitude=22.447&...
```

You pass parameters in the URL query string. The API returns a JSON object.
You can test any request directly in a browser — paste the URL with parameters
and you will see the raw JSON response.

### Try it now
Paste this into your browser address bar and press Enter:

```
https://archive-api.open-meteo.com/v1/archive?latitude=37.636&longitude=22.447&start_date=2019-01-01&end_date=2019-01-07&daily=precipitation_sum&timezone=Europe%2FAthens
```

You should see a JSON response like:

```json
{
  "latitude": 37.6,
  "longitude": 22.5,
  "daily": {
    "time": ["2019-01-01","2019-01-02","2019-01-03","2019-01-04","2019-01-05","2019-01-06","2019-01-07"],
    "precipitation_sum": [12.4, 0.0, 3.1, 0.0, 8.7, 22.1, 4.3]
  }
}
```

Notice: the returned coordinates (`37.6`, `22.5`) are snapped to the nearest
ERA5-Land grid cell, not exactly your input — this is normal at ~9 km resolution.

---

## Step 2 — Understand the variables

The script requests two daily variables:

### `precipitation_sum`
Total daily precipitation in millimetres, aggregated from hourly ERA5-Land
values. Includes rain and snow-water equivalent.

### `et0_fao_evapotranspiration`
**Reference evapotranspiration (ET₀)** in millimetres per day, computed by
Open-Meteo using the FAO Penman-Monteith equation from ERA5-Land meteorological
fields (temperature, humidity, wind speed, solar radiation).

ET₀ represents the water demand of the atmosphere — how much water a standardised
short green grass surface would evaporate under the current weather conditions.
It is a proxy for how quickly a wet surface loses water between rain events.

### Why net balance?

`balance = rain - ET₀` is a rough but useful metric: it tells you how much
water the atmosphere is delivering to the surface *net of losses*. A season
with high rainfall but also high ET₀ (warm, sunny) may produce less flooding
than a season with moderate rainfall and low ET₀ (cool, overcast).

This matters for Argon Pedion because it is a shallow, permeable karst basin —
water drains quickly, so sustained wet conditions matter more than single events.

---

## Step 3 — Walk through the code

### Paths

```python
ROOT = pathlib.Path(__file__).parent.parent
OUT  = ROOT / 'Outputs' / 'CSV' / 'argon_pedion_rainfall.csv'
```

`__file__` is the full path of the running script (e.g. `G:\...\Python\fetch_rainfall.py`).
`.parent` steps up one directory level. `.parent.parent` steps up two levels to the
project root. The `/` operator on `Path` objects builds sub-paths — it is equivalent
to `os.path.join()` but more readable.

This approach means the script finds the output folder correctly regardless of
which directory your terminal is open in when you run it.

### The flood season loop

```python
for yr in YEARS:
    start = f'{yr - 1}-10-01'
    end   = f'{yr}-04-30'
```

The flood season for "year 2019" is defined as October 2018 through April 2019.
`yr - 1` shifts the start date one calendar year back. This is the hydrological
convention: the flood season is named after the year it ends in.

### The API call

```python
resp = requests.get(URL, params={
    'latitude':   LAT,
    'longitude':  LON,
    'start_date': start,
    'end_date':   end,
    'daily':      'precipitation_sum,et0_fao_evapotranspiration',
    'timezone':   'Europe/Athens',
}, timeout=30)
resp.raise_for_status()
```

`requests.get(url, params=dict)` builds the full query-string URL and sends an
HTTP GET request. The server responds with a JSON body.

`timeout=30` — if the server does not respond within 30 seconds, `requests`
raises a `Timeout` exception instead of hanging forever.

`resp.raise_for_status()` — if the HTTP response code is an error (4xx client
error or 5xx server error), this line raises a Python exception immediately with
a clear error message. Without it, the script would silently continue and crash
later with a confusing error when trying to parse the missing data.

**Why `timezone='Europe/Athens'`?**
ERA5-Land data is stored in UTC. A "daily" sum depends on when the day boundary
falls. Without specifying a timezone, a day running midnight–midnight UTC would
split Greek winter nights across two different calendar days, giving slightly
wrong daily totals. Specifying the local timezone ensures the daily aggregation
matches the calendar your rainfall gauge would use.

### Parsing the response

```python
daily = resp.json()['daily']
rain  = round(sum(v for v in daily['precipitation_sum'] if v is not None), 1)
```

`resp.json()` parses the response body from a JSON string into a Python dictionary.
`['daily']` extracts the `daily` sub-dictionary containing the time series arrays.

`v for v in list if v is not None` is a **generator expression** — a compact loop
that yields each value `v` from the list, but only if `v` is not `None`.
`None` values appear when ERA5-Land data for a date is not yet available
(typically the last few days of a recent query). Passing them to `sum()` would
cause a `TypeError`, so we skip them.

---

## Step 4 — Run the script

From a terminal in the project root:

```
py Python\fetch_rainfall.py
```

You should see output like:

```
2017  rain=487 mm  ET0=271 mm  balance=+216 mm
2018  rain=412 mm  ET0=298 mm  balance=+114 mm
2019  rain=621 mm  ET0=303 mm  balance=+318 mm
...
Saved: G:\PYTHON\ARGON PEDION\Outputs\CSV\argon_pedion_rainfall.csv
```

The `+` or `−` prefix on the balance column tells you at a glance whether the
season was net wet or net dry.

---

## Step 5 — Inspect the output CSV

Open `Outputs/CSV/argon_pedion_rainfall.csv` in Excel or any text editor:

```
year,rain_mm,et0_mm,balance_mm,season
2017,487.3,271.1,216.2,2016-10-01 / 2017-04-30
2018,412.0,297.8,114.2,2017-10-01 / 2018-04-30
2019,621.3,303.1,318.2,2018-10-01 / 2019-04-30
...
```

The `season` column records the exact date range used, which is useful if you
ever change the season definition and need to know which version of the data a
particular CSV was generated with.

---

## Adapting the script for a different location or variable

### Different location
Change `LAT` and `LON` at the top of the script. Decimal degrees, WGS84.
If you have a GeoPackage, extract the centroid using the method in
`GEE/LEVEL 1/LEVEL1_WALKTHROUGH.md` (Part A).

### Different date range
Change the `start` and `end` expressions in the loop. For a full calendar year:

```python
start = f'{yr}-01-01'
end   = f'{yr}-12-31'
```

### Additional variables
Add comma-separated variable names to the `daily` parameter.
Full list of available variables: **https://open-meteo.com/en/docs/historical-weather-api**

Common ones for hydrology:
| Variable name | Meaning |
|---|---|
| `precipitation_sum` | Total daily precipitation (mm) |
| `et0_fao_evapotranspiration` | Reference ET₀ (mm/day) |
| `temperature_2m_mean` | Mean daily air temperature at 2 m (°C) |
| `soil_moisture_0_to_7cm` | Volumetric soil moisture (m³/m³) |
| `snowfall_sum` | Daily snowfall in mm water equivalent |

Example — adding mean temperature:

```python
'daily': 'precipitation_sum,et0_fao_evapotranspiration,temperature_2m_mean',
```

Then in the parsing section:

```python
temp = resp.json()['daily']['temperature_2m_mean']
tmean = round(sum(v for v in temp if v is not None) / len([v for v in temp if v is not None]), 1)
```

---

## Troubleshooting

| Problem | Cause | Fix |
|---|---|---|
| `ConnectionError` | No internet / wrong URL | Check connectivity; verify the URL format |
| `HTTPError 400` | Bad parameter name or date format | Check spelling; dates must be `YYYY-MM-DD` |
| `FileNotFoundError` on save | `Outputs/CSV/` folder missing | Create it: `mkdir -p Outputs/CSV` |
| All values are `None` | Date range too recent | ERA5-Land has a ~5-day lag; reduce `end` date |
| Values look too low | Wrong timezone | Make sure `timezone` matches the study area |

---

## What's next

- `plot_comparison_rainfall.py` — visualises the CSV alongside the GEE flood-extent
  results as a dual-axis chart (bars = flood area, lines = rainfall and net balance)
- Level 2 walkthrough (coming) — the GEE Python API uses a similar requests-based
  pattern but authenticates via a service account and returns image data instead of JSON
