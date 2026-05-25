# Level 1 — Build the GEE Flood Mapping Script from Scratch
## Argon Pedion · Sentinel-1 SAR · Google Earth Engine Code Editor

This guide walks you through **writing** the script yourself, line by line.
Each step shows the code to add and explains exactly what it does and why.
No prior GEE or remote sensing experience assumed.

---

## What you will build

A Google Earth Engine JavaScript script that:
1. Defines your study area from a GeoPackage polygon
2. Loads all Sentinel-1 radar scenes over the plain during flood season (Oct–Apr)
3. Classifies each scene into water / not-water using a backscatter threshold
4. Finds the peak flood scene per year and reports the flooded area in km²
5. Displays one coloured map layer per year and a bar chart of annual maxima
6. Adds a flood-frequency layer showing the most persistently flooded pixels

---

## Before you start

### Create a GEE account
GEE is free for non-commercial use.
1. Go to **https://earthengine.google.com** → **Get started**
2. Register a **Noncommercial** project — you need a Google account
3. Approval is usually immediate for research/education use
4. Once approved, open the Code Editor at **https://code.earthengine.google.com**

### Create a new script
1. In the left panel, click the **Scripts** tab
2. Click **+** next to "Owner" → name the script `argon_pedion_s1` → **OK**
3. A blank editor opens — you will write Earth Engine code here as you go
4. Save at any point with **Ctrl+S**

### The Code Editor layout
```
┌─────────────┬──────────────────────────────┬────────────────────┐
│ Left panel  │   Code editor (centre)        │  Right panel       │
│             │                               │                    │
│ Scripts     │  ← write your JavaScript here │  Inspector         │
│ Docs        │                               │  Console  ← output │
│ Assets      │                               │  Tasks    ← exports│
└─────────────┴──────────────────────────────┴────────────────────┘
│                       Map  (bottom half)                         │
│   Layers panel ──────────────────────────────────────────────►   │
└──────────────────────────────────────────────────────────────────┘
```
The **Run** button is at the top of the code editor.
The **Console** tab (right panel) shows `print()` output.
The **Layers** panel (top-right of the map) toggles each map layer on/off.

---

## Part A — Extract AOI coordinates from your GeoPackage

Your study area polygon lives in `floodplain.gpkg` in GGRS87 (EPSG:2100),
a projected coordinate system used for Greek national mapping.
GEE only understands WGS84 (EPSG:4326) — standard decimal-degree GPS coordinates.
You must reproject the polygon and extract its vertex coordinates before you can
define the AOI in GEE.

### What is a GeoPackage?
A `.gpkg` file is a self-contained SQLite database that stores vector geometries
(points, lines, polygons) and their attribute data. It is the modern replacement
for Shapefiles. Libraries like `geopandas` can read it directly in Python.

### What is GGRS87 vs WGS84?
| System | EPSG code | Unit | Used for |
|---|---|---|---|
| GGRS87 (Greek Grid) | 2100 | metres (Easting, Northing) | Greek national maps, cadastre |
| WGS84 | 4326 | decimal degrees (Lon, Lat) | GPS, Google Maps, GEE |

A GGRS87 coordinate like `(363374, 4166038)` has no meaning in GEE — you need
to convert it to something like `(22.453, 37.634)` first.

### Python script to extract and reproject the coordinates

Create a new file `Python/extract_aoi_coords.py` and run it once:

```python
"""
Reads floodplain.gpkg, reprojects from GGRS87 to WGS84, and prints
the polygon vertices in GEE JavaScript format ready to paste.
"""

import pathlib
import geopandas as gpd

ROOT = pathlib.Path(__file__).parent.parent
GPKG = ROOT / 'floodplain.gpkg'

# Read the GeoPackage — geopandas auto-detects the layer name
gdf = gpd.read_file(GPKG)

print(f'Original CRS : {gdf.crs}')          # should show EPSG:2100 or GGRS87
print(f'Feature count: {len(gdf)}')          # should be 1 (one polygon)

# Reproject to WGS84
gdf_wgs = gdf.to_crs('EPSG:4326')

# Extract the exterior ring of the first (only) polygon
polygon = gdf_wgs.geometry.iloc[0]
coords  = list(polygon.exterior.coords)      # list of (lon, lat) tuples

# Print in GEE format — ready to paste directly into the script
print('\n--- paste this into GEE ---\n')
print('var aoi = ee.Geometry.Polygon([')
print('  [[', end='')
for i, (lon, lat) in enumerate(coords):
    sep = ',\n    ' if i > 0 else ''
    print(f'{sep}[{lon:.15f}, {lat:.15f}]', end='')
print(']]')
print(');')
```

Run it from a terminal (make sure geopandas is installed: `pip install geopandas`):

```
cd "G:\PYTHON\ARGON PEDION"
py Python\extract_aoi_coords.py
```

You will see output like:

```
Original CRS : EPSG:2100
Feature count: 1

--- paste this into GEE ---

var aoi = ee.Geometry.Polygon([
  [[[22.457559040088586, 37.66680064822952],
    [22.461397929816396, 37.66641575289452],
    ...
    [22.457559040088586, 37.66680064822952]]]
]);
```

Copy everything from `var aoi = ...` to the closing `);` and keep it ready.

> **Why is the first coordinate repeated at the end?**
> A polygon's exterior ring must be "closed" — the last point must equal the
> first. `polygon.exterior.coords` includes this closing point automatically.

> **Alternative (no Python): QGIS**
> Open the .gpkg in QGIS → right-click the layer → Export → Save Features As →
> Format: GeoJSON, CRS: EPSG:4326 → open the resulting file in a text editor.
> The `coordinates` array inside the JSON is in `[lon, lat]` order, ready to
> paste into GEE (just wrap it in `ee.Geometry.Polygon([...])`).

---

## Part B — Write the script, section by section

Open your blank `argon_pedion_s1` script in the GEE Code Editor.
Add each block below in order, reading the explanation before you type.

---

### Block 1 — Header comment

Every script should start with a plain-text description of what it does.
`//` marks a single-line comment — the Code Editor ignores everything after it on that line.

> **What language is this?**
> The GEE Code Editor uses the **Earth Engine JavaScript API** (sometimes called EEJS).
> The syntax — `var`, `function`, `for`, `//` comments — is plain JavaScript.
> But the `ee.*` objects (`ee.ImageCollection`, `ee.Filter`, `ee.Reducer`, etc.) are
> GEE's own server-side API. When you write `ee.ImageCollection(...).filter(...)`,
> you are not processing data locally — you are building a description of a computation
> that GEE's servers execute remotely. This lazy, server-side evaluation model is what
> makes GEE different from writing plain JavaScript in a browser.

```javascript
// ================================================================
// Argon Pedion — Annual Maximum Flood Extent
// Sentinel-1 SAR (VV, IW, GRD) · 2014–2025
// ================================================================
// Classifies open water in each flood-season (Oct–Apr) S1 scene,
// finds the peak flood scene per year, and maps it.
// ================================================================
```

---

### Block 2 — Area of Interest (AOI)

Paste the `var aoi = ...` block you generated in Part A.

```javascript
var aoi = ee.Geometry.Polygon(
  [[[22.457559..., 37.666800...],
    ...
    [22.457559..., 37.666800...]]]   // closing point = first point
);
```

**What is happening here?**

- `ee.Geometry.Polygon()` is a GEE constructor — it creates a geometry object
  the rest of your code can use for filtering and clipping.
- The coordinate format is `[longitude, latitude]` — longitude (east-west) first,
  latitude (north-south) second. This is the opposite of how most people say it
  ("latitude, longitude") but it is standard in GeoJSON and GEE.
- The triple nesting `[[[ ... ]]]` is GeoJSON MultiPolygon syntax — the outer
  two brackets allow for polygons with holes (islands inside lakes, for example).
  For a simple polygon with no holes, you just have the one ring.

---

### Block 3 — Parameters

Collect all adjustable settings at the top so you can experiment without
editing the rest of the code.

```javascript
// ── Parameters ──────────────────────────────────────────────────
var SAR_THRESHOLD  = -15;   // dB — pixels below this = water
var SPECKLE_KERNEL = 7;     // size of the noise-smoothing window (pixels)
var YEARS  = [2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025];
var COLORS = ['#a6cee3','#1f78b4','#b2df8a','#33a02c',
              '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
              '#cab2d6','#6a3d9a','#ffff99','#b15928'];
```

**SAR_THRESHOLD = −15 dB**
Radar backscatter over calm open water is typically between −20 and −15 dB.
Any pixel below −15 dB is classified as water. If you later see ploughed fields
or gravel roads being wrongly detected, change this to −18 to be stricter.

**SPECKLE_KERNEL = 7**
The smoothing window is 7×7 pixels. A pixel is replaced by the average of itself
and the 48 neighbours around it. Larger values give cleaner edges but blur detail.

**YEARS**
A plain JavaScript array of the years you want to analyse. Sentinel-1 coverage
in GEE starts April 2014, so 2014 is the earliest full flood season available
(Oct 2013 data does not exist; only Jan–Apr 2014).

**COLORS**
Twelve hex colour codes — one per year. These come from the ColorBrewer "Paired"
palette, designed so adjacent colours are easily distinguishable. At 60% opacity
the colours blend where flood extents overlap, revealing the persistently flooded
core of the plain.

---

### Block 4 — Flood-season filter

```javascript
// ── Flood-season filter (Oct–Apr) ───────────────────────────────
var floodSeasonFilter = ee.Filter.or(
  ee.Filter.calendarRange(10, 12, 'month'),   // October, November, December
  ee.Filter.calendarRange(1,  4,  'month')    // January, February, March, April
);
```

**Why not just use a date range?**
The flood season straddles the calendar year: October–December of year Y plus
January–April of year Y+1. A single date range like `2019-10-01` to `2020-04-30`
would work for one winter but not for the whole multi-year collection.
`ee.Filter.or()` combines two independent month filters, making a reusable rule
that applies to any year in the collection.

**ee.Filter**
A Filter is GEE's way of selecting which images to keep. Here:
- `calendarRange(10, 12, 'month')` — keep images where the acquisition month is 10, 11, or 12
- `calendarRange(1, 4, 'month')` — keep images where the month is 1, 2, 3, or 4
- `ee.Filter.or(A, B)` — keep images that pass A **or** B (either is sufficient)

---

### Block 5 — Speckle filter function

```javascript
// ── Speckle filter ───────────────────────────────────────────────
function speckleFilter(image) {
  // Step 1: convert dB to linear power  →  10^(VV/10)
  var linear   = ee.Image(10).pow(image.select('VV').divide(10));
  // Step 2: smooth with a 7×7 neighbourhood average
  var smoothed = linear.focal_mean(SPECKLE_KERNEL, 'square', 'pixels');
  // Step 3: convert back to dB  →  10 × log10(power)
  return smoothed.log10().multiply(10).rename('VV')
    .copyProperties(image, image.propertyNames());
}
```

**What is speckle?**
Radar images have a "salt-and-pepper" noise called speckle. Each pixel records
the interference of many tiny scatterers within its footprint, so a water pixel
can randomly appear bright and a land pixel can appear dark. This makes simple
threshold classification unreliable on raw images.

**Why convert to linear first?**
Backscatter is stored in dB (a logarithm). Averaging logarithms gives the wrong
answer — it is equivalent to computing the geometric mean instead of the arithmetic
mean. The physically correct procedure is:
1. Convert dB → linear power (un-log it)
2. Average in linear space
3. Convert back to dB

`ee.Image(10).pow(image.select('VV').divide(10))` computes `10^(VV/10)`, which
undoes the `10 × log10(power)` formula used to store radar data in dB.

**`.copyProperties(image, image.propertyNames())`**
After all the maths, the returned image would lose its metadata (acquisition date,
orbit number, etc.). This line copies all properties from the original image onto
the processed image so the collection stays identifiable.

---

### Block 6 — Topographic mask

```javascript
// ── Topographic mask ─────────────────────────────────────────────
var floodableMask = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'))
  .lt(5)
  .clip(aoi);
```

**Why mask steep terrain?**
Two radar artefacts on hillsides produce dark (low-backscatter) pixels that
look like water:
- **SAR shadow**: the far side of a hill receives no radar pulse → near-zero signal
- **Layover**: the top of a slope appears geometrically closer than its base

Argon Pedion is surrounded by limestone hills. Without this mask, shadow pixels
on the hillsides pass the −15 dB threshold and are counted as flooded area.

**How the mask works**
- `ee.Image('USGS/SRTMGL1_003')` loads the SRTM 30m global elevation model
- `ee.Terrain.slope()` computes terrain slope in degrees
- `.lt(5)` — "less than 5°" — produces a binary image: **1** where slope < 5° (flat,
  floodable), **0** where slope ≥ 5° (steep, not floodable)
- `.clip(aoi)` trims the mask to the study area boundary

Later, `.and(floodableMask)` multiplies water-detection results by this mask,
zeroing out any steep-slope pixels before area is summed.

---

### Block 7 — Load Sentinel-1 and compute water area per scene

This is the core of the script. It is longer, so read it carefully.

```javascript
// ── Load Sentinel-1 ──────────────────────────────────────────────
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(aoi)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV')
  .filter(floodSeasonFilter)
  .map(speckleFilter)
  .map(function(image) {
    var date = ee.Date(image.get('system:time_start'));

    var areaM2 = ee.Number(
      image.lt(SAR_THRESHOLD)          // 1 where VV < −15 dB (water), 0 elsewhere
           .and(floodableMask)          // remove steep-slope pixels
           .rename('water')             // name the band 'water'
           .clip(aoi)                   // trim to study area
           .multiply(ee.Image.pixelArea())  // each pixel value becomes its area in m²
           .reduceRegion({
             reducer:   ee.Reducer.sum(),   // sum all pixel areas → total flooded m²
             geometry:  aoi,
             scale:     10,                 // Sentinel-1 GRD native resolution (m)
             maxPixels: 1e10
           })
           .get('water')               // extract the numeric result by band name
    );

    return image
      .set('water_area_m2',  areaM2)
      .set('water_area_km2', areaM2.divide(1e6))
      .set('date_str',       date.format('YYYY-MM-dd'))
      .set('year',           date.get('year'))
      .set('month',          date.get('month'));
  });

print('Total S1 flood-season scenes:', s1.size());
```

**Step-by-step breakdown of the filter chain:**

| Line | What it does |
|---|---|
| `'COPERNICUS/S1_GRD'` | GEE's catalogue ID for the Sentinel-1 Ground Range Detected archive |
| `.filterBounds(aoi)` | Keep only images whose footprint overlaps the AOI — discards global data |
| `.filter(eq 'instrumentMode' 'IW')` | IW = Interferometric Wide Swath, the standard land mode for S1 |
| `.filter(listContains ... 'VV')` | Ensure the image has a VV polarisation band (not all modes do) |
| `.select('VV')` | Keep only the VV band, discard others (saves memory) |
| `.filter(floodSeasonFilter)` | Keep only Oct–Apr images (defined in Block 4) |
| `.map(speckleFilter)` | Run the speckle-filter function on every image in the collection |

**The inner `.map()` — computing flooded area per scene:**

`.map()` in GEE works like a loop: it runs your function once for each image
and returns a new collection of the results. For each image the function:

1. `image.lt(SAR_THRESHOLD)` — compares every pixel to −15 dB.
   Output: **1** where the pixel is dark enough to be water, **0** everywhere else.

2. `.and(floodableMask)` — multiplies by the slope mask.
   Output: **1** only where the pixel is both dark AND on flat ground.

3. `.rename('water')` — the result inherits the band name 'VV' from the original
   image. We rename it 'water' so the next line can retrieve the result by name.

4. `.multiply(ee.Image.pixelArea())` — multiplies each pixel's 0-or-1 value by
   the pixel's area in m². Water pixels (value=1) become ≈100 m² (10 m × 10 m).
   Non-water pixels (value=0) become 0 m². Now the image contains areas, not flags.

5. `.reduceRegion({reducer: sum})` — sums all pixel areas within the AOI.
   The result is a dictionary like `{water: 245000}` — 245 000 m² flooded.
   `.get('water')` extracts just the number.

6. `ee.Number(...)` — wraps the result as an explicit GEE Number. This is a
   safeguard: GEE evaluates code lazily, and sorting a collection by an unwrapped
   property can produce incorrect ordering.

7. `.set(...)` — attaches the computed values as metadata on the image.
   `water_area_m2`, `water_area_km2`, `date_str`, `year`, `month` are now
   queryable properties of each image, like columns in a table.

**`print('Total S1 flood-season scenes:', s1.size())`**
Click Run and check this first. A healthy result is around 150–250.
If it is 0, one of your filters has no matching data — most likely `filterBounds`
(check the AOI coordinates) or `floodSeasonFilter`.

---

### Block 8 — Annual maximum statistics

Now that every scene has a `water_area_m2` property, finding the annual peak is
straightforward: filter to one year, sort by area, take the first.

```javascript
// ── Annual max table & chart ─────────────────────────────────────
var annualStats = ee.FeatureCollection(YEARS.map(function(yr) {
  var subset = s1.filter(ee.Filter.eq('year', yr));
  var maxImg  = ee.Image(subset.sort('water_area_m2', false).first());
  return ee.Feature(null, {
    year:      yr,
    water_km2: subset.aggregate_max('water_area_km2'),
    max_date:  maxImg.get('date_str'),
    count:     subset.size()
  });
}));

print('── Annual max flood extent ──');
print(annualStats.filter(ee.Filter.gt('count', 0))
  .reduceColumns(ee.Reducer.toList(3), ['year', 'max_date', 'water_km2']));
```

**`YEARS.map(function(yr) {...})`**
This is a plain JavaScript `.map()` running in your browser (client-side) —
nothing to do with `ee.ImageCollection.map()`. It iterates over the `YEARS`
array locally, building one GEE `ee.Feature` per year. The distinction matters:
plain JS runs immediately in the browser; `ee.*` calls are sent to GEE's servers.

**Inside the function, for each year `yr`:**
- `s1.filter(ee.Filter.eq('year', yr))` — narrow the collection to that year
- `.sort('water_area_m2', false)` — sort descending (largest first)
- `.first()` — grab the top image → the year's peak flood scene
- `subset.aggregate_max('water_area_km2')` — find the maximum km² value across
  all scenes that year (a robust alternative to reading it from the sorted image)
- `ee.Feature(null, {...})` — a GEE Feature with no geometry but with named
  properties: year, area, date, scene count. Think of it as one row in a table.

**`ee.FeatureCollection([...])`**
Wraps the array of Features into a table that GEE can print, chart, and export.

**`print(annualStats.filter(...).reduceColumns(...))`**
`filter(ee.Filter.gt('count', 0))` drops years with zero S1 scenes.
`reduceColumns(ee.Reducer.toList(3), ['year','max_date','water_km2'])` extracts
those three columns as a list. Click the expand arrow in the Console to see
the full table.

---

### Block 9 — Bar chart

```javascript
print(
  ui.Chart.feature.byFeature(
    annualStats.filter(ee.Filter.gt('count', 0)), 'year', 'water_km2'
  )
  .setChartType('ColumnChart')
  .setOptions({
    title:  'Argon Pedion — Annual Maximum Flood Extent (Sentinel-1 SAR)',
    hAxis:  {title: 'Year', format: '####', gridlines: {count: YEARS.length}},
    vAxis:  {title: 'Max Flooded Area (km²)', minValue: 0},
    colors: ['#1565c0'],
    legend: {position: 'none'}
  })
);
```

`ui.Chart.feature.byFeature(collection, xProperty, yProperty)` produces a chart
with one bar per feature, x-axis = `year`, y-axis = `water_km2`.
Click the chart thumbnail in the Console to open it full-size; use the download
icon there to save as PNG or CSV.

---

### Block 10 — Map setup and dry-season basemap

```javascript
// ── Map ──────────────────────────────────────────────────────────
Map.centerObject(aoi, 13);

var s2Dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .median()
  .clip(aoi);
Map.addLayer(s2Dry, {bands: ['B4','B3','B2'], min: 0, max: 2500}, 'S2 dry basemap');
Map.addLayer(aoi, {color: 'FFFFFF'}, 'AOI');
```

**`Map.centerObject(aoi, 13)`**
Zooms and pans the map so the AOI fills the view. `13` is the zoom level
(1 = whole world, 20 = individual buildings). Adjust up or down to taste.

**Why Sentinel-2 as the basemap?**
Sentinel-2 is an optical (camera) satellite — it produces natural-looking
true-colour images, much easier to interpret as a background than radar grey.
We use summer scenes (Jun–Sep) when the plain is dry, to contrast with the flood
layers on top.

**`ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15)`**
Discards scenes where more than 15% of pixels are flagged as cloud.
This is a coarse pre-filter; `.median()` then removes any remaining cloud
pixels (they are bright outliers; the median is robust to outliers).

**`.median()`**
Takes the pixel-wise median across all qualifying summer scenes. This is a
standard cloud-compositing technique: clouds appear in different locations on
different dates, so the median value at each pixel usually comes from a clear day.

**`{bands: ['B4','B3','B2'], min: 0, max: 2500}`**
B4=Red, B3=Green, B2=Blue → true colour. `min`/`max` define the stretch: values
below 0 appear black, above 2500 appear white. Surface Reflectance values for a
summer Mediterranean landscape typically fall in this range.

---

### Block 11 — Annual flood layers

```javascript
for (var i = 0; i < YEARS.length; i++) {
  var yr    = YEARS[i];
  var color = COLORS[i];

  var yearMax = ee.Image(
    s1.filter(ee.Filter.eq('year', yr))
      .sort('water_area_m2', false)
      .first()
  );

  var water = yearMax
    .lt(SAR_THRESHOLD)
    .and(floodableMask)
    .selfMask()
    .clip(aoi);

  Map.addLayer(water, {palette: [color]}, yr + ' MAX', true, 0.6);
}
```

**What the loop does**
Iterates over each year, retrieves its peak-flood scene, creates a binary water
image, and adds it as a named map layer.

**`.selfMask()`**
After thresholding and slope masking, the image contains **1** (water) and **0**
(not water). `.selfMask()` masks out all 0-value pixels, making them transparent.
Without this, non-water pixels would render as a solid colour and cover the basemap.

**`Map.addLayer(image, visParams, name, shown, opacity)`**
- `image` — what to display
- `{palette: [color]}` — colourise 1-values with this hex colour
- `yr + ' MAX'` — layer name shown in the Layers panel
- `true` — layer is visible by default (set to `false` to start hidden)
- `0.6` — 60% opacity; overlapping years' colours blend together

---

### Block 12 — Legend

```javascript
// ── Legend ───────────────────────────────────────────────────────
var legend = ui.Panel({
  style: {position: 'bottom-right', padding: '8px',
          backgroundColor: 'rgba(255,255,255,0.92)'}
});
legend.add(ui.Label('Annual Max Flood Extent',
  {fontWeight: 'bold', fontSize: '12px', margin: '0 0 5px 0'}));

for (var j = 0; j < YEARS.length; j++) {
  legend.add(ui.Panel([
    ui.Label('', {backgroundColor: COLORS[j],
                  padding: '5px 12px', margin: '2px 5px 2px 0'}),
    ui.Label(YEARS[j] + ' MAX', {fontSize: '11px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal')));
}
Map.add(legend);
```

GEE provides a simple UI library (`ui.Panel`, `ui.Label`) for adding widgets to
the map. A `ui.Panel` is a rectangular container. `ui.Label` creates a text
element — here we abuse it with no text and a coloured background to make a
colour swatch. `ui.Panel.Layout.flow('horizontal')` arranges swatch and text
side by side in a row.

---

### Block 13 — Flood-frequency layer

```javascript
// ── Flood frequency ──────────────────────────────────────────────
var floodFreq = s1.map(function(img) {
  return img.lt(SAR_THRESHOLD).and(floodableMask).rename('water');
}).sum()
  .divide(s1.size())
  .multiply(100)
  .rename('flood_pct')
  .clip(aoi);

Map.addLayer(floodFreq,
  {min: 0, max: 50, palette: ['#ffffff','#add8e6','#4fc3f7','#1565c0','#003087']},
  'Flood frequency % (all years)', false);
```

**What this shows**
For each pixel: out of all Sentinel-1 acquisitions in the flood season (all years
combined), what percentage showed that pixel as water?

- **0%** — never detected as water
- **50%** — water in about half of all acquisitions
- **> 50%** — frequently or permanently inundated

**How it is computed**
1. `.map(...)` — threshold and mask every scene → a collection of 0/1 images
2. `.sum()` — pixel-wise sum across the collection → count of wet scenes per pixel
3. `.divide(s1.size())` — divide by total scene count → fraction (0.0–1.0)
4. `.multiply(100)` — convert to percentage

The layer is off by default (`false` at the end of `Map.addLayer`).
Turn it on in the Layers panel to see which part of the plain forms the
permanent core versus the seasonal fringe.

---

## Part C — Run the complete script

Click **Run**. Expect:

1. The map centres on Argon Pedion immediately
2. The Console shows scene count and the table (may take 15–30 s)
3. Map layers appear as GEE finishes computing each one
4. The chart appears in the Console as a thumbnail — click it to expand

> GEE uses *lazy evaluation*: it only computes a layer when you pan or zoom
> to it, not all at once. Layers may appear blank for a moment before loading.

---

## Part D — Experiments to try

| Change | Effect |
|---|---|
| `SAR_THRESHOLD = -18` | Stricter — smaller flood extents, fewer false positives |
| `SAR_THRESHOLD = -12` | Looser — larger extents, more false positives from moist soil |
| `SPECKLE_KERNEL = 3` | Less smoothing — noisier water boundary |
| `SPECKLE_KERNEL = 11` | More smoothing — cleaner boundary, blurs narrow channels |
| Add `2026` to `YEARS` | Include the current season (partial data) |
| Flood-frequency `max: 25` | Stretch the colour scale — shows low-frequency flooding |

After each change, click **Run** to recompute.

---

## Part E — Export results to CSV

Add this block at the end of the script, then click **Run**:

```javascript
Export.table.toDrive({
  collection:     annualStats,
  description:    'argon_pedion_s1_annual_max',
  folder:         'GEE_exports',
  fileNamePrefix: 'argon_pedion_s1_annual_max',
  fileFormat:     'CSV',
  selectors:      ['year', 'max_date', 'water_km2', 'count']
});
```

Go to the **Tasks** tab (right panel) → click **RUN** next to the task.
The CSV appears in `My Drive/GEE_exports/` when the job finishes (1–5 minutes).

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Scene count = 0 | Wrong AOI coords or filter | Check the polygon coordinates; print `s1.size()` after each filter |
| "List index out of range" | A year has zero scenes | Add `.filter(ee.Filter.gt('count', 0))` before printing |
| Map layers never load | GEE servers busy | Wait, then zoom/pan to trigger rendering |
| Flood extent looks huge | Threshold too loose | Raise SAR_THRESHOLD to −18 |
| Hills misclassified as water | Slope mask not applied | Check the `.and(floodableMask)` line is present |
| `null` in the Console | Year has no S1 coverage | Expected for years before 2014 |

---

## Key concepts

| Term | Meaning |
|---|---|
| SAR | Synthetic Aperture Radar — active sensor, works through clouds |
| Backscatter | Fraction of radar energy reflected back to the satellite |
| dB | Decibels — log scale; open water ≈ −20 to −15 dB |
| Speckle | Salt-and-pepper radar noise; removed with a neighbourhood average |
| AOI | Area of Interest — the polygon bounding your study area |
| `ee.ImageCollection` | GEE's catalogue of images from one sensor |
| `.filter()` | Selects images matching a condition |
| `.map()` | Runs a function on every image in a collection |
| `.reduceRegion()` | Computes a statistic (e.g. sum) over all pixels in a region |
| `.selfMask()` | Makes 0-value pixels transparent on the map |
| Flood frequency | Percentage of acquisitions where a pixel was classified as water |

---

## What's next

| Level | Platform | Script |
|---|---|---|
| Level 1 S2 | GEE Code Editor | `gee_s2_argon_pedion.js` — same analysis with Sentinel-2 optical |
| Level 2 | GEE Python API | Reproduce S1 annual max using `geemap` in a Jupyter notebook |
| Level 3 | openEO + CDSE | Reproduce S1 annual max using the openEO Python client |
