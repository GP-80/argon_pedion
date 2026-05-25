// ================================================================
// Argon Pedion — Annual Maximum Flood Extent
// Sentinel-1 SAR (VV, IW, GRD) · 2014–2025
// ================================================================
//
// WHAT THIS SCRIPT DOES
// ─────────────────────
// For each year 2014–2025 this script:
//   1. Collects all Sentinel-1 radar images taken during the flood
//      season (October–April) over the Argon Pedion plain.
//   2. Identifies the single scene with the largest flooded area
//      for that year.
//   3. Displays each year's peak flood extent as a coloured map layer.
//   4. Prints a summary table and bar chart of flooded area per year.
//   5. Adds a flood-frequency layer showing how often each pixel
//      was flooded across all years combined.
//
// WHY SENTINEL-1 (RADAR / SAR)?
// ─────────────────────────────
// Sentinel-1 is a radar satellite — it transmits its own microwave
// pulses and records the reflected signal (backscatter). Unlike an
// optical camera, it does not depend on sunlight and its signal
// passes through clouds. This matters enormously for flood mapping:
// floods happen during rainy, cloudy weather, exactly when optical
// satellites (e.g. Sentinel-2) cannot see the ground.
//
// HOW RADAR DETECTS OPEN WATER
// ────────────────────────────
// Flat water acts like a mirror for radar: almost all energy reflects
// away from the satellite → the sensor records a very low (dark) value.
// Rough surfaces (soil, crops, buildings) scatter energy back strongly
// → bright values. We exploit this contrast with a simple rule:
//   pixels with VV backscatter below –15 dB  →  water
//   pixels at or above –15 dB               →  land
//
// WHAT IS dB?
// ───────────
// Backscatter is measured in decibels (dB), a logarithmic unit.
// –15 dB means the returned power is 10^(–15/10) ≈ 0.032 of the
// transmitted power — i.e. 97% of the energy was NOT reflected back.
// Water is typically –20 to –15 dB; vegetated land is –10 to –5 dB.
//
// LANGUAGE NOTE
// ─────────────
// This script uses the Earth Engine JavaScript API (EEJS).
// The syntax (var, function, for-loops, // comments) is JavaScript,
// but the ee.* objects are GEE's own server-side API. Code after the
// ee. prefix describes a computation; GEE's servers execute it remotely
// when the result is needed (lazy evaluation).
//
// SCRIPT STRUCTURE
// ────────────────
//   §1  Area of Interest (AOI polygon)
//   §2  Parameters (threshold, kernel size, years, colours)
//   §3  Flood-season filter (Oct–Apr)
//   §4  Speckle filter function (radar noise reduction)
//   §5  Topographic mask (exclude steep hillsides)
//   §6  Load S1, compute water area per scene, build annual stats
//   §7  Annual max table & bar chart
//   §8  Map layers (dry basemap + annual max flood extents)
//   §9  Map legend
//   §10 Flood-frequency layer (fraction of scenes where pixel = water)
// ================================================================


// ── 1. Area of Interest ─────────────────────────────────────────
// ee.Geometry.Polygon() defines the study area boundary as a closed
// polygon. GEE uses WGS84 (EPSG:4326) — standard GPS decimal degrees
// in [longitude, latitude] order (note: longitude FIRST, then latitude).
//
// This polygon was traced manually around the visible floodplain in the
// GEE map editor. Every calculation below is clipped to this boundary,
// so nothing runs on the whole planet.
//
// The null and false at the end are optional arguments for projection
// and geodesic mode — the defaults are fine for a small study area.
var aoi = ee.Geometry.Polygon(
  [[[22.457559040088586, 37.66680064822952],
    [22.461397929816396, 37.66641575289452],
    [22.459810479896845, 37.66378477768644],
    [22.46062134638269,  37.66162015205006],
    [22.46141061800174,  37.660499340515784],
    [22.461676126099565, 37.658240581461314],
    [22.459624232294654, 37.65690867188218],
    [22.462616015270022, 37.65512052712688],
    [22.466118896365803, 37.65507912606556],
    [22.46310744567408,  37.652516672358125],
    [22.461710303768218, 37.65128033749649],
    [22.463821225515364, 37.649741692833224],
    [22.464647927653346, 37.64679415969568],
    [22.462404397805575, 37.64415467235887],
    [22.460073186168394, 37.6404698924646],
    [22.459899310180518, 37.63829240205909],
    [22.461177863513583, 37.63730848606322],
    [22.46448612066055,  37.63339269001365],
    [22.466628643866706, 37.63028823499969],
    [22.46742283811807,  37.6289064107416],
    [22.464414156762285, 37.62625701576659],
    [22.465917429223182, 37.62497143710055],
    [22.463748857009442, 37.621375846280586],
    [22.464876442727242, 37.61973735010361],
    [22.466428889010377, 37.618713429657696],
    [22.466441397675492, 37.61810452803893],
    [22.464810018827247, 37.61764826616711],
    [22.461516848769005, 37.6182144338599],
    [22.45988011636385,  37.61801906155693],
    [22.457603443595588, 37.61703222057105],
    [22.45743326496716,  37.61468074966045],
    [22.457476395478373, 37.61259309407664],
    [22.457116150779264, 37.611500769358294],
    [22.456203805362744, 37.61066224859232],
    [22.454484916183645, 37.60916059033295],
    [22.452871830778207, 37.607834312974354],
    [22.451683527694826, 37.60712265271683],
    [22.449323244416842, 37.604916428606764],
    [22.44631873242996,  37.60740017085614],
    [22.4499820722686,   37.6100586247031],
    [22.45371655931099,  37.611934816308555],
    [22.45237557010393,  37.61330938037847],
    [22.449769792305567, 37.615015353320565],
    [22.44828091902933,  37.6156048446798],
    [22.44436243705835,  37.61468315476615],
    [22.44034561101157,  37.61323798162851],
    [22.434083743128888, 37.611936916432306],
    [22.432562054793983, 37.61148166666719],
    [22.43075697462862,  37.6115011775964],
    [22.428229261386157, 37.61207659426142],
    [22.433431167640897, 37.61697476876623],
    [22.433730022559402, 37.61837087793408],
    [22.43283328050836,  37.619403075564016],
    [22.43416190191555,  37.62124790136044],
    [22.43319584253664,  37.62297525110243],
    [22.437232202958675, 37.62611759890499],
    [22.439647452965808, 37.62832478754249],
    [22.439482603290404, 37.6309763734186],
    [22.43975524367894,  37.6336337450014],
    [22.44211176985768,  37.63866788804759],
    [22.446329163056745, 37.643683003985956],
    [22.45138155271777,  37.64540260860606],
    [22.456634973000334, 37.647994715571016],
    [22.456805165764695, 37.650346175866865],
    [22.45617530883862,  37.654340328861515],
    [22.456973152739295, 37.658092137221246],
    [22.456363041328437, 37.661129450964175],
    [22.454251639523036, 37.66266794871616],
    [22.452473834675537, 37.66394976022294],
    [22.457559040088586, 37.66680064822952]]],
  null, false
);


// ── 2. Parameters ───────────────────────────────────────────────
// All tunable settings live here so you can experiment without
// hunting through the rest of the code.

// Water threshold in dB.
// Pixels with VV backscatter BELOW this value are classified as water.
// –15 dB is a standard value for open water in C-band radar (5.4 GHz).
// If dry bare soil or ploughed fields are being wrongly detected as water,
// lower this to –18 dB to make the classifier stricter (fewer detections).
var SAR_THRESHOLD  = -15;

// Speckle-filter kernel size in pixels.
// A 7×7 window means each pixel is averaged with the 48 pixels around it.
// Increase to 9 or 11 for smoother results; decrease to 5 for more detail.
var SPECKLE_KERNEL = 7;

// Years to analyse. Sentinel-1 data in GEE starts April 2014.
var YEARS  = [2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025];

// One distinct colour per year for the map layers (ColorBrewer Paired palette).
// At 60 % opacity the colours blend where extents overlap, revealing which
// part of the plain floods every year vs. only in exceptional years.
var COLORS = ['#a6cee3','#1f78b4','#b2df8a','#33a02c',
              '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
              '#cab2d6','#6a3d9a','#ffff99','#b15928'];


// ── 3. Flood-season filter (Oct–Apr) ────────────────────────────
// Argon Pedion floods during winter rainfall — October through April.
// We cannot use a simple date range because the season straddles the
// calendar year boundary: Oct–Dec belong to one year, Jan–Apr to the next.
//
// ee.Filter.or() combines two separate month-range filters with OR logic:
//   calendarRange(10, 12, 'month')  →  October, November, December
//   calendarRange( 1,  4, 'month')  →  January, February, March, April
// Any image whose acquisition month falls in either range passes through.
var floodSeasonFilter = ee.Filter.or(
  ee.Filter.calendarRange(10, 12, 'month'),  // Oct–Dec
  ee.Filter.calendarRange(1,  4,  'month')   // Jan–Apr
);


// ── 4. Speckle filter (linear-space focal mean) ──────────────────
// WHAT IS SPECKLE?
//   Radar images have a characteristic "salt-and-pepper" noise called
//   speckle. Each pixel records the combined interference of many tiny
//   scatterers within its footprint. A single water pixel can appear
//   randomly bright, and a land pixel can appear randomly dark — even
//   when the surrounding area is clearly one type. Without smoothing,
//   simple threshold classification produces a noisy, dotted result.
//
// THE APPROACH — focal mean in LINEAR power space:
//   We replace each pixel's value with the average of its 7×7 neighbourhood.
//   Crucially, this averaging must be done in LINEAR power, NOT in dB.
//   Why? Because dB is a logarithmic scale: the average of log values ≠
//   the log of the average. Averaging in dB systematically underestimates
//   backscatter. The correct steps are:
//
//     Step 1 — convert dB → linear power:
//              power = 10 ^ (dB / 10)
//              → ee.Image(10).pow(image.select('VV').divide(10))
//
//     Step 2 — apply 7×7 square neighbourhood average:
//              → .focal_mean(7, 'square', 'pixels')
//
//     Step 3 — convert linear power → dB:
//              dB = 10 × log10(power)
//              → .log10().multiply(10)
//
//     Step 4 — rename the band back to 'VV' and copy the original
//              image metadata (date, orbit number, etc.) so the
//              smoothed image is still identifiable.
function speckleFilter(image) {
  var linear   = ee.Image(10).pow(image.select('VV').divide(10)); // dB → linear
  var smoothed = linear.focal_mean(SPECKLE_KERNEL, 'square', 'pixels'); // smooth
  return smoothed.log10().multiply(10).rename('VV')    // linear → dB, rename
    .copyProperties(image, image.propertyNames());     // keep original metadata
}


// ── 5. Topographic mask ──────────────────────────────────────────
// WHY MASK STEEP TERRAIN?
//   Steep hillsides cause two radar artefacts that mimic water:
//
//   • SAR shadow: the far side of a steep hill receives no radar pulse at
//     all (it's in the "shadow" of the hill). The sensor records near-zero
//     backscatter — very dark, exactly like open water. Without masking,
//     these shadow pixels pass the –15 dB threshold and are counted as flood.
//
//   • Layover: the top of a slope appears geometrically closer than its base,
//     creating distorted, unreliable values.
//
//   Argon Pedion is surrounded by limestone hills. This mask removes all
//   pixels on slopes steeper than 5°, keeping only the flat, genuinely
//   floodable valley floor.
//
// HOW IT WORKS:
//   ee.Image('USGS/SRTMGL1_003') loads the SRTM 30m elevation model.
//   ee.Terrain.slope() computes slope in degrees for every pixel.
//   .lt(5) creates a binary image: 1 where slope < 5°, 0 where slope ≥ 5°.
//   .clip(aoi) trims the mask to the study area boundary.
//   Later, .and(floodableMask) multiplies any water-detection result by this
//   mask, zeroing out steep-slope pixels before area is summed.
var floodableMask = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'))
  .lt(5)       // 1 = floodable (flat), 0 = non-floodable (steep)
  .clip(aoi);


// ── 6. Load Sentinel-1 and compute water area per scene ─────────
// GEE stores satellite data as ImageCollections — catalogues of
// thousands of individual images. We filter this catalogue down to
// exactly the images we need, step by step:
//
//   .filterBounds(aoi)
//       Keep only images that spatially overlap our study area.
//       Without this, GEE would process global coverage.
//
//   .filter(ee.Filter.eq('instrumentMode', 'IW'))
//       IW = Interferometric Wide Swath — the standard S1 mode over land.
//       Other modes (EW, SM, WV) exist for ice/ocean; we don't want them.
//
//   .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
//       Each S1 image contains polarisation channels (VV, VH, or both).
//       VV (Vertical–Vertical) is most sensitive to open water surfaces.
//       This filter ensures every image in the collection has a VV band.
//
//   .select('VV')
//       Discard all bands except VV to reduce memory usage.
//
//   .filter(floodSeasonFilter)
//       Keep only images from October–April (flood season, defined in §3).
//
//   .map(speckleFilter)
//       Apply our speckle-smoothing function (§4) to every image.
//       .map() in GEE is like a for-loop: it runs the function once for
//       each image in the collection and returns a new collection of the
//       processed images.
//
// The second .map() adds flood-area statistics to each image as metadata
// (called "properties" in GEE). Here is what happens for each image:
//
//   image.lt(SAR_THRESHOLD)
//       Creates a binary (0/1) image: 1 where VV < –15 dB (likely water).
//
//   .and(floodableMask)
//       Multiplies by the slope mask (§5): removes steep-slope pixels.
//       After this step, 1 = water on flat ground, 0 = everything else.
//
//   .rename('water')
//       Renames the band from 'VV' to 'water' so the next step can
//       retrieve the result by the name 'water'.
//
//   .multiply(ee.Image.pixelArea())
//       Multiplies each pixel (value 0 or 1) by its area in m².
//       Sentinel-1 GRD is processed at 10 m resolution → each pixel ≈ 100 m².
//       Water pixels become ~100 m²; non-water pixels become 0.
//
//   .reduceRegion({reducer: ee.Reducer.sum(), ...})
//       Sums all pixel areas within the AOI → total flooded area in m².
//       Returns a dictionary; .get('water') extracts the numeric value.
//
//   ee.Number(...)
//       Wraps the result as an explicit GEE Number so sorting later
//       works reliably (GEE lazy-evaluation safeguard).
//
//   .set(...)
//       Attaches the computed values as metadata properties to the image
//       so we can sort and filter the collection by flooded area later.
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(aoi)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV')
  .filter(floodSeasonFilter)
  .map(speckleFilter)
  .map(function(image) {
    var date   = ee.Date(image.get('system:time_start')); // acquisition datetime
    var areaM2 = ee.Number(
      image.lt(SAR_THRESHOLD).and(floodableMask).rename('water').clip(aoi)
        .multiply(ee.Image.pixelArea())        // pixel value → pixel area in m²
        .reduceRegion({
          reducer:   ee.Reducer.sum(),          // sum all pixel areas
          geometry:  aoi,
          scale:     10,                        // process at 10 m (native S1 res)
          maxPixels: 1e10                       // raise pixel limit for large areas
        })
        .get('water')                          // extract the 'water' band total
    );
    return image
      .set('water_area_m2',  areaM2)                    // flooded area in m²
      .set('water_area_km2', areaM2.divide(1e6))        // same in km²
      .set('date_str',       date.format('YYYY-MM-dd')) // date as text
      .set('year',           date.get('year'))           // numeric year
      .set('month',          date.get('month'));          // numeric month
  });

// Print to Console as a quick sanity check — if this is 0 a filter is wrong.
print('Total S1 flood-season scenes:', s1.size());


// ── 7. Annual max table & chart ──────────────────────────────────
// For each year we want: which scene had the largest flood, on what date,
// and how large was it?
//
// YEARS.map(function(yr) {...}) iterates over our list of years in
// JavaScript (client-side), running the inner function once per year.
// Inside, for each year:
//
//   s1.filter(ee.Filter.eq('year', yr))
//       Narrow the collection to flood-season scenes from that year only.
//
//   subset.sort('water_area_m2', false).first()
//       Sort descending by flooded area (false = largest first),
//       then take the first image → the annual maximum-flood scene.
//
//   subset.aggregate_max('water_area_km2')
//       Alternative way to get the maximum value across all scenes —
//       used here so the chart value comes from the collection aggregate
//       rather than a single scene property (more robust).
//
//   ee.Feature(null, {...})
//       Creates a GEE Feature with no geometry but with attribute table
//       columns: year, water_km2, max_date, count.
//       Wrapping everything in ee.FeatureCollection() gives us a table
//       we can print, chart, and export.
var annualStats = ee.FeatureCollection(YEARS.map(function(yr) {
  var subset = s1.filter(ee.Filter.eq('year', yr));
  var maxImg  = ee.Image(subset.sort('water_area_m2', false).first());
  return ee.Feature(null, {
    year:      yr,
    water_km2: subset.aggregate_max('water_area_km2'),  // km² of peak flood
    max_date:  maxImg.get('date_str'),                   // date of peak scene
    count:     subset.size()                             // number of S1 scenes
  });
}));

// Print the table to the Console (years with no scenes are filtered out).
// Click the list icon next to the result to expand the full table.
print('── Annual max flood extent ──');
print(annualStats.filter(ee.Filter.gt('count', 0))
  .reduceColumns(ee.Reducer.toList(3), ['year', 'max_date', 'water_km2']));

// Bar chart — one bar per year, height = peak flooded area in km².
// Click the chart in the Console to open it in a larger window where
// you can hover over bars to see exact values and download as CSV/PNG.
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


// ── 8. Map ───────────────────────────────────────────────────────
// Map.centerObject() zooms the map to the AOI. Zoom level 13 shows
// the plain at roughly 1:10 000 scale — adjust up (14, 15) to zoom in.
Map.centerObject(aoi, 13);

// DRY-SEASON TRUE-COLOUR BASEMAP (Sentinel-2)
// A summer (Jun–Sep) Sentinel-2 median composite serves as the background.
// Using Sentinel-2 (optical) for the basemap because it produces a natural-
// looking true-colour image — much easier to interpret than a grey radar image.
// .median() takes the pixel-wise median across all summer scenes, effectively
// removing remaining clouds (cloud pixels are high outliers, median ignores them).
// Bands B4=Red, B3=Green, B2=Blue → standard true colour.
// min/max of 0–2500 is a typical stretch for S2 Surface Reflectance (SR).
var s2Dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15)) // pre-filter cloudy scenes
  .filter(ee.Filter.calendarRange(6, 9, 'month'))       // summer only
  .median().clip(aoi);
Map.addLayer(s2Dry, {bands: ['B4','B3','B2'], min: 0, max: 2500}, 'S2 dry basemap');
Map.addLayer(aoi,   {color: 'FFFFFF'}, 'AOI'); // white outline of the study area

// ANNUAL MAXIMUM FLOOD LAYERS — one coloured layer per year.
// The loop adds 12 layers (one per year) to the map.
// For each year:
//   • Find the scene with the largest water area (.sort descending, .first())
//   • Threshold at SAR_THRESHOLD → binary water/non-water image
//   • .and(floodableMask) → zero out steep hillside pixels
//   • .selfMask() → make 0-value pixels transparent (only water pixels render)
//   • Map.addLayer(..., 0.6) → display at 60% opacity
//
// At 60% opacity the colours blend wherever multiple years overlap, showing
// the permanently or near-permanently flooded core of the plain in a mixed hue.
// Toggle individual years on/off in the Layers panel (top-right of the map).
for (var i = 0; i < YEARS.length; i++) {
  var yr    = YEARS[i];
  var color = COLORS[i];
  var yearMax = ee.Image(
    s1.filter(ee.Filter.eq('year', yr)).sort('water_area_m2', false).first()
  );
  var water = yearMax.lt(SAR_THRESHOLD).and(floodableMask).selfMask().clip(aoi);
  Map.addLayer(water, {palette: [color]}, yr + ' MAX', true, 0.6);
}


// ── 9. Legend ────────────────────────────────────────────────────
// GEE's UI library lets you build custom HTML-style panels on the map.
// ui.Panel() creates a container; ui.Label() creates a text element.
// The outer loop adds one colour swatch + label row per year.
// Map.add() places the finished panel on the map canvas.
var legend = ui.Panel({
  style: {
    position:        'bottom-right',       // anchor to bottom-right corner
    padding:         '8px',
    backgroundColor: 'rgba(255,255,255,0.92)' // semi-transparent white
  }
});

legend.add(ui.Label('Annual Max Flood Extent', {
  fontWeight: 'bold', fontSize: '12px', margin: '0 0 5px 0'
}));

for (var j = 0; j < YEARS.length; j++) {
  legend.add(ui.Panel([
    ui.Label('', {                          // empty label styled as colour swatch
      backgroundColor: COLORS[j],
      padding:         '5px 12px',
      margin:          '2px 5px 2px 0'
    }),
    ui.Label(YEARS[j] + ' MAX', {fontSize: '11px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal')));  // arrange swatch + text side by side
}

Map.add(legend);


// ── 10. Flood frequency (all flood-season scenes) ─────────────────
// This layer answers: "what fraction of all flood-season radar acquisitions
// showed a pixel as water?" — a proxy for how often and how persistently
// each part of the plain floods.
//
// HOW IT IS COMPUTED:
//   For every scene in the collection, threshold at SAR_THRESHOLD and apply
//   the slope mask → produces a 0/1 image per scene (1 = water that day).
//   .sum() adds all these 0/1 images pixel by pixel → count of wet scenes.
//   .divide(s1.size()) divides by total number of scenes → fraction (0–1).
//   .multiply(100) converts to percentage (0–100%).
//
// HOW TO READ THE MAP:
//   Dark blue  → flooded in nearly every scene → permanent or near-permanent water
//   Light blue → flooded occasionally → seasonal margins
//   White      → rarely or never flooded
//
// The layer is off by default (false at the end of Map.addLayer).
// Turn it on in the Layers panel to compare with the annual max layers above.
var floodFreq = s1.map(function(img) {
  return img.lt(SAR_THRESHOLD).and(floodableMask).rename('water'); // 1=wet, 0=dry
}).sum()                        // count of scenes where pixel was wet
  .divide(s1.size())            // ÷ total scenes = fraction
  .multiply(100)                // → percentage
  .rename('flood_pct')
  .clip(aoi);

Map.addLayer(floodFreq,
  {min: 0, max: 50, palette: ['#ffffff','#add8e6','#4fc3f7','#1565c0','#003087']},
  'Flood frequency % (all years)', false); // false = layer off by default


