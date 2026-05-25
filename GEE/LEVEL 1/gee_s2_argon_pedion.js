// ================================================================
// Argon Pedion — Annual Maximum Flood Extent (Script 2)
// Sentinel-2 Monthly Max-MNDWI Composites · 2017–2026
// ================================================================
// Water index: MNDWI = (B3 Green − B11 SWIR1) / (B3 + B11)
//   MNDWI is the improved NDWI for water mapping: replacing NIR
//   with SWIR1 suppresses false positives from vegetation and
//   built-up areas, which is important at flood plain edges.
//   NDWI = (B3−B8)/(B3+B8) — swap INDEX parameter below if needed.
//
// Compositing strategy:
//   For each month, take the PIXEL-WISE MAX of all cloud-free
//   MNDWI values. If the flood persisted for even a few days,
//   at least one Sentinel-2 pass (5-day revisit with A+B) will
//   have been cloud-free over that pixel, and the max composite
//   captures it. Remaining cloud-masked pixels are filled with -1
//   (permanently treated as non-water).
// ================================================================


// ── 1. Area of Interest ─────────────────────────────────────────
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
var INDEX           = 'MNDWI';  // 'MNDWI' (recommended) or 'NDWI'
var WATER_THRESHOLD = 0.0;      // fallback threshold if Otsu fails
var THRESH_MIN      = -0.1;     // Otsu result clamped to this lower bound
var THRESH_MAX      =  0.3;     // Otsu result clamped to this upper bound
                                // (sensitivity range — Sharma & Saharia 2025)
var CS_THRESHOLD    = 0.65;     // Cloud Score Plus clarity score: 1=clear, 0=opaque
                                // replaces QA60 — catches cloud shadows that QA60 misses

// S2 SR data starts 2017-03-28 → first usable flood season Oct 2017
// 2017 has Oct–Dec only (no Jan–Apr S2 data); 2026 has Jan–Apr only
var YEARS  = [2017,2018,2019,2020,2021,2022,2023,2024,2025,2026];
var COLORS = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
              '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'];

var FLOOD_MONTHS = [10, 11, 12, 1, 2, 3, 4];


// ── 3. Topographic mask ──────────────────────────────────────────
// Pixels with slope > 5° are excluded as non-floodable.
// Removes SAR/optical artefacts on surrounding limestone hillsides.
// Source: Sharma & Saharia 2025.
var floodableMask = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'))
  .lt(5).clip(aoi);


// ── 4. Cloud masking ─────────────────────────────────────────────
// QA60 kept for the dry-season basemap only (simple, no join needed).
function maskS2Clouds(image) {
  var qa = image.select('QA60');
  return image.updateMask(
    qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0))
  );
}

// Cloud Score Plus masking for flood composites — catches cloud shadows
// that QA60 misses (shadow pixels have low reflectance, similar to water,
// causing false MNDWI positives). Source: Sharma & Saharia 2025.
function applyCSPlus(image) {
  return image.updateMask(
    ee.Image(image.get('cs_plus')).select('cs').gte(CS_THRESHOLD)
  );
}


// ── 5. Otsu automatic threshold ──────────────────────────────────
// Finds the WI value that maximises between-class variance (water vs land).
// Adapts to inter-annual variation in turbidity and soil moisture that
// shifts the water/land spectral boundary. Source: Kordelas et al. 2018.
function otsu(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means  = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size   = means.length().get([0]);
  var total  = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum    = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean   = sum.divide(total);
  var bss = ee.List.sequence(1, size).map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount  = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMean   = means.slice(0, 0, i).multiply(aCounts)
                    .reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount);
    var bCount  = total.subtract(aCount);
    var bMean   = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2))
           .add(bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  return ee.Number(means.sort(bss).get([-1]));
}


// ── 6. Water index ───────────────────────────────────────────────
function addWaterIndex(image) {
  var wi = (INDEX === 'NDWI')
    ? image.normalizedDifference(['B3', 'B8'])
    : image.normalizedDifference(['B3', 'B11']);
  return image.addBands(wi.rename('WI'));
}


// ── 7. Monthly max composite for one year-month ──────────────────
// CS+ join masks cloud shadows that QA60 misses (shadows mimic water in WI).
// Otsu threshold adapts to inter-annual turbidity/soil-moisture variation.
// Slope mask excludes non-floodable terrain. Water band stored for map reuse.
function monthlyMaxComposite(yr, m) {
  var start = ee.Date.fromYMD(yr, m, 1);
  var end   = start.advance(1, 'month');

  var s2col  = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi).filterDate(start, end);
  var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    .filterBounds(aoi).filterDate(start, end);

  var joined = ee.Join.saveFirst('cs_plus').apply({
    primary:   s2col,
    secondary: csPlus,
    condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
  });

  var masked     = ee.ImageCollection(joined).map(applyCSPlus).map(addWaterIndex);
  var sceneCount = masked.size();

  // Pixel-wise max WI (fallback -1 guards against empty-month crash)
  var wiMax = masked.select('WI')
    .merge(ee.ImageCollection([ee.Image.constant(-1).rename('WI')]))
    .max().clip(aoi);

  // Otsu threshold from real pixels only; clamp and fall back if no scenes
  var histogram = masked.select('WI').max().clip(aoi)
    .reduceRegion({
      reducer: ee.Reducer.histogram({maxBuckets: 200, minBucketWidth: 0.01}),
      geometry: aoi, scale: 20, maxPixels: 1e10
    }).get('WI');

  var thresh = ee.Number(ee.Algorithms.If(
    sceneCount.gt(0),
    ee.Number(otsu(histogram)).max(THRESH_MIN).min(THRESH_MAX),
    WATER_THRESHOLD
  ));

  var waterMask = wiMax.gte(thresh).and(floodableMask);

  // ee.Number() wrap avoids lazy-eval sort bug in batch exports
  var areaM2 = ee.Number(
    waterMask.rename('water')
      .multiply(ee.Image.pixelArea())
      .reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10})
      .get('water')
  );

  return wiMax.addBands(waterMask.rename('water'))
    .set('year',          yr)
    .set('month',         m)
    .set('date_str',      start.format('YYYY-MM'))
    .set('scene_count',   sceneCount)
    .set('water_area_m2', areaM2)
    .set('water_area_km2', areaM2.divide(1e6));
}


// ── 8. Build annual max composites ───────────────────────────────
// For each year: build 7 monthly composites, pick the one with the
// largest water area. That composite IS the annual max flood scene.
var yearlyMaxImages = YEARS.map(function(yr) {
  var monthly = ee.ImageCollection(
    FLOOD_MONTHS.map(function(m) { return monthlyMaxComposite(yr, m); })
  );
  return ee.Image(monthly.sort('water_area_m2', false).first())
    .set('year', yr);
});


// ── 9. Console: annual table ─────────────────────────────────────
var annualStats = ee.FeatureCollection(yearlyMaxImages.map(function(img) {
  return ee.Feature(null, {
    year:        img.get('year'),
    peak_month:  img.get('date_str'),
    water_km2:   img.get('water_area_km2'),
    scene_count: img.get('scene_count')
  });
}));

print('Index used: ' + INDEX);
print('── Annual max flood extent ──');
print('Columns: year | peak_month | water_km2 | S2_scenes_that_month');
print(annualStats.reduceColumns(
  ee.Reducer.toList(4),
  ['year', 'peak_month', 'water_km2', 'scene_count']
));


// ── 10. Chart ────────────────────────────────────────────────────
print(
  ui.Chart.feature.byFeature(annualStats, 'year', 'water_km2')
  .setChartType('ColumnChart')
  .setOptions({
    title:  'Argon Pedion — Annual Max Flood Extent (Sentinel-2 ' + INDEX + ')',
    hAxis:  {title: 'Year', format: '####', gridlines: {count: YEARS.length}},
    vAxis:  {title: 'Max Flooded Area (km²)', minValue: 0},
    colors: ['#2e7d32'],
    legend: {position: 'none'}
  })
);


// ── 11. Map ──────────────────────────────────────────────────────
Map.centerObject(aoi, 13);

// Dry-season true-colour basemap
var s2Dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .map(maskS2Clouds)
  .median().clip(aoi);
Map.addLayer(s2Dry, {bands: ['B4','B3','B2'], min: 0, max: 3000}, 'S2 dry basemap');
Map.addLayer(aoi,   {color: 'FFFFFF'}, 'AOI');

// Annual max layer per year (60% opacity so overlapping years remain visible)
var legend = ui.Panel({
  style: {
    position:        'bottom-right',
    padding:         '8px',
    backgroundColor: 'rgba(255,255,255,0.92)'
  }
});
legend.add(ui.Label('Annual Max Flood (' + INDEX + ')', {
  fontWeight: 'bold', fontSize: '12px', margin: '0 0 5px 0'
}));

for (var i = 0; i < YEARS.length; i++) {
  var yr    = YEARS[i];
  var color = COLORS[i];
  var img   = ee.Image(yearlyMaxImages[i]);
  // 'water' band = pre-computed combined mask (open water + flooded veg + slope filter)
  var water = img.select('water').selfMask().clip(aoi);
  Map.addLayer(water, {palette: [color]}, yr + ' MAX', true, 0.6);

  legend.add(ui.Panel([
    ui.Label('', {backgroundColor: color, padding: '5px 12px', margin: '1px 5px 1px 0'}),
    ui.Label(yr + ' MAX', {fontSize: '11px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal')));
}

Map.add(legend);


// ── Notes ────────────────────────────────────────────────────────
// • scene_count: CS+-cleared S2 acquisitions available for the peak
//   month. Low counts (1–2) mean sparse coverage — treat with caution.
//
// • Cloud masking: Cloud Score Plus (CS+) replaces QA60 in flood
//   composites. CS+ catches cloud shadows, which are spectrally
//   similar to water and cause MNDWI false positives (QA60 does not
//   mask shadows). QA60 is retained for the dry-season basemap only.
//   Source: Sharma & Saharia 2025.
//
// • Adaptive threshold: each month uses Otsu's method on the WI
//   histogram, clamped to [THRESH_MIN, THRESH_MAX]. This adapts to
//   inter-annual variation in turbidity and soil moisture that shifts
//   the water/land spectral boundary. Source: Kordelas et al. 2018.
//
// • Slope mask: pixels with terrain slope > 5° are excluded as
//   non-floodable (removes limestone hillside artefacts).
//   Source: Sharma & Saharia 2025.
//
// • 2017: no Jan–Apr S2 SR data (collection starts Mar 2017);
//   annual max drawn from Oct–Dec 2017 only.
//   2026: no Oct–Dec data yet; max drawn from Jan–Apr 2026.
//
// • To switch to plain NDWI: change INDEX = 'NDWI' at the top.
//   MNDWI is recommended — SWIR1 is less sensitive to vegetation
//   shade than NIR, reducing false positives at floodplain edges.
