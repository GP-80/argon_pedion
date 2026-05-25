// ================================================================
// Argon Pedion — Sentinel-1 vs Sentinel-2 Comparison Script
// ================================================================
// Compares the two flood detection approaches:
//   S1: VV backscatter threshold on individual scenes (cloud-free)
//   S2: MNDWI monthly max composite (flood-persistence assumption)
//
// Outputs
// ───────
//   Tasks tab → Google Drive (CSV):
//     1. Annual comparison table — year | s1_km2 | s2_km2 | s1_date | s2_peak_month
//        (use in Sheets / Excel / Python to plot the S1 vs S2 bar chart)
//   Map:
//     2. Spatial agreement for COMPARE_YEAR:
//          Purple = both sensors detected water
//          Blue   = S1 only
//          Green  = S2 only
//     3. Jaccard similarity index + area stats in results panel
// ================================================================


// ── 1. AOI ──────────────────────────────────────────────────────
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
var SAR_THRESHOLD      = -15;   // dB — S1 VV open-water threshold
var MNDWI_THRESHOLD    =  0.0;  // S2 fallback threshold when Otsu has no data
var THRESH_MIN         = -0.1;  // Otsu result clamped lower bound
var THRESH_MAX         =  0.3;  // Otsu result clamped upper bound
var CS_THRESHOLD   =  0.65; // Cloud Score Plus clarity — replaces QA60 in S2
var SPECKLE_KERNEL =  7;

// Change this to any year in COMPARE_YEARS to inspect a different year
var COMPARE_YEAR  = 2019;

// Overlapping period: S2 SR starts 2017, S1 starts 2014
var COMPARE_YEARS = [2017,2018,2019,2020,2021,2022,2023,2024,2025,2026];
var FLOOD_MONTHS  = [10, 11, 12, 1, 2, 3, 4];


// ══════════════════════════════════════════════════════════════════
// SHARED — topographic mask (used by both S1 and S2)
// ══════════════════════════════════════════════════════════════════
// Slope > 5° excluded as non-floodable (limestone hillside artefacts).
// Source: Sharma & Saharia 2025.
var floodableMask = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'))
  .lt(5).clip(aoi);


// ══════════════════════════════════════════════════════════════════
// SENTINEL-1 PROCESSING
// ══════════════════════════════════════════════════════════════════

function speckleFilter(image) {
  var linear   = ee.Image(10).pow(image.select('VV').divide(10));
  var smoothed = linear.focal_mean(SPECKLE_KERNEL, 'square', 'pixels');
  return smoothed.log10().multiply(10).rename('VV')
    .copyProperties(image, image.propertyNames());
}

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(aoi)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV')
  .filter(ee.Filter.or(
    ee.Filter.calendarRange(10, 12, 'month'),
    ee.Filter.calendarRange(1,  4,  'month')
  ))
  .map(speckleFilter)
  .map(function(img) {
    var date   = ee.Date(img.get('system:time_start'));
    // Apply slope mask; wrap in ee.Number immediately to avoid lazy-eval sort bug
    var areaM2 = ee.Number(
      img.lt(SAR_THRESHOLD).and(floodableMask).rename('water').clip(aoi)
        .multiply(ee.Image.pixelArea())
        .reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale: 10, maxPixels: 1e10})
        .get('water')
    );
    return img
      .set('water_area_m2',  areaM2)
      .set('water_area_km2', areaM2.divide(1e6))
      .set('date_str',       date.format('YYYY-MM-dd'))
      .set('year',           date.get('year'));
  });

// Annual max S1 scene per year
function s1AnnualMax(yr) {
  return ee.Image(
    s1.filter(ee.Filter.eq('year', yr)).sort('water_area_m2', false).first()
  ).set('year', yr);
}


// ══════════════════════════════════════════════════════════════════
// SENTINEL-2 PROCESSING
// ══════════════════════════════════════════════════════════════════

// QA60 kept for the dry-season basemap only (no join needed)
function maskS2Clouds(image) {
  var qa = image.select('QA60');
  return image.updateMask(
    qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0))
  );
}

// Cloud Score Plus — catches cloud shadows that QA60 misses
// (shadows have low reflectance, mimic water in MNDWI)
function applyCSPlus(image) {
  return image.updateMask(
    ee.Image(image.get('cs_plus')).select('cs').gte(CS_THRESHOLD)
  );
}

// Otsu threshold: maximises between-class variance (water vs land)
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

function addWaterIndex(image) {
  return image.addBands(
    image.normalizedDifference(['B3', 'B11']).rename('WI')
  );
}

function s2MonthlyMax(yr, m) {
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

  var wiMax = masked.select('WI')
    .merge(ee.ImageCollection([ee.Image.constant(-1).rename('WI')]))
    .max().clip(aoi);

  var histogram = masked.select('WI').max().clip(aoi)
    .reduceRegion({
      reducer: ee.Reducer.histogram({maxBuckets: 200, minBucketWidth: 0.01}),
      geometry: aoi, scale: 20, maxPixels: 1e10
    }).get('WI');

  var thresh = ee.Number(ee.Algorithms.If(
    sceneCount.gt(0),
    ee.Number(otsu(histogram)).max(THRESH_MIN).min(THRESH_MAX),
    MNDWI_THRESHOLD
  ));

  var waterMask = wiMax.gte(thresh).and(floodableMask);

  var areaM2 = ee.Number(
    waterMask.rename('water')
      .multiply(ee.Image.pixelArea())
      .reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10})
      .get('water')
  );

  return wiMax.addBands(waterMask.rename('water'))
    .set('year',          yr).set('month', m)
    .set('date_str',      start.format('YYYY-MM'))
    .set('scene_count',   sceneCount)
    .set('water_area_m2', areaM2)
    .set('water_area_km2', areaM2.divide(1e6));
}

// Annual max S2 composite per year
function s2AnnualMax(yr) {
  var monthly = ee.ImageCollection(
    FLOOD_MONTHS.map(function(m) { return s2MonthlyMax(yr, m); })
  );
  return ee.Image(monthly.sort('water_area_m2', false).first()).set('year', yr);
}


// ══════════════════════════════════════════════════════════════════
// ANNUAL COMPARISON TABLE
// ══════════════════════════════════════════════════════════════════

var annualComparison = ee.FeatureCollection(
  COMPARE_YEARS.map(function(yr) {
    var s1Max = s1AnnualMax(yr);
    var s2Max = s2AnnualMax(yr);
    // Recompute S2 area from the pre-computed 'water' band (open water + flooded veg
    // + slope mask) rather than reading the stored property — guards against lazy-eval
    // caching in batch export.
    var s2AreaKm2 = ee.Number(
      s2Max.select('water')
        .multiply(ee.Image.pixelArea())
        .reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10})
        .get('water')
    ).divide(1e6);
    return ee.Feature(null, {
      year:     yr,
      s1_km2:   s1Max.get('water_area_km2'),
      s2_km2:   s2AreaKm2,
      s1_date:  s1Max.get('date_str'),
      s2_month: s2Max.get('date_str')
    });
  })
);

// ══════════════════════════════════════════════════════════════════
// EXPORT — annual S1 vs S2 comparison table → Google Drive (CSV)
// ══════════════════════════════════════════════════════════════════
// Printing annualComparison as a Console chart triggers ~70 concurrent
// reduceRegion calls (7 flood months × 10 years × S2) and hits the GEE
// session memory cap. Export.table.toDrive runs as a batch job instead,
// which has no session memory limit.
//
// How to run:
//   1. Click Run — the task appears in the Tasks tab (right panel).
//   2. Click RUN next to 'argon_pedion_s1_vs_s2'.
//   3. When the job completes, the CSV lands in My Drive › GEE_exports.
//   4. Open in Google Sheets / Excel / Python to plot the bar chart.
//
// CSV columns: year | s1_km2 | s2_km2 | s1_date | s2_peak_month
Export.table.toDrive({
  collection:     annualComparison,
  description:    'argon_pedion_s1_vs_s2',
  folder:         'GEE_exports',
  fileNamePrefix: 'argon_pedion_s1_vs_s2',
  fileFormat:     'CSV',
  selectors:      ['year', 's1_km2', 's2_km2', 's1_date', 's2_month']
});


var s1MaxScene     = s1AnnualMax(COMPARE_YEAR);
var s2MaxComposite = s2AnnualMax(COMPARE_YEAR);


// ══════════════════════════════════════════════════════════════════
// MAP — spatial agreement for COMPARE_YEAR
// ══════════════════════════════════════════════════════════════════
Map.centerObject(aoi, 13);

// Dry basemap
var s2Dry = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .map(maskS2Clouds).median().clip(aoi);
Map.addLayer(s2Dry, {bands: ['B4','B3','B2'], min: 0, max: 3000}, 'S2 dry basemap');
Map.addLayer(aoi, {color: 'FFFFFF'}, 'AOI');

// Binary water masks for COMPARE_YEAR (both include slope filter)
var s1Water = s1MaxScene.lt(SAR_THRESHOLD).and(floodableMask).clip(aoi);
var s2Water = s2MaxComposite.select('water').clip(aoi);

// Agreement categories
var bothWater = s1Water.and(s2Water).selfMask();
var s1Only    = s1Water.and(s2Water.not()).selfMask();
var s2Only    = s2Water.and(s1Water.not()).selfMask();

Map.addLayer(s2Only,    {palette: ['2e7d32']}, 'S2 only (COMPARE_YEAR)');
Map.addLayer(s1Only,    {palette: ['1565c0']}, 'S1 only (COMPARE_YEAR)');
Map.addLayer(bothWater, {palette: ['7b1fa2']}, 'Both sensors (COMPARE_YEAR)');

// Raw images for context (off by default)
Map.addLayer(s1MaxScene.clip(aoi),
  {min: -25, max: 0, palette: ['000000','FFFFFF']},
  'S1 VV backscatter — ' + COMPARE_YEAR, false);

Map.addLayer(s2MaxComposite.select('WI').clip(aoi),
  {min: -0.5, max: 0.5, palette: ['#795548','#FFFFFF','#0D47A1']},
  'S2 WI composite — ' + COMPARE_YEAR, false);

// Map legend
var legend = ui.Panel({
  style: {position: 'bottom-right', padding: '8px',
          backgroundColor: 'rgba(255,255,255,0.92)'}
});
legend.add(ui.Label('Spatial Agreement — ' + COMPARE_YEAR,
  {fontWeight:'bold', fontSize:'12px', margin:'0 0 5px 0'}));

function legendRow(color, label) {
  return ui.Panel([
    ui.Label('', {backgroundColor: color, padding:'5px 12px', margin:'1px 5px 1px 0'}),
    ui.Label(label, {fontSize:'11px', margin:'2px 0'})
  ], ui.Panel.Layout.flow('horizontal'));
}
legend.add(legendRow('#7b1fa2', 'Both detected water'));
legend.add(legendRow('#1565c0', 'S1 only'));
legend.add(legendRow('#2e7d32', 'S2 only'));
Map.add(legend);


// ── Results panel (map): Jaccard + area stats for COMPARE_YEAR ───
// All three area reductions run at scale=20 so S1 and S2 are compared
// on a consistent grid (S1 native is 10 m but 20 m is fine for ratios).
var scale20 = {reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10};

// After .and()/.or() the output band name comes from the first image.
// bothWater = s1Water.and(s2Water) → band 'VV' (from s1Water)
// Rename everything to 'water' to make the .get() key unambiguous.
var intersectArea = bothWater.rename('water')
  .multiply(ee.Image.pixelArea()).reduceRegion(scale20).get('water');

var unionAreaVal = s1Water.or(s2Water).rename('water')
  .multiply(ee.Image.pixelArea()).reduceRegion(scale20).get('water');

var jaccard = ee.Number(intersectArea).divide(ee.Number(unionAreaVal));

var statsPanel = ui.Panel({
  style: {position: 'bottom-left', padding:'10px',
          backgroundColor: 'rgba(255,255,255,0.92)', width:'260px'}
});
statsPanel.add(ui.Label('Agreement Stats — ' + COMPARE_YEAR,
  {fontWeight:'bold', fontSize:'13px', margin:'0 0 6px 0'}));

var sLabel = ui.Label('Computing…', {fontSize:'12px'});
statsPanel.add(sLabel);

ee.Dictionary({
  s1_km2:       s1MaxScene.get('water_area_km2'),
  s1_date:      s1MaxScene.get('date_str'),
  s2_km2:       s2MaxComposite.get('water_area_km2'),
  s2_month:     s2MaxComposite.get('date_str'),
  intersect_m2: intersectArea,
  union_m2:     unionAreaVal,
  jaccard:      jaccard
}).evaluate(function(d) {
  var jPct = (d.jaccard * 100).toFixed(1);
  sLabel.setValue(
    'S1:  ' + d.s1_km2.toFixed(2) + ' km²  (' + d.s1_date + ')\n' +
    'S2:  ' + d.s2_km2.toFixed(2) + ' km²  (' + d.s2_month + ')\n' +
    '────────────────────────\n' +
    'Jaccard index:  ' + jPct + ' %\n' +
    '(overlap / union — 100% = perfect)'
  );
});

statsPanel.add(ui.Label(
  'Jaccard < 50%: sensors disagree on location\n' +
  'Jaccard > 75%: sensors agree well',
  {fontSize:'10px', color:'#666', margin:'6px 0 0 0'}
));

Map.add(statsPanel);
