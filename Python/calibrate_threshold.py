"""
Tests multiple SAR thresholds against dry-season false positives and
JRC max water extent to find the best calibrated value for this site.
"""
import ee

ee.Initialize(project='geo-tut')

aoi = ee.Geometry.Polygon(
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
    None, False)

SPECKLE_KERNEL = 7

floodable_mask = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003')).lt(5).clip(aoi)

def speckle_filter(image):
    linear   = ee.Image(10).pow(image.select('VV').divide(10))
    smoothed = linear.focal_mean(SPECKLE_KERNEL, 'square', 'pixels')
    return (smoothed.log10().multiply(10).rename('VV')
            .copyProperties(image, image.propertyNames()))

# Dry-season median (Jun-Aug, 2019-2023 — multiple years for robustness)
dry_s1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
          .filterBounds(aoi)
          .filter(ee.Filter.eq('instrumentMode', 'IW'))
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
          .select('VV')
          .filter(ee.Filter.calendarRange(6, 8, 'month'))
          .filter(ee.Filter.calendarRange(2019, 2023, 'year'))
          .map(speckle_filter)
          .median().clip(aoi))

# Best wet year: February 2019 (highest S1 in 2019, +318mm water balance)
wet_s1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
          .filterBounds(aoi)
          .filter(ee.Filter.eq('instrumentMode', 'IW'))
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
          .select('VV')
          .filterDate('2019-02-15', '2019-02-22')
          .map(speckle_filter)
          .median().clip(aoi))

# JRC max water extent (ever-flooded pixels, occurrence >= 5%)
jrc_max_km2 = (ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
               .select('max_extent')
               .clip(aoi)
               .rename('w')
               .multiply(ee.Image.pixelArea())
               .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                             scale=30, maxPixels=1e10)
               .get('w'))
jrc_km2 = ee.Number(jrc_max_km2).divide(1e6).getInfo()
print(f'JRC max water extent (reference ceiling): {jrc_km2:.3f} km2\n')

print(f'{"Threshold":>10}  {"Dry FP (km2)":>13}  {"Wet detect (km2)":>17}  {"Dry FP %":>9}')
print('-' * 58)

for thresh in [-15, -16, -17, -18, -19, -20]:
    dry_fp = (dry_s1.lt(thresh).And(floodable_mask).rename('w')
              .multiply(ee.Image.pixelArea())
              .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                            scale=10, maxPixels=1e10).get('w'))
    wet_det = (wet_s1.lt(thresh).And(floodable_mask).rename('w')
               .multiply(ee.Image.pixelArea())
               .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                             scale=10, maxPixels=1e10).get('w'))
    result = ee.Dictionary({'dry': dry_fp, 'wet': wet_det}).getInfo()
    dry_km2 = result['dry'] / 1e6
    wet_km2 = result['wet'] / 1e6
    pct     = dry_km2 / 9.37 * 100
    marker  = '  <-- JRC ceiling' if abs(wet_km2 - jrc_km2) < 0.3 else ''
    print(f'{thresh:>10} dB  {dry_km2:>10.3f} km2  {wet_km2:>14.3f} km2  {pct:>8.1f}%{marker}')
