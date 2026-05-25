"""
Runs the Sentinel-1 vs Sentinel-2 annual flood comparison for Argon Pedion
via the Earth Engine Python API, saving results to CSV without needing
Google Drive or the Code Editor.
"""

import ee
import pandas as pd
import pathlib

ee.Initialize(project='geo-tut')

ROOT    = pathlib.Path(__file__).parent.parent
CSV_OUT = ROOT / 'Outputs' / 'CSV' / 'argon_pedion_s1_vs_s2.csv'
CSV_OUT.parent.mkdir(parents=True, exist_ok=True)

# ── AOI ─────────────────────────────────────────────────────────────────────
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

# ── Parameters ───────────────────────────────────────────────────────────────
SAR_THRESHOLD   = -16  # calibrated: -15 dB gave 3.7% dry-season false positives on this site
MNDWI_THRESHOLD =  0.0
THRESH_MIN      = -0.1
THRESH_MAX      =  0.3
CS_THRESHOLD    =  0.65
SPECKLE_KERNEL  =  7
COMPARE_YEARS   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026]
FLOOD_MONTHS    = [10, 11, 12, 1, 2, 3, 4]

# ── Topographic mask ──────────────────────────────────────────────────────────
floodable_mask = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003')).lt(5).clip(aoi)


# ── Sentinel-1 ───────────────────────────────────────────────────────────────
def speckle_filter(image):
    linear   = ee.Image(10).pow(image.select('VV').divide(10))
    smoothed = linear.focal_mean(SPECKLE_KERNEL, 'square', 'pixels')
    return (smoothed.log10().multiply(10).rename('VV')
            .copyProperties(image, image.propertyNames()))

def add_s1_stats(img):
    date    = ee.Date(img.get('system:time_start'))
    area_m2 = ee.Number(
        img.lt(SAR_THRESHOLD).And(floodable_mask).rename('water').clip(aoi)
           .multiply(ee.Image.pixelArea())
           .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                         scale=10, maxPixels=1e10)
           .get('water'))
    return img.set({
        'water_area_m2':  area_m2,
        'water_area_km2': area_m2.divide(1e6),
        'date_str':       date.format('YYYY-MM-dd'),
        'year':           date.get('year'),
    })

s1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
      .filterBounds(aoi)
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
      .select('VV')
      .filter(ee.Filter.Or(
          ee.Filter.calendarRange(10, 12, 'month'),
          ee.Filter.calendarRange(1,  4,  'month')))
      .map(speckle_filter)
      .map(add_s1_stats))

def s1_annual_max(yr):
    return (ee.Image(s1.filter(ee.Filter.eq('year', yr))
                       .sort('water_area_m2', False).first())
            .set('year', yr))


# ── Sentinel-2 ───────────────────────────────────────────────────────────────
def apply_cs_plus(image):
    cs = ee.Image(image.get('cs_plus'))
    return image.updateMask(cs.select('cs').gte(CS_THRESHOLD))

def add_water_index(image):
    return image.addBands(image.normalizedDifference(['B3', 'B11']).rename('WI'))

def otsu(histogram):
    counts  = ee.Array(ee.Dictionary(histogram).get('histogram'))
    means   = ee.Array(ee.Dictionary(histogram).get('bucketMeans'))
    size    = means.length().get([0])
    total   = counts.reduce(ee.Reducer.sum(), [0]).get([0])
    sum_val = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
    mean    = sum_val.divide(total)

    def bss(i):
        i       = ee.Number(i)
        ac      = counts.slice(0, 0, i)
        a_count = ac.reduce(ee.Reducer.sum(), [0]).get([0])
        a_mean  = (means.slice(0, 0, i).multiply(ac)
                   .reduce(ee.Reducer.sum(), [0]).get([0]).divide(a_count))
        b_count = total.subtract(a_count)
        b_mean  = sum_val.subtract(a_count.multiply(a_mean)).divide(b_count)
        return (a_count.multiply(a_mean.subtract(mean).pow(2))
                .add(b_count.multiply(b_mean.subtract(mean).pow(2))))

    bss_list = ee.List.sequence(1, size).map(bss)
    return ee.Number(means.sort(bss_list).get([-1]))

def s2_monthly_max(yr, m):
    start = ee.Date.fromYMD(yr, m, 1)
    end   = start.advance(1, 'month')

    s2col  = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
              .filterBounds(aoi).filterDate(start, end))
    csplus = (ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
              .filterBounds(aoi).filterDate(start, end))

    joined = ee.Join.saveFirst('cs_plus').apply(
        primary=s2col, secondary=csplus,
        condition=ee.Filter.equals(leftField='system:index', rightField='system:index'))

    masked      = ee.ImageCollection(joined).map(apply_cs_plus).map(add_water_index)
    scene_count = masked.size()

    wi_max = (masked.select('WI')
              .merge(ee.ImageCollection([ee.Image.constant(-1).rename('WI')]))
              .max().clip(aoi))

    histogram = (masked.select('WI').max().clip(aoi)
                 .reduceRegion(
                     reducer=ee.Reducer.histogram(maxBuckets=200, minBucketWidth=0.01),
                     geometry=aoi, scale=20, maxPixels=1e10)
                 .get('WI'))

    thresh = ee.Number(ee.Algorithms.If(
        scene_count.gt(0),
        ee.Number(otsu(histogram)).max(THRESH_MIN).min(THRESH_MAX),
        MNDWI_THRESHOLD))

    water_mask = wi_max.gte(thresh).And(floodable_mask)

    area_m2 = ee.Number(
        water_mask.rename('water')
                  .multiply(ee.Image.pixelArea())
                  .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                                scale=20, maxPixels=1e10)
                  .get('water'))

    return (wi_max.addBands(water_mask.rename('water'))
            .set('year', yr).set('month', m)
            .set('date_str',       start.format('YYYY-MM'))
            .set('scene_count',    scene_count)
            .set('water_area_m2',  area_m2)
            .set('water_area_km2', area_m2.divide(1e6)))

def s2_annual_max(yr):
    monthly = ee.ImageCollection([s2_monthly_max(yr, m) for m in FLOOD_MONTHS])
    return ee.Image(monthly.sort('water_area_m2', False).first()).set('year', yr)


# ── Build annual comparison ───────────────────────────────────────────────────
print('Processing year by year...')
rows = []

for yr in COMPARE_YEARS:
    print(f'  {yr}...', end=' ', flush=True)
    s1_max = s1_annual_max(yr)
    s2_max = s2_annual_max(yr)
    s2_area_km2 = ee.Number(
        s2_max.select('water')
              .multiply(ee.Image.pixelArea())
              .reduceRegion(reducer=ee.Reducer.sum(), geometry=aoi,
                            scale=20, maxPixels=1e10)
              .get('water')).divide(1e6)

    feat = ee.Feature(None, {
        'year':     yr,
        's1_km2':   s1_max.get('water_area_km2'),
        's2_km2':   s2_area_km2,
        's1_date':  s1_max.get('date_str'),
        's2_month': s2_max.get('date_str'),
    })

    p = feat.getInfo()['properties']
    s1_val = round(float(p['s1_km2']), 4) if p.get('s1_km2') is not None else None
    s2_val = round(float(p['s2_km2']), 4) if p.get('s2_km2') is not None else None
    rows.append({
        'year':     int(p['year']),
        's1_km2':   s1_val,
        's2_km2':   s2_val,
        's1_date':  p.get('s1_date'),
        's2_month': p.get('s2_month'),
    })
    print(f"S1={s1_val} km²  S2={s2_val} km²  "
          f"(S1: {p.get('s1_date')}  S2: {p.get('s2_month')})")

df = pd.DataFrame(rows).sort_values('year').reset_index(drop=True)
df.to_csv(CSV_OUT, index=False)
print(f'\nSaved → {CSV_OUT}')
