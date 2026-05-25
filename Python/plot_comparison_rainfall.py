import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

ROOT    = pathlib.Path(__file__).parent.parent
CSV_FLD = ROOT / 'Outputs' / 'CSV' / 'argon_pedion_s1_vs_s2.csv'
CSV_RN  = ROOT / 'Outputs' / 'CSV' / 'argon_pedion_rainfall.csv'
OUT     = ROOT / 'Outputs' / 'PNG' / 'argon_pedion_s1_vs_s2_rainfall.png'

flood = pd.read_csv(CSV_FLD).sort_values('year').reset_index(drop=True)
rain  = pd.read_csv(CSV_RN).sort_values('year').reset_index(drop=True)
df    = flood.merge(rain, on='year')

years = df['year'].astype(int)
x     = np.arange(len(years))
width = 0.38

fig, ax = plt.subplots(figsize=(13, 5.5))

bars_s1 = ax.bar(x - width / 2, df['s1_km2'], width,
                 color='#1565c0', label='Sentinel-1 SAR (VV, -16 dB)')
bars_s2 = ax.bar(x + width / 2, df['s2_km2'], width,
                 color='#2e7d32', label='Sentinel-2 MNDWI (> 0.0)')

ax.set_xlabel('Year', fontsize=11)
ax.set_ylabel('Max Flooded Area (km²)', fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(years)
ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
ax.set_ylim(0, df[['s1_km2', 's2_km2']].max().max() * 1.25)
ax.grid(axis='y', linewidth=0.5, alpha=0.5)
ax.grid(axis='y', which='minor', linewidth=0.3, alpha=0.3)
ax.spines[['top', 'right']].set_visible(False)

# Value labels on bars
for bar in list(bars_s1) + list(bars_s2):
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 0.06,
            f'{h:.2f}', ha='center', va='bottom', fontsize=7, color='#333')

# Secondary axis — rainfall and net balance lines
ax2 = ax.twinx()

line_rain = ax2.plot(x, df['rain_mm'], color='#e65100', linewidth=1.5,
                     linestyle='--', marker='o', markersize=5, zorder=4,
                     label='Rainfall Oct-Apr (mm)')
line_bal  = ax2.plot(x, df['balance_mm'], color='#b71c1c', linewidth=2,
                     linestyle='-', marker='s', markersize=6, zorder=5,
                     label='Net balance = rainfall - ET0 (mm)')

ax2.axhline(0, color='#b71c1c', linewidth=0.7, linestyle=':', alpha=0.5)
ax2.set_ylabel('mm', fontsize=11, color='#7f0000')
ax2.tick_params(axis='y', colors='#7f0000')
ax2.spines['right'].set_color('#7f0000')
ax2.spines[['top', 'left']].set_visible(False)
ymax = max(df['rain_mm'].max(), abs(df['balance_mm']).max()) * 1.4
ax2.set_ylim(df['balance_mm'].min() * 1.5, ymax)

# Value labels on balance line only (less cluttered)
for xi, val in zip(x, df['balance_mm']):
    offset = 12 if val >= 0 else -18
    ax2.text(xi, val + offset, f'{val:.0f}',
             ha='center', va='bottom', fontsize=7, color='#b71c1c')

# Combined legend
bar_handles = [bars_s1, bars_s2]
line_handles = line_rain + line_bal
handles = bar_handles + line_handles
labels  = [h.get_label() for h in handles]
ax.legend(handles, labels, framealpha=0.9, fontsize=9.5,
          loc='upper left', ncol=2)

ax.set_title(
    'Argon Pedion — Annual Max Flood Extent vs Flood-Season Water Balance\n'
    'rainfall & ET0 = Oct(Y-1) to Apr(Y)  |  ERA5-Land via Open-Meteo',
    fontsize=12, pad=10
)

fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches='tight')
print(f'Saved: {OUT}')
