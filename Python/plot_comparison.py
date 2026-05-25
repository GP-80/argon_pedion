import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

ROOT = pathlib.Path(__file__).parent.parent
CSV  = ROOT / 'Outputs' / 'CSV' / 'argon_pedion_s1_vs_s2.csv'
OUT  = ROOT / 'Outputs' / 'PNG' / 'argon_pedion_s1_vs_s2.png'

df = pd.read_csv(CSV)
df = df.sort_values('year').reset_index(drop=True)

years  = df['year'].astype(int)
x      = np.arange(len(years))
width  = 0.38

fig, ax = plt.subplots(figsize=(11, 5))

bars_s1 = ax.bar(x - width / 2, df['s1_km2'], width,
                 color='#1565c0', label='Sentinel-1 SAR (VV, -16 dB)')
bars_s2 = ax.bar(x + width / 2, df['s2_km2'], width,
                 color='#2e7d32', label='Sentinel-2 MNDWI (> 0.0)')

ax.set_title('Argon Pedion — Annual Max Flood Extent', fontsize=13, pad=12)
ax.set_xlabel('Year', fontsize=11)
ax.set_ylabel('Max Flooded Area (km²)', fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(years)
ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
ax.set_ylim(0, df[['s1_km2', 's2_km2']].max().max() * 1.18)
ax.legend(framealpha=0.9, fontsize=10)
ax.grid(axis='y', linewidth=0.5, alpha=0.6)
ax.grid(axis='y', which='minor', linewidth=0.3, alpha=0.4)
ax.spines[['top', 'right']].set_visible(False)

# value labels above each bar
for bar in list(bars_s1) + list(bars_s2):
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 0.05,
            f'{h:.2f}', ha='center', va='bottom', fontsize=7.5, color='#333')

fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches='tight')
print(f'Saved: {OUT}')
