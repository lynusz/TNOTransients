import TNODatabase as tno
from astropy.coordinates import SkyCoord
from astropy import units as u
import easyaccess as ea
import pandas as pd
import time
import numpy as np
import matplotlib.pyplot as plt

zeropoints = pd.read_csv('fgcm_zeropoints_v2_0.csv')
zeropoints_Y4 = pd.read_csv('Y4N_zeropoints_03.09.2017.csv')
all_exps = pd.read_csv('exposures.csv')


db = tno.Connect()
df_accepted = db.execute("SELECT * FROM LZULLO.TNOBS INNER JOIN TNOLINK ON TNOBS.OBJID = TNOLINK.OBJID " +
                         "WHERE SPREADERR_MODEL IS NOT NULL AND SPREAD_MODEL IS NOT NULL AND MAG " +
                         "IS NOT NULL AND ID NOT LIKE 'S%'")
print df_accepted.head()
df_accepted['starlike'] = df_accepted.SPREAD_MODEL + 3 * df_accepted.SPREADERR_MODEL

plt.plot(df_accepted['MAG'], df_accepted['starlike'], linestyle='None', marker='o')
plt.xlabel('Magnitude')
plt.ylabel(r'$\mathrm{SM} + 3\mathrm{SEM}$')
plt.title('Starlike Quantity vs. Mag for Submitted Observations')
plt.savefig('accepted_starlike.png')

plt.figure()
plt.hist(df_accepted['starlike'], 50, range=(0, 0.1))
plt.savefig('accepted_starlike_hist.png')

df_cand = db.execute("SELECT * FROM LZULLO.TNOBS WHERE SPREADERR_MODEL IS NOT NULL " +
                     "AND SPREAD_MODEL IS NOT NULL AND MAG IS NOT NULL")

df_cand['starlike'] = df_cand.SPREAD_MODEL + 3 * df_cand.SPREADERR_MODEL
df_cand_res = df_cand[(df_cand.MAG > 20)]

plt.figure()
plt.plot(df_cand['MAG'], df_cand['starlike'], linestyle='None', marker=',')
plt.xlabel('Magnitude')
plt.ylabel(r'$\mathrm{SM} + 3\mathrm{SEM}$')
plt.title('Starlike Quantity vs. Mag for Candidate Observations')
plt.savefig('cand_starlike.png')

plt.figure()
plt.hist(df_cand['starlike'], 50, range=(0, 0.1))
plt.savefig('cand_starlike_hist.png')

plt.figure()
plt.plot(df_cand_res['MAG'], df_cand_res['starlike'], linestyle='None', marker=',')
plt.xlabel('Magnitude')
plt.ylabel(r'$\mathrm{SM} + 3\mathrm{SEM}$')
plt.title('Starlike Quantity vs. Mag for Candidate Observations')
plt.savefig('cand_starlike_res.png')

plt.figure()
plt.plot(df_cand['MAG'], df_cand['FLUX'] / df_cand['FLUX_ERR'], linestyle='None', marker=',')
plt.xlabel('Magnitude')
plt.ylabel('Flux/FluxErr')
plt.title('Signal to Noise vs. Mag for Candidate Observations')
plt.savefig('cand_sig2noise.png')

plt.figure()
plt.plot(df_cand_res['MAG'], df_cand_res['FLUX'] / df_cand_res['FLUX_ERR'], linestyle='None', marker=',')
plt.xlabel('Magnitude')
plt.ylabel('Flux/FluxErr')
plt.title('Signal to Noise vs. Mag for Candidate Observations')
plt.savefig('cand_sig2noise_res.png')

plt.figure()
plt.plot(df_accepted['MAG'], df_accepted['FLUX'] / df_accepted['FLUX_ERR'], linestyle='None', marker='o')
plt.xlabel('Magnitude')
plt.ylabel('Flux/FluxErr')
plt.title('Signal to Noise vs. Mag for Submitted Observations')
plt.savefig('accepted_sig2noise.png')

df_mag = db.execute('SELECT * FROM TNOBS WHERE MAG IS NOT NULL')

plt.figure()
plt.hist(df_mag['MAG'], 50, range=(15, 28))
plt.savefig('cand_mag_hist.png')
