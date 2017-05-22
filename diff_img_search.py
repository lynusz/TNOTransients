
from __future__ import division

import easyaccess as ea
import ephem
import numpy as np
import pandas as pd
import os
import extendOrbit as ext

#Generic testing



def searchBox(season_df, ra_deg, dec_deg):
    cut = ((season_df['ra'] >= ra_deg) & (season_df['ra'] <= season_df['ra'] + ra_deg) & (season_df['dec'] >= season_df['dec'] - dec_deg) & (season_df['dec'] <= season_df['dec'] + dec_deg))
    season_df = season_df[cut]
    return season_df

#test = pd.read_csv('/Volumes/lsa-gerdes/wsdiff_catalogs/season242/nofakes/wsdiff_season242_Y3_griz_nofakes.csv')
mybox = pd.read_csv('/Volumes/lsa-gerdes/wsdiff_catalogs/season241/nofakes/wsdiff_season241_Y4_griz_nofakes.csv')
#y4_dir = ('/Volumes/lsa-gerdes/wsdiff_catalogs/season242/nofakes/')
print os.getcwd()
#os.chdir('/Volumes/lsa-gerdes/wsdiff_catalogs/season242/nofakes/')
print os.getcwd()
#mybox = pd.read_csv('wsdiff_season242_Y4_griz_nofakes.csv')


mybox = searchBox(mybox, 0.01, 0.01)

print "done"