import pandas as pd
import numpy as np
from Orbit import Orbit
import ephem
import matplotlib.pyplot as plt
import os
import glob

#y4_fakes = pd.read_csv('/Volumes/lsa-gerdes/wsdiff_catalogs/season200/fakesonly/wsdiff_season200_Y4_griz_fakesonly.csv')


class Triplets:

    def fit_orbit(self, df_obs):

        #df_obs = df_obs.ix[['#' not in row['date'] for ind, row in df_obs.iterrows()]]  # filter comment lines
        nobs = len(df_obs)
        ralist = [ephem.hours(r) for r in df_obs['ra'].values]
        declist = [ephem.degrees(r) for r in df_obs['dec'].values]
        datelist = [ephem.date(d) for d in df_obs['date'].values]
        obscode = np.ones(nobs, dtype=int) * 807
        orbit = Orbit(dates=datelist, ra=ralist, dec=declist, obscode=obscode, err=0.15)
        # print orbit.chisq
        return orbit

    def groupObs(self, csv):
        lrg_df_obs = pd.read_csv(csv)
        list_unique = lrg_df_obs.fakeid.unique()
        baryc_df = pd.DataFrame(columns=('fakid', 'baryc_dist', 'inc'))
        for fakeid in list_unique:
            indiv_df_obs = lrg_df_obs.loc[lrg_df_obs['fakeid'] == fakeid]
            if len(indiv_df_obs) > 2:
                indiv_orb = self.fit_orbit(indiv_df_obs)
                indiv_orbitdf = indiv_orb.orbit2df()
                data = pd.DataFrame({"fakid": fakeid, "baryc_dist": indiv_orb.barycentric_distance()[0], "inc": indiv_orbitdf['inc']})
                baryc_df = baryc_df.append(data)
        return baryc_df

    def indivCSV(self, csv, year):
        '''
        
        :param csv: Must have at least 4 observations
        :return: 
        '''
        indiv_df_obs = pd.read_csv(csv)
        indiv_orb = self.fit_orbit(indiv_df_obs)
        indiv_orbitdf = indiv_orb.orbit2df()
        baryc_df = pd.DataFrame(columns=('baryc_dist', 'inc', 'year'))
        data = pd.DataFrame(
            {"baryc_dist": indiv_orb.barycentric_distance()[0], "inc": indiv_orbitdf['inc'], "year": "Y" + str(year)})
        baryc_df = baryc_df.append(data)

        return baryc_df

    def extend(self):
        pass



triplet = Triplets()

all_orbs_df = pd.DataFrame()

for i in range(1,5):
    os.chdir('/Volumes/lsa-gerdes/wsdiff_catalogs/season200/nofakes/cands_Y' + str(i))
    nitelist = glob.glob('4nite*')
    plotdf = pd.DataFrame()
    for csv in nitelist:
        gold_df = triplet.indivCSV(str(csv), i)
        plotdf = plotdf.append(gold_df)
        #plt.plot(gold_df['inc'], gold_df['baryc_dist'], 'ro')

    if not plotdf.empty:
        all_orbs_df = all_orbs_df.append(plotdf)
        os.chdir('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search')
        plt.plot(plotdf['inc'], plotdf['baryc_dist'], marker = 'o', linestyle='none', alpha= 0.5)



plt.ylabel('baryc_dist')
plt.xlabel ('inclination')
plt.xlim([0,37])
plt.title("Barycentric Dist vs Inclination")
plt.savefig('TRUEbarycvsinc.png')


result = all_orbs_df.sort_values(['baryc_dist', 'inc'], ascending=[0, 1])


print "done"
# gold_df = triplet.groupObs('/Volumes/lsa-gerdes/wsdiff_catalogs/season200/fakesonly/wsdiff_season200_Y4_griz_fakesonly.csv')
#
# plt.plot(gold_df['inc'], gold_df['baryc_dist'], 'ro')
# plt.ylabel('baryc_dist')
# plt.xlabel ('inclination')
# plt.title("Barycentric Dist vs Inclination")
# plt.xlim([0,190])
# plt.ylim([0,400])
# plt.savefig('barycvsinc.png')

#baryc = orb.barycentric_distance()
#baryc_dist = str(baryc[0])
#baryc_error = str(baryc[1])

