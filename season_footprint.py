import pandas as pd
import matplotlib.pyplot as plt


def diffimg(season):

    csv_file = ('/Volumes/lsa-gerdes/wsdiff_catalogs/season' + str(season) +
                '/nofakes/wsdiff_season' + str(season) + '_Y4_griz_nofakes.csv')
    raw_diff_img = pd.read_csv(csv_file)
    return raw_diff_img

for szn in [240, 241, 242, 250]:

    diff_df = diffimg(szn)
    plt.figure()
    plt.plot(diff_df['ra'], diff_df['dec'], linestyle='None', marker=',')
    plt.title("Season " + str(szn))
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.savefig("Season" + str(szn) + ".png")
    print "Finished Season " + str(szn)