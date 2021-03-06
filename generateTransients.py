from __future__ import division

import easyaccess as ea
import ephem
import numpy as np
import pandas as pd
from linkmap import build_kdtree
import matplotlib.pyplot as plt
import os
import Find_imgs
import time
from ExposureCheck import MPChecker
from MPCRecord import MPCToTNO
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

rows_so_far = 0.
zeropoints = pd.read_csv('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/fgcm_zeropoints_v2_0.csv')
zeropoints_Y4 = pd.read_csv('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/Y4N_zeropoints_03.09.2017.csv')
all_exps = pd.read_csv('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/exposures.csv')


def get_coadd_cutout(connection, ra, dec, box):
    '''
    Gets Y3A1_COADD_OBJECTs in the given search box.
    :param connection: easyaccess database connection to DESSCI
    :param ra: RA of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param dec: DEC of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param box: width of search region (+/- from center) in arcsec (i.e. the units returned by predict_pos)
    :return: dataframe of Y3A1_COADD_OBJECTs.
    '''

    query_Y3A1_coadd = r""" SELECT o.ra, o.dec, o.mag_auto_g, o.mag_auto_r, o.mag_auto_i, o.mag_auto_z, 
    o.nepochs_g, o.nepochs_r, o.nepochs_i, o.nepochs_z, o.nepochs_Y,
    o.spread_model_g, o.spread_model_r, o.spread_model_i, o.spread_model_z, 
    o.spreaderr_model_g, o.spreaderr_model_r, o.spreaderr_model_i, o.spreaderr_model_z,
    o.flags_g, o.flags_r, o.flags_i, o.flags_z
    FROM Y3A1_COADD_OBJECT_SUMMARY o WHERE 
    """

    ra_deg, dec_deg = np.degrees(ra), np.degrees(dec)
    box_deg = box / 3600
    query = query_Y3A1_coadd + "o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg) + "+" + str(
        box_deg) \
            + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and " + str(dec_deg) + "+" + str(box_deg)
    result = connection.query_to_pandas(query)
    result.columns = map(str.lower, result.columns)  # convert column names to lowercase
    result.rename(columns={'date_obs': 'date'}, inplace=True)
    return result


def get_SE_detections(connection, ra=None, dec=None, box=None, expnum=None):
    '''
    Gets single-epoch detections in all exposures from year for within a search box.

    :param connection: easyaccess database connection. Currently should be to DESOPER if year='Y4' and to DESSCI for earlier years.
    :param expnum: exposure number to search
    :param ra: RA of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param dec: DEC of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param box: width of search region (+/- from center) in arcsec (i.e. the units returned by predict_pos)
    :return: dataframe of single-epoch detections.
    '''

    query = r"""
    with x as (
        select qa.expnum as expnum,
            max(qa.lastchanged_time) as timestamp
        from prod.firstcut_eval qa
        where qa.analyst!='SNQUALITY'

        group by qa.expnum
    )
    select c.expnum, c.filename,
        cast(o.band as VARCHAR(1)) as BAND,
        c.ccdnum as CCD,
        o.object_number, o.ra, o.dec,
        o.flux_auto, o.fluxerr_auto,
        o.spread_model, o.spreaderr_model, o.flags_model
    from prod.se_object o, prod.catalog c, prod.proctag t, prod.firstcut_eval f, x
    where t.tag='Y4N_FIRSTCUT'
        and t.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_firstcut'
        and c.expnum=f.expnum
        and f.expnum=x.expnum
        and f.lastchanged_time=x.timestamp
        and f.accepted='True'
        and o.filename=c.filename
    """
    if expnum:
        query += " and c.expnum = " + str(expnum)

    if ra or dec or box:
        ra_deg, dec_deg = np.degrees(ra), np.degrees(dec)
        box_deg = box / 3600
        query += (" and o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg)
                 + "+" + str(box_deg) + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and "
                 + str(dec_deg) + "+" + str(box_deg))

    print query
    result = connection.query_to_pandas(query)
    result.columns = map(str.lower, result.columns)  # convert column names to lowercase
    result.rename(columns={'date_obs': 'date'}, inplace=True)
    merged = result.merge(all_exps[['expnum', 'nite', 'band', 'date', 't_eff', 'fwhm_asec', 'exptime']], how='inner',
                          on=['expnum', 'band'])
    merged = merged.loc[merged['exptime'] > 80.0]
    return merged


def get_transient_detections(df_SE, df_coadd, threshold):
    '''
    Given dataframes of single-epoch and coadd detections, removes SE detections within 1 arcsec of a coadd detection 
    and returns the remaining SE detections as possible transients. This is a catalog-level attempt at difference-imaging.
    :param df_SE: dataframe of SE detections (the search sample)
    :param df_coadd: dataframe of coadd detections (the template sample)
    :return: dataframe of unmatched SE detections
    '''
    if len(df_SE) and len(df_coadd):
        df_coadd['ra'] = np.radians(df_coadd['ra'])
        df_coadd['dec'] = np.radians(df_coadd['dec'])
        tree_coadd = build_kdtree(df_coadd)

        df_SE['ra'] = np.radians(df_SE['ra'])
        df_SE['dec'] = np.radians(df_SE['dec'])
        tree_SE = build_kdtree(df_SE)

        near_list = tree_SE.query_ball_tree(tree_coadd, threshold * np.pi / 648000)  # require 1 arcsec match to a coadd object
        transients = [i for i in range(len(near_list)) if near_list[i] == []]
        df_coadd['ra'] = np.degrees(df_coadd['ra'])
        df_coadd['dec'] = np.degrees(df_coadd['dec'])
        df_SE['ra'] = np.degrees(df_SE['ra'])
        df_SE['dec'] = np.degrees(df_SE['dec'])

        return df_SE.ix[transients]
    else:
        return pd.DataFrame()


def get_diffimg_cutout(ra, dec, box, season, path=None):
    if not path:
        csv_file = '/Volumes/lsa-gerdes/wsdiff_catalogs'
    else:
        csv_file = path

    csv_file += '/season' + str(season) + '/nofakes/wsdiff_season'\
                + str(season) + '_Y4_griz_nofakes.csv'

    raw_diff_img = pd.read_csv(csv_file)
    ra_deg, dec_deg = np.degrees(ra), np.degrees(dec)

    box_deg = box / 3600
    cut = (((float(ra_deg) - box_deg <= raw_diff_img['ra']) & (raw_diff_img['ra'] <= float(ra_deg) + box_deg)) &
           ((float(dec_deg) - box_deg <= raw_diff_img['dec']) & (raw_diff_img['dec'] <= float(dec_deg) + box_deg)))
    return raw_diff_img[cut]


def get_diffimg_season(season, path=None):
    if not path:
        csv_file = '/Volumes/lsa-gerdes/wsdiff_catalogs'
    else:
        csv_file = path

    csv_file += '/season' + str(season) + '/nofakes/wsdiff_season' \
                + str(season) + '_Y4_griz_nofakes.csv'

    return pd.read_csv(csv_file, engine='python')


def overlap(df1, df2=None, datematch=True, dropOverlap=True, anti_datematch=False, threshold=1):
    """
    Calculates the overlap between two dataframes on ra and dec and puts those rows in a third dataframe.
    Only puts one row per overlap. Removes those rows from df1 and df2.
    :param df1: 
    :type df1: pd.DataFrame
    :param df2: 
    :type df2: pd.DataFrame
    :return: 
    """
    df1['ra'] = np.radians(df1['ra'])
    df1['dec'] = np.radians(df1['dec'])
    if type(df2) == pd.DataFrame:
        df_other = df2
        df2['ra'] = np.radians(df2['ra'])
        df2['dec'] = np.radians(df2['dec'])
    else:
        df_other = df1

    tree1 = build_kdtree(df1)
    tree2 = build_kdtree(df_other)

    near1 = tree1.query_ball_tree(tree2, threshold * np.pi / 648000)

    overlap_index1 = []
    overlap_index2 = []

    for i in range(len(near1)):
        if near1[i]:
            appended = False
            for j in near1[i]:
                if datematch:
                    if df1.iloc[i]['date'] == df_other.iloc[j]['date']:
                        overlap_index2.append(j)
                        if not appended:
                            overlap_index1.append(i)
                            appended = True
                elif anti_datematch:
                    if df1.iloc[i]['nite'] != df_other.iloc[j]['nite']:
                        overlap_index2.append(j)
                        if not appended:
                            overlap_index1.append(i)
                            appended = True

                else:
                    overlap_index2.append(j)
                    if not appended:
                        overlap_index1.append(i)
                        appended = True

    df1['ra'] = np.degrees(df1['ra'])
    df1['dec'] = np.degrees(df1['dec'])
    if type(df2) == pd.DataFrame:
        df2['ra'] = np.degrees(df2['ra'])
        df2['dec'] = np.degrees(df2['dec'])

    overlap_df = df1.iloc[overlap_index1]

    if dropOverlap:
        df1.drop(df1.index[overlap_index1], inplace=True)
        if type(df2) == pd.DataFrame:
            df2.drop(df2.index[overlap_index2], inplace=True)

    return overlap_df

def self_overlap(df, dropOverlap=True):
    return overlap(df, datematch=False, dropOverlap=dropOverlap, anti_datematch=True)


def mag_calib(row):
    '''
    Calculates the calibrated magnitude from the flux.
    :param row: object data (as a pandas Series)
    :return: calibrated magnitude
    '''
    global rows_so_far
    rows_so_far += 1
    try:
        if row['expnum']<552000: # Y3A1
            zpt = zeropoints.ix[(zeropoints['expnum']==row['expnum']) & (zeropoints['ccd']==row['ccd'])]
            mag = zpt['fgcm_zpt'] - 2.5*np.log10(row['flux_auto'])
        elif 552000<row['expnum']<622400:   # Y4
            try:
                if row['ccd'] == -99:
                    mag = 99
                else:
                    zpt = zeropoints_Y4.ix[(zeropoints_Y4['EXPNUM']==row['expnum']) & (zeropoints_Y4['CCDNUM']==row['ccd'])]
                    mag = -zpt['NewZP'].iloc[0] - 2.5 * np.log10(row['flux_auto'])
            except Exception as oops:
                print oops
                print row
    except KeyError as e:
        print 'Error, ', e
#        print row
        mag = 99
    return mag

def drawT_eff(cataloguedf, diffimgdf):
    print "Test"

    cataloguedf['freq'] = cataloguedf.groupby('expnum')['expnum'].transform('count')
    diffimgdf['freq'] = diffimgdf.groupby('expnum')['expnum'].transform('count')

    cat_teff = cataloguedf
    diff_teff = diffimgdf

    cat_teff = cat_teff[~cat_teff.expnum.duplicated(keep='first')]
    diff_teff = diff_teff[~diff_teff.expnum.duplicated(keep='first')]

    t_eff_diff = []
    cat_teff_diff = []

    for index, row in cat_teff.iterrows():
        cat_teff_diff.append(float((all_exps.loc[all_exps['expnum'] == row['expnum']]['t_eff']).to_string().split()[1]))

    cat_teff['t_eff'] = cat_teff_diff

    for index, row in diff_teff.iterrows():
        t_eff_diff.append(float((all_exps.loc[all_exps['expnum'] == row['expnum']]['t_eff']).to_string().split()[1]))

    diff_teff['t_eff'] = t_eff_diff

    t_eff_cat = []
    for index, row in cat_teff.iterrows():
        for i in range(row['freq']):
            t_eff_cat.append(row['t_eff'])

    t_eff_diff = []
    for index, row in diff_teff.iterrows():
        for i in range(row['freq']):
            t_eff_diff.append(row['t_eff'])

    plt.hist([t_eff_cat, t_eff_diff])
    # plt.bar(diff_teff['t_eff'], diff_teff['freq'], color='r')
    plt.savefig('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/catalog_teff.png')


    print "done"

def df_intersection(df1, df2, threshold):

    overlap = pd.DataFrame()
    for index, row in df1.iterrows():
        est_ra = row['ra']
        est_dec = row['dec']
        for index, row in df2.iterrows():
            if (abs(est_ra - row['ra']) < threshold and abs(est_dec - row['dec']) < threshold):
                overlap = overlap.append(row)

    return overlap



def pickle(coadd=None, se=None, catalog=None, diff_img=None):
    if type(coadd) == pd.DataFrame:
        coadd.to_pickle('coadd.pickle')
    if type(se) == pd.DataFrame:
        se.to_pickle('se.pickle')
    if type(catalog) == pd.DataFrame:
        catalog.to_pickle('catalog_df.pickle')
    if type(diff_img) == pd.DataFrame:
        diff_img.to_pickle('diff_img_df.pickle')


def unpickle(coadd=False, se=False, catalog=False, diff_img=False):
    """

    :param coadd:
    :param se:
    :param catalog:
    :param diff_img:
    :return: Tuple of format (coadd, se, catalog, diff_img).
    """
    coadd_df, se_df, catalog_df, diff_img_df = None, None, None, None
    if coadd:
        coadd_df = pd.read_pickle('coadd.pickle')
    if se:
        se_df = pd.read_pickle('se.pickle')
    if catalog:
        catalog_df = pd.read_pickle('catalog_df.pickle')
    if diff_img:
        diff_img_df = pd.read_pickle('diff_img_df.pickle')

    return coadd_df, se_df, catalog_df, diff_img_df


def grab_bleedtrails(db, expnum_ccd_band_set):
    """

    :param expnum_ccd_band_list: Set of (expnum, ccd, band) triplets
    :type db: ea.connect
    :return: dataframe of all bleed trails in those expnum, ccd pairs
    """
    query = "SELECT * FROM PROD.BLEEDTRAIL WHERE "

    for expnum, ccd, band in expnum_ccd_band_set:
        query += "FILENAME LIKE 'D00" + str(expnum) + "_" + str(band) + "_" + str(ccd).zfill(2) + "%' OR "

    query = query[:-4]
    return db.query_to_pandas(query)


def main():
    desoper = ea.connect(section='desoper')
    dessci = ea.connect(section='dessci')

    #y4_exps = pd.read_csv('/Volumes/lsa-gerdes/wsdiff_catalogs/season240/nofakes/wsdiff_season240_Y4_griz_nofakes.csv')

    #y4_exps.to_pickle('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/y4_pickle')

    #s240_df = pd.read_pickle('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/y4_pickle')
    #s240_expnums = s240_df['expnum'].drop_duplicates()

    ra_deg = 315
    dec_deg = -46

    catalogue = 0
    diffimg = 0
    astcount = 0
    diff_asteroids = 0
    cat_asteroids = 0
    catdiffast = 0

    #for exp in s240_expnums:

    ra = ephem.degrees(ra_deg * ephem.pi / 180)
    dec = ephem.degrees(dec_deg * ephem.pi / 180)
    box = 10000.  # arcsec
    season = 240
    expnum = 572126

    exposure_info = desoper.query_to_pandas('SELECT RADEG, DECDEG, DATE_OBS FROM PROD.EXPOSURE WHERE ' +
                                            'EXPNUM = ' + str(expnum))
    ra_exp, dec_exp = np.radians(exposure_info.RADEG[0]), np.radians(exposure_info.DECDEG[0])
    ra_exp_deg, dec_exp_deg = np.degrees(ra), np.degrees(dec)
    date = exposure_info.DATE_OBS[0]  # type: str
    date = date.replace('-', '/')
    date = date.replace('T', ' ')
    date = date.split('.')[0]

    c = SkyCoord(ra=ra_exp * u.rad, dec=dec_exp * u.rad, frame='icrs')

    asteroids = MPChecker(str(c.ra.to_string(unit=u.hour, sep=':')),
                          str(c.dec.to_string(unit=u.degree, sep=':')),
                          date)#, box=box // 60 + 1)

    diff_img_df = get_diffimg_season(season)
    diff_img_df = diff_img_df[diff_img_df['expnum'] == expnum]


    asteroids = asteroids.rename(columns={'mpcRA': 'ra', 'mpcDec': 'dec'})
    asteroids.ra = asteroids.ra.apply(MPCToTNO.ra_convert)
    asteroids.dec = asteroids.dec.apply(MPCToTNO.dec_convert)

    asteroids.ra = asteroids.ra.apply(lambda ra: float(Angle(ra + ' hours').to_string(unit=u.degree, decimal=True)))
    asteroids.dec = asteroids.dec.apply(lambda dec: float(Angle(dec + ' degrees').to_string(unit=u.degree, decimal=True)))

    assert not asteroids.empty
    print "Asteroids Found!"

    se_df = get_SE_detections(desoper, expnum=expnum)
    assert not se_df.empty
    print "Exposure Not Empty!"
    # se_df = get_SE_detections(desoper, ra, dec, box, expnum=expnum)
    coadd_df = get_coadd_cutout(dessci, ra_exp, dec_exp, box)

    catalog_df = get_transient_detections(se_df, coadd_df, 1)
    catalog_df['date'] = catalog_df['date'].apply(lambda date: str(ephem.date(date)))
    # diff_img_df = get_diffimg_cutout(ra, dec, box, season, "../wsdiff_catalogs")
    # diff_img_df = get_diffimg_season(season, "../wsdiff_catalogs")

    pickle(coadd=coadd_df, catalog=catalog_df, diff_img=diff_img_df)
    # _, _, catalog_df, diff_img_df = unpickle(catalog=True, diff_img=True)

    # diff_img_df.to_pickle('diff_img_df.pickle')

    # catalog_df = pd.read_pickle('catalog_df_LARGE.pickle')
    # diff_img_df = pd.read_pickle('diff_img_df_LARGE.pickle')

    # drawT_eff(catalog_df, diff_img_df)

    # diff_img_df['star_like'] = diff_img_df['spread_model'] + 3 * diff_img_df['spreaderr_model']
    #
    # diff_img_df = diff_img_df[(diff_img_df['star_like'] < 0.1) & (diff_img_df['star_like'] > -0.1)]
    # plt.plot(diff_img_df['mag'], diff_img_df['star_like'], linestyle='None', marker='o')
    # plt.xlabel('Magnitude')
    # plt.ylabel("SM + 3 * SEM")
    # plt.title("Starlike Diff_img Detections over Magnitude")
    # plt.tight_layout()
    # plt.savefig("Mag_vs_Starlike_Catalog.png")

    catalog_df['mag'] = catalog_df.apply(lambda row: mag_calib(row), axis=1)

    catalog_df['star_like'] = catalog_df['spread_model'] + 3 * catalog_df['spreaderr_model']
    # catalog_df.to_pickle('catalog_df_LARGE_mag.pickle')
    # catalog_df = pd.read_pickle('catalog_df_LARGE_mag.pickle')
    print 'star_like column created'
    # print catalog_df.head(10)
    catalog_df = catalog_df[(catalog_df['star_like'] > -0.05) & (catalog_df['star_like'] < 0.05)]
    print 'cuts made'

    # catalog_df = catalog_df[(catalog_df['star_like'] > -0.1) & (catalog_df['star_like'] < 0.1)]
    # plt.plot(catalog_df['mag'], catalog_df['star_like'], linestyle='None', marker='o')
    # plt.xlabel('Magnitude')
    # plt.ylabel("SM + 3 * SEM")
    # plt.title("Starlike Catalog Detections over Magnitude")
    # plt.tight_layout()
    # plt.savefig("Mag_vs_Starlike_Catalog.png")

    # overlap_coadd = overlap(diff_img_df, coadd_df, datematch=False, dropOverlap=True)
    selfoverlap_df = self_overlap(catalog_df)

    # overlap_df = overlap(diff_img_df, catalog_df)
    try:
        asteroid_catalog = overlap(catalog_df, asteroids, dropOverlap=False, datematch=False, threshold=4)
        asteroid_diffimg = overlap(diff_img_df, asteroids, dropOverlap=False, datematch=False, threshold=4)
        astcat_astdiff = df_intersection(asteroid_catalog, asteroid_diffimg, 1 / 3600)
        # = overlap(asteroid_catalog, asteroid_diffimg, dropOverlap=True, datematch=False,
        #                        threshold=1)
    except ValueError:  # raised if `y` is empty.
        pass

    print "Catalog Only: ", len(catalog_df)
    print "Diff Img Only: ", len(diff_img_df)
    # print "Both: ", len(overlap_df)
    # print "Coadd Matches: ", len(overlap_coadd)
    print "Self Overlap: ", len(selfoverlap_df)
    print "Asteroids: ", len(asteroids)
    print "Asteroids Found via Cat: ", len(asteroid_catalog)
    print "Asteroids Found via DiffImg: ", len(asteroid_diffimg)

    # expnum_ccd_band_set = set()
    #
    # for index, row in catalog_df.iterrows():
    #     expnum_ccd_band_set.add((row['expnum'], row['ccd'], row['band']))
    #
    # bleed_df = grab_bleedtrails(desoper, expnum_ccd_band_set)
    # print bleed_df.head()
    # print len(bleed_df)

    cat_detect = plt.scatter(catalog_df['ra'], catalog_df['dec'], linestyle='None', color='g', marker='.', s=0.5)
    diff_detect = plt.scatter(diff_img_df['ra'], diff_img_df['dec'], linestyle='None',color='c', marker='.', s=0.5) #see if any fall in coadd already
    asteroid_plt = plt.scatter(asteroids['ra'], asteroids['dec'], linestyle='None', color='r', marker='o', alpha=0.4, label='asteroids')
    cat_asteroid = plt.scatter(asteroid_catalog['ra'], asteroid_catalog['dec'], linestyle='None', color='blue', marker='o', label = 'catalogue')
    diff_asteroid = plt.scatter(asteroid_diffimg['ra'], asteroid_diffimg['dec'], linestyle='None', color='black', marker='o', label = 'diffimg')
    cat_diff_asteroid = plt.scatter(astcat_astdiff['ra'], astcat_astdiff['dec'], linestyle='None', color='purple', marker='o', label = 'catdiff')

    plt.legend((cat_detect, diff_detect, asteroid_plt, cat_asteroid, diff_asteroid, cat_diff_asteroid),
    ('Catalogue', 'DiffImg', 'Asteroids', 'Catalogue Detections', 'Diffimg Detections', 'Catalogue and Diffimg'),
        scatterpoints=1,
        loc='lower left',
        ncol=3,
        fontsize=8)

    title = "Exposure " + str(expnum)
    plt.title(title)
    plt.savefig('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/' + title)
    print "Saved as " + title
    # plt.savefig('detectionsLarge_0.05_cut_Overlap_Cat_Removed.png')
    #
    # plt.figure()
    # plt.hist(catalog_df['star_like'], 50)
    # plt.savefig('hist_starlike_Large_0.05_cut_Overlap_Cat_Removed.png')
    #
    # plt.figure()
    # plt.hist(catalog_df['mag'], 50, range=(20, 30))
    # plt.savefig('hist_mag_Large_0.05_cut_Overlap_Cat_Removed.png')
    #
    # plt.savefig('detectionsLargeCoadd_DI_Removed_Overlap_Cat_Removed_SpreadCut.png')
    # plt.savefig("exposure" + str(expnum) + "_zoom.png")
    expnumCCD_list = []
    catalog_list = []
    diff_img_list = []
    overlap_list = []
    coadd_list = []
    asteroid_list = []

    for index, row in catalog_df.iterrows():
        expnumCCD_list.append((row['expnum'], row['ccd']))
        catalog_list.append((row['ra'], row['dec'], row['expnum'], row['ccd']))

    for index, row in diff_img_df.iterrows():
        expnumCCD_list.append((row['expnum'], row['ccd']))
        diff_img_list.append((row['ra'], row['dec'], row['expnum'], row['ccd']))

    for index, row in asteroid_catalog.iterrows():
        expnumCCD_list.append((row['expnum'], row['ccd']))
        overlap_list.append((row['ra'], row['dec'], row['expnum'], row['ccd']))

    # for index, row in coadd_df.iterrows():
    #     coadd_list.append((row['ra'], row['dec']))

    for index, row in asteroids.iterrows():
        asteroid_list.append((row['ra'], row['dec']))
    os.chdir('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search')
    Find_imgs.findImgs(expnumCCD_list, catalog_list, diff_img_list, overlap_list, ra_deg,
                       dec_deg, box * 2, coadd_list=coadd_list, asteroid_list=asteroid_list,
                       keep_fits=False, only_asteroids=True)


if __name__ == '__main__':
    main()


# def codeArchive():

    # fig_cat, ax_cat = plt.subplots(2, 2, sharex='col', sharey='row')
    # fig_diff, ax_diff = plt.subplots(2, 2, sharex='col', sharey='row')
    # fig_combined, ax_combined = plt.subplots(2, 2, sharex='col', sharey='row')
    #
    # catalog_df = get_transient_detections(se_df, coadd_df, 1)
    # catalog_df['date'] = catalog_df['date'].apply(lambda date: str(ephem.date(date)))
    #
    # diff_img_df = get_diffimg_cutout(ra, dec, box, season)
    #
    # catalog_df.to_pickle('catalog_df.pickle')
    # diff_img_df.to_pickle('diff_img_df.pickle')
    # diff_img_df.to_csv('diff_img_df.csv')
    #
    # overlap_df = overlap(diff_img_df, catalog_df)
    #
    # for i in [1, 5, 10, 15]:
    #     catalog_df = get_transient_detections(se_df, coadd_df, i)
    #     catalog_df['date'] = catalog_df['date'].apply(lambda date: str(ephem.date(date)))
    #
    #     diff_img_df = get_diffimg_cutout(ra, dec, box, season)
    #
    #     catalog_df.to_pickle('catalog_df.pickle')
    #     diff_img_df.to_pickle('diff_img_df.pickle')
    #
    #     overlap_df = overlap(diff_img_df, catalog_df)
    #
    #     print "Catalog Only: ", len(catalog_df)
    #     print "Diff Img Only: ", len(diff_img_df)
    #     print "Both: ", len(overlap_df)
    #
    #     print diff_img_df
    #
    #     if i == 1:
    #         j = 0
    #         k = 0
    #     elif i == 5:
    #         j = 0
    #         k = 1
    #     elif i == 10:
    #         j = 1
    #         k = 0
    #     else:
    #         j = 1
    #         k = 1
    #
    #     l_cat, = ax_cat[j, k].plot(catalog_df['ra'], catalog_df['dec'], linestyle='None', marker='o', color='g')
    #     l_c_overlap, = ax_cat[j, k].plot(overlap_df['ra'], overlap_df['dec'], linestyle='None', marker='o', color='b')
    #     ax_cat[j, k].ticklabel_format(useOffset=False)
    #     ax_cat[j, k].grid()
    #     ax_cat[j, k].set_title('Threshold = ' + str(i) + '"')
    #     if k == 0:
    #         ax_cat[j, k].set_ylabel('DEC\n(deg)')
    #     if j == 1:
    #         ax_cat[j, k].set_xlabel('RA\n(deg)')
    #     # plt.savefig('catalog' + str(i) + '.png')
    #
    #     l_diff, = ax_diff[j, k].plot(diff_img_df['ra'], diff_img_df['dec'], linestyle='None', marker='o', color='r')
    #     l_d_overlap, = ax_diff[j, k].plot(overlap_df['ra'], overlap_df['dec'], linestyle='None', marker='o', color='b')
    #     # ax_diff[j, k].ticklabel_format(useOffset=False)
    #     ax_diff[j, k].grid()
    #     ax_diff[j, k].set_title('Threshold = ' + str(i) + '"')
    #     if k == 0:
    #         ax_diff[j, k].set_ylabel('DEC\n(deg)')
    #     if j == 1:
    #         ax_diff[j, k].set_xlabel('RA\n(deg)')
    #     # plt.savefig('combined.png')
    #
    #     l_combined_cat, = ax_combined[j, k].plot(catalog_df['ra'], catalog_df['dec'],
    #                                             linestyle='None', marker='o', color='g')
    #     l_combined_diff, = ax_combined[j, k].plot(overlap_df['ra'], overlap_df['dec'],
    #                                              linestyle='None', marker='o', color='b', alpha=0.5)
    #     l_combined_overlap, = ax_combined[j, k].plot(diff_img_df['ra'], diff_img_df['dec'],
    #                                                 linestyle='None', marker='o', color='r')
    #     # ax_combined[j, k].ticklabel_format(useOffset=False)
    #     ax_combined[j, k].grid()
    #     ax_combined[j, k].set_title('Threshold = ' + str(i) + '"')
    #     if k == 0:
    #         ax_combined[j, k].set_ylabel('DEC\n(deg)')
    #     if j == 1:
    #         ax_combined[j, k].set_xlabel('RA\n(deg)')
    #
    #     # plt.savefig('diff_img' + str(i) + '.png')
    #     # print ra
    #     # print dec
    #
    # fig_cat.legend((l_cat, l_c_overlap), ('Catalog Only', 'Both'), 'best')
    # fig_diff.legend((l_diff, l_d_overlap), ('Diff Img Only', 'Both'), 'best')
    # fig_combined.legend((l_combined_cat, l_combined_diff, l_combined_overlap),
    #                     ('Catalog Only', 'Diff Img Only', 'Both'), 'best')
    #
    # fig_cat.tight_layout()
    # fig_diff.tight_layout()
    # fig_combined.tight_layout()
    #
    # fig_cat.savefig('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/catalog.png')
    # fig_diff.savefig('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/diff_img.png')
    # fig_combined.savefig('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/combined.png')