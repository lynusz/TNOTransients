from __future__ import division

import easyaccess as ea
import ephem
import numpy as np
import pandas as pd
from KBO import compute_chip
from extendOrbit import getOrbit
from linkmap import build_kdtree
import matplotlib.pyplot as plt

#zeropoints = pd.read_csv('fgcm_zeropoints_v2_0.csv')
#zeropoints_Y4 = pd.read_csv('Y4N_zeropoints_03.09.2017.csv')
all_exps = pd.read_csv('exposures.csv')


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


def get_SE_detections(connection, ra, dec, box):
    '''
    Gets single-epoch detections in all exposures from year for within a search box.

    :param connection: easyaccess database connection. Currently should be to DESOPER if year='Y4' and to DESSCI for earlier years.
    :param expnum: exposure number to search
    :param ra: RA of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param dec: DEC of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param box: width of search region (+/- from center) in arcsec (i.e. the units returned by predict_pos)
    :return: dataframe of single-epoch detections.
    '''

    query_Y4 = r"""
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
        and
    """

    ra_deg, dec_deg = np.degrees(ra), np.degrees(dec)
    box_deg = box / 3600
    query = (query_Y4 + "o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg)
             + "+" + str(box_deg) + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and "
             + str(dec_deg) + "+" + str(box_deg))

    result = connection.query_to_pandas(query)
    result.columns = map(str.lower, result.columns)  # convert column names to lowercase
    result.rename(columns={'date_obs': 'date'}, inplace=True)
    merged = result.merge(all_exps[['expnum', 'nite', 'band', 'date', 't_eff', 'fwhm_asec', 'exptime']], how='inner',
                          on=['expnum', 'band'])
    merged = merged.loc[merged['exptime'] > 80.0]
    return merged


def get_transient_detections(df_SE, df_coadd):
    '''
    Given dataframes of single-epoch and coadd detections, removes SE detections within 1 arcsec of a coadd detection 
    and returns the remaining SE detections as possible transients. This is a catalog-level attempt at difference-imaging.
    :param df_SE: dataframe of SE detections (the search sample)
    :param df_coadd: dataframe of coadd detections (the template sample)
    :return: dataframe of unmatched SE detections
    '''
    if len(df_SE) and len(df_coadd):
        tree_coadd = build_kdtree(df_coadd)
        tree_SE = build_kdtree(df_SE)
        near_list = tree_SE.query_ball_tree(tree_coadd, 1 / 3600)  # require 1 arcsec match to a coadd object
        transients = [i for i in range(len(near_list)) if near_list[i] == []]
        return df_SE.ix[transients]
    else:
        return pd.DataFrame()

def get_diffimg_cutout(ra, dec, box, season):
    csv_file = ('/Volumes/lsa-gerdes/wsdiff_catalogs/season' + str(season) +
                '/nofakes/wsdiff_season' + str(season) + '_Y4_griz_nofakes.csv')
    raw_diff_img = pd.read_csv(csv_file)
    ra_deg, dec_deg = np.degrees(ra), np.degrees(dec)

    box_deg = box / 3600
    cut = (((float(ra_deg) - box_deg <= raw_diff_img['ra']) & (raw_diff_img['ra'] <= float(ra_deg) + box_deg)) &
           ((float(dec_deg) - box_deg <= raw_diff_img['dec']) & (raw_diff_img['dec'] <= float(dec_deg) + box_deg)))
    return raw_diff_img[cut]


def main():
    desoper = ea.connect(section='desoper')
    dessci = ea.connect(section='dessci')

    ra = ephem.degrees(310.0 * ephem.pi / 180)
    dec = ephem.degrees(-51.0 * ephem.pi / 180)
    box = 1000.
    season = 240

    result = get_transient_detections(get_SE_detections(desoper, ra, dec, box),
                                      get_coadd_cutout(dessci, ra, dec, box))

    diff_img = get_diffimg_cutout(ra, dec, box, season)

    print len(result)
    print len(diff_img)

    print result.head()
    print diff_img.head()

    fig, ax = plt.subplots()
    ax.plot(result['dec'], result['ra'], linestyle='None', marker='o')
    ax.ticklabel_format(useOffset=False)
    plt.savefig('catalog.png')

    fig, ax = plt.subplots()
    ax.plot(diff_img['dec'], diff_img['ra'], linestyle='None', marker='o')
    ax.ticklabel_format(useOffset=False)
    plt.savefig('diff_img.png')

if __name__ == '__main__':
    main()