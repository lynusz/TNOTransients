from __future__ import division

import easyaccess as ea
import ephem
import numpy as np
import pandas as pd
from linkmap import build_kdtree
import matplotlib.pyplot as plt

#zeropoints = pd.read_csv('fgcm_zeropoints_v2_0.csv')
#zeropoints_Y4 = pd.read_csv('Y4N_zeropoints_03.09.2017.csv')
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
        and T_EFF > 0.3
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


def get_transient_detections(df_SE, df_coadd, threshold):
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
        near_list = tree_SE.query_ball_tree(tree_coadd, threshold / 3600)  # require 1 arcsec match to a coadd object
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


def overlap(df1, df2):
    """
    Calculates the overlap between two dataframes on ra and dec and puts those rows in a third dataframe.
    Only puts one row per overlap. Removes those rows from df1 and df2.
    :param df1: 
    :type df1: pd.DataFrame
    :param df2: 
    :type df2: pd.DataFrame
    :return: 
    """
    if len(df1) > len(df2):
        df_long = df1
        df_short = df2
    else:
        df_long = df2
        df_short = df1

    df_hash = {}
    for index, row in df_short.iterrows():
        df_hash[(round(row['ra'], 2), round(row['dec'], 2), row['date'])] = index

    overlap_index_short = []
    overlap_index_long = []

    diff_img_buds = pd.DataFrame()

    for index, row in df_long.iterrows():
        if (round(row['ra'], 2), round(row['dec'], 2), row['date']) in df_hash:
            overlap_index_short.append(df_hash[(round(row['ra'], 2), round(row['dec'], 2), row['date'])])
            overlap_index_long.append(index)
            if row['expnum'] == 581186:
                diff_img_buds = diff_img_buds.append(row)

    diff_img_buds.to_csv("BuddiesOf581186.csv")

    df_overlap = pd.DataFrame(df_short.ix[overlap_index_short])
    df_short.drop(overlap_index_short, inplace=True)
    df_long.drop(overlap_index_long, inplace=True)
    return df_overlap


def main():
    desoper = ea.connect(section='desoper')
    dessci = ea.connect(section='dessci')

    ra = ephem.degrees(310.0 * ephem.pi / 180)
    dec = ephem.degrees(-51.0 * ephem.pi / 180)
    box = 100.
    season = 240
    se_df = get_SE_detections(desoper, ra, dec, box)
    coadd_df = get_coadd_cutout(dessci, ra, dec, box)

    catalog_df = get_transient_detections(se_df, coadd_df, 1)
    catalog_df['date'] = catalog_df['date'].apply(lambda date: str(ephem.date(date)))

    diff_img_df = get_diffimg_cutout(ra, dec, box, season)

    catalog_df.to_pickle('catalog_df.pickle')
    diff_img_df.to_pickle('diff_img_df.pickle')

    fig_cat, ax_cat = plt.subplots(2, 2, sharex='col', sharey='row')
    fig_diff, ax_diff = plt.subplots(2, 2, sharex='col', sharey='row')
    fig_combined, ax_combined = plt.subplots(2, 2, sharex='col', sharey='row')

    catalog_df = get_transient_detections(se_df, coadd_df, 1)
    catalog_df['date'] = catalog_df['date'].apply(lambda date: str(ephem.date(date)))

    diff_img_df = get_diffimg_cutout(ra, dec, box, season)

    catalog_df.to_pickle('catalog_df.pickle')
    diff_img_df.to_pickle('diff_img_df.pickle')
    diff_img_df.to_csv('diff_img_df.csv')

    overlap_df = overlap(diff_img_df, catalog_df)

    for i in [1, 5, 10, 15]:
        catalog_df = get_transient_detections(se_df, coadd_df, i)
        catalog_df['date'] = catalog_df['date'].apply(lambda date: str(ephem.date(date)))

        diff_img_df = get_diffimg_cutout(ra, dec, box, season)

        catalog_df.to_pickle('catalog_df.pickle')
        diff_img_df.to_pickle('diff_img_df.pickle')

        overlap_df = overlap(diff_img_df, catalog_df)

        print "Catalog Only: ", len(catalog_df)
        print "Diff Img Only: ", len(diff_img_df)
        print "Both: ", len(overlap_df)

        print diff_img_df

        if i == 1:
            j = 0
            k = 0
        elif i == 5:
            j = 0
            k = 1
        elif i == 10:
            j = 1
            k = 0
        else:
            j = 1
            k = 1

        l_cat, = ax_cat[j, k].plot(catalog_df['ra'], catalog_df['dec'], linestyle='None', marker='o', color='g')
        l_c_overlap, = ax_cat[j, k].plot(overlap_df['ra'], overlap_df['dec'], linestyle='None', marker='o', color='b')
        ax_cat[j, k].ticklabel_format(useOffset=False)
        ax_cat[j, k].grid()
        ax_cat[j, k].set_title('Threshold = ' + str(i) + '"')
        if k == 0:
            ax_cat[j, k].set_ylabel('DEC\n(deg)')
        if j == 1:
            ax_cat[j, k].set_xlabel('RA\n(deg)')
        # plt.savefig('catalog' + str(i) + '.png')

        l_diff, = ax_diff[j, k].plot(diff_img_df['ra'], diff_img_df['dec'], linestyle='None', marker='o', color='r')
        l_d_overlap, = ax_diff[j, k].plot(overlap_df['ra'], overlap_df['dec'], linestyle='None', marker='o', color='b')
        # ax_diff[j, k].ticklabel_format(useOffset=False)
        ax_diff[j, k].grid()
        ax_diff[j, k].set_title('Threshold = ' + str(i) + '"')
        if k == 0:
            ax_diff[j, k].set_ylabel('DEC\n(deg)')
        if j == 1:
            ax_diff[j, k].set_xlabel('RA\n(deg)')
        # plt.savefig('combined.png')

        l_combined_cat, = ax_combined[j, k].plot(catalog_df['ra'], catalog_df['dec'],
                                                linestyle='None', marker='o', color='g')
        l_combined_diff, = ax_combined[j, k].plot(overlap_df['ra'], overlap_df['dec'],
                                                 linestyle='None', marker='o', color='b', alpha=0.5)
        l_combined_overlap, = ax_combined[j, k].plot(diff_img_df['ra'], diff_img_df['dec'],
                                                    linestyle='None', marker='o', color='r')
        # ax_combined[j, k].ticklabel_format(useOffset=False)
        ax_combined[j, k].grid()
        ax_combined[j, k].set_title('Threshold = ' + str(i) + '"')
        if k == 0:
            ax_combined[j, k].set_ylabel('DEC\n(deg)')
        if j == 1:
            ax_combined[j, k].set_xlabel('RA\n(deg)')

        # plt.savefig('diff_img' + str(i) + '.png')
        # print ra
        # print dec

    fig_cat.legend((l_cat, l_c_overlap), ('Catalog Only', 'Both'), 'best')
    fig_diff.legend((l_diff, l_d_overlap), ('Diff Img Only', 'Both'), 'best')
    fig_combined.legend((l_combined_cat, l_combined_diff, l_combined_overlap),
                        ('Catalog Only', 'Diff Img Only', 'Both'), 'best')

    fig_cat.tight_layout()
    fig_diff.tight_layout()
    fig_combined.tight_layout()

    fig_cat.savefig('catalog.png')
    fig_diff.savefig('diff_img.png')
    fig_combined.savefig('combined.png')

if __name__ == '__main__':
    main()
