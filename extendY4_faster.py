#!/usr/bin/env/python

from __future__ import division

import easyaccess as ea
import ephem
import numpy as np
import pandas as pd
from KBO import compute_chip
from extendOrbit import getOrbit
from linkmap import build_kdtree

zeropoints = pd.read_csv('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/fgcm_zeropoints_v2_0.csv')
zeropoints_Y4 = pd.read_csv('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/Y4N_zeropoints_03.09.2017.csv')
all_exps = pd.read_csv('/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/exposures.csv')


def possibleExposures(orbit, year='Y4'):
    '''
    Returns a list of exposures near the target orbit in a given year. Assumes the object moves slowly enough that
    it won't wander more than 10 degrees from its location on 2014/01/01. Returns a  preliminary exposure 
    list from which we do a more careful search. 
    :param orbit: orbit from pyOrbfit
    :param year: DES observing year to search, e.g. 'Y4'
    :return: dataframe of exposures to be searched.
    '''
    pos_ref = orbit.predict_pos(ephem.date('2014/01/01'))
    ra_ref, dec_ref = pos_ref['ra'], pos_ref['dec']
    assert year in ['all', 'SV', 'Y1', 'Y2', 'Y3', 'Y4', 'Y5']
    if year == 'all':
        search_exps = all_exps
    elif year == 'SV':
        search_exps = all_exps.loc[
            all_exps.apply(lambda row: 20121102 < row['nite'] < 20130220 and row['band'] != 'Y', axis=1)]
    elif year == 'Y1':
        search_exps = all_exps.loc[
            all_exps.apply(lambda row: 20130810 < row['nite'] < 20140220 and row['band'] != 'Y', axis=1)]
    elif year == 'Y2':
        search_exps = all_exps.loc[
            all_exps.apply(lambda row: 20140810 < row['nite'] < 20150220 and row['band'] != 'Y', axis=1)]
    elif year == 'Y3':
        search_exps = all_exps.loc[
            all_exps.apply(lambda row: 20150810 < row['nite'] < 20160220 and row['band'] != 'Y', axis=1)]
    elif year == 'Y4':
        search_exps = all_exps.loc[
            all_exps.apply(lambda row: 20160810 < row['nite'] < 20170220 and row['band'] != 'Y', axis=1)]
    elif year == 'Y5':
        search_exps = all_exps.loc[
            all_exps.apply(lambda row: 20170810 < row['nite'] < 20180220 and row['band'] != 'Y', axis=1)]
    else:
        print 'Oops'
        search_exps = None
    try:
        near_exps = search_exps.loc[search_exps.apply(lambda row:
                                                      ephem.separation(
                                                          (ephem.hours(row['ra']), ephem.degrees(row['dec'])),
                                                          (ra_ref, dec_ref)) < ephem.degrees('10:00:00'), axis=1)]
    except:
        near_exps = None
    in_exposures = []
    in_ccd = []
    if near_exps is not None:
        for ind, row in near_exps.iterrows():
            pos = orbit.predict_pos(row['date'])
            ra, dec = ephem.hours(pos['ra']), ephem.degrees(pos['dec'])
            ccdname, ccdnum = compute_chip(ra, dec, row['ra'], row['dec'])
            if ccdnum > -99 and ccdnum not in [2, 31, 61]:
                in_exposures.append(row['expnum'])
                in_ccd.append(ccdnum)
        near_exps = near_exps.loc[near_exps['expnum'].isin(in_exposures)]
        near_exps['ccd'] = pd.Series(in_ccd, index=near_exps.index)
    return near_exps


def hasDetection(expnum, cand_df):
    return True if expnum in cand_df['expnum'].values else False


def search_exposures(cand_df, orbit, year='Y4', teff_min=0.3):
    '''
    Get list of exposures to be searched from the desired opposition. These are exposures that the object may have appeared in but which
    do not already have a detection. A minimum t_effective cut is imposed.
    param cand_df: dataframe of observations of the candidate
    param orbit: orbit (from pyOrbfit) corresponding to the observations in cand_df
    param year: DES observing year to search, e.g. 'Y4'
    param teff_min: minimum t_eff of exposures to consider
    '''
    search_exps = possibleExposures(orbit, year=year)
    print len(search_exps)
    try:
        search_exps['hasDetection'] = pd.Series([hasDetection(e, cand_df) for e in search_exps['expnum'].values],
                                                index=search_exps.index)
        return search_exps.loc[(search_exps['hasDetection'] == False) & (search_exps['t_eff'] > teff_min)]
    except:
        return pd.DataFrame()


'''
def get_coadd_cutout(connection, df):


    #Gets Y3A1_COADD_OBJECTs in the given search box.
    #:param connection: easyaccess database connection to DESSCI
    #:param ra: RA of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    #:param dec: DEC of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    #:param box: width of search region (+/- from center) in arcsec (i.e. the units returned by predict_pos)
    #:return: dataframe of Y3A1_COADD_OBJECTs.
    

    query_Y3A1_coadd = r""" select o.ra, o.dec, o.mag_auto_g, o.mag_auto_r, o.mag_auto_i, o.mag_auto_z, 
    o.nepochs_g, o.nepochs_r, o.nepochs_i, o.nepochs_z, o.nepochs_Y,
    o.spread_model_g, o.spread_model_r, o.spread_model_i, o.spread_model_z, 
    o.spreaderr_model_g, o.spreaderr_model_r, o.spreaderr_model_i, o.spreaderr_model_z,
    o.flags_g, o.flags_r, o.flags_i, o.flags_z
    from Y3A1_COADD_OBJECT_SUMMARY o where (
    """

    query = ""
    for i, (index, row) in enumerate(df.iterrows()):

        ra_deg, dec_deg = np.degrees(row['ra_pred']), np.degrees(row['dec_pred'])
        box_deg = row['err_a'] / 3600
        query += "( o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg) + "+" + str(
            box_deg) \
                + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and " + str(dec_deg) + "+" + str(box_deg) + ')'

        if i != len(df) - 1:
            query += ' or '
        else:
            pass

    query = query_Y3A1_coadd + query + ')'

    result = connection.query_to_pandas(query)
    result.columns = map(str.lower, result.columns)  # convert column names to lowercase
    result.rename(columns={'date_obs': 'date'}, inplace=True)
    return result
'''

def get_coadd_cutout(connection, ra, dec, box):

    '''
    Gets Y3A1_COADD_OBJECTs in the given search box.
    :param connection: easyaccess database connection to DESSCI
    :param ra: RA of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param dec: DEC of search center. Should be an ephem.angle object (i.e. expressed in *** radians ***)
    :param box: width of search region (+/- from center) in arcsec (i.e. the units returned by predict_pos)
    :return: dataframe of Y3A1_COADD_OBJECTs.
    '''

    query_Y3A1_coadd = r""" select o.ra, o.dec, o.mag_auto_g, o.mag_auto_r, o.mag_auto_i, o.mag_auto_z, 
    o.nepochs_g, o.nepochs_r, o.nepochs_i, o.nepochs_z, o.nepochs_Y,
    o.spread_model_g, o.spread_model_r, o.spread_model_i, o.spread_model_z, 
    o.spreaderr_model_g, o.spreaderr_model_r, o.spreaderr_model_i, o.spreaderr_model_z,
    o.flags_g, o.flags_r, o.flags_i, o.flags_z
    from Y3A1_COADD_OBJECT_SUMMARY o where 
    """

    ra_deg, dec_deg = np.degrees(ra), np.degrees(dec)
    box_deg = box / 3600
    query = query_Y3A1_coadd + "o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg) + "+" + str(
        box_deg) \
            + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and " + str(dec_deg) + "+" + str(box_deg)

    print query
    result = connection.query_to_pandas(query)
    result.columns = map(str.lower, result.columns)  # convert column names to lowercase
    result.rename(columns={'date_obs': 'date'}, inplace=True)
    return result


def get_SE_detections(connection, df):
    '''
    Gets single-epoch detections in the given exposure number and search box.
    
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
        and (
    """
    query_Y3A1 = r"""
    select o.expnum, o.filename, cast(o.band as VARCHAR(1)) as BAND,
    o.object_number, o.ra, o.dec, o.flux_auto, o.fluxerr_auto, o.spread_model, o.spreaderr_model, o.imaflags_iso
    from Y3A1_FINALCUT_OBJECT o where o.band !='Y' 
    and ( 
    """


    query_y4_appended = "("
    query_Y3A1_appended = "("

    y4_obs = 0
    y3a1_obs = 0

    for i, (index, row) in enumerate(df.iterrows()):

        box_deg = row['err_a'] / 3600

        ra_deg, dec_deg = np.degrees(row['ra_pred']), np.degrees(row['dec_pred'])

        if row['expnum'] > 552000:  # Y4 exposure, use Y4N_FIRSTCUT. This will need to be updated when Y5 happens.
            y4_obs = 1
            query_y4_appended += "(o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg) + "+" + str(
                box_deg) \
                    + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and " + str(dec_deg) + "+" + str(
                box_deg) \
                    + " and c.expnum=" + str(row['expnum']) + ')'

            if i != len(df) - 1:
                query_y4_appended += ' or '
            else:
                query_y4_appended += ')'

        else:  # Y1, Y2, or Y3... use Y3A1_FINALCUT_OBJECT
            y3a1_obs = 1
            query_Y3A1_appended = "(o.ra between " + str(ra_deg) + "-" + str(box_deg) + " and " + str(ra_deg) + "+" + str(
                box_deg) \
                    + " and o.dec between " + str(dec_deg) + "-" + str(box_deg) + " and " + str(dec_deg) + "+" + str(
                box_deg) \
                    + "and o.expnum=" + str(row['expnum']) + ')'

            if index != len(df) - 1:
                query_Y3A1_appended += ' or '
            else:
                query_Y3A1_appended += ')'

    query_Y4 = query_Y4 + query_y4_appended + ')'
    query_Y3A1 = query_Y3A1 + query_Y3A1_appended + ')'

    #initialize to empty df
    result = result_Y3A1 = result_y4 = pd.DataFrame()

    if y4_obs:
        result_y4 = connection.query_to_pandas(query_Y4)
        result_y4.columns = map(str.lower, result_y4.columns)  # convert column names to lowercase
        result_y4.rename(columns={'date_obs': 'date'}, inplace=True)
    if y3a1_obs:
        result_Y3A1 = connection.query_to_pandas(query_Y3A1)
        result_Y3A1.columns = map(str.lower, result_Y3A1.columns)  # convert column names to lowercase
        result_Y3A1.rename(columns={'date_obs': 'date'}, inplace=True)

    bool1 = result_y4.empty
    bool2 = result_Y3A1.empty
    #if not (result_y4.empty or result_Y3A1.empty):
    result = result_y4.append(result_Y3A1, ignore_index = True)
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
        near_list = tree_SE.query_ball_tree(tree_coadd, 1/3600)   # require 1 arcsec match to a coadd object
        transients = [i for i in range(len(near_list)) if near_list[i]==[]]
        return df_SE.ix[transients]
    else:
        return pd.DataFrame()


def point_in_ellipse(x, y, x0, y0, a, b, phi):
    '''
    Tests to see if the point (x,y) is inside the ellipse with center at (x0,y0), major and minor axes a,b, and 
    rotation angle phi as measured from the positive x-axis. 

    x,y,x0,y0 are assumed to be in DECIMAL DEGREES
    a,b are assumed to be in ARCSEC (as returned by orbit.predict_pos)
    phi is assumed to be in DEGREES

    Note 1: The position angle (PA) returned by orbit.predict_pos() is in degrees east of north. Use phi=90-PA
    to get the correct orientation here. 

    Note 2: This will misbehave if x and x0 (the RA coordinate) are on opposite sides of x=0 (or 360). Need to 
    fix this.
    '''
    phi = np.radians(phi)  # degrees to radians
    a = a / 3600  # arcsec to degrees
    b = b / 3600
    dist_norm = (np.cos(phi) * (x - x0) + np.sin(phi) * (y - y0)) ** 2 / a ** 2 + (np.sin(phi) * (x - x0) - np.cos(
        phi) * (y - y0)) ** 2 / b ** 2
    return True if dist_norm <= 1 else False

def mag_calib(row):
    '''
    Calculates the calibrated magnitude from the flux.
    :param row: object data (as a pandas Series)
    :return: calibrated magnitude
    '''
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
                    mag = -zpt['NewZP'] - 2.5*np.log10(row['flux_auto'])
            except Exception as oops:
                print oops
                print row
    except KeyError as e:
        print 'Error, ', e
#        print row
        mag = 99
    return mag


def get_ccd(row, exp_df):
    '''
    Returns the ccd containing a given SE observations.
    :param row: pandas Series containing observation data
    :param exp_df: dataframe of exposure info.
    :return: ccd number
    '''
    exp_ra = exp_df[exp_df['expnum']==row['expnum']]['ra'].values[0]
    exp_dec = exp_df[exp_df['expnum']==row['expnum']]['dec'].values[0]
    ccdname, ccdnum = compute_chip(np.radians(row['ra']), np.radians(row['dec']),exp_ra, exp_dec)
    return ccdnum


def find_observations(search_exps, cand_df, orbit, year='Y4', nsigma=3):
    '''
    Finds SE transients in the search exposures from the given year that fall inside the error ellipse of the candidate.
    :param search_exps: dataframe of exposures to be searched
    :param cand_df: dataframe of observations of the candidate
    :param orbit: the candidate's corresponding orbit
    :param year: DES observing year of the search exposures, e.g. 'Y4'
    :param nsigma: multiplier for error ellipse dimensions
    :return: dataframe of transients
    '''
    conn_desoper = ea.connect(section='desoper')
    conn_dessci = ea.connect(section='dessci')
    matched_obs = pd.DataFrame()
    transients_inside = pd.DataFrame()
#    search_exps = search_exps[search_exps['expnum']!=374572]  # to exclude particular exposures

    df = pd.DataFrame(columns=['expnum', 'ra_pred', 'dec_pred', 'err_a', 'err_b', 'err_PA'])
    full_df_coadd = pd.DataFrame()

    for ind, exp in search_exps.iterrows():
        print 'Searching exposure... ', exp['expnum'], ephem.date(exp['date']), exp['band'], exp['t_eff'],'...',
        pos_pred = orbit.predict_pos(exp['date'])
        ra_pred, dec_pred, err_a, err_b, err_PA = pos_pred['ra'], pos_pred['dec'], pos_pred['err']['a'],\
            pos_pred['err']['b'], pos_pred['err']['PA']
        print 'RA_PRED, DEC_PRED, ERR_A, ERR_B, PA = ', ra_pred, dec_pred, err_a, err_b, err_PA
        if err_a<15.0: err_a=15.0
        if err_b<15.0: err_b=15.0
        coadd_obs = get_coadd_cutout(conn_dessci, ra_pred, dec_pred, nsigma * err_a)
        full_df_coadd = full_df_coadd.append(coadd_obs)

        df.loc[ind] = [exp['expnum'], ra_pred, dec_pred, nsigma*err_a, nsigma*err_b, err_PA]

    if year=='Y4':
        se_obs = get_SE_detections(conn_desoper, df)
    else:
        se_obs = get_SE_detections(conn_dessci, df)

    #coadd_obs = get_coadd_cutout(conn_dessci, df)
    transients = get_transient_detections(se_obs, full_df_coadd)



    if len(transients):
        transients['delta_chisq'] = transients.apply(lambda row: getOrbit(cand_df, row).chisq-orbit.chisq,axis=1)
#            transients['mag'] = good_matches.apply(lambda row: mag_calib(row), axis=1)
        transients['mag'] = 99.    # temp workaround
        transients['ml_score']=10.0   # placeholder
        transients['fakeid']=0
        transients.rename(columns={'object_number': 'objid'}, inplace=True)
# These transients fall inside the error ellipse

        transients_inside = transients.loc[transients.apply(lambda row:
                                                            point_in_ellipse(row['ra'],
                                                                             row['dec'],
                                                                             np.degrees(df.loc[df['expnum'] == row['expnum'], 'ra_pred'].iloc[0]),
                                                                             np.degrees(df.loc[df['expnum'] == row['expnum'], 'dec_pred'].iloc[0]),
                                                                             df.loc[df['expnum'] == row['expnum'], 'err_a'].iloc[0],
                                                                             df.loc[df['expnum'] == row['expnum'], 'err_b'].iloc[0],
                                                                             90 - df.loc[df['expnum'] == row['expnum'], 'err_PA'].iloc[0]),
                                                                             axis=1)]

            #transients_inside = transients.loc[transients.apply(lambda row: point_in_ellipse(row['ra'], row['dec'], df), axis=1)]
#            for ind2, row2 in transients_inside.iterrows():
#                print ephem.date(exp['date']), str(ephem.hours(np.radians(row2['ra']))), str(ephem.degrees(np.radians(row2['dec']))), exp['expnum'],exp['band']
        matched_obs = pd.concat([matched_obs, transients_inside])
    return matched_obs


def extend(target, year='Y4', chisq_cut=10, teff_min=0.3):
    '''
    Identifies matches to the target object in the SE catalog from the given year and returns them as a dataframe.
    :param target: dataframe of observations
    :param year: DES observing year to be searched, e.g. 'Y4'
    :param chisq_cut: delta_chisq cut for matches, default=10
    :param teff_min: minimum t_eff of exposures to consider
    :return: dataframe of matches
    '''
    matched_good = pd.DataFrame()
    orbit = getOrbit(target)
    search_exps = search_exposures(target, orbit, year=year, teff_min=teff_min)
    if len(search_exps):
        print 'Number of exposures to be searched: ', len(search_exps)
        print (search_exps[['expnum', 'date','band','ccd','nite','t_eff','fwhm_asec']]).head(20)
        matched_obs = find_observations(search_exps, target, orbit, year=year, nsigma=5)
    else:
        matched_obs = []
        print 'No exposures to search'

    if len(matched_obs):
        matched_obs['ccd'] = matched_obs.apply(lambda row: get_ccd(row, search_exps), axis=1)
        matched_obs['ccd'] = matched_obs['ccd'].astype(int)
        matched_good = matched_obs[matched_obs['delta_chisq'] < chisq_cut]
        print 'Number of good matches: ', len(matched_good)
        if len(matched_good) > 0:
            matched_good['date'] = matched_good.apply(lambda row: str(ephem.date(row['date'])), axis=1)
            matched_good['ra'] = matched_good.apply(lambda row: str(ephem.hours(np.radians(row['ra']))), axis=1)
            matched_good['dec'] = matched_good.apply(lambda row: str(ephem.degrees(np.radians(row['dec']))), axis=1)
            try:
                matched_good['mag'] = matched_good.apply(lambda row: round(mag_calib(row), 2), axis=1)
            except TypeError as e:
                print e

    return matched_good


def main():
    infile = '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/DeeDee-noY4.csv'
    target = pd.read_csv(infile, comment='#')
    matches = extend(target, year='Y4', chisq_cut=5)
    if len(matches):
        matches[['date', 'ra', 'dec', 'expnum', 'exptime', 'band', 'ccd', 'mag', 'delta_chisq', 'objid']].to_csv(
            '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/tmp_good.csv', index=None)
        print matches.head(len(matches))
    else:
        print 'No good matches found.'
    print 'Done!'


if __name__=='__main__':
    main()