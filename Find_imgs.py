from __future__ import division

import easyaccess as ea
import numpy as np
import pandas as pd
import urllib
import time
import os
import subprocess
import requests
# from KBO import compute_chip
# from extendOrbit import getOrbit
# from linkmap import build_kdtree
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from glob import glob
from astropy.wcs import WCS

ccdWidth = .298  # width of a CCD in degrees
ccdHeight = .149349  # height of a CCD in degrees


def findImgs(expnumCCD_list, cat_list, diff_img_list, both_list, ra, dec, side=5,
             coadd_list=None, asteroid_list=None, keep_fits=False, only_asteroids=False):
    """
    
    :param expnumCCDs: List of (expnum, ccd) pairs
    :return: 
    """

    if not expnumCCD_list:
        return pd.DataFrame()

    expnumCCDSet = set(expnumCCD_list)
    desoper = ea.connect(section='desoper')

    query = "SELECT FILENAME, NITE, EXPNUM, CCDNUM FROM PROD.IMAGE WHERE FILENAME LIKE '%immasked.fits' AND ("
    for expnum, ccd in expnumCCDSet:
        query += "(EXPNUM = " + str(expnum) + " AND CCDNUM = " + str(ccd) + ") OR "

    query = query[0:-4]
    query += ')'
    print query
    imglist = desoper.query_to_pandas(query)

    pathlist = []

    for index, row in imglist.iterrows():
        # print row
        query = "SELECT PATH FROM PROD.FILE_ARCHIVE_INFO WHERE FILENAME = '" + str(row['FILENAME']) + "'"
        path = desoper.query_to_pandas(query)
        new_path = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + path['PATH'][0] + "/" + str(
            row['FILENAME'])  # + '.fz'
        pathlist.append((new_path, row['EXPNUM'], row['CCDNUM']))

    # https://desar2.cosmology.illinois.edu/DESFiles/desarchive/
    dir = '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/Exposures'

    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)

    for elt in pathlist:
        # url = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + str(row['PATH'])
        fits_filename = download_file(elt[0], elt[1], elt[2])
        # cut_fits(fits_filename)

        ds9cut(fits_filename, elt[1], elt[2], ra, dec, cat_list,
               diff_img_list, both_list, coadd_list=coadd_list, asteroid_list=asteroid_list, side=side,
               only_asteroids=only_asteroids)
    if not keep_fits:
        cleanDir()

    return pathlist


# def findImgs(ra_in, dec_in):
#
#     conn_desoper = ea.connect(section='desoper')
#
#     right_ra = ra_in + ccdWidth
#     left_ra = ra_in - ccdWidth
#
#     upper_dec = dec_in + ccdHeight
#     lower_dec = dec_in - ccdHeight
#
#     query = ("SELECT FILENAME, NITE FROM PROD.IMAGE WHERE RA_CENT BETWEEN " + str(left_ra) + " AND " + str(right_ra) +
#              " AND DEC_CENT BETWEEN "+ str(lower_dec) + " AND " + str(upper_dec) +
#              " AND EXPNUM > 552000 AND FILENAME LIKE '%immasked.fits'")  # Y4 -> expnum > 552000
#
#     imglist = conn_desoper.query_to_pandas(query)
#
#     pathlist = pd.DataFrame()
#
#     for index, row in imglist.iterrows():
#         #print row
#         query = "SELECT PATH FROM PROD.FILE_ARCHIVE_INFO WHERE FILENAME = '"+ str(row['FILENAME']) + "'"
#         path = conn_desoper.query_to_pandas(query)
#         path['PATH'] = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + path['PATH'] + "/" + str(row['FILENAME']) #+ '.fz'
#         pathlist = pathlist.append(path)
#
#     #https://desar2.cosmology.illinois.edu/DESFiles/desarchive/
#     dir = 'FitsFiles'
#
#     if not os.path.exists(dir):
#         os.mkdir(dir)
#     os.chdir(dir)
#
#     for index, row in pathlist.iterrows():
#         #url = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + str(row['PATH'])
#         fits_filename = download_file(row['PATH'])
#         ds9cut(fits_filename, ra_in, dec_in)
#
#     cleanDir()
#
#     return pathlist


# TRUE https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/firstcut/Y4N/20170216-r2876/D00621526/p01/red/immask/D00621526_i_c01_r2876p01_immasked.fits.fz
# MINE https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/finalcut/Y2A1/Y3-2379/20160109/D00509722/p01/red/immaskD00509722_i_c36_r2379p01_immasked.fits.fz

def download_file(url, expnum, ccd):
    local_filename = "Exp_" + str(expnum) + "_ccd_" + str(ccd) + ".fits.fz"
    # NOTE the stream=True parameter
    try:
        r = requests.get(url, stream=True, auth=('lzullo', 'lzu70chips'))
        if r.status_code == 404:
            url = url + '.fz'
            r = requests.get(url, stream=True, auth=('lzullo', 'lzu70chips'))
    except:
        print "THIS IS AN ERROR YO " + str(url)

    print url
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
    return local_filename


def descut_convert(ra, dec):
    fits_filenames = glob('*.fits')
    for file in fits_filenames:
        cmdstr = ""
        #cmdstr += 'rm -f /tmp/.X1-lock; export DISPLAY=:1; Xvfb :1 -screen 0 1024x768x16 & '
        cmdstr += 'ds9x ' + str(file) + ' -scale zscale -scale squared'
        cmdstr += (' -colorbar no')

        cmdstr += (' -regions command "ICRS;circle(' + ra.to_string(unit=u.hour, sep=':', alwayssign=True)
                   + ',' + dec.to_string(unit=u.degree, sep=':', alwayssign=True) + ',10i)#color=yellow"')
        cmdstr += ' -zoom to fit -saveimage ' + str(file) + '.png'
        cmdstr += ' -exit'
        print cmdstr
        subprocess.check_call(cmdstr, shell=True)



def fetchFiles(jobid, token):

    request = 'http://descut.cosmology.illinois.edu/api/jobs/?token=' + token + '&jobid=' + jobid
    req = requests.get(request)

    # Checks to see if request is successful
    while req.status_code != 200 or req.json()['job_status'] != "SUCCESS":
        time.sleep(0.1)
        req = requests.get(request)

    text_file = open('file_list_' + jobid + '.txt', "w")

    json = req.json()

    for link in req.json()['links']:
        text_file.write(link + '\n')

    text_file.close()
    os.system('wget -i ' + 'file_list_' + jobid + '.txt')


def grab_descut_imgs(ra, dec):
    req = requests.post('http://descut.cosmology.illinois.edu/api/token/',
                        data={'username': 'kfranson', 'password': 'kfr70chips'})
    print(req)
    print(req.text)

    ra_rad, dec_rad = ra.degree, dec.degree
    token = str(req.json()['token'])
    # create body of request
    body = {
        'token': token,  # required
        'ra': ra_rad,  # required
        'dec': dec_rad,  # required
        'band': 'g,r,i,z',
        'job_type': 'single',
        'list_only': 'true'
    }
    req = requests.post('http://descut.cosmology.illinois.edu/api/jobs/', data=body)
    print(req)
    print(req.text)
    fetchFiles(req.json()['job'], token)
    descut_convert(ra, dec)


def ds9cut(fits_filename, expnum, ccd, ra, dec, cat_list, diff_img_list, both_list,
           asteroid_list=None, coadd_list=None, side=5, only_asteroids=False, comparison_imgs=True):
    # (ra, dec, expnum)
    angles = []
    if not only_asteroids:
        angles.append(SkyCoord(ra=Angle(str(ra) + 'd'), dec=Angle(str(dec) + 'd')))
    else:
        side = 60
        hdu_list = fits.open(fits_filename)
        header = hdu_list['SCI'].header
        for ra_ast, dec_ast in asteroid_list:
            racmin = float(header['RACMIN'])
            racmax = float(header['RACMAX'])
            deccmin = float(header['DECCMIN'])
            deccmax = float(header['DECCMAX'])

            if (float(header['RACMIN']) < ra_ast < float(header['RACMAX']) and
                    float(header['DECCMIN']) < dec_ast < float(header['DECCMAX'])):
                angles.append(SkyCoord(ra=Angle(str(ra_ast) + 'd'), dec=Angle(str(dec_ast) + 'd')))

    if not angles:
        print "No asteroid in ccd " + str(ccd)
        return

    for ang in angles:
        ra_ang = ang.ra
        dec_ang = ang.dec
        print ra_ang.to_string(unit=u.hour, sep=':', alwayssign=True)
        print dec_ang.to_string(unit=u.degree, sep=':', alwayssign=True)
        png_name = fits_filename.split('.')[0]

        cmdstr = ""
        #cmdstr += 'rm -f /tmp/.X1-lock; export DISPLAY=:1; Xvfb :1 -screen 0 1024x768x16 & '
        cmdstr += 'ds9x ' + str(fits_filename) + ' -scale zscale -scale squared'
        cmdstr += (' -crop ' + ra_ang.to_string(unit=u.hour, sep=':', alwayssign=True) + ' ' +
                   dec_ang.to_string(unit=u.degree, sep=':', alwayssign=True) + ' ' + str(side) + ' ' + str(side) +
                   ' wcs icrs arcsec')
        cmdstr += (' -colorbar no -grid yes -grid type publication -grid system wcs -grid axes type interior' +
                   ' -grid axes style 1 -grid format1 d.2 -grid format2 d.2 -grid numerics yes -grid numerics' +
                   ' vertical no')

        for elt in cat_list:
            if elt[2] == expnum and elt[3] == ccd:
                ra_elt = Angle(str(elt[0]) + 'd')
                dec_elt = Angle(str(elt[1]) + 'd')

                cmdstr += (' -regions command "ICRS;circle(' + ra_elt.to_string(unit=u.hour, sep=':', alwayssign=True)
                           + ',' + dec_elt.to_string(unit=u.degree, sep=':', alwayssign=True) + ',10i)#color=green"')
        for elt in diff_img_list:
            if elt[2] == expnum and elt[3] == ccd:
                ra_elt = Angle(str(elt[0]) + 'd')
                dec_elt = Angle(str(elt[1]) + 'd')

                cmdstr += (' -regions command "ICRS;circle(' + ra_elt.to_string(unit=u.hour, sep=':', alwayssign=True)
                           + ',' + dec_elt.to_string(unit=u.degree, sep=':', alwayssign=True) + ',10i)#color=red"')

        for elt in both_list:
            if elt[2] == expnum and elt[3] == ccd:
                ra_elt = Angle(str(elt[0]) + 'd')
                dec_elt = Angle(str(elt[1]) + 'd')

                cmdstr += (' -regions command "ICRS;circle(' + ra_elt.to_string(unit=u.hour, sep=':', alwayssign=True)
                           + ',' + dec_elt.to_string(unit=u.degree, sep=':', alwayssign=True) + ',10i)#color=blue"')

        if coadd_list:
            for elt in coadd_list:
                ra_elt = Angle(str(elt[0]) + 'd')
                dec_elt = Angle(str(elt[1]) + 'd')

                cmdstr += (' -regions command "ICRS;circle(' + ra_elt.to_string(unit=u.hour, sep=':', alwayssign=True)
                           + ',' + dec_elt.to_string(unit=u.degree, sep=':', alwayssign=True) + ',10i)#color=pink"')

        if asteroid_list:
            for elt in asteroid_list:
                ra_elt = Angle(str(elt[0]) + 'd')
                dec_elt = Angle(str(elt[1]) + 'd')

                cmdstr += (' -regions command "ICRS;circle(' + ra_elt.to_string(unit=u.hour, sep=':', alwayssign=True)
                           + ',' + dec_elt.to_string(unit=u.degree, sep=':', alwayssign=True) + ',10i)#color=yellow"')

        cmdstr += ' -zoom to fit -saveimage ' + str(png_name) + '.png'
        cmdstr += ' -exit'
        print cmdstr
        subprocess.check_call(cmdstr, shell=True)
        if comparison_imgs:
            folder_name = 'ra_' + str(ra_ang.degree) + '_dec_' + str(dec_ang.degree)
            if not os.path.exists(folder_name):
                os.mkdir(folder_name)
            os.system('mv ' + png_name + '.png ' + folder_name + '/' + png_name + '.png')
            os.chdir(folder_name)
            grab_descut_imgs(ra_ang, dec_ang)
            os.chdir('..')
        # os.system("rm -f /tmp/.X1-lock")


def cleanDir():
    filelist = glob("*.fz")

    for file in filelist:
        os.system('rm -r ' + file)

    print "Fits files deleted"


if __name__ == '__main__':
    df = findImgs(75.7, -24.)
    print "done"

# dir = '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/FitsFiles'
# os.chdir(dir)
# filename = download_file('https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/firstcut/Y4N/20170129-r2843/D00614390/p01/red/immask/D00614390_z_c01_r2843p01_immasked.fits.fz')
# ds9cut(filename)


'''

def cut_fits(fits_filename):#, ra_min, ra_max, dec_min, dec_max):
    """

    :param fits_filename:
    :param ra:
    :param dec:
    :param box: length in arcmin of side of box
    :return:
    """
    hdu_list = fits.open(fits_filename)

    header = hdu_list[0].header
    wmap = WCS(hdu_list[0].header)

    ra_img_min = header['RACMIN']
    ra_img_max = header['RACMAX']
    dec_img_min = header['DECCMIN']
    dec_img_max = header['DECCMAX']

    # in arcsec/pixel
    ra_scale = header['PIXSCAL2']
    dec_scale = header['PIXSCAL1']

    # in deg/pixel
    ra_scale /= 3600.
    dec_scale /= 3600.

    box_deg = box_size / 60.

    ra_min = ra - box_deg / 2.
    ra_max = ra + box_deg / 2.

    dec_min = dec - box_deg / 2.
    dec_max = dec + box_deg / 2.

    x_min = max(int((ra_min - ra_img_min) / ra_scale), 0)
    x_max = min(int((ra_max - ra_img_min) / ra_scale), header['NAXIS2'] - 1)

    y_min = max(int((dec_min - dec_img_min) / dec_scale), 0)
    y_max = min(int((dec_max - dec_img_min) / dec_scale), header['NAXIS2'] - 1)

    hdu_list[0].data = hdu_list[0].data[x_min:(x_max + 1)][y_min:(y_max + 1)]
    if x_min:
        header['RACMIN'] = ra_min
        header['RAC1'] = ra_min

    if y_min:
        header['DECCMIN'] = dec_min
'''
