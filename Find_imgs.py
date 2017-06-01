from __future__ import division

import easyaccess as ea
import numpy as np
import pandas as pd
import urllib

import os
import requests
#from KBO import compute_chip
#from extendOrbit import getOrbit
#from linkmap import build_kdtree
from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord


ccdWidth = .298   #width of a CCD in degrees
ccdHeight = .149349   #height of a CCD in degrees


def findImgs(ra_in, dec_in):

    conn_desoper = ea.connect(section='desoper')

    right_ra = ra_in + ccdWidth
    left_ra = ra_in - ccdWidth

    upper_dec = dec_in + ccdHeight
    lower_dec = dec_in - ccdHeight

    query = ("SELECT FILENAME, NITE FROM PROD.IMAGE WHERE RA_CENT BETWEEN " + str(left_ra) + " AND " + str(right_ra) +
             " AND DEC_CENT BETWEEN "+ str(lower_dec) + " AND " + str(upper_dec) +
             " AND EXPNUM > 552000 AND FILENAME LIKE '%immasked.fits'")  # Y4 -> expnum > 552000

    imglist = conn_desoper.query_to_pandas(query)

    pathlist = pd.DataFrame()

    for index, row in imglist.iterrows():
        #print row
        query = "SELECT PATH FROM PROD.FILE_ARCHIVE_INFO WHERE FILENAME = '"+ str(row['FILENAME']) + "'"
        path = conn_desoper.query_to_pandas(query)
        path['PATH'] = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + path['PATH'] + "/" + str(row['FILENAME']) #+ '.fz'
        pathlist = pathlist.append(path)

    #https://desar2.cosmology.illinois.edu/DESFiles/desarchive/
    dir = '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/FitsFiles'

    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)

    for index, row in pathlist.iterrows():
        #url = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + str(row['PATH'])
        fits_filename = download_file(row['PATH'])
        ds9cut(fits_filename, ra_in, dec_in)

    return pathlist


# TRUE https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/firstcut/Y4N/20170216-r2876/D00621526/p01/red/immask/D00621526_i_c01_r2876p01_immasked.fits.fz
# MINE https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/finalcut/Y2A1/Y3-2379/20160109/D00509722/p01/red/immaskD00509722_i_c36_r2379p01_immasked.fits.fz

def download_file(url):

    local_filename = url.split('/')[-5] + ".fits.fz"
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
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
    return local_filename



def ds9cut(fits_filename, ra, dec, side=5):

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    print c.ra.hms
    print c.dec.dms

    ra_ang = Angle(str(ra) + 'd')
    dec_ang = Angle(str(dec) + 'd')

    ra_hms = ra_ang.hms
    dec_dms = dec_ang.dms

    ra_str = str(ra_hms[0]) + ':' + str(ra_hms[1]) + ':' + str(ra_hms[2])
    dec_str = str(dec_dms[0]) + ':' + str(dec_dms[1]) + ':' + str(dec_dms[2])
    png_name = fits_filename.split('.')[0]
    cmdstr = "ds9x " + str(fits_filename) + ' -scale zscale -scale squared -crop ' + ra_str + ' ' + dec_str + ' ' + str(side) + ' ' + str(side) + 'wcs fk5 -colorbar no -saveimage ' + str(png_name) + '.png -exit'
    #cmdstr = "ds9x " + str(fits_filename) + ' -scale zscale -scale squared -colorbar no -saveimage ' + str(png_name) + '.png -exit'
    print cmdstr
    os.system(cmdstr)

df = findImgs(58.54, -27.6)
print "done"

#dir = '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/FitsFiles'
#os.chdir(dir)
#filename = download_file('https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/firstcut/Y4N/20170129-r2843/D00614390/p01/red/immask/D00614390_z_c01_r2843p01_immasked.fits.fz')
#ds9cut(filename)




# def cut_fits(fits_filename, ra, dec, box_size):
#     """
#
#     :param fits_filename:
#     :param ra:
#     :param dec:
#     :param box: length in arcmin of side of box
#     :return:
#     """
#     hdu_list = fits.open(fits_filename)
#
#     header = hdu_list[0].header
#
#     ra_img_min = header['RACMIN']
#     ra_img_max = header['RACMAX']
#     dec_img_min = header['DECCMIN']
#     dec_img_max = header['DECCMAX']
#
#     # in arcsec/pixel
#     ra_scale = header['PIXSCAL2']
#     dec_scale = header['PIXSCAL1']
#
#     # in deg/pixel
#     ra_scale /= 3600.
#     dec_scale /= 3600.
#
#     box_deg = box_size / 60.
#
#     ra_min = ra - box_deg / 2.
#     ra_max = ra + box_deg / 2.
#
#     dec_min = dec - box_deg / 2.
#     dec_max = dec + box_deg / 2.
#
#     x_min = max(int((ra_min - ra_img_min) / ra_scale), 0)
#     x_max = min(int((ra_max - ra_img_min) / ra_scale), header['NAXIS2'] - 1)
#
#     y_min = max(int((dec_min - dec_img_min) / dec_scale), 0)
#     y_max = min(int((dec_max - dec_img_min) / dec_scale), header['NAXIS2'] - 1)
#
#     hdu_list[0].data = hdu_list[0].data[x_min:(x_max + 1)][y_min:(y_max + 1)]
#     if x_min:
#         header['RACMIN'] = ra_min
#         header['RAC1'] = ra_min
#
#
#     if y_min:
#         header['DECCMIN'] = dec_min
