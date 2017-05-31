from __future__ import division

import easyaccess as ea
import numpy as np
import pandas as pd
import urllib

import wget
import os
import requests


ccdWidth = .298   #width of a CCD in degrees
ccdHeight = .149349   #height of a CCD in degrees

def findImgs(ra_in, dec_in):

    conn_desoper = ea.connect(section='desoper')

    right_ra = ra_in + ccdWidth
    left_ra = ra_in - ccdWidth

    upper_dec = dec_in + ccdHeight
    lower_dec = dec_in - ccdHeight

    query = "SELECT FILENAME, NITE FROM PROD.IMAGE WHERE RA_CENT BETWEEN " + str(left_ra) + " AND " + str(right_ra) + " AND DEC_CENT BETWEEN"+ str(lower_dec) + " AND " + str(upper_dec) + " AND FILENAME LIKE '%immasked.fits'"

    imglist = conn_desoper.query_to_pandas(query)

    pathlist = pd.DataFrame()

    for index, row in imglist.iterrows():
        #print row
        query = "SELECT PATH FROM PROD.FILE_ARCHIVE_INFO WHERE FILENAME = '"+ str(row['FILENAME']) + "'"
        path = conn_desoper.query_to_pandas(query)
        path['PATH'] = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + path['PATH'] + "/" + str(row['FILENAME']) + '.fz'
        pathlist = pathlist.append(path)

    #https://desar2.cosmology.illinois.edu/DESFiles/desarchive/
    dir = '/Users/lynuszullo/pyOrbfit/Y4_Transient_Search/FitsFiles'

    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)

    for index, row in pathlist.iterrows():
        #url = 'https://desar2.cosmology.illinois.edu/DESFiles/desarchive/' + str(row['PATH'])
        download_file(row['PATH'])

    return pathlist


# TRUE https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/firstcut/Y4N/20170216-r2876/D00621526/p01/red/immask/D00621526_i_c01_r2876p01_immasked.fits.fz
# MINE https://desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/finalcut/Y2A1/Y3-2379/20160109/D00509722/p01/red/immaskD00509722_i_c36_r2379p01_immasked.fits.fz

def download_file(url):
    print url
    local_filename = url.split('/')[-4]
    # NOTE the stream=True parameter
    r = requests.get(url, stream=True, auth=('lzullo', 'lzu70chips'))
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
                #f.flush() commented by recommendation from J.F.Sebastian
    return local_filename

df = findImgs(58.54, -27.6)
print "done"