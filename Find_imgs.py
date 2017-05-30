from __future__ import division

import easyaccess as ea
import ephem
import numpy as np
import pandas as pd
from KBO import compute_chip
from extendOrbit import getOrbit
from linkmap import build_kdtree


ccdWidth = .3/2   #width of a CCD in degrees
ccdHeight = .3   #height of a CCD in degrees

def findImgs(ra_in, dec_in):

    conn_desoper = ea.connect(section='desoper')

    right_ra = ra_in + ccdWidth
    left_ra = ra_in - ccdWidth

    upper_dec = dec_in + ccdHeight
    lower_dec = dec_in - ccdHeight

    query = "SELECT FILENAME FROM PROD.IMAGE WHERE RA_CENT <= " + str(right_ra) + " AND RA_CENT >= "+ str(left_ra) + " AND DEC_CENT >="+ str(lower_dec)+ " AND DEC_CENT <=" + str(upper_dec)

    imglist = conn_desoper.query_to_pandas(query)

    return imglist
