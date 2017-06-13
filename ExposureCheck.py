
from contextlib import contextmanager
import sys, os

import pandas as pd
import MPCRecord
import mechanize
from bs4 import BeautifulSoup
from optparse import OptionParser


def MPChecker(ra, dec, date, box=100):

    def parseResponses(MPSoup):
        '''

        :param MPSoup: Beautiful Soup object with output of MPChecker
        :type MPSoup: BeautifulSoup
        :return: A data frame with all the matched observations
        '''
        outFrame = pd.DataFrame(columns=['tmpid', 'designation', 'mpcRA', 'mpcDec', 'v', 'offsetRA', 'offsetDec', 'motion_ra', 'motion_dec', 'orbit', 'comment'])
        responseSet = MPSoup.find_all('pre')

        if responseSet.__len__() == 0:
            return outFrame

        lines = responseSet[0].text.split('\n')  # type: List[str]

        for i in range(4, len(lines)):
            if lines[i] == "" or lines[i] == " ":
                continue

            object = pd.Series(index=['tmpid', 'designation', 'mpcRA', 'mpcDec', 'v', 'offsetRA',
                                      'offsetDec', 'motion_ra', 'motion_dec', 'orbit', 'comment'])  # type: List[str]

            object['tmpid'] = lines[i][0:9]
            object['designation'] = lines[i][9:25]
            object['mpcRA'] = lines[i][25:36]
            object['mpcDec'] = lines[i][36:47]
            object['v'] = lines[i][47:53]
            object['offsetRA'] = lines[i][53:60]
            object['offsetDec'] = lines[i][60:68]
            object['motion_ra'] = lines[i][68:72]
            object['motion_dec'] = lines[i][75:81]
            object['orbit'] = lines[i][81:87]
            object['comment'] = lines[i][87:]

            outFrame.loc[i - 4] = object

        return outFrame

    def buildmatch(matches, recid):
        matchframe = pd.DataFrame({'matchtext': matches})

        matchframe['tmpid'] = recid
        matchframe['designation'] = matchframe.matchtext.str[0:25]
        matchframe['mpcRA'] = matchframe.matchtext.str[25:36]
        matchframe['mpcDec'] = matchframe.matchtext.str[36:47]
        matchframe['v'] = matchframe.matchtext.str[47:54]
        matchframe['offsetRA'] = matchframe.matchtext.str[54:61]
        matchframe['offsetDec'] = matchframe.matchtext.str[61:69]
        matchframe['motion_ra'] = matchframe.matchtext.str[69:76]
        matchframe['motion_dec'] = matchframe.matchtext.str[76:81]
        matchframe['orbit'] = matchframe.matchtext.str[81:87]
        matchframe['comment'] = matchframe.matchtext.str[87:]

        del matchframe['matchtext']
        return matchframe

    br = mechanize.Browser()
    br.set_handle_robots(False)
    br.open("http://www.minorplanetcenter.net/cgi-bin/checkmp.cgi")
    br.form = list(br.forms())[1]

    br.form.find_control(name="which").value = ["pos"]

    br['ra'] = MPCRecord.MPC_Obs.setRA(ra)
    br['decl'] = MPCRecord.MPC_Obs.setDec(dec)
    mpc_date = MPCRecord.MPC_Obs.setDate(date)
    year = mpc_date[0:4]
    month = mpc_date[5:7]
    day = mpc_date[8:]
    day = "%.2f" % float(day)

    br['year'] = year
    br['month'] = month
    br['day'] = day

    br['oc'] = 'W84'
    br['radius'] = str(min(box, 300))
    br['limit'] = '24'
    response = br.submit()
    resp = response.read()
    soup = BeautifulSoup(resp)

    responses = parseResponses(soup)
    return responses
