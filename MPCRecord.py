from __future__ import print_function
import ephem
import pandas as pd
import numpy as np
import sys
import traceback


# designed to convert an object from MPC format to our format
class MPCToTNO:
    def __init__(self, observations):
        """

        Args:
            observations (str): Text of observations in MPC format
        """
        obs_list = observations.split('\n')
        observations = []
        for obs in obs_list:
            if obs:
                observations.append(self.mpc_to_dict(obs))

        self.observations = pd.DataFrame(observations)
        self.designation = self.observations['designation'][0]

    @classmethod
    def mpc_to_dict(cls, obs):
        """

        Args:
            obs (str):

        Returns:
            dict
        """

        prov_des = obs[5:12]
        date_mpc = obs[15:32]
        ra_mpc = obs[32:44]
        dec_mpc = obs[44:56]
        mag_mpc = obs[65:70]
        band_mpc = obs[70]
        obs_code = obs[77:80]
        observation = dict()
        observation['designation'] = cls.__unpack_designation(prov_des)
        observation['date'] = cls.date_convert(date_mpc)
        observation['ra'] = cls.ra_convert(ra_mpc)
        observation['dec'] = cls.dec_convert(dec_mpc)
        observation['mag'] = mag_mpc.strip()
        observation['band'] = band_mpc
        observation['obs_code'] = obs_code
        return observation

    @staticmethod
    def __unpack_designation(designation):
        packed_yr = designation[0]

        # 'I' -> 18
        century = ord(packed_yr) - ord('A') + 10
        year = str(century) + designation[1:3]
        month = designation[3] + designation[6]

        cycle = designation[4:6]

        letters = dict()

        if ord('0') <= ord(cycle[0]) <= ord('9'):
            pass
        else:
            # gets alphabet
            for i in range(65, 123):
                if i not in range(91, 97):
                    letters[chr(i)] = i - 55
            cycle = str(letters[cycle[0]]) + cycle[1]

        return year + " " + month + cycle

    @staticmethod
    def date_convert(mpc_date):
        """

        Args:
            mpc_date (str):
        Returns:
            str
        """
        mpc_date = mpc_date.strip()
        year = mpc_date[0:4]
        month = mpc_date[5:7]
        day = mpc_date[8:]

        return str(ephem.date((int(year), int(month), float(day))))

    @staticmethod
    def ra_convert(mpc_ra):
        mpc_ra = mpc_ra.strip()
        hours = mpc_ra[0:2]
        minutes = mpc_ra[3:5]
        seconds = mpc_ra[6:]
        return hours + ":" + minutes + ":" + seconds

    @staticmethod
    def dec_convert(mpc_dec):
        """

        Args:
            mpc_dec (str):

        Returns:
            str
        """
        mpc_dec = mpc_dec.strip()
        sign = mpc_dec[0]
        deg = mpc_dec[1:3]
        minutes = mpc_dec[4:6]
        seconds = mpc_dec[7:]

        if sign == '-':
            return sign + deg + ":" + minutes + ":" + seconds
        else:
            return deg + ":" + minutes + ":" + seconds


class MPCObj(object):
    # REQUIRES: observations is pandas dataframe with each row a single observation
    #           Must have columns DATE_OBS, RA, DEC, Mag, Band
    def __init__(self, observations, newObj=False, unpackedDesignation='       ',
                 observers='', measurers='D. Gerdes', tel_details='4.0-m f/2.87 reflector + CCD',
                 ob_code='W84', comments='',
                 contact='D. Gerdes, Dept. of Physics, University of Michigan, Ann Arbor, MI 48109',
                 con_email='gerdes@umich.edu', ack_message='', ack_email='kfranson@umich.edu'):
        """
        :type observations: pandas.DataFrame
        :param observations: DataFrame with columns DATE_OBS, RA, DEC, Mag, Band representing all observations
        :return:
        """

        self.obsList = np.array([], dtype=MPC_Obs)

        observations.sort(columns='DATE_OBS')

        for obs in observations.iterrows():
            discovery = ' '
            if obs[0] == 0:
                discovery = '*'

            # obs[1] is data series, obs[0] is index
            self.obsList = np.append(self.obsList,
                                     [MPC_Obs(obs[1], provNum=unpackedDesignation,
                                              discovery=discovery, observCode=ob_code, newObj=newObj)])

        self.observers = observers
        self.measurers = measurers
        self.tel_details = tel_details
        self.ob_code = ob_code
        self.comments = comments
        self.contact = contact
        self.con_email = con_email
        self.ack_message = ack_message
        self.ack_email = ack_email

    def get_submission(self):
        """

        :return: Submission string for MPC
        """

        cod = "COD " + self.ob_code
        assert (self.valid_line(cod))
        output = cod + '\n'

        con1 = "CON " + self.contact
        assert (self.valid_line(con1))
        output += con1 + '\n'

        con2 = "CON [" + self.con_email + "]"
        assert (self.valid_line(con2))
        output += con2 + "\n"

        output += self.__pack_observers()
        output += self.__pack_measurers()

        tel = "TEL " + self.tel_details
        assert (self.valid_line(tel))
        output += tel + '\n'

        ack = "ACK " + self.ack_message
        assert (self.valid_line(ack))
        output += ack + '\n'

        ac2 = "AC2 " + self.ack_email
        assert (self.valid_line(ac2))
        output += ac2 + '\n'

        for obs in self.obsList:
            output += obs.getRecord() + '\n'

        return output

    def __pack_observers(self):
        return self.__pack_flag(self.observers, "OBS")

    def __pack_measurers(self):
        return self.__pack_flag(self.measurers, "MEA")

    @classmethod
    def __pack_flag(cls, names, flag):

        obs_lines = [flag + " " + names]

        while not cls.valid_line(obs_lines[-1]):
            i = obs_lines[-1].rindex(',', 0, 79)

            new_line = obs_lines[-1][i:]
            obs_lines[-1] = obs_lines[-1][0:i]
            obs_lines.append(flag + " " + new_line[3:])

        output_lines = ""
        for line in obs_lines:
            assert (cls.valid_line(line))
            output_lines += line + "\n"

        return output_lines

    @staticmethod
    def valid_line(line):
        """

        :rtype: bool
        """
        return len(line) <= 80

    def get_observation(self, n):
        """
        :param n: Index of observation
        :return: 80 character observation string in MPC format
        """
        return self.obsList[n].getRecord()


# defines a minor planet optical observation in MPC format
# See http://www.minorplanetcenter.net/iau/info/OpticalObs.html for more info
class MPC_Obs(object):
    # REQUIRES: obs is pandas series that represents single observation
    #           & has columns DATE_OBS, RA, DEC, Mag, Band

    def __init__(self, obs, obsNum='     ', provNum='       ', discovery=' ', note1=' ', note2='C', observCode='W84',
                 newObj=True):
        """
        :type newObj: bool
        :type provNum: str
        :type obs: pandas.Series
        """

        self.__dFrame_check(obs)
        self.newObj = newObj
        self.obsNum = self.__setObsNum(obsNum)
        self.provNum = self.__setProvNum(provNum)
        self.discovery = self.__setDiscovery(discovery)

        self.__parse(obs)

        if len(note1) != 1:
            self.__error("Invalid note1")
        else:
            self.note1 = note1

        if len(note2) != 1:
            self.__error("Invalid note2")
        else:
            self.note2 = note2

        self.date = self.setDate(self.date)
        self.ra = self.setRA(self.ra)
        self.dec = self.setDec(self.dec)
        self.magBand = self.__setMagAndBand(self.mag, self.band)
        self.observCode = observCode

    def getRecord(self):
        """
        Prints 80 character string that represents single observation in MPC format
        :rtype: str
        """
        record = self.obsNum + self.provNum + self.discovery + self.note1 + self.note2 + \
                 self.date + self.ra + self.dec + '         ' + self.magBand + '      ' + self.observCode
        assert len(record) == 80

        return record

    def __dFrame_check(self, obs):
        """
        Checks dataFrame for correct formatting
        :type obs: pandas.Series
        """

        if ('DATE_OBS' not in obs.index or 'RA' not in obs.index or 'DEC' not in obs.index or
                    'MAG' not in obs.index or 'BAND' not in obs.index):
            indices = ""
            if 'DATE_OBS' not in obs.index:
                indices += 'DATE_OBS'
            if 'RA' not in obs.index:
                if indices != "":
                    indices += ', RA'
                else:
                    indices += 'RA'
            if 'DEC' not in obs.index:
                if indices != "":
                    indices += ', DEC'
                else:
                    indices += 'DEC'

            if 'MAG' not in obs.index:
                if indices != "":
                    indices += ', MAG'
                else:
                    indices += 'MAG'

            if 'BAND' not in obs.index:
                if indices != "":
                    indices += ', BAND'
                else:
                    indices += 'BAND'

            self.__error("obs failed to have the following indices: " + indices)

    # MODIFIES: num
    # RETURNS: modified version of num
    def __setObsNum(self, num):
        """

        :rtype: str
        """

        if self.newObj:
            num = '     '
        else:
            num.zfill(5)

            if len(num) > 5:
                num = self.__pack(num)

        assert len(num) == 5
        return num

    def __setProvNum(self, num):
        if self.newObj:
            if len(num) > 7 or (len(num) == 7 and num[0] != ' '):
                self.__error("Temporary Designation must be less than 7 characters in length")

            while len(num) < 7:
                num = ' ' + num

        else:
            if ' ' in num:
                # assume num is in unpacked prov designation (i.e. 2008 AA360)
                num = self.__pack(num)

            elif len(num) != 7:
                self.__error("Unrecognized format for provisional designation")

        assert len(num) == 7
        return num

    @classmethod
    def setDate(cls, date_in):
        """

        :param date_in: string that represents date of observation as 'YYYY/MM/DD HH:MM:SS'
        :type date_in: date
        :return: string that represents date in MPC format ('YYYY MM DD.dddddd')
        :rtype: str
        """
        try:
            date = ephem.Date(date_in)
            (year, month, day) = date.triple()
            day = "%.6f" % day

            day_str = str(day)
            month_str = str(month)
            year_str = str(year)

            if len(day_str) == 8:
                day_str = "0" + day_str
            if len(month_str) == 1:
                month_str = "0" + month_str

            if len(day_str) != 9 or len(month_str) != 2 or len(year_str) != 4:
                raise ValueError

        except ValueError:
            cls.__error("Date of observation must be in format YYYY/MM/DD HH:MM:SS")
        else:
            return year_str + ' ' + month_str + ' ' + day_str

    @classmethod
    def setRA(cls, ra_in):
        """

        :param ra_in: string that represents RA as hour:minute:second.second
        :type ra_in: str
        :return: string that represents RA in MPC format ('HH MM SS.dd')
        :rtype: str
        """
        try:
            if type(ra_in) == ephem.Angle:
                hours, minutes, seconds = str(ra_in).split(":")
            else:
                hours, minutes, seconds = ra_in.split(":")

            seconds = "%.2f" % float(seconds)
            if len(hours) == 1:
                hours = "0" + hours

            if len(minutes) == 1:
                minutes = "0" + minutes

            if len(seconds) == 4:
                seconds = "0" + seconds
            seconds += ' '
            out = hours + ' ' + minutes + ' ' + seconds

            if len(out) != 12:
                raise ValueError

        except ValueError:
            cls.__error("RA of observation must be in format hour:minute:second.second")

        else:
            return out

    @classmethod
    def setDec(cls, dec_in):
        """
        :param dec_in: string that represents Dec as hour:minute:second.second or -hour:minute:second.second
        :type dec_in: str
        :return: string that represents RA in MPC format ('HH MM SS.dd')
        :rtype: str
        """
        try:
            if type(dec_in) == ephem.Angle:
                hours, minutes, seconds = str(dec_in).split(":")
            else:
                hours, minutes, seconds = dec_in.split(":")

            seconds = "%.1f" % float(seconds)
            if hours[0] == "-":
                sign = "-"
                hours = hours[1:]
            elif hours[0] == "+":
                sign = "+"
                hours = hours[1:]
            else:
                sign = "+"

            if len(hours) == 1:
                hours = "0" + hours

            if len(minutes) == 1:
                minutes = "0" + minutes

            if len(seconds) == 3:
                seconds = "0" + seconds

            # as we're only going to a precision of one decimal place
            seconds += ' '

            out = sign + hours + ' ' + minutes + ' ' + seconds

            if len(out) != 12:
                raise ValueError

        except ValueError:
            cls.__error("Dec of observation must be in format hour:minute:second.second")

        else:
            return out

    def __setMagAndBand(self, mag_in, band_in):
        """

        :param mag_in: magnitude to set
        :param band_in: band of magnitude measurement
        :return: string that represents magnitude and band in MPC format
        """

        # round magnitude to one decimal place
        if mag_in < 10 or mag_in > 30:
            mag = '     '
            # talk to Gerdes about this
            assert 1 == 0
        else:
            mag = "%4.1f" % mag_in + ' '

        if len(mag) != 5:
            self.__error("Invalid magnitude")

        if len(band_in) != 1 or band_in not in ['g', 'r', 'i', 'z', 'Y']:
            self.__error("Invalid band")

        return mag + band_in

    def __pack(self, num):
        """
        :param num: unpacked provisional designation (i.e. "2008 AA360")
        :type num: str
        :return: packed provisional designation
        :rtype: str
        """
        out = ""

        # split up year and other info
        year, code = num.split(' ')

        hundredYr = year[:2]

        if hundredYr == "18":
            out = "I"
        elif hundredYr == "19":
            out = "J"
        elif hundredYr == "20":
            out = "K"
        else:
            self.__error("Invalid century for prov designation")

        if len(code) < 2:
            self.__error("Invalid half-month letters for prov designation")

        # adds decade and year digits of year to designation
        restOfYr = year[2:]
        out += restOfYr

        # adds half-month letter 1
        out += code[0]

        # adds half-month repeat info
        if len(code) == 2:
            out += "00"
        elif len(code) == 3:
            out += "0" + code[2]
        elif len(code) == 4:
            out += code[2:]
        elif len(code) == 5:
            out += self.__condenseNum(code[2:4])
        else:
            self.__error("Invalid half-month letters for prov designation")

        # adds half-month letter2
        out += code[5]

        return out

    def __condenseNum(self, numStr):
        """
        Takes two digit number between 10 and 61 inclusive and represents it as a letter

        :param numStr: String that represents 2 digit number <= 61
        :return: length 1 string that represents letter
        :rtype: str
        """
        assert len(numStr) == 2 and numStr[0] != '0'
        if int(numStr) not in range(10, 62):
            self.__error("Invalid two digit number to string conversion")

        letters = ""

        # gets alphabet
        for i in range(65, 123):
            if i not in range(91, 97):
                letters += chr(i)

        numStr = letters[int(numStr) - 10]
        return numStr

    def __setDiscovery(self, disc):
        if disc == ' ' or disc == '*':
            return disc
        else:
            self.__error("Discovery value not valid")

    @staticmethod
    def __error(mssg):
        """
        Prints error message to stderr and exits program
        :param mssg: Error message that you want printed. Follows 'ERROR: ' in message
        :type mssg: str

        """

        print("ERROR: " + mssg, file=sys.stderr)
        traceback.print_stack()
        exit(1)

    # gets data from obsSeries (data series for single observation)
    #
    def __parse(self, obsSeries):
        """

        :param obsSeries: data series for single observation
        :type obsSeries: pandas.Series
        """

        self.date = obsSeries["DATE_OBS"]
        self.ra = obsSeries["RA"]
        self.dec = obsSeries["DEC"]
        self.mag = obsSeries["MAG"]
        self.band = obsSeries["BAND"]
        return
