""" This script will do calibration testing for ETSI."""
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from config import Configuration
import os
from astropy.io import fits
import time
import numpy as np


class Calibration:

    @staticmethod
    def mk_color_flats(colors):
        """ This is the wrapper script to make each color flat, and then generate comparison files.
        :parameter colors - A list of strings which denote each color of the flat in nm.
        """

        for color in colors:
            st = time.time()  # clock started

            # get the file list
            Utils.log("Getting file list...", "info")
            if color == '500':
                blu = Preprocessing.mk_flat(Configuration.FLAT_DIRECTORY + Configuration.DATE + "_" + color, color)
            if color == '600':
                grn = Preprocessing.mk_flat(Configuration.FLAT_DIRECTORY + Configuration.DATE + "_" + color, color)
            if color == '763':
                red = Preprocessing.mk_flat(Configuration.FLAT_DIRECTORY + Configuration.DATE + "_" + color, color)
            fn = time.time()  # clock stopped
            Utils.log("Imaging cleaning complete in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

        bmr = blu - red
        fits.writeto(Configuration.FLAT_DIRECTORY + "500_m_763.fits", bmr, overwrite=True)
        bdr = blu / red
        fits.writeto(Configuration.FLAT_DIRECTORY + "500_d_763.fits", bdr, overwrite=True)
        Utils.log("The mean value of 500 minus 763 is: " + str(np.around(np.mean(bmr), decimals=4)) + "+/-" +
                  str(np.around(np.std(bmr), decimals=4)), "info")
        Utils.log("The mean value of 500 divided by 763 is: " + str(np.around(np.mean(bdr), decimals=4)) + "+/-" +
                  str(np.around(np.std(bdr), decimals=4)), "info")

        bmg = blu - grn
        fits.writeto(Configuration.FLAT_DIRECTORY + "500_m_600.fits", bmg, overwrite=True)
        bdg = blu / grn
        fits.writeto(Configuration.FLAT_DIRECTORY + "500_d_600.fits", bdg, overwrite=True)
        Utils.log("The mean value of 500 minus 600 is: " + str(np.around(np.mean(bmg), decimals=4)) + "+/-" +
                  str(np.around(np.std(bmg), decimals=4)), "info")
        Utils.log("The mean value of 500 divided by 600 is: " + str(np.around(np.mean(bdg), decimals=4)) + "+/-" +
                  str(np.around(np.std(bdg), decimals=4)), "info")

        gmr = grn - red
        fits.writeto(Configuration.FLAT_DIRECTORY + "600_m_763.fits", gmr, overwrite=True)
        gdr = grn / red
        fits.writeto(Configuration.FLAT_DIRECTORY + "600_d_763.fits", gdr, overwrite=True)
        Utils.log("The mean value of 600 minus 763 is: " + str(np.around(np.mean(gmr), decimals=4)) + "+/-" +
                  str(np.around(np.std(gmr), decimals=4)), "info")
        Utils.log("The mean value of 600 divided by 763 is: " + str(np.around(np.mean(gdr), decimals=4)) + "+/-" +
                  str(np.around(np.std(gdr), decimals=4)), "info")
        print('hold')

