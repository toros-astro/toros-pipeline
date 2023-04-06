""" This class of functions is primarily for image manipulation or coordinate information"""
import numpy as np
import pandas as pd
from astropy.wcs import WCS
from astropy.io import fits


class Images:

    @staticmethod
    def get_coords(mast_head):
        """ This function will generate the angular distance from the center of the master frame, to the edge of
        the ccd, so the TIC can be searched.

        :parameter mast_head - The header of the image

        :return coords - The center coordinates of the frame, and the distance to the edge"""

        # get the relative positions of the edges of the frame, and the center of the frame
        # SHOULD THIS BE UPDATED TO NOT BE SPECIFIC TO TESS?
        x = [0., 0., 1024., 2048., 2048.]
        y = [0., 2048., 1024., 0., 2048.]

        # convert the x & y positions to ra & dec positions based on the image header
        w = WCS(mast_head)
        ra, dec = w.wcs_pix2world(x, y, 0)

        dist = np.zeros(5)
        for zz in range(0, 4):
            dist[zz] = np.arccos(np.sin(dec[2] * np.pi / 180.) * np.sin(dec[zz] * np.pi / 180.) +
                                 np.cos(dec[2] * np.pi / 180.) * np.cos(dec[zz] * np.pi / 180.) *
                                 np.cos(ra[2] * np.pi / 180. - ra[zz] * np.pi / 180.)) * 180. / np.pi
        mx_dist = np.nanmax(dist)

        return mx_dist, ra[2], dec[2]

    @staticmethod
    def get_fits(path):
        """ This function will read in a TESS FFI fits file into a header and data file.

        :parameter path - A string with the path to the file
        :parameter clip - An optional string to clip the over read section of the file

        :return header - The header is returned for the fits file
        :return data - The data is returned as a data frame
        """

        # get the hdu of the data
        hdu = fits.open(path)

        # read in the data of the fits file
        dt = hdu[1].data
        df1 = pd.DataFrame(dt)
        df = df1.apply(lambda x: x.values.byteswap().newbyteorder())

        # close the hdu of the file
        hdu.close()

        return df
