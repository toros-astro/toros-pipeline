""" This class is used for basic functions to be used by the broker: S/N calculations, field generation, etc."""

from config import Configuration
import numpy as np
import pandas as pd
import healpy as hp
from libraries.utils import Utils
from astropy.coordinates import SkyCoord, AltAz, HADec
from astroplan import Observer
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.time import Time
from datetime import timedelta
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from ligo.skymap.postprocess import crossmatch


class Priority:

    @staticmethod
    def sort_toros_fields_skymap(toros_fields, skymap):
        """ This function will generate a prioritized list of toros fields based on the provided SkyMap when no NED
        result is available.

        :parameter toros_fields - This is the standard set of TOROS fields
        :parameter skymap - The skymap output directly from the NED alert

        :return selected_fields - A list of fields within the skymap banana sorted by declination and priority
        """

        # get the resolution of the healpy map
        npix = len(skymap)
        nside = hp.npix2nside(npix)

        # copy the data frame so we don't ruin the OG data
        selected_fields = toros_fields.copy().reset_index(drop=True)

        # get the coordinate values for each toros field in the appropriate units
        field_coords = SkyCoord(selected_fields.ra, selected_fields.dec, unit=u.deg)

        # now crossmatch the TOROS fields with the sky map probabilities
        cross_match = crossmatch(skymap, field_coords)

        # now move through each TOROS field and get the integrated area within the field area
        selected_fields['prob'] = cross_match.probdensity

        # sort all TOROS fields based on the probability strip
        selected_fields = selected_fields.sort_values(by='prob', ascending=False).copy().reset_index(drop=True)

        return selected_fields

    @staticmethod
    def ephemeris_observable(field_ephemeris, field_period, field_tolerance, time_of_exp):

        """ This function will use a fields ephemeris date to determine whether or not the field is observable. If it
        is, it will return a 1.

        :parameter field_ephemeris - The ephemeris in JD of the target of opportunity
        :parameter field_period - The period of the target of opportunity (for exact date, the period should be 1)
        :parameter field_tolerance - The tolerance in phase, that the observation is allowed to occur
        :parameter time_of_exp - The start of the exposure in JD

        :return observe - either 0/1 if the field can be observed right now
        """

        # does the field require an observation at anytime on this date?
        if field_period == -1:
            # check to see if we are on the current observation date
            if int(field_ephemeris) == int(time_of_exp):
                observe = 1

        # does the field require a phased observation?
        elif field_period > 0:
            # get the current phase on the field epehmeris
            phase = ((time_of_exp - field_ephemeris) / field_period) % 1

            # determine if the current phase is within the suggested tolerance
            if np.abs(phase) <= field_tolerance:
                observe = 1
            else:
                observe = 0

        # the field requires not ephemeris observation (i.e. normal science or survey field)
        else:
            observe = 0

        return observe

    @staticmethod
    def return_toros_fields(target_list, toros_fields):
        """ This function will return the unique toros fields for a given set of coordiantes

        :parameter target_list - A pandas data frame with the target list
        :parameter toros_fields - A pandas data frame with the toros fields

        :return toros - A data frame with the fields to observe in the set of coordinates
        """

        # set up initial values for each field
        target_list['toros_field'] = 'N/A'
        target_list['toros_ra'] = 0
        target_list['toros_dec'] = 90

        Utils.log("Aligning targets with TOROS fields...", "info")

        # iterate through each target and make the target selection
        for idx, row in target_list.iterrows():

            if idx % 10 == 0:
                Utils.log(str(len(target_list) - idx) + " object remain. Finding next 10...", "info")

            # find the right toros field
            field = Priority.find_toros_field(toros_fields, row.ra, row.dec)

            if isinstance(field, str) is False:
                # update the target list
                target_list.loc[idx, 'toros_field'] = field['toros_field_id']
                target_list.loc[idx, 'toros_ra'] = field['ra']
                target_list.loc[idx, 'toros_dec'] = field['dec']
            else:
                Utils.log("No TOROS field found for " + str(row['objname']), "info")

        # make sure there are no duplicates
        grouped_targets = target_list.groupby(['toros_field',
                                               'toros_ra',
                                               'toros_dec'])['objname'].apply(list).reset_index(name='galaxies')

        # pull out the maximum probability for each field
        priority_fields = target_list.groupby('toros_field')['dP_dV'].max()

        # merge the two lists to generate the final list
        final_list = pd.merge(grouped_targets,
                              priority_fields,
                              on='toros_field').reset_index(drop=True).sort_values(by='dP_dV', ascending=False)

        Utils.log("Finished aligning targets to TOROS fields.", "info")

        return final_list

    @staticmethod
    def angular_distance(ra1, dec1, ra2, dec2):
        """ This function will determine the angular distance between two points.

        :parameter ra1 - The right ascension of point 1
        :parameter dec1 - The declination of point 1
        :parameter ra2 - The right ascension of point 2
        :parameter dec2 - the declinatino of point 2

        return - The angular distance of the two points in degrees

        """

        # convert to radians
        ra1_rad = np.deg2rad(ra1)
        dec1_rad = np.deg2rad(dec1)
        ra2_rad = np.deg2rad(ra2)
        dec2_rad = np.deg2rad(dec2)

        # the angular distance in radians
        ang_dist_rad = np.arccos((np.sin(dec1_rad) * np.sin(dec2_rad)) +
                                 (np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra1_rad - ra2_rad)))

        # the angular distance in degrees
        ang_dist_deg = np.rad2deg(ang_dist_rad)

        return ang_dist_deg

    @staticmethod
    def find_toros_field(toros_fields, ra, dec):
        """ This function will select the appropriate TOROS field based on the RA/Dec of the target you are interested
        in observing.

        :parameter toros_fields - The set of TOROS survey fields
        :parameter ra - The right ascension of the target in degrees
        :parameter dec - The declination of the target in degrees

        return field - The TOROS field name to be observed for the specific science
        """

        # get the angular distance between the given position and the toros fields
        ang_dist = toros_fields.apply(lambda x: Priority.angular_distance(x.ra, x.dec, ra, dec), axis=1)

        if np.min(ang_dist) <= 2:
            # pull out the field row
            field = toros_fields.loc[np.argmin(ang_dist)].copy()
        else:
            field = 'No field found.'

        return field

    @staticmethod
    def toros_field_selector(survey_fields):
        """ This function will select the best field for observation at the current UTC, moon, ligo alert etc.

        :parameter survey_fields - The full list of survey fields with current observing history.

        :return priority_field - The TOROS field with the best priority
        """

        # set up the Macon Parameters
        observatory = EarthLocation(lat=Configuration.TOROS_LATITUDE * u.deg,
                                    lon=Configuration.TOROS_LONGITUDE * u.deg,
                                    height=Configuration.TOROS_ELEVATION * u.m)

        # for the observing position to be at the toros location
        toros = Observer(location=observatory, name='TOROS', timezone='America/Argentina/Salta')

        # get the field positions
        field_c = SkyCoord(ra=survey_fields['ra'].to_list(), dec=survey_fields['dec'].to_list(),
                           frame='icrs', unit='deg')

        # get the current time
        time_now = Time.now()

        # get the moon information
        if toros.moon_illumination(time_now) >= 0.5:
            moon_bright = 1
        else:
            moon_bright = 0

        # the priority order for field types
        field_priority = ['lvc', 'main_survey']

        # get the current moon
        moon = toros.moon_altaz(time=time_now)
        obj_altaz = field_c.transform_to(AltAz(obstime=time_now, location=observatory))
        obj_hadec = field_c.transform_to(HADec(obstime=time_now, location=observatory))
        targets_exp = toros.altaz(time_now, field_c)

        # make sure that the moon is above the horizon before we settle on dark time
        if (moon.alt.value > 0) & (moon_bright == 1):
            moon_phase = 1
        else:
            moon_phase = 0

        # update the hour angles for the expected exposure time
        survey_fields['field_ha_exp'] = np.abs(obj_hadec.ha.value)
        survey_fields['airmass_exp'] = obj_altaz.secz.value
        survey_fields['altitude_exp'] = obj_altaz.alt.value
        survey_fields['azimuth_exp'] = obj_altaz.az.value
        survey_fields['moon_dist'] = survey_fields.apply(lambda x: Priority.angular_distance(
            x.azimuth_exp, x.altitude_exp, moon.az.value, moon.alt.value), axis=1)

        # get the possible fields based on observing cadences
        possible_fields = survey_fields[((survey_fields['field_ha_exp'] < 5) |
                                        ((survey_fields['airmass_exp'] >= 1) &
                                         (survey_fields['airmass_exp'] <= 2))) &
                                        (survey_fields['moon_phase'] == moon_phase) &
                                        (survey_fields['altitude_exp'] > 0) &
                                        (survey_fields['observations'] ==
                                         survey_fields['observations'].min())].sort_values(by=['dec', 'moon_dist'],
                                                                                           ascending=[True, True])

        # look for ligo alerts, then science fields, then survey fields, then high cadence fields
        for prio in field_priority:
            field_chk = possible_fields[possible_fields['program'] == prio]

            if len(field_chk) > 1:
                field = field_chk.iloc[0]
                break
        return field

    @staticmethod
    def toros_night_field_selector(survey_fields, science_fields, hc_fields):
        """ This function will select fields for tonight based on the current UTC, moon, ligo alert etc.

        :parameter survey_fields - The full set of survey fields with current observing history.
        :parameter science_fields - The full set of science fields with priority
        :parameter hc_fields - The set of high cadence fields with priority & cadence requirements
        :parameter ligo_fields - The set of alert ligo fields

        :return obs_list - The list of fields to observe tonight (in order)

        """

        # merge the lists together to make a giant list selection
        toros_fields = survey_fields.append([science_fields, hc_fields]).reset_index(drop=True)

        # set up the Macon Parameters
        location = EarthLocation.from_geodetic(Configuration.TOROS_LONGITUDE * u.deg,
                                               Configuration.TOROS_LATITUDE * u.deg,
                                               Configuration.TOROS_ELEVATION * u.m)

        # for the observing position to be at the toros location
        toros = Observer(location=location, name='TOROS', timezone='America/Argentina/Salta')
        field_c = SkyCoord(ra=toros_fields.ra.to_list(), dec=toros_fields.dec.to_list(), frame='icrs', unit='deg')

        # get the time for the evening hours
        jd_start = toros.tonight(horizon=-6 * u.deg)[0]
        jd_end = toros.tonight(horizon=-6 * u.deg)[1]

        # get the moon information
        if toros.moon_illumination(jd_start) >= 0.5:
            moon_bright = 1
        else:
            moon_bright = 0

        time_of_exp = jd_start
        field_priority = ['ligo', 'science', 'survey', 'high_cadence']

        while time_of_exp < jd_end:
            moon = toros.moon_altaz(time=time_of_exp)
            targets_exp = toros.altaz(time_of_exp, field_c)

            # make sure that the moon is above the horizon before we settle on dark time
            if (moon.alt.value > 0) & (moon_bright == 1):
                moon_phase = 1
            else:
                moon_phase = 0

            # update the hour angles for the expected exposure time
            toros_fields['field_ha_exp'] = toros.target_hour_angle(time_of_exp, field_c)
            toros_fields['airmass_exp'] = targets_exp.secz
            toros_fields['altitude_exp'] = targets_exp.alt
            toros_fields['azimuth_exp'] = targets_exp.az
            toros_fields['horizon_exp'] = toros.target_is_up(time_of_exp, field_c)

            # get the possible fields based on observing cadences
            possible_fields = toros_fields[((toros_fields.field_ha_exp < 5) |
                                            ((toros_fields.airmass_exp >= 1) &
                                             (toros_fields.airmass_exp <= 2))) &
                                           (toros_fields.moon_phase == moon_phase) &
                                           (toros_fields.observations == toros_fields.observations.min())].sort_values(by=['dec'],
                                                                                                ascending=True)

            # look for ligo alerts, then science fields, then survey fields, then high cadence fields
            for prio in field_priority:
                field_chk = possible_fields[possible_fields['field_type'] == prio]

                if len(field_chk) > 1:

                    # get the alt and az of the positions of the fields, assuming the moon is up
                    #if moon_phase == 1:
                    #    field_chk['moon_dist'] = field_chk.apply(lambda x: Priority.angular_distance(
                    #        x.azimuth_exp, x.altitude_exp, moon.az.value, moon.alt.value), axis=1)
                    #else:
                    #    field_chk['moon_dist'] = Configuration.MOON_DISTANCE + 10.

                    #if len(field_chk[field_chk['moon_dist'] > Configuration.MOON_DISTANCE]) > 0:

                    field = field_chk.index[0]
                    toros_fields.loc[field, 'observations'] += 1

                    next_exp = time_of_exp + timedelta(days=toros_fields.loc[field, 'cadence'])
                    time_of_exp = next_exp
                    break

        plt.figure(figsize=(32, 24))
        plt.subplot(111, projection="aitoff")
        toros_fields['ra_wrap'] = np.where(toros_fields.ra > 180, toros_fields.ra - 360, toros_fields.ra)
        plt.scatter(toros_fields[toros_fields.moon_phase == 0].ra_wrap * np.pi / 180,
                    toros_fields[toros_fields.moon_phase == 0].dec * np.pi / 180, s=80, edgecolors='b', facecolors='none', label='Dark Time', alpha=0.2)
        plt.scatter(toros_fields[toros_fields.moon_phase == 1].ra_wrap * np.pi / 180,
                    toros_fields[toros_fields.moon_phase == 1].dec * np.pi / 180,  s=80, edgecolors='g', facecolors='none',label='Bright Time', alpha=0.2)
        plt.scatter(toros_fields[toros_fields.field_type == 'science'].ra_wrap * np.pi / 180,
                    toros_fields[toros_fields.field_type == 'science'].dec * np.pi / 180,  s=80, edgecolors='r', facecolors='none',label='Science Fields', alpha=0.2)
        plt.scatter(toros_fields[toros_fields.field_type == 'high_cadence'].ra_wrap * np.pi / 180,
                    toros_fields[toros_fields.field_type == 'high_cadence'].dec * np.pi / 180,  s=80, edgecolors='k', facecolors='none',label='High Cadence', alpha=0.2)
        plt.grid(True)

        plt.scatter(toros_fields[(toros_fields.observations >= 1) & (toros_fields.field_type == 'science')].ra_wrap * np.pi / 180,
                    toros_fields[(toros_fields.observations >= 1) & (toros_fields.field_type == 'science')].dec * np.pi / 180.,
                    c='r')
        plt.scatter(toros_fields[(toros_fields.observations >= 1) & (toros_fields.field_type == 'high_cadence')].ra_wrap * np.pi / 180,
                    toros_fields[(toros_fields.observations >= 1) & (toros_fields.field_type == 'high_cadence')].dec * np.pi / 180.,
                    c='k')
        plt.scatter(toros_fields[(toros_fields.moon_phase == 1) & (toros_fields.observations >= 1) & (toros_fields.field_type == 'survey')].ra_wrap * np.pi / 180,
                    toros_fields[(toros_fields.moon_phase == 1) & (toros_fields.observations >= 1) & (toros_fields.field_type == 'survey')].dec * np.pi / 180.,
                    c='g')
        plt.scatter(toros_fields[(toros_fields.moon_phase == 0) & (toros_fields.observations >= 1) & (toros_fields.field_type == 'survey')].ra_wrap * np.pi / 180,
                    toros_fields[(toros_fields.moon_phase == 0) & (toros_fields.observations >= 1) & (toros_fields.field_type == 'survey')].dec * np.pi / 180.,
                    c='b')
        plt.xlabel('Right Ascension [deg]')
        plt.ylabel('Declination [deg]')
        plt.legend()
        plt.show()
        plt.close()

        return

    @staticmethod
    def toros_survey_simulator(toros_fields, ligo_fields):
        """ This function will select the fields to observe based on the current UTC.

        :parameter toros_fields - The set of TOROS survey fields
        :parameter ligo_fields - The set of LIGO fields from NED to be used for observing

        :return fields_to_observe - The selection of fields to observe, in priority order
        """

        toros_fields['observations'] = 0
        toros_fields['moon_phase'] = np.where((toros_fields.b < 15) & (toros_fields.b > -15), 1, 0)

        # set up the Macon Parameters
        location = EarthLocation.from_geodetic(Configuration.TOROS_LONGITUDE * u.deg,
                                               Configuration.TOROS_LATITUDE * u.deg,
                                               Configuration.TOROS_ELEVATION * u.m)

        # for the observing position to be at the toros location
        toros = Observer(location=location, name='TOROS', timezone='America/Argentina/Salta')
        field_c = SkyCoord(ra=toros_fields.ra.to_list(), dec=toros_fields.dec.to_list(), frame='icrs', unit='deg')

        jd_st = toros.tonight()[0]
        for ngt in range(0, 365):

            dte = jd_st+timedelta(days=ngt)

            # get the time for the evening hours
            jd_start = toros.tonight(time=dte, horizon=-6 * u.deg)[0]
            jd_end = toros.tonight(time=dte, horizon=-6 * u.deg)[1]

            # get the moon information
            if toros.moon_illumination(jd_start) >= 0.5:
                moon_bright = 1
            else:
                moon_bright = 0

            time_of_exp = jd_start
            total_exp_time = Configuration.TOTAL_EXPOSURE / 60. / 60. / 24.

            n_iters = 0

            while time_of_exp < jd_end:
                next_exp = time_of_exp + timedelta(days=total_exp_time)
                moon = toros.moon_altaz(time=time_of_exp)
                targets_exp = toros.altaz(time_of_exp, field_c)
                targets_next_exp = toros.altaz(next_exp, field_c)

                # make sure that the moon is above the horizon before we settle on dark time
                if (moon.alt.value > 0) & (moon_bright == 1):
                    moon_phase = 1
                else:
                    moon_phase = 0

                # update the hour angles for the expected exposure time
                toros_fields['field_ha_exp'] = toros.target_hour_angle(time_of_exp, field_c)
                toros_fields['airmass_exp'] = targets_exp.secz
                toros_fields['altitude_exp'] = targets_exp.alt
                toros_fields['azimuth_exp'] = targets_exp.az
                toros_fields['horizon_exp'] = toros.target_is_up(time_of_exp, field_c)

                # update the hour angles for the next expected exposure time
                toros_fields['field_ha_next_exp'] = toros.target_hour_angle(next_exp, field_c)
                toros_fields['airmass_next_exp'] = targets_next_exp.secz
                toros_fields['altitude_next_exp'] = targets_next_exp.alt
                toros_fields['azimuth_next_exp'] = targets_next_exp.az
                toros_fields['horizon_next_exp'] = toros.target_is_up(next_exp, field_c)

                # get the possible fields based on observing cadences
                possible_fields = toros_fields[(((toros_fields.field_ha_exp < 5) |
                                                 (toros_fields.airmass_exp > 21)) &
                                                (toros_fields.airmass_exp < 2)) &
                                               (((toros_fields.field_ha_next_exp < 5) |
                                                (toros_fields.airmass_next_exp > 21)) &
                                               (toros_fields.airmass_next_exp < 2)) &
                                               (toros_fields.observations == toros_fields.observations.min()) &
                                               (toros_fields.moon_phase == moon_phase)].sort_values(by=['dec'],
                                                                                                    ascending=True)

                # get the alt and az of the positions of the fields, assuming the moon is up
                if moon_phase == 1:
                    possible_fields['moon_dist'] = possible_fields.apply(lambda x:
                                                                         Priority.angular_distance(
                                                                             x.azimuth_exp,
                                                                             x.altitude_exp,
                                                                             moon.az.value,
                                                                             moon.alt.value),
                                                                         axis=1)
                else:
                    possible_fields['moon_dist'] = Configuration.MOON_DISTANCE + 10.

                if len(possible_fields[possible_fields['moon_dist'] > Configuration.MOON_DISTANCE]) > 0:

                    field = possible_fields[possible_fields['moon_dist'] > Configuration.MOON_DISTANCE].index[0]
                    toros_fields.loc[field, 'observations'] += 1

                    time_of_exp = next_exp

                else:
                    time_of_exp = next_exp

            if (ngt % 30 == 0) & (ngt > 0):
                plt.figure()
                plt.subplot(111, projection="aitoff")
                toros_fields['ra_wrap'] = np.where(toros_fields.ra > 180, toros_fields.ra - 360, toros_fields.ra)
                plt.scatter(toros_fields[toros_fields.moon_phase == 0].ra_wrap * np.pi / 180,
                            toros_fields[toros_fields.moon_phase == 0].dec * np.pi / 180, c='b', label='Dark Time')
                plt.scatter(toros_fields[toros_fields.moon_phase == 1].ra_wrap * np.pi / 180,
                            toros_fields[toros_fields.moon_phase == 1].dec * np.pi / 180, c='g', label='Bright Time')
                plt.grid(True)

                plt.scatter(toros_fields[toros_fields.observations >= 1].ra_wrap * np.pi / 180,
                            toros_fields[toros_fields.observations >= 1].dec * np.pi / 180.,
                            c='r', label='Observed Field')
                plt.xlabel('Right Ascension [deg]')
                plt.ylabel('Declination [deg]')
                plt.legend()
                plt.show()
                plt.close()
        return

    @staticmethod
    def toros_field_generator(field_size):
        """ This function will generate the TOROS fields based on the provided separation of each field.

        :parameter field_size - The desired field size, in degrees, required for each field center

        :return field_list - The field list is both printed and returned as a pandas data frame

        """

        if os.path.isfile(Configuration.ANALYSIS_DIRECTORY + 'toros_fields.dat') is False or \
                Configuration.FEILD_GENERATION == 'Y':
            Utils.log("Now generating TOROS fields.", "info")

            # set up the TOROS fields data frames and start with the first field
            c = SkyCoord(ra=0.e0, dec=-90e0, frame='icrs', unit='deg')
            toros_fields = pd.DataFrame(data=[['00.000', 0e0, -90e0,
                                               c.galactic.l.to_value(), c.galactic.b.to_value(),
                                               'main_survey', 300., 1, 0, 0]],
                                        columns=['toros_field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time',
                                                 'cadence', 'ephemeris', 'period'])

            # set up the TOROS field of view
            toros_fov = Configuration.PIXEL_SIZE * Configuration.NUM_PIXELS / 3600.
            toros_field_number = int(np.ceil((90 + Configuration.TOROS_DEC_LIMIT) / field_size))

            # set up variables necessary for geometry
            deg_to_rad = np.pi / 180.  # degree to radians conversion
            field_sep = toros_fov * np.sqrt(2.2 ** 2 + 3.73 ** 2) / 5.0

            # separation of the fields in declination
            declination_strips = -90 + np.arange(0, toros_field_number) * field_sep

            # now loop through and generate the fields
            eo = 0
            for idx in range(1, toros_field_number):
                nfr = np.ceil(360. * np.cos(declination_strips[idx] * deg_to_rad) / field_size)
                ra_sep = 360. / nfr

                if eo == 1:
                    ra_off = 0.5 * ra_sep
                    eo = 0
                else:
                    ra_off = 0.
                    eo = 1

                for idy in range(0, int(nfr)):
                    # set up the field name with the first hex field
                    if len(hex(idx).split('x')[1]) == 1:
                        field_1 = '0' + hex(idx).split('x')[1]
                    else:
                        field_1 = hex(idx).split('x')[1]

                    # set up the second hex field
                    if len(hex(idy).split('x')[1]) == 1:
                        field_2 = '00' + hex(idy).split('x')[1]
                    elif len(hex(idy).split('x')[1]) == 2:
                        field_2 = '0' + hex(idy).split('x')[1]
                    else:
                        field_2 = hex(idy).split('x')[1]

                    # get the galactic coordinates of the field
                    c_idy = SkyCoord(ra=idy * ra_sep + ra_off, dec=declination_strips[idx], frame='icrs', unit='deg')
                    if (c_idy.galactic.b.to_value() < 15) & (c_idy.galactic.b.to_value() > -15):
                        moon_phase = 1
                    else:
                        moon_phase = 0

                    # set up the series for appending
                    field = pd.Series(data=[field_1 + '.' + field_2,
                                            idy * ra_sep + ra_off,
                                            declination_strips[idx],
                                            c_idy.galactic.l.to_value(),
                                            c_idy.galactic.b.to_value(),
                                            'main_survey', 300., 1, 0, 0, 0, moon_phase],
                                      index=['toros_field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time',
                                             'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

                    # append the series
                    toros_fields = toros_fields.append(field, ignore_index=True)

            toros_fields.to_csv(Configuration.ANALYSIS_DIRECTORY + 'toros_fields.dat',
                                sep=' ', header=True, index=False, float_format='%.3f')

        else:
            # if the file exists already, then just read the field list in
            toros_fields = pd.read_csv(Configuration.ANALYSIS_DIRECTORY + 'toros_fields.dat', header=0, sep=' ')

        return toros_fields

    @staticmethod
    def signal_to_noise(exp_time, ab_brightness, number_of_reads):
        """ This function will calculate the TOROS signal to noise ratio for the given exposure time and
        ab_brightness of the object.

        :parameter exp_time - The exposure time in seconds of the observation
        :parameter ab_brightness - The AB brightness of the object as a scalar or array
        :parameter number_of_reads - The number of reads for the exposure (default = 1).

        :return snr_total - The total SNR of the exposure is returned
        """

        # calculate the telscopes area
        telescope_area = Configuration.TOROS_MIRROR_R ** 2 * np.pi
        telescope_aperture = np.pi * ((2.04 * Configuration.SEEING) / 2.0) ** 2

        # calculate flux from a zero magnitude star
        bandpass_hz = (299792458000000000. / (np.array(Configuration.CWL) - (np.array(Configuration.BP) / 2.))) -\
                      (299792458000000000. / (np.array(Configuration.CWL) + (np.array(Configuration.BP) / 2.)))
        hc_lambda = 6.626E-34 * 299800000000000000.0 / np.array(Configuration.CWL)
        zero_mag = 3.631E-23 * bandpass_hz / hc_lambda  # in photon/sec/m^2
        zero_mag_photons = zero_mag * telescope_area * exp_time
        zero_mag = zero_mag * 1000.0 / np.array(Configuration.BP)  # in photon/sec/m^2/micron

        # calculate the sky surface brightness (# joules/sec/m^2/micron/arcsec^2)
        # sky_surface_brightness = (10.0 ** (np.array(Configuration.SKY) / (-2.5))) * zero_mag * hc_lambda
        # sky_surface_brightness = sky_surface_brightness * 42545250225  # W/m^2/micron/sr
        # sky_surface_brightness = sky_surface_brightness / 10000.  # W/cm^2/micron/sr

        # calculate the signal from various sources
        signal_0mag_source = zero_mag_photons * Configuration.TOTAL_THROUGHPUT * np.array(Configuration.ATMOSPHERE_QE)
        # signal_20mag_source = signal_0mag_source * (10.0 ** (20.0/(-2.5)))
        signal_specific_source = signal_0mag_source * (10.0 ** (ab_brightness/(-2.5)))
        signal_sky = (10 ** (np.array(Configuration.SKY) / (-2.5)) * telescope_aperture * signal_0mag_source)
        readnoise_aperture = Configuration.READOUT_NOISE ** 2 * (telescope_aperture / (0.468 ** 2))

        # calculate the snr
        snr_per_band = signal_specific_source / ((signal_specific_source +
                                                  signal_sky +
                                                  number_of_reads * readnoise_aperture) ** (1.0 / 2.0))

        # total signal to noise
        snr = np.sum(signal_specific_source) / ((np.sum(signal_specific_source) + np.sum(signal_sky) +
                                                 number_of_reads * readnoise_aperture) ** (1.0 / 2.0))

        return snr_per_band, snr
