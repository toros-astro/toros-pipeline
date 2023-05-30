""" This wrapper script will be the script which executes operations for each night for the TOROS hardware. This
includes connections to the telescope, slewing, homing, exposing, updating the observations table etc. """
from astroplan import Observer
from astropy.coordinates import EarthLocation
from config import Configuration
from libraries.utils import Utils
from libraries.dbaccess import DBaccess
from libraries.alerts import Alerts
from libraries.telescope import Telescope
from libraries.pwi4_client import PWI4
from libraries.priority import Priority
import time
import astropy.units as u
from astropy.time import Time
import subprocess


# connect to PWI4
Utils.log("Connecting to PWI4...", "info")
# subprocess.Popen("C:\\Program Files (x86)\\PlaneWave Instruments\\PlaneWave Interface 4\\PWI4.exe")
# subprocess.Popen("C:\\Program Files (x86)\\Common Files\\ASCOM\\DeviceHub\\ASCOM.DeviceHub.exe")
# time.sleep(5)
pwi4 = PWI4()
Utils.log("Connected to PWI4.", "info")

# connect to the telescope mount
Telescope.connect_to_mount()

# connect to the dome

# connect to the camera

# set up observatory positions
# check the current nighttime status
# set up the Macon Parameters
location = EarthLocation.from_geodetic(Configuration.TOROS_LONGITUDE * u.deg,
                                       Configuration.TOROS_LATITUDE * u.deg,
                                       Configuration.TOROS_ELEVATION * u.m)

# for the observing position to be at the toros location
toros = Observer(location=location, name='TOROS', timezone='America/Argentina/Salta')

# The operations script will run at all times with a sleep of 5 minutes to check if the Sun is low enough to observe
while True:

    # make sure the time is before civil twilight
    if Time.now().jd > toros.twilight_evening_civil(Time.now()).value:
        Utils.log("Operations beginning for tonight.", "info")

        while True:
            # get the target to observe
            try:
                # pull the appropriate query from the SQL file
                obs_query = DBaccess.get_query(Configuration.QUERIES_DIRECTORY + 'obs_list.sql')

                # pull the target from the database
                obs_target = DBaccess.query_torosdb(obs_query, Configuration.MACHINE)
            except ConnectionError:
                # send an alert that connection has been lost
                Alerts.alert_toros_team("DATABASE")
                Utils.log("No connection to TOROSdb! Pulling from legacy file.", "info")

                # make the toros field list
                survey_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)

                # select the best field to observe
                obs_target = Priority.toros_field_selector(survey_fields)

            # slew to the target
            pwi4.mount_goto_ra_dec_j2000(obs_target['ra'], obs_target['dec'])
            Utils.log("Starting slew to field " + str(obs_target['toros_field_id']) + ". Forcing mount to stop.", "info")
            while True:
                # check the mount status
                s = pwi4.status()
                # if the mount is not slewing then break and sleep
                if not s.mount.is_slewing:
                    break
                time.sleep(0.2)

            Utils.log("Slew complete to field " + str(obs_target['toros_field_id']) + ". Forcing mount to stop.", "info")
            pwi4.mount_stop()

            # expose on the target

            # break if past morning civil twilight
            if Time.now.jd > toros.twilight_morning_civil(Time.now()).value:
                Utils.log("It is now past Morning Civil Twilight. Shutting down TOROS operations.", "info")

                # park telescope
                Telescope.find_home()

                # close dome

                # stop operations
                break

    # wait 5 minutes before checking if the sun has set
    Utils.log("The Sun is up. TOROS operations have not begun.", "info")
    time.sleep(300)
