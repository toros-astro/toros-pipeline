""" This set of library functions will control the telescope primarily. These include connecting/disconnecting to the
mount, slewing the mount to objects, and finding the home position. """
from libraries.utils import Utils
from libraries.pwi4_client import PWI4
import time


class Telescope:

    @staticmethod
    def connect_to_mount():
        """ This function will connect to the mount and wait if not connecting properly. No inputs are needed."""
        s = PWI4.status()
        if not s.mount.is_connected:
            Utils.log('Connecting to the mount...', "info")
            PWI4.mount_connect()
            while not PWI4.status().mount.is_connected:
                time.sleep(1)
        Utils.log("Mount now connected.", "info")

    @staticmethod
    def enable_motors():
        """ This function will enable the required motors on the mount."""
        Utils.log("Enabling mount motors...", "info")
        PWI4.mount_enable(0)
        PWI4.mount_enable(1)
        while True:
            status = PWI4.status()
            if status.mount.axis0.is_enabled and status.mount.axis1.is_enabled:
                break
            time.sleep(1)
        Utils.log("Mount moters now enabled.", "info")

    @staticmethod
    def find_home():
        Utils.log("Sending telescope to home positon.", "info")
        PWI4.mount_find_home()
        last_axis0_pos_degs = -99999
        last_axis1_pos_degs = -99999
        while True:
            status = PWI4.status()
            delta_axis0_pos_degs = status.mount.axis0.position_degs - last_axis0_pos_degs
            delta_axis1_pos_degs = status.mount.axis1.position_degs - last_axis1_pos_degs

            if abs(delta_axis0_pos_degs) < 0.001 and abs(delta_axis1_pos_degs) < 0.001:
                break

            last_axis0_pos_degs = status.mount.axis0.position_degs
            last_axis1_pos_degs = status.mount.axis1.position_degs

            time.sleep(1)
        Utils.log("Telescope now at the home position.", "info")
