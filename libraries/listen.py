""" This is the listener function. This will listen for GCN alerts, generate the field list,
and alert the TOROS team."""
from gcn_kafka import Consumer
from libraries.utils import Utils
from libraries.alerts import Alerts
from libraries.priority import Priority
from config import Configuration
import time
import os
import logging
from ligo.skymap.io import read_sky_map
import healpy as hp
from astropy.io import fits
logging.getLogger("healpy").setLevel(logging.WARNING)
logging.getLogger("gcn").setLevel(logging.WARNING)

# credentials for NASA (B4rn d00r)

# config = {'group.id': 'TOROS',
#           'auto.offset.reset': 'earliest',
#           'enable.auto.commit': False}

#consumer = Consumer(config=config,
consumer = Consumer(client_id='276gq6h6e75qnuoapm1uvgeo9d',
                    client_secret='1qt4th4b0damb9rd1tn2nvukujbnb5mv508thgkregaiur2h9ni7')

# Subscribe to topics and receive alerts
consumer.subscribe(['gcn.classic.text.LVC_INITIAL',
                    'gcn.classic.text.LVC_PRELIMINARY',
                    'gcn.classic.text.LVC_RETRACTION'])
while True:
    for message in consumer.consume(timeout=1):
        # convert the message to UNICODE
        message_contents = message.value().decode('utf-8').split('\n')

        keys = []
        items = []

        # split the message into keys
        for content in message_contents:
            if len(content.split(': ')) > 1:
                keys.append(content.split(': ')[0].strip())
                items.append(content.split(': ')[1].strip())

        params = dict(zip(keys, items))
        alert_type = params['NOTICE_TYPE']
        event_name = params['TRIGGER_NUM']

        # make sure the alert is not a retraction
        if params['NOTICE_TYPE'] == 'LVC Retraction':
            Utils.log("Alert for " + event_name + " was retracted. Deleting sky map files.", "info")
            os.system('rm ' + Configuration.ANALYSIS_DIRECTORY + event_name + '_SkyMap_toros_fields.txt')

            Utils.log("Alerting team of retraction.", "info")
            Alerts.alert_toros_team(alert_type=alert_type, event_name=event_name)
            continue

        Utils.log("Alert for " + event_name + " was heard! Alerting TOROS team.", "info")

        # alert the TOROS team
        Alerts.alert_toros_team(alert_type=alert_type, event_name=event_name)

        # send in a request to NED to get the galaxy list
        ned_result = Alerts.query_ned(event_name)

        if ned_result is None:
            Utils.log("No NED galaxy priority yet. Extracting TOROS fields from SkyMap in the meantime.", "info")
            wait_iter = 0
            while wait_iter < Configuration.LISTEN_NED_WAIT:
                if 'SKYMAP_FITS_URL' in params:
                    # check to see if field list from skymap exists
                    if os.path.isfile(Configuration.ANALYSIS_DIRECTORY +
                                      event_name + '_SkyMap_toros_fields.txt') is False:

                        # Read the HEALPix sky map and the FITS header.
                        skymap = read_sky_map(params['SKYMAP_FITS_URL'], moc=True)

                        # pull in the toros fields
                        toros_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)

                        # generate a ranked list of toros fields within the skymap
                        toros_fields_prio = Priority.sort_toros_fields_skymap(toros_fields, skymap, event_name)
                        toros_fields_prio.to_csv(Configuration.ALERTS_DIRECTORY +
                                                 event_name +
                                                 '_SkyMap_toros_fields.txt',
                                                 sep=' ', index=False)

                        Utils.log("TOROS field list generated based on SkyMap from GCN. "
                                  "Waiting 5 minutes to re-query NED.", "info")

                        # increase the wait iterations
                        wait_iter += 1
                        time.sleep(300)

                        # requery NED
                        ned_result = Alerts.query_ned(event_name)

                        if ned_result is None:
                            Utils.log("NED queried again after " + str(wait_iter * 5) +
                                      " minutes. Still no priority list.", "info")
                        else:
                            Utils.log("NED result obtained. Identifying necessary TOROS fields.", "info")
                            ligo_fields = Priority.return_toros_fields(ned_result, toros_fields)
                            ligo_fields.to_csv(Configuration.ALERTS_DIRECTORY +
                                               event_name +
                                               '_NED_toros_fields.txt',
                                               sep=' ', index=False)
                            continue
                    else:
                        Utils.log("TOROS field list already generated from SkyMap. Skipping for now and "
                                  "waiting 5 minutes to re-query NED.", "info")
                        wait_iter += 1
                        time.sleep(300)
                        ned_result = Alerts.query_ned(event_name)

                        if ned_result is None:
                            Utils.log("NED queried again after " + str(wait_iter * 5) +
                                      " minutes. Still no priority list.", "info")
                        else:
                            Utils.log("NED result obtained. Identifying necessary TOROS fields.", "info")
                            toros_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)
                            ligo_fields = Priority.return_toros_fields(ned_result, toros_fields)
                            ligo_fields.to_csv(Configuration.ALERTS_DIRECTORY +
                                               event_name + '_NED_toros_fields.txt',
                                               sep=' ', index=False)
                            continue

                else:
                    Utils.log("SkyMap not found. Check manually.", "info")
                    continue

        else:
            Utils.log("NED result obtained. Identifying necessary TOROS fields.", "info")

            # pull in the toros fields
            toros_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)
            ligo_fields = Priority.return_toros_fields(ned_result, toros_fields)

            ligo_fields.to_csv(Configuration.ALERTS_DIRECTORY + event_name + '_NED_toros_fields.txt',
                               sep=' ', index=False)
            continue

        Utils.log("Maximum iterations for NED reached.", "info")
        continue
