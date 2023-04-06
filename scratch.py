# Define your custom handler here.
# @gcn.include_notice_types(
#    gcn.notice_types.LVC_EARLY_WARNING,
#    gcn.notice_types.LVC_PRELIMINARY,
#    gcn.notice_types.LVC_INITIAL,
#    gcn.notice_types.LVC_UPDATE,
#    gcn.notice_types.LVC_RETRACTION)

def handler(payload, root):

    Utils.log("New " + root.attrib['role'] + ' event detected!', "info")
    if root.attrib['role'] == 'test':
        Utils.log("Test event found. No need to alert team.", "info")
        #return

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']: elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    alert_type = params['AlertType']
    event_name = params['GraceID']

    # make sure the alert is not a retraciton
    if params['AlertType'] == 'Retraction':
        Utils.log("Alert for " + event_name + " was retracted. Deleting sky map files.", "info")
        os.system('rm ' + Configuration.ANALYSIS_DIRECTORY + event_name + '_SkyMap_toros_fields.txt')

        Utils.log("Alerting team of retraction.", "info")
        Alerts.alert_toros_team(alert_type=alert_type, event_name=event_name)
        return

    Utils.log("Alerting TOROS team.", "info")

    # alert the TOROS team
    Alerts.alert_toros_team(alert_type=root.attrib['role'],  event_name=event_name)

    # send in a request to NED to get the galaxy list
    ned_result = Alerts.query_ned(event_name)

    if ned_result is None:
        Utils.log("No NED galaxy priority yet. Waiting to re-query and extracting TOROS fields from SkyMap "
                  "in the meantime.", "info")
        wait_iter = 0
        while wait_iter < Configuration.LISTEN_NED_WAIT:
            if 'skymap_fits' in params:
                # check to see if field list from skymap exists
                if os.path.isfile(Configuration.ANALYSIS_DIRECTORY +
                                  event_name + '_SkyMap_toros_fields.txt') is False:
                    # Read the HEALPix sky map and the FITS header.
                    skymap, header = hp.read_map(params['skymap_fits'], h=True)
                    header = dict(header)

                    # pull in the toros fields
                    toros_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)

                    # generate a ranked list of toros fields within the skymap
                    toros_fields_prio = Priority.sort_toros_fields_skymap(toros_fields, skymap)
                    toros_fields_prio.to_csv(Configuration.ANALYSIS_DIRECTORY +
                                             event_name +
                                             '_SkyMap_toros_fields.txt',
                                             sep=',')

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
                        ligo_fields.to_csv(Configuration.ANALYSIS_DIRECTORY +
                                           event_name +
                                           '_NED_toros_fields.txt',
                                           sep=',')
                        return ned_result
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
                        ligo_fields.to_csv(Configuration.ANALYSIS_DIRECTORY +
                                           event_name + '_NED_toros_fields.txt', sep=',')
                        return
            else:
                Utils.log("SkyMap not found. Check manually.", "info")
                return

        return
    else:
        Utils.log("NED result obtained. Identifying necessary TOROS fields.", "info")
        # pull in the toros fields
        toros_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)
        ligo_fields = Priority.return_toros_fields(ned_result, toros_fields)

        ligo_fields.to_csv(Configuration.ANALYSIS_DIRECTORY + event_name + '_NED_toros_fields.txt', sep=',')

        return

    Utils.log("Maximum iterations for NED reached.", "info")
    return

# Listen for VOEvents until killed with Control-C.
gcn.listen(handler=handler)

from libraries.utils import Utils
from config import Configuration
import pandas as pd
import numpy as np
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# do necessary prep work such as making output directories
files = Utils.get_file_list(Configuration.LIGHTCURVE_DIRECTORY, '.lc')

nfiles = len(files)
mags = np.zeros(nfiles)
errs = np.zeros(nfiles)

for idx, file in enumerate(files):
    lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + '/' + file, header=0, sep=' ')
    mmag, dmag, smag = sigma_clipped_stats(lc.cln, sigma=2.5)
    mags[idx] = dmag
    errs[idx] = smag
    if idx == 105:
        plt.scatter(lc.jd - 2459587, lc.cln-0.88, c='k', marker='.')
        plt.ylabel('Instrumental Magnitude [+offset]', size=15)
        plt.xlabel('JD - 2459587 [days]', size=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylim([16.9, 16.5])
        plt.xlim([0.6, 0.9])
        plt.show()

plt.scatter(mags-0.88, errs, c='k')
plt.yscale('log')
plt.ylabel('rms')
plt.xlabel('Instrumental Magnitude [+offset]')
plt.show()

nms = ['s1', 's2', 's3', 's4', 's5']
corr = pd.DataFrame(columns=nms)
for idx, file in enumerate(files):

    dmag = np.abs(mags-mags[idx])
    ntrs = len(np.argwhere((dmag < 0.1) & (dmag > 0)))
    lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + '/' + file, header=0, sep=' ')
    trd = np.zeros((ntrs, len(lc)))
    idz = 0
    for idy in range(0, nfiles):
        if (dmag[idy] < 0.1) & (dmag[idy] > 0):
            ss = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + '/' + files[idy], header=0, sep=' ')
            trd[idz, :] = ss.cln.to_numpy()
            idz += 1
    tr = np.median(trd, axis=0)

    if idx == 105:
        plt.scatter(lc.jd - 2459587, lc.cln-0.85, c='k', marker='.')
        plt.ylabel('Instrumental Magnitude [+offset]', size=15)
        plt.xlabel('JD - 2459587 [days]', size=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylim([16.9, 16.5])
        plt.show()

    mmag, dmag, smag = sigma_clipped_stats(lc.mag - (tr-np.median(tr)), sigma=2.5)
    mags[idx] = dmag
    errs[idx] = smag

plt.scatter(mags-0.85, errs, c='k')
#plt.yscale('log')
plt.ylabel('rms')
plt.xlabel('Instrumental Magnitude [+offset]')
plt.show()
