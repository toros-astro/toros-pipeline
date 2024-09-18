from libraries.priority import Priority
from config import Configuration
import numpy as np
import pandas as pd

# calculate the field separations
sep = (Configuration.FOV * np.sqrt(2.2 ** 2 + 3.75 ** 2)) / 5.0
field_list = Priority.toros_field_generator(sep)
# Priority.toros_survey_simulator(field_list)

# determine the toros fields for Alessi 13 and Blanco 1
al_13_ra = 52.0500000
al_13_dec = -35.9000000
al_13_ext = 0.7

bl_1_ra = 1.0291667
bl_1_dec = -29.8333333
bl_1_ext = 0.7

tuc_47_ra = 6.0223292
tuc_47_dec = -72.0814444
tuc_47_ext = 44. / 60.

tuc_47_list, tuc_47_ang = Priority.find_toros_field(field_list, tuc_47_ra, tuc_47_dec,
                                                    ang_extend=tuc_47_ext, program='commissioning')
bl_1_list, bl_1_ang = Priority.find_toros_field(field_list, bl_1_ra, bl_1_dec,
                                                ang_extend=bl_1_ext, program='commissioning')
al_13_list, al_13_ang = Priority.find_toros_field(field_list, al_13_ra, al_13_dec,
                                                  ang_extend=al_13_ext, program='commissioning')

comms_list = pd.concat([tuc_47_list, bl_1_list, al_13_list]).reset_index(drop=True)
Priority.toros_night_field_selector(field_list, comms_list)
print('hold')