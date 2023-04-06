UPDATE torosdb.obs_queue
SET toros_field_id = '%(toros_field_id)s', ra = %(ra)s, dec = %(dec)s, program = '%(program)s',
exposure_time = %(exposure_time)s, cadence = %(cadence)s, ephemeris = %(ephemeris)s, period = %(period)s
WHERE queue = 1
