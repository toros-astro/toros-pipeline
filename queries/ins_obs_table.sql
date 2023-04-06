INSERT INTO torosdb.observations (toros_field_id, ra, dec, program, exposure_time, cadence, ephemeris, period, obs_date)
VALUES ('%(toros_field_id)s', %(ra)s, %(dec)s, '%(program)s', %(exposure_time)s, %(cadence)s,
%(ephemeris)s, %(period)s, '%(obs_date)s');
