UPDATE torosdb.toros_fields
SET observations = observations + 1
WHERE toros_field_id = '%(toros_field_id)s'
