DROP TABLE IF EXISTS builds_temp;
CREATE TEMP TABLE builds_temp(dbkey unique, label);
.mode csv
.separator '	'
.import $__app__.config.builds_file_path builds_temp
INSERT OR REPLACE INTO builds_temp (dbkey, label)
SELECT (r.rep_accnum || '_' || r.rep_version)                                           AS dbkey,
    (REPLACE(r.definition, ', complete genome.', '') || ' [' || r.rep_accnum || '.' || r.rep_version ||
     ']')                                                                            AS label
FROM genomeproject
      JOIN replicon r ON genomeproject.gpv_id = r.gpv_id AND r.rep_type = 'chromosome'
WHERE genomeproject.version_id = (SELECT MAX(version_id) FROM version)
    AND genomeproject.file_types IS NOT NULL
    AND genomeproject.file_types LIKE '%.fna%'
    AND r.rep_type = 'chromosome';
.mode list
.separator '	'
.once $__app__.config.builds_file_path
SELECT dbkey, label FROM builds_temp ORDER BY dbkey;