SELECT json_group_array(json_object('value', unique_build_id, 'dbkey', dbkey, 'name', label, 'path', file_path))
FROM (
 SELECT (r.rep_accnum || '.' || r.rep_version)                                           AS unique_build_id,
        (r.rep_accnum || '_' || r.rep_version)                                           AS dbkey,
        (REPLACE(r.definition, ', complete genome.', '') || ' [' || r.rep_accnum || '.' || r.rep_version ||
         ']')                                                                            AS label,
        (:basepath || '/' || genomeproject.gpv_directory || '/' || genomeproject.filename || '_genomic.fna') AS file_path
 FROM genomeproject
          JOIN replicon r ON genomeproject.gpv_id = r.gpv_id AND r.rep_type = 'chromosome'
 WHERE genomeproject.version_id = (SELECT version_id FROM version WHERE is_current == 1)
   AND genomeproject.file_types IS NOT NULL
   AND genomeproject.file_types LIKE '%.fna%'
   AND r.rep_type = 'chromosome'
)