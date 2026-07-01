-- Every orthogroup must have exactly one row per species.
-- Any orthogroup whose row count != number of species is a broken matrix.
-- A passing test returns zero rows.



select
    orthogroup_id,
    count(*) as row_count
from {{ ref('int_orthogroup_presence') }}
group by orthogroup_id
having count(*) <> (select count(*) from {{ ref('stg_species') }})