{{ config(materialized='table') }}

with orthogroups as (
    select distinct orthogroup_id
    from {{ ref('stg_orthogroup_membership') }}
),

species as (
    select species_label, is_calcifying
    from {{ ref('stg_species') }}
),

grid as (
    select o.orthogroup_id, s.species_label, s.is_calcifying
    from orthogroups o
    cross join species s
),

counts as (
    select orthogroup_id, species_label, count(*) as gene_count
    from {{ ref('stg_orthogroup_membership') }}
    group by 1, 2
)

select
    g.orthogroup_id,
    g.species_label,
    g.is_calcifying,
    coalesce(c.gene_count, 0)      as gene_count,
    coalesce(c.gene_count, 0) > 0  as is_present
from grid g
left join counts c
    on  c.orthogroup_id = g.orthogroup_id
    and c.species_label = g.species_label