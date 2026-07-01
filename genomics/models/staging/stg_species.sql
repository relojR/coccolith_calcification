select
    species_label,
    species_name,
    cast(is_calcifying as boolean) as is_calcifying,
    clade
from {{ ref('species_metadata') }}