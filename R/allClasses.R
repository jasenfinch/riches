
setClass('EnrichmentParameters',
         slots = list(
           organism = 'character',
           functional = 'list',
           structural = 'list'
         )
)

setClass('FunctionalEnrichment',
         slots = list(
           network = 'FELLA.DATA',
           hits = 'tbl_df',
           explanatory = 'tbl_df',
           results = 'list'
         )
)

setClass('StructuralEnrichment')

setClass('Enrichment',
         slots = list(
           workflow = 'Workflow',
           functional = 'FunctionalEnrichment',
           structural = 'StructuralEnrichment'
         )
)