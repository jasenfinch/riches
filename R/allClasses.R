
setClass('EnrichmentParameters',
         slots = list(
           functional = 'list',
           structural = 'list'
         )
)

setClass('FunctionalEnrichment',
         slots = list(
           network = 'FELLA.DATA',
           hits = 'tbl_df',
           explanatory = 'character',
           results = 'FELLA.USER'
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