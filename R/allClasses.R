
setClass('EnrichmentParameters',
         slots = list(
           organism = 'character',
           features = 'list',
           functional = 'list',
           structural = 'list'
         ),
         prototype = list(
           functional = list(methods = c('hypergeom',
                                         'diffusion',
                                         'pagerank'))
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

# setClass('Enrichment',
#          slots = list(
#            workflow = 'Workflow',
#            functional = 'FunctionalEnrichment',
#            structural = 'StructuralEnrichment'
#          )
# )