
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
