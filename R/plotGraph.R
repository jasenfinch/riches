
setMethod('plotGraph',signature = 'FunctionalEnrichment',
          function(x,pairwise,type = 'diffusion'){
            en <- x@results[[pairwise]]
            
            pg <- generateResultsGraph(method = type,object = en,data = FELLA,nlimit = 1000,LabelLengthAtPlot = 100) %>%
              as_tbl_graph() %>%
              activate(nodes)
            
            types <- pg %>%
              nodes() %>%
              .$com
            
            types[types == 1] <- 'Pathway'
            types[types == 2] <- 'Module'
            types[types == 3] <- 'Enzyme'
            types[types == 4] <- 'Reaction'
            types[types == 5] <- 'Compound'
            
            pg <- pg %>%
              mutate(Type = factor(types,levels = c('Compound','Reaction','Enzyme','Module','Pathway')),
                     label = str_split_fixed(label,coll(' - '),2)[,1])
            
            pageLayout <- pg %>%
              create_layout('nicely') 
            
            pagePathways <- pageLayout %>%
              filter(Type == 'Pathway')
            
            pageCompounds <- pageLayout %>%
              filter(Type == 'Compound')
            
            nudgey <- ({pageLayout$y %>% max()} - {pageLayout$y %>% min()}) / 25
            
            colours <- tibble(Type = c('Compound','Reaction','Enzyme','Module','Pathway'),
                              Colour = ptol_pal()(5))
            
            ggraph(pageLayout) +
              geom_edge_link(alpha = 0.4) +
              geom_node_point(aes(colour = Type),size = 3,alpha = 0.5) +
              geom_node_text(data = pagePathways,aes(label = label),nudge_y = nudgey,size = 3,colour = 'black') +
              geom_node_text(data = pageCompounds,aes(label = label),nudge_y = nudgey,size = 3,colour = 'black') +
              scale_colour_manual(values = colours %>% filter(Type %in% unique(pg %>% nodes() %>% .$Type)) %>% .$Colour) +
              theme_graph(base_size = 11,base_family = '') +
              coord_fixed() +
              labs(title = pairwise) 
          }
)