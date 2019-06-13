
#' @export

setMethod('functionalEnrichment',signature = 'Workflow',
          function(x,parameters){
            FELLA <- organismNetwork(organism)
            
            oc <- organismCompounds(FELLA)
            
            suppressMessages(oc <- metabolites %>%
                               filter(ACCESSION_ID %in% oc$name) %>%
                               {metaboliteDB(.,descriptors = descriptors(.))})
            
            mfs<- x %>%
              resultsAnnotation() %>%
              assignments() %>%
              select(Name,MF,Adduct) %>%
              distinct()
            
            MFhits <- mfs %>%
              split(1:nrow(.)) %>%
              map(~{
                m <- .
                oc %>%
                  filterMF(m$MF) %>%
                  filterIP(Adducts$Rule[Adducts$Name == m$Adduct]) %>%
                  getAccessions() %>%
                  mutate(Name = m$Name,MF = m$MF,Adduct = m$Adduct)
              }) %>%
              bind_rows()
            
            bc <- MFhits$ACCESSION_ID %>%
              unique()
            
            explanFeat <- x %>%
              resultsAnalysis() %>%
              featureSelectionResults() %>%
              filter(Method == parameters@features$method,
                     Pvalue < parameters@features$threshold)
            
            pairwises <- explanFeat$Pairwise %>%
              unique()
            
            enrichRes <- pairwises %>%
              map(~{
                p <- .
                feat <- explanFeat %>%
                  filter(Pairwise == p)
                ec <- MFhits %>%
                  filter(Name %in% feat$Feature) %>%
                  .$ACCESSION_ID %>%
                  unique()
                defineCompounds(
                  compounds = ec,
                  compoundsBackground = bc,
                  data = FELLA) %>%
                  runHypergeom(data = FELLA) %>%
                  runDiffusion(data = FELLA) %>%
                  runPagerank(data = FELLA)
              }) %>%
              set_names(pairwises)
            
            new('FunctionalEnrichment',
                network = FELLA,
                hits = MFhits,
                explanatory = explanFeat,
                results = enrichRes
                )
          }
)