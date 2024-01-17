## Key post-processing functions used for calculating scores and ALCS
## Yuyao Song
## ysong@ebi.ac.uk
##

# all_metrics is the concatenated file across all batch_metrics/cross-species/**_scIB_metrics.csv
# I did some str replace to get the integration_method and homology_method from the file names
# min max scaling and average, without LISI scores as in final paper

## Batch correction metrics and rank

batch_metrics = all_metrics %>% select(PCR, iLISI, bASW, GC, kBET, integration_method, homology_method) %>%
mutate(type  = paste(integration_method, homology_method, sep = " ")) %>% 
mutate(bASW_scaled = (bASW - min(bASW, na.rm = TRUE)) / (max(bASW, na.rm = TRUE) - min(bASW, na.rm = TRUE))) %>% 
mutate(GC_scaled = (GC - min(GC)) / (max(GC) - min(GC)))  %>% 
mutate(kBET_scaled = (kBET - min(kBET)) / (max(kBET) - min(kBET))) %>% 
mutate(PCR_scaled = (max(PCR, na.rm = TRUE) - PCR) / (max(PCR, na.rm = TRUE) - min(PCR, na.rm = TRUE))) %>% 
mutate(avg_score = (bASW_scaled + GC_scaled + kBET_scaled + PCR_scaled ) / 4) %>% 
arrange(desc(avg_score)) %>% 
mutate(GC_scaled_rank = dense_rank(desc(GC_scaled)))%>% 
mutate(PCR_scaled_rank = dense_rank(desc(PCR_scaled)))%>% 
mutate(bASW_scaled_rank = dense_rank(desc(bASW_scaled))) %>% 
mutate(kBET_scaled_rank = dense_rank(desc(kBET_scaled))) %>% 
mutate(avg_score_scaled_rank = dense_rank(desc(avg_score)))

## Bio conservation metrics and rank

bio_metrics = all_metrics %>% select(cASW, cLISI, ARI, NMI, iso_F1, integration_method, homology_method) %>%
mutate(type  = paste(integration_method, homology_method, sep = " ")) %>% 
mutate(cASW_scaled = (cASW - min(cASW, na.rm = TRUE)) / (max(cASW, na.rm = TRUE) - min(cASW, na.rm = TRUE))) %>% 
mutate(ARI_scaled = (ARI - min(ARI)) / (max(ARI) - min(ARI)))  %>% 
mutate(NMI_scaled = (NMI - min(NMI)) / (max(NMI) - min(NMI))) %>% 
mutate(iso_F1_scaled = (iso_F1 - min(iso_F1)) / (max(iso_F1) - min(iso_F1))) %>% 
mutate(avg_score = (cASW_scaled + ARI_scaled + NMI_scaled + iso_F1_scaled) / 4) %>% 
arrange(desc(avg_score)) %>% 
mutate(NMI_scaled_rank = dense_rank(desc(NMI_scaled)))%>% 
mutate(ARI_scaled_rank = dense_rank(desc(ARI_scaled)))%>% 
mutate(cASW_scaled_rank = dense_rank(desc(cASW_scaled))) %>% 
mutate(iso_F1_scaled_rank = dense_rank(desc(iso_F1_scaled))) %>% 
mutate(avg_score_scaled_rank = dense_rank(desc(avg_score))) 

## Integrated score
integrated_metrics = merge(batch_metrics, bio_metrics, by = 'type', suffixes = c("_batch", "_bio")) %>% 
mutate(integrated_score = avg_score_batch*0.4 + avg_score_bio*0.6) 


### ALCS: compare per-species before integration and after integration


# sccaf is the concatenated file across all <method>/cross-species/SCCAF_projection/**SCCAF_accuracy_summary.csv
# per_species is the concatenated file across all species in per-species/**SCCAF_accuracy_summary_<species>.csv

sccaf_acc = sccaf %>% 
select(test_acc, type_label, from_species, integration_method, input_file) %>% unique() %>% 
filter(type_label == 'original') %>%
group_by(input_file) %>% 
mutate(mean_acc = mean(test_acc)) %>%  ## mean accuracy of two species
select(-c(from_species, test_acc)) %>% unique() %>% 
ungroup() %>% 
mutate(type = paste(integration_method, homology_method))

per_species_acc = per_species %>% select(test_acc, from_species) %>% unique()

colnames(per_species_acc) = c("acc_prior", "from_species")

sccaf_results = sccaf_acc %>% 
merge(per_species_acc, by = 'from_species') %>% 
mutate(acc_loss = acc_prior - test_acc) %>% 
arrange(desc(acc_loss)) 
