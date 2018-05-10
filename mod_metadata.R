library(tidyverse)

mgtstMetadata <-read_tsv("data/mgtstMetadata.tsv") %>%
    mutate(r = str_extract(pos,"[:ALPHA:]"),  c = str_sub(pos, 2)) %>% 
    rowwise() %>% 
    mutate(biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 1, "E01JH0011", biosample_id),
           biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 2, "E01JH0016", biosample_id),
           biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 3, "E01JH0017", biosample_id),
           biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 4, "E01JH0004", biosample_id)) %>% 
    select(-r, -c)

write_tsv(mgtstMetadata, "data/mgtstMetadata.tsv")
