library(tidyverse)

mgtstMetadata <-read_tsv("data/mgtstMetadata.tsv") %>%
    mutate(r = str_extract(pos,"[:ALPHA:]"),  c = str_sub(pos, 2)) %>% 
    rowwise() %>% 
    mutate(biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 1, "E01JH0011", biosample_id),
           biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 2, "E01JH0016", biosample_id),
           biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 3, "E01JH0017", biosample_id),
           biosample_id = if_else(pcr_16S_plate == 1 && seq_lab == "nist" && c == 4, "E01JH0004", biosample_id)) %>% 
    select(-r, -c)

check_samplesheet <- function(sample_df, fix = TRUE){
    ## Add columns with correct biosample (b_check) and titration (t_check)
    t20_pos <- paste0("A", c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11))
    t0_pos <- paste0(rep(c("B", "C", "E", "F", "G"), each = 2), 
                     rep(c(6, 12), 5))
    ## Setting titration factor as character to prevent comparison errors
    sample_df <- sample_df %>% 
        mutate(t_fctr = as.character(t_fctr),
               biosample_id = as.character(biosample_id))
    
    sam_check_df <- sample_df %>%
        mutate(t_fctr = as.character(t_fctr),
               t_check = case_when(pos %in% t20_pos ~ "20",
                                   pos %in% t0_pos ~ "0",
                                   TRUE ~ t_fctr
               ),
               b_check = case_when(
                   pos %in% c("B6", "B12") ~ "E01JH0004",
                   pos %in% c("C6", "C12") ~ "E01JH0011",
                   pos %in% c("E6", "E12") ~ "E01JH0016",
                   pos %in% c("F6", "F12") ~ "E01JH0017",
                   pos %in% c("G6", "G12") ~ "E01JH0038",
                   TRUE ~ biosample_id
               )
        )
    if (all(sam_check_df$t_fctr == sam_check_df$t_check,na.rm = TRUE) |
        all(sam_check_df$biosample_id == sam_check_df$b_check, na.rm = TRUE)) {
        if (fix) {
            msg <- str_c(
                "Titration or biosample metadata is not corrected. ",
                "Returning data frame with corrected metadata"
            )
            message(msg)
            
            corrected_sample_df <- sam_check_df %>%
                mutate(t_fctr = t_check,
                       biosample_id = b_check) %>%
                select(-t_check, -b_check)
            
            return(corrected_sample_df)
        } else {
            message(str_c(
                "Titration or biosample metadata is not correct ",
                "and has not been corrected. Returning input data frame."
            ))
            return(sample_df)
        }
    }
    
    message(str_c("Titration and biosample metadata is correct", 
                  ", returning input data frame"))
    return(sample_df)
}

## Correcting unmixed Post biosample ids
mgtstMetadata <- check_samplesheet(mgtstMetadata, fix = TRUE)

write_tsv(mgtstMetadata, "data/mgtstMetadata.tsv")
