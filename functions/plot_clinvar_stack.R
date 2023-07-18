plot_clinvar_stack <- function(gr_clinvar_full, clinvar_default_results, execute_analysis_and_plot = F) {
    if (execute_analysis_and_plot == T) {
      library(ggplot2)
      library(Biostrings)
    
      clinvar_default_results = unlist(clinvar_default_results)
      clinvar_default_results$res_aenmd = cbind(clinvar_default_results$res_aenmd, 
                                          nmd_escaping = apply(clinvar_default_results$res_aenmd[,c(2,3,4,5,6)], 1, any))
      
      # only keep is ptc
      clinvar_annotated_ptc = clinvar_default_results[clinvar_default_results$res_aenmd$is_ptc]
      
      
      ###BENIGN:
      #How many total benign variants in clinvar?
      clinvar_benign <- gr_clinvar_full[stringr::str_detect(gr_clinvar_full$clinsig, "(?i)benign")]
      #How many benign PTC?
      clinvar_annotated_ptc_benign  <- clinvar_annotated_ptc[stringr::str_detect(clinvar_annotated_ptc$clinsig, "(?i)benign")]
      
      # view each variant as a different key not a position; so one genomic position can have more than 1 key due to alternative allele
      # go through each key, count the number of escaping, not escaping and conflict
      outcome_benign = lapply(unique(clinvar_annotated_ptc_benign$key), function(x){
        a = clinvar_annotated_ptc_benign[clinvar_annotated_ptc_benign$key == x]
        if(length(unique(a$res_aenmd$nmd_escaping))>1){
          df = data.frame(key = x, nmd_escaping =  "conflicting")
        }else{
          df = data.frame(key = x, nmd_escaping =  unique(a$res_aenmd$nmd_escaping))
        }
        return(df)
      })
      
      outcome_benign_total = do.call("rbind", outcome_benign)
      b1 = length(unique(clinvar_benign$id))
      b2 = length(unique(clinvar_annotated_ptc_benign))
      b3 = table(outcome_benign_total$nmd_escaping)
      
      
      ###Pathogenic:
      #How many total pathogenic variants in clinvar?
      clinvar_patho = clinvar_gr[stringr::str_detect(clinvar_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(clinvar_gr$clinsig, "(?i)conflicting")]
      #How many pathogenic PTC?
      clinvar_annotated_ptc_patho = clinvar_annotated_ptc[stringr::str_detect(clinvar_annotated_ptc$clinsig, "(?i)pathogenic") & !stringr::str_detect(clinvar_annotated_ptc$clinsig, "(?i)conflicting")]
      
      # view each variant as a different key not a position; so one genomic position can have more than 1 key due to alternative allele
      # go through each key, count the number of escaping, not escaping and conflict
      outcome_path = lapply(unique(clinvar_annotated_ptc_patho$key), function(x){
        a = clinvar_annotated_ptc_patho[clinvar_annotated_ptc_patho$key == x]
        if(length(unique(a$res_aenmd$nmd_escaping))>1){
          df = data.frame(key = x, nmd_escaping =  "conflicting")
        }else{
          df = data.frame(key = x, nmd_escaping =  unique(a$res_aenmd$nmd_escaping))
        }
        return(df)
      })
      
      outcome_path_total = do.call("rbind", outcome_path)
      p1=length(unique(clinvar_patho$id))
      p2=length(unique(clinvar_annotated_ptc_patho))
      p3=table(outcome_path_total$nmd_escaping)
      
      
      ###Conflicting / Uncertain:
      #How many total Conflicting / Uncertain variants in clinvar?
      clinvar_uncertain = clinvar_gr[stringr::str_detect(clinvar_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(clinvar_gr$clinsig, "(?i)uncertain")]
      #How many Conflicting / Uncertain PTC?
      clinvar_annotated_ptc_uncertain = clinvar_annotated_ptc[stringr::str_detect(clinvar_annotated_ptc$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(clinvar_annotated_ptc$clinsig, "(?i)uncertain")]
      
      # view each variant as a different key not a position; so one genomic position can have more than 1 key due to alternative allele
      # go through each key, count the number of escaping, not escaping and conflict
      outcome_uncrt = lapply(unique(clinvar_annotated_ptc_uncertain$key), function(x){
        a = clinvar_annotated_ptc_uncertain[clinvar_annotated_ptc_uncertain$key == x]
        if(length(unique(a$res_aenmd$nmd_escaping))>1){
          df = data.frame(key = x, nmd_escaping =  "conflicting")
        }else{
          df = data.frame(key = x, nmd_escaping =  unique(a$res_aenmd$nmd_escaping))
        }
        return(df)
      })
      
      outcome_uncrt_total = do.call("rbind", outcome_uncrt)
      u1=length(unique(clinvar_uncertain$id))
      u2=length(unique(clinvar_annotated_ptc_uncertain))
      u3=table(outcome_uncrt_total$nmd_escaping)
      
      
      df = data.frame(class = rep(c("patho/likely patho", "benign/likely benign", "uncertain/conflicting"),each = 5),
                      nmd = rep(c("non-PTC", "PTC", "nmd-escaping", "nmd-triggering", "transcript dependent"), 3),
                      number = c(p1-p2, p2, p3["TRUE"], p3["FALSE"], p3["conflicting"],
                                 b1-b2, b2, b3["TRUE"], b3["FALSE"], b3["conflicting"],
                                 u1-u2, u2, u3["TRUE"], u3["FALSE"], u3["conflicting"]))
      
      
      
      df = df[!df$nmd %in% "PTC",]
      df$class = factor(df$class, levels = c("patho/likely patho",  "uncertain/conflicting", "benign/likely benign"))
      df$nmd   = factor(df$nmd, levels = c("transcript dependent",  "nmd-escaping", "nmd-triggering", "non-PTC"))
      first_plot = ggplot(data=df, aes(x=class, y=number, fill=nmd)) +
        geom_bar(stat="identity") + 
        # scale_fill_brewer(palette="YlGnBu", name = "Number of  \ndiseases associated") + 
        scale_fill_manual(values=c("#d1b2ff",
                                   "#33FFFF",
                                   "red",
                                   "darkblue")) + 
        theme_classic() +
        labs(title="", 
             x="", y = "number of variants", 
             fill = "aenmd annotation") + 
        ggtitle("Clinvar variants") + 
        theme(plot.title = element_text(hjust = 0.5, size = 20))
      
      first_plot
      df2 = df[!df$nmd %in% "non-PTC",]
      
      cowplot::save_plot(filename="./plot/clinvar_nmd_outcome.pdf",first_plot, base_height = 5.2, base_width = 6.3)
      
      
      df2 = df2 %>% dplyr::group_by(class) %>% dplyr::mutate(percentage = number/sum(number))
      df2 <- df2 %>%
        dplyr::group_by(class) %>%
        dplyr::arrange(class, desc(nmd)) %>%
        dplyr::mutate(lab_ypos = cumsum(percentage) - 0.5 * percentage) 
      
      second_plot = ggplot(data=df2, aes(x=class, y=percentage, fill=nmd)) +
        geom_bar(stat="identity") + 
        # scale_fill_brewer(palette="YlGnBu", name = "Number of  \ndiseases associated") + 
        scale_fill_manual(values=c("#d1b2ff",
                                   "#33FFFF",
                                   "red")) + 
        geom_text(aes(x = class, y = lab_ypos, label =number, group = nmd), color = "black") +
        theme_classic() +
        labs(title="", 
             x="", y = "percentage", fill = "aenmd annotation",) +
        ggtitle("Clinvar PTC variants") + 
        theme(plot.title = element_text(hjust = 0.5, size = 20))
      second_plot
      cowplot::save_plot(filename="./plot/clinvar_nmd_outcome_percentage.pdf",second_plot, base_height = 5.2, base_width = 6.3)
    
      #put the plots in to an object to return, naming them in the process
      saved_plots <-list(first_plot,second_plot)
      names(saved_plots) <- c("clinvar_nmd_outcome", "clinvar_nmd_outcome_percentage")
      
      return(saved_plots)
  } else {message("execute_analysis_and_plot variable is set to FALSE, change to TRUE to execute.")}
}