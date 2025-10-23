# Load required libraries
library(tidyverse); library(anytime); library(incidence2); library(truncdist);

## Old Notes from Initial Processing 
### ID #54 is missing - why is this?
### Missing "LOCATION" entries all have "ZONE_PEAL" filled in - can I just fill in "LOCATION" to be "NORTE" given they're part of PEAL?
### Latitude and Longitude have had an autofill issue - also there are 5 missing - any idea why?
### Some entries in "Date_notification" have "No_collection" written in them (3 of them) and one has "Nao notificado" there.
###    presumably they still received a notification even if no collection was made? And if no notification, how was it discovered?
###    Was it just randomly and so didn't go through the standard notification -> collection route.
###    Can I treat the 4 above as the same i.e. that they were identified through non-notification routes and hence no collection was made?
### What does "Não tem ou Positivo" mean for IHQ and PCR results? Same with Positivo ou Amostra inadequada - sem fígado.
###   ==> I ask this because monkeys with ID 66 and 67 have "Não tem ou Positivo" for IHQ and PCR but have been assigned positive - what should I do here?
### Some with the ambiguous dates = could be either of them (thinking the ones that are like "24th October or 24th December")

## Data dictionary
## ID - individual sample ID
## Species - monkey species
## LOCATION - municipality of Sao Paulo where found
## Zone Peal - if found in Horto/Peak, which zone 
## Latitude and Longitude - self-explanatory, but autofill has ruined the data - ask Nuno
## Date_death - their estimated date of death based on the decomposition
## Date_notification - date which the notification of the dead monkey was received
## Date_of_identification_or_collection - if notified, when the team reached there and collected the sample. If not notified, when the monkey was identified.
## Collection - whether a collection was undertaken or not (though there are some that are "Nao" for collection but then PCR positive - how does that work?)
##              also none on the non-PEAL have this filled in - why is that?
## Stage_decomposition - how decomposed the body was
## PCR_YFV_Result - PCR result for the monkey (Não tem means "it does not have")
## IHQ_YFV_Result - immunohistochemistry result for the monkey
## Final_YFV_Result - whether the monkey was assigned positive/negative based on IHQ/PCR results 
##   (check this column and others with Nuno as there's heterogeneity which doesn't enable me to figure out whether absence of result mean not tested
##    or tested and inconclusive result)
## STATUS - appears to describe whether or not the monkey's YFV was sequenced (I count only 20 in the spreadsheet for Horto, is that wrong?)

# Loading data
df <- read.csv(file = "data/Febre Amarela Horto.csv", fileEncoding="latin1")[1:159, ] %>%
  dplyr::select(ID, STATUS, Species, LOCATION, ZONE_PEAL, Latitude, Longitude, 
         Date_notification, Date_of_identification_or_collection, Collection, Stage_decomposition, Date_death,
         PCR_YFV_Result, IHQ_YFV_Result, FINAL_YFV_Result) %>%
  rename(sao_paulo_municipality = LOCATION, estimated_date_death = Date_death) %>%
  rename_with(tolower) %>%
  mutate(status = ifelse(status == "", "Not_Sequenced", status)) %>%
  mutate(species = case_when(species == "Alouatta guariba clamitans" ~ "Alouatta guariba clamitans",
                             species == "Callithrix aurita" | species == "Callithrix sp" | species == "Callithrix sp." ~ "Callithrix sp.",
                             species == "Sapajus sp" ~ "Sapajus sp.")) %>%
  mutate(sao_paulo_municipality = ifelse(sao_paulo_municipality == "", "Norte", sao_paulo_municipality)) %>%
  mutate(identified_through_notification = ifelse(date_notification == "NÃO NOTIFICADO" | date_notification == "There_was_no_collection", 
                                                  "not_idenfitied_through_notification", "identified_through_notification")) %>%
  mutate(date_notification = ifelse(date_notification == "NÃO NOTIFICADO" | date_notification == "There_was_no_collection", NA, date_notification)) %>%
  mutate(stage_decomposition = case_when(stage_decomposition == "" ~ "not_determined",
                                         stage_decomposition %in% c("Avançado", "Avançado (+5 a 10 dias)") ~ "3_advanced",
                                         stage_decomposition %in% c("Médio", "Médio-Avançado", "Médio (3 a 5 dias)") ~ "2_medium",
                                         stage_decomposition %in% c("íntegro", "Íntegro", "Íntegro (vivo _ 0 a 2 dias)", 
                                                                    "Íntegro (vivo _0 a 2 dias)", "Íntegro (vivo + 1 a 2 dias)",
                                                                    "Íntegro (vivo_0 a 2 dias)", "Íntegro (Vivo_0 a 2 dias)") ~ "1_intact",
                                         stage_decomposition == "Vivo" ~ "0_alive",
                                         TRUE ~ stage_decomposition)) %>%
  mutate(estimated_date_death = ifelse(estimated_date_death %in% c("desconhecida", "desconhecido"), "unknown", estimated_date_death)) %>%
  mutate(pcr_yfv_result = case_when(pcr_yfv_result %in% c("", "Não tem", "Não Tem", "Não tem ou Positivo", "Sample_not_sent") ~ "not_available",
                                    pcr_yfv_result %in% c("Negativo", "Negativo (fígado)", "Negativo (encéfalo), Negativo (fígado)", "Negativo (fígado), negativo (encéfalo)") ~ "negative",
                                    pcr_yfv_result %in% c("Negativo (fígado) / Positivo (encéfalo)", "Positivo") ~ "positive",
                                    TRUE ~ pcr_yfv_result)) %>%
  mutate(ihq_yfv_result = case_when(ihq_yfv_result %in% c("Amostra inadequada", "Amostra inadequada- sem amostra figado",
                                                          "Amostra inadequada, autólise", "Inadequado", "Inadequado (autólise)",
                                                          "Inadequado/SemFígado", "HP não compativel com FA", "Não realizada (HP ñ compatível com FA)",
                                                          "Não recebido", "Não tem", "Não Tem", "Não tem ou Positivo", "inconclusiva",
                                                          "Positivo ou Amostra inadequada - sem fígado", "Sample_not_sent") ~ "not_tested_or_inconclusive",
                                    ihq_yfv_result == "Negativo" ~ "negative", 
                                    ihq_yfv_result == "Positivo" ~ "positive")) %>%
  mutate(final_yfv_result = case_when(final_yfv_result == "NÃO DETERMINADO" ~ "undetermined",
                                      final_yfv_result == "NEGATIVO" ~ "negative", 
                                      final_yfv_result == "POSITIVO" ~ "positive")) %>%
  mutate(date_collection = substr(date_of_identification_or_collection, start = 1, stop = 10)) %>%
  mutate(date_collection = as.Date(date_collection, format = "%d/%m/%Y"))
saveRDS(object = df, file = "analyses/1_IBM_estimationR0/data/processed_HortoData.rds")
