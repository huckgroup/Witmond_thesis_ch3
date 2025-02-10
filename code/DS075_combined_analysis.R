library(tidyverse)

plateA <- read_delim("plateA_counts_and_percentages.csv", delim = ";", locale = locale(decimal_mark = ","))
plateB <- read_delim("plateB_counts_and_percentages.csv", delim = ";", locale = locale(decimal_mark = ","))
plateB <- subset(plateB, Plate == "Plate_B")
meta <- read_delim("meta_antibody_barcodes_samples_new.csv", delim = ",")

bc_data <- full_join(plateA, plateB)
bc_data <- left_join(bc_data, subset(meta, select = c(expected_Ab_BC_correct, Plate, Row, Column, well_nr)))

write_csv(bc_data, "combined_counts_and_percentages_corrected.csv")

# Select smaller dataset and order
dataset <- bc_data %>%
    select(Plate, Row, Column, well_nr, 
           expected_Ab_BC_correct, FeatureBarcode_nr, Barcodename, 
           nCount, counts, Percentage) %>%
    rename(plate = Plate, 
           row = Row, 
           column = Column, 
           well = well_nr, 
           expected_BC = expected_Ab_BC_correct, 
           measured_BC = FeatureBarcode_nr, 
           measured_BC_long = Barcodename, 
           count_total = nCount, 
           count_BC = counts, 
           count_percent = Percentage)  %>%
    filter(count_total > 100) %>%
    arrange(plate, well)

# Filter contamination on percentage contaminated: 1%
cat("\n\n Contamination threshold: 1% \n")
dataset_1_count <- dataset %>%
    filter(count_percent >= 1) %>%
    group_by(plate, well, expected_BC) %>%
    count() %>%
    rename(numberof_BC = n) %>%
    filter(numberof_BC > 1) %>%
    arrange(expected_BC)

dataset_1_list <- dataset %>%
    filter(count_percent >= 1) %>%
    group_by(plate, well, expected_BC) %>%
    filter(length(measured_BC) > 1) %>%
    summarise(allpresent_BC = paste(sort(unique(measured_BC)), collapse = ", ")) %>%
    arrange(expected_BC)

dataset_1 <- full_join(dataset_1_count, dataset_1_list)
dataset_1 %>% head(5)
cat("Number of contaminated BC:", nrow(dataset_1), "\n")
dataset_1_Ab <- dataset_1 %>%
    filter(expected_BC <= 64)
cat("Number of contaminated antibodies:", nrow(dataset_1_Ab), "\n")

# Filter contamination on percentage contaminated: 2.5%
cat("\n\n Contamination threshold: 2.5% \n")
dataset_2.5_count <- dataset %>%
    filter(count_percent >= 2.5) %>%
    group_by(plate, well, expected_BC) %>%
    count() %>%
    rename(numberof_BC = n) %>%
    filter(numberof_BC > 1) %>%
    arrange(expected_BC)

dataset_2.5_list <- dataset %>%
    filter(count_percent >= 2.5) %>%
    group_by(plate, well, expected_BC) %>%
    filter(length(measured_BC) > 1) %>%
    summarise(allpresent_BC = paste(sort(unique(measured_BC)), collapse = ", ")) %>%
    arrange(expected_BC)

dataset_2.5 <- full_join(dataset_2.5_count, dataset_2.5_list) %>% filter(!is.na(expected_BC))
dataset_2.5 %>% head(5)
cat("Number of contaminated BC:", nrow(dataset_2.5), "\n")
dataset_2.5_Ab <- dataset_2.5 %>%
    filter(expected_BC <= 64)
cat("Number of contaminated antibodies:", nrow(dataset_2.5_Ab), "\n")
dataset_2.5_free <- dataset_2.5 %>%
    filter(expected_BC > 64 & expected_BC <= 100)
cat("Number of free clean BC:", (100 - 64 - nrow(dataset_2.5_free)), "\n")

write_csv(dataset_2.5, "contaminated_BC_2.5_percent.csv")

# Filter contamination on percentage contaminated: 3%
cat("\n\n Contamination threshold: 3% \n")
dataset_3_count <- dataset %>%
    filter(count_percent >= 3) %>%
    group_by(plate, well, expected_BC) %>%
    count() %>%
    rename(numberof_BC = n) %>%
    filter(numberof_BC > 1) %>%
    arrange(expected_BC)

dataset_3_list <- dataset %>%
    filter(count_percent >= 3) %>%
    group_by(plate, well, expected_BC) %>%
    filter(length(measured_BC) > 1) %>%
    summarise(allpresent_BC = paste(sort(unique(measured_BC)), collapse = ", ")) %>%
    arrange(expected_BC)

dataset_3 <- full_join(dataset_3_count, dataset_3_list)
dataset_3 %>% head(5)
cat("Number of contaminated BC:", nrow(dataset_3), "\n")
dataset_3_Ab <- dataset_3 %>%
    filter(expected_BC <= 64)
cat("Number of contaminated antibodies:", nrow(dataset_3_Ab), "\n")

# Filter contamination on percentage contaminated: 5%
cat("\n\n Contamination threshold: 5% \n")
dataset_5_count <- dataset %>%
    filter(count_percent >= 5) %>%
    group_by(plate, well, expected_BC) %>%
    count() %>%
    rename(numberof_BC = n) %>%
    filter(numberof_BC > 1) %>%
    arrange(expected_BC)

dataset_5_list <- dataset %>%
    filter(count_percent >= 5) %>%
    group_by(plate, well, expected_BC) %>%
    filter(length(measured_BC) > 1) %>%
    summarise(allpresent_BC = paste(sort(unique(measured_BC)), collapse = ", ")) %>%
    arrange(expected_BC)

dataset_5 <- full_join(dataset_5_count, dataset_5_list)
dataset_5 %>% head(5)
cat("Number of contaminated BC:", nrow(dataset_5), "\n")
dataset_5_Ab <- dataset_5 %>%
    filter(expected_BC <= 64)
cat("Number of contaminated antibodies:", nrow(dataset_5_Ab), "\n")