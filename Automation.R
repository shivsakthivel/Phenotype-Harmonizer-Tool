# Code adapted from UW-GAC TOPMed Phenotype Harmonization Example Process Files
# Citation: Stilp AM, Emery LS, et al. A System for Phenotype Harmonization in the NHLBI Trans-Omics for Precision Medicine (TOPMed) Program. Am J Epidemiol. 2021 Apr 16:kwab115. doi: 10.1093/aje/kwab115. Epub ahead of print. PMID: 33861317.

# Note: For variables not requiring any previously harmonized DCC variables, the script can be run directly. When harmonized variables are required, the script needs to be run in proper succession. For example, in the set of VTE variables, followup_start_age, prior_history and case_status need to be run in that particular order.

# Input for Terminal: > Rscript [Script Name] [phenotype_concept + concept_variant]
# Output: Harmonized dataset as a tab-separated values file (tsv) and a pdf of descriptive plots for appraisal

# Current Tool Functions (August 2021): It is able to successfully concatenate datasets for the 63 phenotype variables with documentation provided by the UW-GAC, within the cohort of TOPMed Studies available on the SRTI Servers. The JSON documentation for harmonized variables needs to be present in the working directory. The companion plots file provides simple visualizations of distributions for the transformed dataset by individual study, and an overall distribution for the harmonized set.

# Required Libraries
library(readr)
library(plyr)
library(dplyr)
library(stringr)
library(jsonlite)
library(xml2)
library(R.utils)
library(ggplot2)
library(Rmisc)

# Source File for dbgaptools and constants. Instead of importing and running dbGaP functions directly from the library, the required functions and constants for this project were manually collected in this file. Credits for the functions go the code available in the dbGaP GitHub repository. This simplifies the import process for the very limited subset of external functions required.
source("phen_harm.R")

# Command line arguments for required harmonized variable - Must be entered with correct concept variant as well. This information is available with the documentation and also, the paper describing the harmonization framework.
arguments = commandArgs(trailingOnly=TRUE)
if (length(arguments)==0) {
  stop("At least one argument must be supplied (input phenotype).n", call.=FALSE)
}
variable <- arguments[1]

# All components from the json file
json_file <- paste(variable, "json", sep = ".")
json <- fromJSON(json_file, simplifyDataFrame = FALSE)
age_variable <- json$has_age_variable

# Variable Type
var_type <- NULL
if(length(json$encoded_values) == 0){
    var_type <- "continuous"
} else {
    var_type <- "discrete"
}

# Names
full_name <- paste(json$phenotype_concept, json$concept_variant, sep = "_")
age_name <- paste("age", "at", full_name, sep = "_")

# Empty Df (for concatenation through loop) and Plots list for graphics
all_studies <- NULL
all_plots <- list()

# Innermost loop
for (i in 1:length(json$harmonization_units)){
    
    # Variables constructed from solely DCC harmonized variables don't need to go through the loop. This applies to BMI and LDL
    if(variable == 'bmi_baseline_1'){
        table1 <- read.table('height_baseline_1_harmonized.tsv', header = TRUE)
        table2 <- read.table('weight_baseline_1_harmonized.tsv', header = TRUE)
        phen_list <- list()
        phen_list$harmonized_data$height_baseline_1 <- table1
        phen_list$harmonized_data$weight_baseline_1 <- table2
        overall <- json$harmonization_units[[1]]
        harmonize <- eval(parse(text = overall$harmonization_function))
        out <- harmonize(phen_list)
        names(out)[names(out) == json$phenotype_concept] = full_name
        if(age_variable){
            names(out)[names(out) == "age"] = age_name
        }
        all_studies <- out
        break
    }
    
    if(variable == 'ldl_1'){
        table1 <- read.table('total_cholesterol_1_harmonized.tsv', header = TRUE)
        table2 <- read.table('triglycerides_1_harmonized.tsv', header = TRUE)
        table3 <- read.table('hdl_1_harmonized.tsv', header = TRUE)
        phen_list <- list()
        phen_list$harmonized_data$total_cholesterol_1 <- table1
        phen_list$harmonized_data$triglycerides_1 <- table2
        phen_list$harmonized_data$hdl_1 <- table3
        overall <- json$harmonization_units[[1]]
        harmonize <- eval(parse(text = overall$harmonization_function))
        out <- harmonize(phen_list)
        names(out)[names(out) == json$phenotype_concept] = full_name
        if(age_variable){
            names(out)[names(out) == "age"] = age_name
        }
        all_studies <- out
        break
    }
    
    
    # Information unit by unit from the json file
    unit <- json$harmonization_units[[i]]
    component_variables <- unit$component_study_variables
    matches <- stringr::str_match(component_variables, REGEX_DBGAP)
    components <- as.data.frame(matches, stringsAsFactors = FALSE) %>%
    transmute(
        accession_string = V1,
        phs = as.integer(V2),
        phs_version = as.integer(V3),
        pht = as.integer(V4),
        pht_version = as.integer(V5),
        phv = as.integer(V6),
        phv_version = as.integer(V7)
      )
    phs <- unique(components$phs)
    phs_version <- unique(components$phs_version)
    
    # Study Name from the json file and mapping may not be identical. Study number is the appropriate identifier.
    phs_mapping_file <- read_tsv("phs-mapping.tsv")
    study_number <- 0
    studies <- as.integer(phs_mapping_file[['phs']])
    for(j in 1:length(studies)){
        if(studies[[j]] == phs){
            study_number <- studies[[j]]
            break
        }
    }
    phs_mapping <- NULL
    
    if(study_number != 0){
       phs_mapping <- phs_mapping_file %>% filter(phs == study_number) 
    } else {
        break
    }
    
    print(study_number)
    
    # Built-in Caveats - Manually fixed Abnormal Harmonizations
    # Approach: Manually build dataset for that study and add it during its turn
    
    # For HVH on the smoker variables, the primary key was being dropped inexplicably during the automated harmonization, and hence this fix was required.
    if (phs == 1013 && variable == "current_smoker_baseline_1") {
        csb_HVH <- read.table(file = 'csb_HVH_harmonized.tsv', sep = '\t', header = TRUE) 
        csb_HVH <- csb_HVH %>% mutate(Individual_ID = as.character(Individual_ID))
        if(var_type == "continuous"){
            all_plots[[i]] <- ggplot(data = csb_HVH, aes_string(full_name)) + geom_histogram(aes(y = ..ncount..)) + geom_density(aes(y = ..scaled..)) + labs(title = sprintf("Distribution of %s in %s", full_name, unit$name))
        } else {
            all_plots[[i]] <- ggplot(data = csb_HVH) + geom_bar(aes_string(full_name)) + labs(title = sprintf("Distribution of %s in %s", full_name, unit$name))
        }
        all_studies <- bind_rows(all_studies, csb_HVH)
        next
    }
    
    if (phs == 1013 && variable == "ever_smoker_baseline_1") {
        esb_HVH <- read.table(file = 'esb_HVH_harmonized.tsv', sep = '\t', header = TRUE) 
        esb_HVH <- esb_HVH %>% mutate(Individual_ID = as.character(Individual_ID))
        if(var_type == "continuous"){
            all_plots[[i]] <- ggplot(data = esb_HVH, aes_string(full_name)) + geom_histogram(aes(y = ..ncount..)) + geom_density(aes(y = ..scaled..)) + labs(title = sprintf("Distribution of %s in %s", full_name, unit$name))
        } else {
            all_plots[[i]] <- ggplot(data = esb_HVH) + geom_bar(aes_string(full_name)) + labs(title = sprintf("Distribution of %s in %s", full_name, unit$name))
        }
        all_studies <- bind_rows(all_studies, esb_HVH)
        next
    }
    
        
    # Now, using the location of the study, traverse directories to get the files required for harmonization
    server <- as.character(phs_mapping[["location"]])
    common <- NULL
    if (server == "n"){
        next
    } else if (server == "stsi3") {
        common <- common_path_stsi3
    } else {
        common <- common_path_p1
    }
    
    # Depth: These will be the different consent groups, and all of these need to be identified and concatenated.
    args <- list.files(path = common, pattern = sprintf('phs%06d', phs))
    # Caveat: There are two identical MGH directories on the same level, so one needs to be skipped
    if(phs == 1001){
        args <- args[-2]
    }
    directories <- c()
    change <- NULL
    for (j in 1:length(args)){
        version <- list.files(path = paste(common, args[j], sep = "/"), pattern = sprintf('phs%06d.v%d', phs, phs_version))
        if (length(version) == 0){
            folders <- list.files(path = paste(common, args[j], sep = "/"))
            matches_2 <- stringr::str_match(folders, REGEX_PHS)
            organized <- as.data.frame(matches_2, stringsAsFactors = FALSE) %>% transmute(
                accession_string = V1,
                phs = as.integer(V2),
                phs_version = as.integer(V3))
            change <- (organized$phs_version)[1]
            version <- folders[[1]]
        } 
        # Only irregularity - There is an extra inner directory. All other directories for all the other studies is good to go!
        if (phs == 956){
            version <- paste(version, "PhenoGenotypeFiles", sep = "/")
        }
        root_study <- list.files(path = paste(common, args[j], version, sep = "/"), pattern = "^Root*")
        phenotype_files <- list.files(path = paste(common, args[j], version, root_study, sep = "/"), pattern = "^Phenotype*")
        full_path = paste(common, args[j], version, root_study, phenotype_files, sep = "/")
        directories <- append(directories, full_path)
    }
    
    # Indication of whether different version is attempted for use - The versions used by harmonization in the DCC are often older versions than those available on SRTI servers. So, the closest version is found and used instead with the expectation that the required variable has not been phased out. Currently, no such instances were identified so using the closest version has been completely viable.
    if (is.null(change)){
    } else {
    phs_version <- change
    }

    
    # This gives the innermost common directory paths for the consent groups to obtain the phenotype files later on.
    # Now here is the subject file. This is common across all consent groups for a single study, so it is just obtained from the first argument.
    if(length(directories) == 0){
        next
    }
    files_1 <- list.files(directories[1])
    subj_file <- files_1[grepl("Subject\\.MULTI", files_1)]
    
    # Only the identifier columns are required from the subject file. dbGaP id is the common column name across ALL studies
    # But, the secondary identifiers are different across different studies
    # To automate finding the right variable, the number of unique entries are compared to dbGaP id. When equal, that is the correct identifier.
    ext <- tools::file_ext(subj_file)
    subj <- NULL
    if (ext == "gz"){
        subj <- read_ds_file(gunzip(paste(directories[1], subj_file, sep = "/")))
    } else {
        subj <- read_ds_file(paste(directories[1], subj_file, sep = "/"))
    }
    columns <- colnames(subj)
    entries <- length(subj[['dbGaP_Subject_ID']])
    identifier <- NULL
    for (j in 1:length(columns)){
        if(columns[j] == 'dbGaP_Subject_ID'){
            next
        } else {
            current <- length(subj[[columns[j]]])
            if(current == entries){
                identifier <- columns[j]
                break
            }
        }
    }
    
    # Now, we'll only keep the necessary columns and mutate as necessary. Individual_ID regulates colnames
    subj <- subj %>% select(dbGaP_Subject_ID, as.character(identifier))
    if(identifier == 'Individual_ID'){
        subj <- subj %>% mutate(unique_subject_key = sprintf("%s_%s", unit$name, Individual_ID))
    } else {
        subj['Individual_ID'] = subj[as.character(identifier)]
        subj <- subj %>%
            mutate(unique_subject_key = sprintf("%s_%s", unit$name, Individual_ID)) %>%
            select(-c(2))
    }
    study_integer <- phs_mapping$unique_id
    subj <- subj %>%
        mutate(topmed_subject_id = 1:n() + study_integer)
    
    # The subject file is ready. Now, there is some other information required about phv and pht numbers possibly.
    pht <- components$pht
    pht_version <- components$pht_version
    regex_patterns <- c()
    for (j in 1:length(pht)){
        regex_patterns[j] <- sprintf('^phs%06d\\.v%d.pht%06d\\.', phs, phs_version, pht[j])
    }
    regex_patterns <- unique(regex_patterns)
    pht <- unique(pht)
    
    # File Containers! Each element of the outer list is a list of files for that pht number.
    file_container <- list()
    for (j in 1:length(regex_patterns)){
        current <- list()
        for (k in 1:length(directories)){
            current[[k]] <- list.files(directories[k], pattern = regex_patterns[j], full.names = TRUE) 
        }
        current <- unlist(current)
        file_container[[j]] <- current
    }
    
    # Source data containing lists of dataframes
    source_data <- list()
    for (j in 1:length(file_container)){
        current <- file_container[[j]]
        dd_file <- current[endsWith(current, 'data_dict.xml')][1]
        dd <- read_dd_xml(dd_file) %>% filter(phv %in% components$phv)
        phen_files <- current[endsWith(current, '.txt')]
        
        tmp_list <- NULL
        if (length(phen_files) == 0){
            phen_files <- current[endsWith(current, '.gz')]
            # Very weird fix - Suggestions Welcome - Issue: Unzipping without dynamic variable creation in list. This works because none of the studies in the SRTI servers have more than 4 consent groups. It will break or be improperly harmonized if that is exceeded.
            # This was required because unzipping and concatenating the dataset simultaneously is necessary but not possible by iteration.
            if(length(phen_files) == 1){
              tmp_list <- list(read_ds_file(gunzip(phen_files[1])))
            } else if (length(phen_files) == 2){
                tmp_list <- list(read_ds_file(gunzip(phen_files[1])), read_ds_file(gunzip(phen_files[2])))
            } else if (length(phen_files) == 3){
                tmp_list <- list(read_ds_file(gunzip(phen_files[1])), read_ds_file(gunzip(phen_files[2])), read_ds_file(gunzip(phen_files[3])))
            } else if (length(phen_files) == 4){
                tmp_list <- list(read_ds_file(gunzip(phen_files[1])), read_ds_file(gunzip(phen_files[2])), read_ds_file(gunzip(phen_files[3])), read_ds_file(gunzip(phen_files[4])))
            } 
        } else {
            tmp_list <- lapply(phen_files, read_ds_file)
        }
        phen <- bind_rows(tmp_list) %>% left_join(subj, by = 'dbGaP_Subject_ID')
        component_variable_names <- as.list(dd$name)
        phen <- phen %>% select(topmed_subject_id, !!! component_variable_names)
        source_data[[sprintf("pht%06d", pht)[j]]] <- phen
    }
    phen_list <- list(source_data = source_data)
    
    # Caveat Alert: FHS, WHI on VTE Variables
    # For these VTE variables, information from previously harmonized variables are required. So, those are read in manually and joined.
    if (phs == 7 && variable == "vte_prior_history_1"){
        followup <-  read.table('vte_followup_start_age_1_harmonized.tsv', header = TRUE)
        just_FHS <- followup[grep("FHS", followup$unique_subject_key), ]
        phen_list$harmonized_data$vte_followup_start_age_1 <- just_FHS
    }
    
    if (phs == 7 && variable == "vte_case_status_1"){
        followup <-  read.table('vte_prior_history_1_harmonized.tsv', header = TRUE)
        just_FHS <- followup[grep("FHS", followup$unique_subject_key), ]
        phen_list$harmonized_data$vte_prior_history_1 <- just_FHS
    }
    
    if (phs == 200 && variable == "vte_case_status_1"){
        followup <-  read.table('vte_prior_history_1_harmonized.tsv', header = TRUE)
        just_WHI <- followup[grep("WHI", followup$unique_subject_key), ]
        phen_list$harmonized_data$vte_prior_history_1 <- just_WHI
    }
        
    # Harmonization!
    harmonize <- eval(parse(text = unit$harmonization_function))
    out <- harmonize(phen_list)
    
    # Variable Name Changes
    names(out)[names(out) == json$phenotype_concept] = full_name
    if(age_variable){
        names(out)[names(out) == "age"] = age_name
    }
    
    # Caveat Alert: CRA unit conversion for Height (Inconsistent between cm and m in the Root Phenotype Files - This is an outright inaccuracy in those files)
    if (phs == 988 && variable == "height_baseline_1") {
        out$height_baseline_1 <- ifelse(out$height_baseline_1 < 3, 100 * out$height_baseline_1, out$height_baseline_1)
    }
    
    if(age_variable){
        harmonized_study <- out %>% 
            left_join(subj, by = "topmed_subject_id") %>%
            select(
            Individual_ID,
            unique_subject_key,
            topmed_subject_id,
            as.character(full_name),
            as.character(age_name)
            ) 
        #%>% mutate_if(str_detect(colnames(.), "age"), as.double)
        } else {
           harmonized_study <- out %>% 
            left_join(subj, by = "topmed_subject_id") %>%
            select(
            Individual_ID,
            unique_subject_key,
            topmed_subject_id,
            as.character(full_name))
        }
    
    # Caveat: Variable Type Mismatch Resolution
    if(variable == "sleep_duration_1"){
        harmonized_study$age_at_sleep_duration_1 = as.double(as.character(harmonized_study$age_at_sleep_duration_1))
    }
    
    if(var_type == "continuous"){
        all_plots[[i]] <- ggplot(data = harmonized_study, aes_string(full_name)) + geom_histogram(aes(y = ..ncount..)) + geom_density(aes(y = ..scaled..)) + labs(title = sprintf("Distribution of %s in %s", full_name, unit$name))
    } else {
        all_plots[[i]] <- ggplot(data = harmonized_study) + geom_bar(aes_string(full_name)) + labs(title = sprintf("Distribution of %s in %s", full_name, unit$name))
    }
    
    all_studies <- bind_rows(all_studies, harmonized_study)
     
}

# Plot Output File
pdf_file <- paste(variable, "pdf", sep = ".")
pdf(pdf_file)

if(var_type == "continuous"){
    ggplot(data = all_studies, aes_string(full_name)) + geom_histogram(aes(y = ..ncount..)) + geom_density(aes(y = ..scaled..)) + labs(title = sprintf("Overall Distribution of %s", full_name))
} else {
    ggplot(data = all_studies) + geom_bar(aes_string(full_name)) + labs(title = sprintf("Overall Distribution of %s", full_name))
}

multiplot(all_plots)

dev.off()

print(nrow(all_studies))

file_name <- paste(full_name, "harmonized", sep = "_")
file_name <- paste(file_name, "tsv", sep = ".")
readr::write_tsv(all_studies, file_name)