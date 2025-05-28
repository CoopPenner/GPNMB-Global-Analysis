library(dplyr)
library(survival)
library(ggplot2)
library(readxl) 
library(lubridate)  # for time calculations
library(lme4)

# ---- Initializing ----


# Read in Excel file from INDD
dataTable <- read_excel("/Volumes/PC60/InqueryDatasets/FinalizedSets/ALS_all_Clinical_data.xlsx")

# Extract relevant columns
ptID <- dataTable$INDDID                          # individual patient IDs
visitDate <- dataTable$VisitDate                 # visit dates
Weight <- dataTable$Weight                       # weight
FRS <- dataTable$FRSTotal                        # ALSFRS score
ageAtDeath <- dataTable$AgeatDeath               # age at death
ageAtOnset <- dataTable$AgeatOnset               # age at Onset
ageAtDiagnosis <- dataTable$AgeatDiagnosed       # age at Diagnosis
bioSex <- dataTable$Sex                          # biological sex
snp <- dataTable$rs199347                        # rs199347 stats
geneStat <- dataTable$Mutation_Summary           # all mutations

# Specifically filter for c9 
c9Pos <- grepl("C9orf72", geneStat)             
# Pull birthdate and test date 
testDate <- dataTable$VisitDate
birthDay <- dataTable$DOB

# Date of tracheostomy
trachDate <- dataTable$TracheostomyDate

# Diagnostic category
dx <- dataTable$ElEscorialVisit

# Site of symptom onset
onsetSite <- dataTable$ALSSymptomOnsetSite

# Calculate age at each test (in years), accounting for leap years
ageAtTest <- as.numeric(difftime(testDate, birthDay, units = "days")) / 365.25



#I have some trouble converting date formats, again prob me being unfamiliar with R syntax
safe_as_date <- function(x) {
  suppressWarnings(as.Date(x, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d-%b-%Y")))
}

# Convert date columns to Date class if not already
dataTable <- dataTable %>%
  mutate(across(ends_with("Date"), safe_as_date))

# ---- populating tables ----
#output summary data
summaryData <- dataTable %>%
  group_by(INDDID) %>%
  summarize(
    snpStatus = unique(rs199347)[1],
    sex = unique(Sex)[1],
    ageAtDiag = unique(AgeatDiagnosed)[1],
    ageAtDeath = unique(AgeatDeath)[1],
    trachDate = unique(TracheostomyDate)[1],
    deathDate = unique(DOD)[1],
    onsetDate = unique(ALSSymptomOnsetDate)[1],
    diagDate = unique(DiagnosisDate)[1],
    FRS = first(na.omit(FRSTotal)),
    firstVisitDate = min(VisitDate, na.rm = TRUE),
    numVisits = n(),
    .groups = "drop"
  )


summaryData <- summaryData %>% #populating my table in full
  mutate(
    startDate = coalesce(diagDate, onsetDate), #as in matlab if there's no diagnosis date we take onset date
    endDate = pmin(trachDate, deathDate, na.rm = TRUE), # we take whichever comes first tracheostomy or passing away
    timeToEvent = as.numeric(difftime(endDate, startDate, units = "days")) / 30.44,
    sexStat = case_when(
      sex == "Male" ~ 0,
      sex == "Female" ~ 1,
      TRUE ~ NA_real_
    ),
    snpStat = case_when(
      snpStatus == "CC" ~ 0, #converting to numeric
      snpStatus == "CT" ~ 1,
      snpStatus == "TT" ~ 2,
      TRUE ~ NA_real_
    ),
    estAgeAtDeath = ageAtDiag + timeToEvent / 12, #outputting age at Death (this is to avoid patients who have a nanned out age at death )
    avgAge = (ageAtDiag + estAgeAtDeath) / 2
  )


survData <- summaryData %>% # filtering 
  filter(
    !is.na(timeToEvent), #obvi remove nans
    !is.na(snpStat),  
    timeToEvent < 360, # I remove people who survived more than 30 years post diag
    ageAtDiag >= 18 # and people who were less than 18 when diagnosed initially
  ) %>%
  mutate(
    ageAtPlay = (estAgeAtDeath)  # z-score
  )

#fitting model using age and sex as covariates
cox_model <- coxph(Surv(timeToEvent) ~ snpStat + ageAtPlay + sexStat, data = survData)
summary(cox_model)

# outputting Hazard ratios
exp(coef(cox_model)) 


# as in the matlab code I am plotting raw survival data
survData$SNPgroup <- factor(survData$snpStat, levels = 0:2, labels = c("CC", "CT", "TT"))

cox_full <- coxph(Surv(timeToEvent, rep(1, nrow(survData))) ~ SNPgroup, data = survData)

fit <- survfit(Surv(timeToEvent, rep(1, nrow(survData))) ~ SNPgroup, data = survData)

# Create a new survival object, here we are just outputting raw survivals for plotting
surv_no_censor <- Surv(time = survData$timeToEvent, event = rep(1, nrow(survData)))

# Make sure SNP group is still defined correctly (had a problem w/ this early on, harder for me to debug)
survData$SNPgroup <- factor(
  survData$snpStat,
  levels = 0:2,
  labels = c("CC", "CT", "TT")
)


#plotting as we did in matlab code

colors <- c("cyan", "magenta", "green3")
snp_levels <- c("CC", "CT", "TT")

plot(NULL, xlim = c(0, 120), ylim = c(0, 1),
     xlab = "Time (Months)", ylab = "Survival Probability",
     main = " Survival Probability stratified by rs199347 status")

for (i in seq_along(snp_levels)) {
  group_data <- survData[survData$SNPgroup == snp_levels[i], ]
  print(group_data)
  # Cox model with no covariates
  cox_fit <- coxph(Surv(timeToEvent, rep(1, nrow(group_data))) ~ 1, data = group_data)
  
  # Extract cumulative hazard and convert to survival
  base_surv <- survfit(cox_fit)
  
  lines(base_surv$time, exp(-base_surv$cumhaz), col = colors[i], lwd = 2)
  
  print(exp(-base_surv$cumhaz))
  print(base_surv$time)
  
}

legend("topright", legend = snp_levels, col = colors, lwd = 2)








