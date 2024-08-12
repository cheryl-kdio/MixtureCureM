library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(rlang)
library(latex2exp)
library(rlang)
source("functions1.R")

############### TABLES per SCENARIO

# Compute Adjusted Results
compute_adjusted_results <- function(scenario_data, confidence_level) {
  adjusted_total <- list()
  adjusted_results <- list()
  ci_lower <- list()
  ci_upper <- list()
  
  z_score <- qnorm(1 - (1 - confidence_level) / 2)
  for (link in names(scenario_data$test_results)) {
    adjusted_total[[link]] <- 10000 - scenario_data$null_results[[link]]
    adjusted <- (1 - scenario_data$test_results[[link]] / adjusted_total[[link]]) * 100
    adjusted_results[[link]] <- adjusted
    std_error <- sqrt((adjusted / 100) * (1 - (adjusted / 100)) / adjusted_total[[link]])
    ci_lower[[link]] <- adjusted - z_score * std_error * 100
    ci_upper[[link]] <- adjusted + z_score * std_error * 100
  }
  
  return(list(adjusted_total = adjusted_total, adjusted_results = adjusted_results, ci_lower = ci_lower, ci_upper = ci_upper))
}

process_all <- function(scenario_data_file, confidence_level = 0.95) {
  # Load data
  scenario_data <- readRDS(scenario_data_file)
  
  adjusted_output <- compute_adjusted_results(scenario_data, confidence_level)
  
  # Format and compile results
  results_list <- list()
  for(model in names(adjusted_output$adjusted_results)) {
    model_results <- adjusted_output$adjusted_results[[model]]
    formatted_results <- sprintf("%0.2f (%0.2f; %0.2f)", model_results, adjusted_output$ci_lower[[model]], adjusted_output$ci_upper[[model]])
    results_table <- data.frame(
      Model = c(names(model_results), "nb_conv"),
      Results_with_CI = c(formatted_results, round(adjusted_output$adjusted_total[[model]]*100 / 10000,2)),
      row.names = NULL
    )
    results_list[[model]] <- results_table
  }
  
  # Transpose and rename columns
  transpose_and_rename <- function(df, model_name) {
    t_df <- as.data.frame(t(df), row.names = F)
    colnames(t_df) <- df$Model
    t_df <- cbind(Model = model_name, t_df)
    return(t_df[-1,])
  }
  
  # Apply the transpose function
  transposed_results <- lapply(names(results_list), function(model) {
    transpose_and_rename(results_list[[model]], model)
  })
  
  # Merge all results into a single data frame
  merged_results <- bind_rows(transposed_results)
  
  return(merged_results)
}

scenario1 <- process_all("resultats_sim/scenario1.rds")
scenario2 <- process_all("resultats_sim/scenario2.rds")
scenario3 <- process_all("resultats_sim/scenario3.rds")
scenario4 <- process_all("resultats_sim/scenario4.rds")

#merge
all_scenarios <- bind_rows(scenario1, scenario2, scenario3, scenario4, .id = "Scenario")
all_scenarios


library(xtable)
print(xtable(all_scenarios), type = "latex", include.rownames = FALSE)

#SCENARIO LOGLOG
scenario1_loglog <- process_all("resultats_sim/scenario1_loglog.rds")
scenario2_loglog <- process_all("resultats_sim/scenario2_loglog.rds")
scenario3_loglog <- process_all("resultats_sim/scenario3_loglog.rds")
scenario4_loglog <- process_all("resultats_sim/scenario4_loglog.rds")

#merge
all_scenarios_loglog <- bind_rows(scenario1_loglog, scenario2_loglog, scenario3_loglog, scenario4_loglog, .id = "Scenario")
all_scenarios_loglog

print(xtable(all_scenarios_loglog), type = "latex", include.rownames = FALSE)

#SCENARIO PROBIT
scenario1_probit <- process_all("resultats_sim/scenario1_probit.rds")
scenario2_probit <- process_all("resultats_sim/scenario2_probit.rds")
scenario3_probit <- process_all("resultats_sim/scenario3_probit.rds")
scenario4_probit <- process_all("resultats_sim/scenario4_probit.rds")

#merge
all_scenarios_probit <- bind_rows(scenario1_probit, scenario2_probit, scenario3_probit, scenario4_probit, .id = "Scenario")
all_scenarios_probit

print(xtable(all_scenarios_probit), type = "latex", include.rownames = FALSE)

########### Comptage du critère défaillant
scenario1_d <- readRDS("resultats_sim/scenario1.rds")
scenario2_d <- readRDS("resultats_sim/scenario2.rds")
scenario3_d <- readRDS("resultats_sim/scenario3.rds")
scenario4_d <- readRDS("resultats_sim/scenario4.rds")

count_default1 <- readRDS("resultats_sim/count_default_1.rds")
count_default2 <- readRDS("resultats_sim/count_default_2.rds")
count_default3 <- readRDS("resultats_sim/count_default_3.rds")
count_default4 <- readRDS("resultats_sim/count_default_4.rds")

compute_adjusted_default <- function(scenario_data,link,count_default) {
  adjusted_total <- 10000 - scenario_data$null_results[[link]]
  adjusted_both <- adjusted_total - scenario_data$test_results[[link]][["both"]]
  adjusted_overall <- adjusted_total - scenario_data$test_results[[link]][["overall"]]
  count_default_atleast <- round(unlist(count_default[[link]][c("nb_default_pi","nb_default_rn")])*100/adjusted_both,2)
  count_default_overall <- round(unlist(count_default[[link]][c("nb_default_pi_overall","nb_default_rn_overall")])*100/adjusted_overall,2)
  
  return(list(count_default_atleast,count_default_overall))
}

for(i in c("logit","loglog","probit")){
  print(compute_adjusted_default(scenario4_d,i,count_default4))
}



########### Visualisation des résultats de simulation

# Load data
# Create a list to store the adjusted results for each scenario
results_adjusted <- list(
  scenario1 = compute_adjusted_results(scenario1_d, 0.95)$adjusted_results,
  scenario2 = compute_adjusted_results(scenario2_d, 0.95)$adjusted_results,
  scenario3 = compute_adjusted_results(scenario3_d, 0.95)$adjusted_results,
  scenario4 = compute_adjusted_results(scenario4_d, 0.95)$adjusted_results
)

# Initialize a list to store plots for each scenario
plot_list <- list()

# Loop through each scenario to create plots
for (scenario in names(results_adjusted)) {
  results_list <- list()
  
  for (model in names(results_adjusted[[scenario]])) {
    model_results <- results_adjusted[[scenario]][[model]]
    results_table <- data.frame(
      model = model,
      conclusion = names(model_results),
      proportion = model_results,
      row.names = NULL
    )
    results_list[[model]] <- results_table
  }
  
  combined_results <- bind_rows(results_list)
  
  combined_results$conclusion <- factor(combined_results$conclusion, 
                                        levels = c("both", "on_avg", "at_least_one", "overall"),
                                        labels = c("Both Arms",  "On Average","At Least One Arm", "Overall"))
  
  plot_list[[scenario]] <- ggplot(combined_results, aes(x = model, y = proportion, fill = conclusion)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("black", "darkgrey", "grey", "lightgrey", "white")) +
    labs(y = "Proportion de rejet de pertinence du modèle",
         fill = "Règles de décision",
         x = "Fonctions de lien") +
    theme_minimal() +
    theme(legend.position = "none") 
}

# Combine the plots using patchwork
combined_plot <- plot_list[["scenario1"]] + plot_list[["scenario2"]] +
  plot_list[["scenario3"]] + plot_list[["scenario4"]] +
  plot_annotation(tag_levels = '1', tag_prefix = 'Scenario ', tag_sep = '.', tag_suffix = ')')

legend <- get_legend(plot_list[["scenario1"]]  + theme(legend.position = "top"))

final_plot <- combined_plot + plot_layout(guides = 'collect',nrow=2) & theme(legend.position = "right")

print(final_plot)



############# VARIATION TABLE

compute_results <- function(variable, values, result_type) {
  results_list <- lapply(1:length(values), function(i) {
    results <- compute_adjusted_results(variable[[i]], 0.95)
    adjusted_results <- results$adjusted_results
    nb_conv <- results$adjusted_total
    
    # Récupération des résultats pour chaque type de modèle
    list(
      logit_both = adjusted_results$logit[["both"]] / 100,
      logit_on_avg = adjusted_results$logit[["on_avg"]] / 100,
      logit_at_least_one = adjusted_results$logit[["at_least_one"]] / 100,
      logit_overall = adjusted_results$logit[["overall"]] / 100,
      logit_nb_conv = nb_conv$logit / 10000,
      
      loglog_both = adjusted_results$loglog[["both"]] / 100,
      loglog_on_avg = adjusted_results$loglog[["on_avg"]] / 100,
      loglog_at_least_one = adjusted_results$loglog[["at_least_one"]] / 100,
      loglog_overall = adjusted_results$loglog[["overall"]] / 100,
      loglog_nb_conv = nb_conv$loglog / 10000,
      
      probit_both = adjusted_results$probit[["both"]]/100,
      probit_on_avg = adjusted_results$probit[["on_avg"]]/100,
      probit_at_least_one = adjusted_results$probit[["at_least_one"]]/100,
      probit_overall = adjusted_results$probit[["overall"]]/100,
      probit_nb_conv = nb_conv$probit / 10000
    )
  })
  
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))
  if(result_type == "pi") {
    results_df$pi <- values
  }else if(result_type == "HR") {
    results_df$HR <- values
  }else if(result_type == "OR") {
    results_df$OR <- values
  }
  
  return(results_df)
}

# cure_var <- readRDS("resultats_sim/cure_variation.rds")
HR_var <- readRDS("resultats_sim/HR_variation.rds")
OR_var <- readRDS("resultats_sim/OR_variation.rds")
scenario5 <- readRDS("resultats_sim/scenario5.rds")
# Calcul des résultats pour pi_values
# pi_values <- seq(0.01, 0.8, 0.01)
# resultat_cure <- compute_results(cure_var, pi_values, "pi")

# Calcul des résultats pour HR_values
HR_values <- seq(0.5, 1.5, 0.05)
resultat_HR <- compute_results(HR_var, HR_values, "HR")

# Calcul des résultats pour OR_values
OR_values <- seq(0.5, 1.5, 0.05)
resultat_OR <- compute_results(OR_var, OR_values, "OR")

pi_values <- seq(0.025,0.8,0.025)
resultat_cure <- compute_results(scenario5, pi_values, "pi")

create_plot <- function(data, var_names) {
  data <- data %>% select(all_of(c(var_names)))
  data_long <- pivot_longer(data, cols = all_of(var_names[-1]), names_to = "Conclusions", values_to = "value")
  #refactor
  data_long$Conclusions <- factor(data_long$Conclusions, 
                                        levels = var_names[-1],
                                        labels = c("Both Arms", "On Average", "At Least One Arm", "Overall"))
  
  
  if (var_names[1]=="pi"){
    x <- TeX(r"( Proportion de guérison $\pi(x=1)$)")}
  else if(var_names[1]=="HR"){
    x <- TeX(r"( $ Hazard~ratio~des~non~guéris$)")}
  else if(var_names[1]=="OR"){
    x <- TeX(r"( $ Odds~ratio~de~la~guérison$)")}
  
  ggplot(data_long, aes(x =  !!sym(var_names[1]), y = value, color = Conclusions, shape = Conclusions)) +
    geom_point(size = 3) +
    geom_line(aes(group = Conclusions), linewidth = 1) +
    labs(#x = TeX(r"( Rapport de risques instantanés $\lambda(x=1)$)"), 
      y = TeX(r"(Proportion de rejet de $H_0$ ($\alpha$))", 
              output = "expression"), 
      x=x,
      colour = "Règles de décision") +
    scale_color_manual(values = c("black", "darkgrey", "grey", "lightgrey")) +
    scale_shape_manual(values = c(1, 2, 3, 4)) +
    theme_minimal() +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    guides(color = guide_legend(title = "Règle de décision"), shape = guide_legend(title = "Règle de décision"))
}


pi_logit_plot <- create_plot(resultat_cure, c("pi","logit_at_least_one","logit_overall"))
HR_logit_plot <- create_plot(resultat_HR, c("HR","logit_both","logit_on_avg","logit_at_least_one","logit_overall"))
OR_logit_plot <- create_plot(resultat_OR, c("OR","logit_both","logit_on_avg","logit_at_least_one","logit_overall"))

pi_loglog_plot <- create_plot(resultat_cure, c("pi","loglog_at_least_one","loglog_overall"))
HR_loglog_plot <- create_plot(resultat_HR, c("HR","loglog_both","loglog_on_avg","loglog_at_least_one","loglog_overall"),  1.45)
OR_loglog_plot <- create_plot(resultat_OR, c("OR","loglog_both","loglog_on_avg","loglog_at_least_one","loglog_overall"))

pi_probit_plot <- create_plot(resultat_cure, c("pi","probit_at_least_one","probit_overall"))
HR_probit_plot <- create_plot(resultat_HR, c("HR","probit_both","probit_on_avg","probit_at_least_one","probit_overall"))
OR_probit_plot <- create_plot(resultat_OR, c("OR","probit_both","probit_on_avg","probit_at_least_one","probit_overall"))


############# PI PLOT
combined_plot <- (pi_logit_plot + pi_loglog_plot+ pi_probit_plot ) +
  plot_annotation(tag_levels = c('A', '1'), tag_prefix = 'Fig. ',
                  tag_sep = '.', tag_suffix = ')')

legend <- get_legend(pi_logit_plot + theme(legend.position = "bottom"))

pi_plot <- combined_plot + plot_layout(guides = 'collect',nrow=3) & theme(legend.position = "right")

print(pi_plot)

############# HR PLOT

combined_plot <- (HR_logit_plot + HR_loglog_plot+ HR_probit_plot ) +
  plot_annotation(tag_levels = c('A', '1'), tag_prefix = 'Fig. ',
                  tag_sep = '.', tag_suffix = ')')

legend <- get_legend(HR_logit_plot + theme(legend.position = "bottom"))

HR_plot <- combined_plot + plot_layout(guides = 'collect',nrow=3) & theme(legend.position = "right")

print(HR_plot)

############# OR PLOT

combined_plot <- (OR_logit_plot + OR_loglog_plot+ OR_probit_plot ) +
  plot_annotation(tag_levels = c('A', '1'), tag_prefix = 'Fig. ',
                  tag_sep = '.', tag_suffix = ')')

legend <- get_legend(OR_logit_plot + theme(legend.position = "bottom"))

OR_plot <- combined_plot + plot_layout(guides = 'collect',nrow=3) & theme(legend.position = "right")

print(OR_plot)
