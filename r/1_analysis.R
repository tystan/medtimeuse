




# ---- libs ----

suppressPackageStartupMessages({
  require("compositions")
  
  require("dplyr")
  require("tidyr")
  require("readr")
  require("forcats")
  library("ggplot2")
  
  library("knitr")
  
  require("lme4")
  require("lmerTest")
  library("optimx")
  library("performance")
})



# ---- funcs_and_consts ----


add_alpha <- function(col, alpha = 1) {
  apply(
    sapply(col, col2rgb) / 255, 2,
    function(x) rgb(x[1], x[2], x[3], alpha = alpha)
  )
}


stage_ins_col <- add_alpha(c("cyan", "magenta"), 0.25)
stage_out_col <-  add_alpha(c("cyan", "magenta"), 0.75)
names(stage_ins_col) <- names(stage_out_col) <- NULL


med_ins_col <- add_alpha(c("orange", "purple"), 0.25)
med_out_col <-  add_alpha(c("orange", "purple"), 0.75)
names(med_ins_col) <- names(med_out_col) <- NULL



pal_use <- "Plasma" # "Temps", "Zissou 1"
plas_pal <- hcl.colors(n = 10, palette = pal_use, rev = FALSE)
sed_ins_col <- add_alpha(plas_pal, 0.25)
sed_out_col <-  add_alpha(plas_pal, 0.75)
names(sed_ins_col) <- names(sed_out_col) <- NULL



pal_use <- "Viridis" 
vir_pal <- hcl.colors(n = 11, palette = pal_use, rev = FALSE)
ach_ins_col <- add_alpha(vir_pal, 0.25)
ach_out_col <-  add_alpha(vir_pal, 0.75)
names(ach_ins_col) <- names(ach_out_col) <- NULL


# ---- read ----



sedach_dat <-
  read_rds("dat/sedach_dat.rds") %>%
  as_tibble(.)


sedach_dat$TrialStage <- fct_infreq(sedach_dat$TrialStage)

# ---- ilr_calcs ----

# these are the time-use compositions
time_use_cols <- paste0("tu_", c("sl", "sed", "lp", "mv"))
tu_dat <- sedach_dat[, time_use_cols]

# make isometric log ratios for compositional analysis of time-use composition
tu_comp <- acomp(tu_dat)
tu_ilrs <- as.data.frame(ilr(tu_comp))
D <- ncol(tu_ilrs)
colnames(tu_ilrs) <- paste0("ilr", 1:D)

# add ilrs to analysis dataset
sedach_dat <- bind_cols(sedach_dat, tu_ilrs)
colnames(sedach_dat)


# ---- exploratory_1 ----


corr_dat <-
  sedach_dat %>%
  group_by(TrialStage) %>%
  summarise(
    sp_cor = cor(x = sed_score, y = ach_score, method =  "spearman"),
    .groups = "drop"
  ) %>%
  mutate(
    cor_lab = sprintf("%s Spearman corr: %1.3f", TrialStage, sp_cor),
    x0 = 3,
    y0 = 9.5
  )

sedach_dat %>%
  ggplot(., aes(x = sed_score, y = ach_score)) +
  geom_jitter(
    aes(fill = TrialStage, col = TrialStage),
    size = 2, shape = 21,
    width = 0.2, height = 0.2
  ) +
  geom_text(
    inherit.aes = FALSE, 
    data = corr_dat, 
    aes(x = x0, y = y0, label = cor_lab)
  ) +
  scale_fill_manual(values = stage_ins_col) +
  scale_color_manual(values = stage_out_col) +
  facet_wrap(~ TrialStage, labeller = label_both) +
  theme_bw() +
  labs(
    fill = "Trial\nstage",
    col = "Trial\nstage",
    x = "Sedative load score",
    y = "Anticholinergic load score"
  )


# ---- exploratory_2 ----


sedach_dat %>% 
  pivot_longer(
    ., 
    cols = c(sed_score, ach_score), 
    names_to = "med", 
    values_to = "val"
  ) %>%
  pivot_wider(
    ., 
    id_cols = c(StudyID, med), 
    names_from = "TrialStage", 
    values_from = "val", 
    values_fill = NA
  )  %>%
  mutate(
    med = 
      case_when(
        med == "ach_score"  ~ "Anticholinergic",
        med == "sed_score" ~ "Sedative",
        TRUE                 ~ "ERROR!!!"
      )
  ) %>%
  na.omit(.) %>%
  ggplot(., aes(x = Baseline, y = `12 Months`)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_jitter(
    aes(fill = med, col = med),
    size = 2, shape = 21,
    width = 0.2, height = 0.2
  ) +
  xlim(-0.75, 10) +
  ylim(-0.75, 10) +
  scale_fill_manual(values = med_ins_col) +
  scale_color_manual(values = med_out_col) +
  facet_wrap(~ med, labeller = label_both) +
  theme_bw() +
  labs(
    fill = "Medicine load\nscore",
    col = "Medicine load\nscore"
  )


# ---- exploratory_3a ----

d_tu_long <- 
  sedach_dat %>%
  dplyr::select(-starts_with("ilr")) %>%
  pivot_longer(
    ., 
    cols = starts_with("tu_"), 
    names_to = "timeuse_cat", 
    values_to = "val"
  ) %>%
  mutate(
    timeuse_cat = gsub("tu_", "", timeuse_cat),
    timeuse_cat = 
      case_when(
        timeuse_cat == "sl"  ~ "Sleep",
        timeuse_cat == "sed" ~ "Sedentary",
        timeuse_cat == "lp"  ~ "Light PA",
        timeuse_cat == "mv"  ~ "Mod-to-vig PA",
        TRUE                 ~ "ERROR!!!"
      ),
    timeuse_cat = fct_inorder(timeuse_cat)
  )


d_tu_long %>%
  ggplot(., aes(x = TrialStage, y = val, group = StudyID)) +
  geom_line(linewidth = 0.5, alpha = 0.25) + 
  geom_jitter(
    aes(fill = factor(sed_score), col = factor(sed_score)),
    size = 2, shape = 21,
    width = 0.05, height = 0
  ) +
  scale_fill_manual(values = sed_ins_col) +
  scale_color_manual(values = sed_out_col) +
  facet_wrap(~ timeuse_cat , scales = "free_y", labeller = label_both, nrow = 1) +
  theme_bw() +
  labs(
    fill = "Sedative\nload score",
    col = "Sedative\nload score",
    y = "Minutes in time-use category"
  )

# ---- exploratory_3b ----

d_tu_long %>%
  ggplot(., aes(x = TrialStage, y = val, group = StudyID)) +
  geom_line(linewidth = 0.5, alpha = 0.25) + 
  geom_jitter(
    aes(fill = factor(ach_score), col = factor(ach_score)),
    size = 2, shape = 21,
    width = 0.05, height = 0
  ) +
  scale_fill_manual(values = ach_ins_col) +
  scale_color_manual(values = ach_out_col) +
  facet_wrap(~ timeuse_cat , scales = "free_y", labeller = label_both, nrow = 1) +
  theme_bw() +
  labs(
    fill = "Anticholinergic\nload score",
    col = "Anticholinergic\nload score",
    y = "Minutes in time-use category"
  )




# ---- exploratory_4a ----

d_ilr_long <- 
  sedach_dat %>%
  dplyr::select(-starts_with("tu_")) %>%
  pivot_longer(
    ., 
    cols = ilr1:ilr3, 
    names_to = "ilr_no", 
    values_to = "val"
  )


d_ilr_long %>%
  ggplot(., aes(x = TrialStage, y = val, group = StudyID)) +
  geom_line(linewidth = 0.5, alpha = 0.25) + # +aes(col = facility_group)
  geom_jitter(
    aes(fill = factor(sed_score), col = factor(sed_score)),
    size = 2, shape = 21,
    width = 0.05, height = 0
  ) +
  scale_fill_manual(values = sed_ins_col) +
  scale_color_manual(values = sed_out_col) +
  facet_wrap(~ ilr_no , scales = "free_y", labeller = label_both) +
  theme_bw() +
  labs(
    fill = "Sedative\nload score",
    col = "Sedative\nload score"
  )

# ---- exploratory_4b ----

d_ilr_long %>%
  ggplot(., aes(x = TrialStage, y = val, group = StudyID)) +
  geom_line(linewidth = 0.5, alpha = 0.25) + # +aes(col = facility_group)
  geom_jitter(
    aes(fill = factor(ach_score), col = factor(ach_score)),
    size = 2, shape = 21,
    width = 0.05, height = 0
  ) +
  scale_fill_manual(values = ach_ins_col) +
  scale_color_manual(values = ach_out_col) +
  facet_wrap(~ ilr_no , scales = "free_y", labeller = label_both) +
  theme_bw() +
  labs(
    fill = "Anticholinergic\nload score",
    col = "Anticholinergic\nload score"
  )



# ---- analysis_data ----

# create stacked dataset for multi-level model, because the dependent variable will be
# the activity composition ILRS which are multivariate (there are 3 of them),
# and lmer can't handle multi-variate dependent variables.


dat_lng <- 
  sedach_dat %>%
  dplyr::select(-starts_with("tu_")) %>%
  pivot_longer(
    ., 
    cols = ilr1:ilr3, 
    names_to = "ilr.no", 
    values_to = "val"
  )




# ---- mod_sed ----

# sedative load

set.seed(123)

mod_sed <- 
  lmer(
    val ~ -1 + 
      ilr.no + 
      ilr.no:TrialStage + ilr.no:sed_score + 
      ilr.no:TrialStage:sed_score +
      (0 + ilr.no | StudyID), 
    data = dat_lng, 
    control = lmerControl(
      optimizer = "Nelder_Mead",
      check.conv.singular = 
        .makeCC(action = "ignore", tol = formals(isSingular)$tol)
    )
  )

summary(mod_sed)



car::Anova(mod_sed, test.statistic = "F", type = "III")


check_model(
  mod_sed,
  check = c("reqq", "qq", "linearity", "homogeneity", "outliers", "pp_check")
)





# ---- mod_ach ----

# Anti-cholinergic load


mod_ach <-
  lmer(
    val ~ 
      -1 + ilr.no + 
      ilr.no:TrialStage + ilr.no:ach_score + 
      TrialStage:ach_score:ilr.no +
      (0 + ilr.no | StudyID),
    data = dat_lng, 
    control = lmerControl(
      optimizer = "bobyqa",
      check.conv.singular = 
        .makeCC(action = "ignore", tol = formals(isSingular)$tol)
    )
  )

summary(mod_ach)

car::Anova(mod_ach, test.statistic = "F", type = "III") 

check_model(
  mod_ach,
  check = c("reqq", "qq", "linearity", "homogeneity", "outliers", "pp_check")
)



# ---- pred_df ----


newd1 <- 
  expand.grid(
    ilr.no = c("ilr1", "ilr2", "ilr3"),
    TrialStage = c("Baseline", "12 Months"),
    score = 4
  ) 
rownames(newd1) <- apply(newd1, 1, paste, collapse = "_")
newd2 <- 
  expand.grid(
    ilr.no = c("ilr1", "ilr2", "ilr3"),
    TrialStage = c("12 Months", "Baseline"),
    score = 2
  )
rownames(newd2) <- apply(newd2, 1, paste, collapse = "_")

newd <- rbind(newd1, newd2)

newd <- 
  newd %>% 
  dplyr::filter(
    (TrialStage == "Baseline" & score == 2) |
      (TrialStage == "12 Months" & score == 4)
  )

newd_sed <- 
  newd %>%
  rename(sed_score = score)

newd_ach <- 
  newd %>%
  rename(ach_score = score)



get_sed_diff <- function(.) {
  pred_val <- predict(., newdata = newd_sed, re.form = NA)
  pred_newd <- cbind.data.frame(pred_val, newd_sed)
  pred_newd_w <- spread(pred_newd, key = ilr.no, value = pred_val)
  ilr_cols <- grepl("ilr", colnames(pred_newd_w))
  time_use <- 1440 * unclass(ilrInv(pred_newd_w[, ilr_cols]))
  colnames(time_use) <- c("sl", "sed", "lp", "mv")
  # return(time_use[2, "sed"] - time_use[1, "sed"])
  return(time_use[2, ] - time_use[1, ])
}

cat(
  "This is expected change in minutes to the time-use composition\n",
  "when going from sed load = 2 to sed load = 4 from baseline to 12 months.\n"
)
get_sed_diff(mod_sed) %>% 
  tibble(time_use_cat = names(.), minutes = .) %>%
  kable(., digits = 1)



get_ach_diff <- function(.) {
  pred_val <- predict(., newdata = newd_ach, re.form = NA)
  pred_newd <- cbind.data.frame(pred_val, newd_ach)
  pred_newd_w <- spread(pred_newd, key = ilr.no, value = pred_val)
  ilr_cols <- grepl("ilr", colnames(pred_newd_w))
  time_use <- 1440 * unclass(ilrInv(pred_newd_w[, ilr_cols]))
  colnames(time_use) <- c("sl", "sed", "lp", "mv")
  # return(time_use[2, "sed"] - time_use[1, "sed"])
  return(time_use[2, ] - time_use[1, ])
}

cat(
  "This is expected change in minutes to the time-use composition\n",
  "when going from anticholinergic load = 2 to sed load = 4\n",
  "from baseline to 12 months.\n"
)
get_ach_diff(mod_ach) %>% 
  tibble(time_use_cat = names(.), minutes = .) %>%
  kable(., digits = 1)



