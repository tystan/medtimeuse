




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

pal_use <- "Classic Tableau" 
ct_pal <- palette.colors(n = 4, palette = pal_use)
timeuse_col <-  add_alpha(ct_pal, 0.75)
names(timeuse_col) <- NULL




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


# ---- exploratory_0 ----

join_dat <-
  sedach_dat %>% 
    select(StudyID, TrialStage, sed_score, ach_score) 

join_dat <-
  inner_join(
    join_dat %>% dplyr::filter(TrialStage == "Baseline"),
    join_dat %>% dplyr::filter(TrialStage == "12 Months"),
    "StudyID"
  ) %>%
    mutate(
      sed_ot = 
        case_when(
          sed_score.x < sed_score.y ~ "(c) sed increase",
          sed_score.x == sed_score.y ~ "(b) sed constant",
          sed_score.x > sed_score.y ~ "(a) sed decrease",
          TRUE ~ as.character(NA)
        ),
      ach_ot = 
        case_when(
          ach_score.x < ach_score.y ~ "(c) ach increase",
          ach_score.x == ach_score.y ~ "(b) ach constant",
          ach_score.x > ach_score.y ~ "(a) ach decrease",
          TRUE ~ as.character(NA)
        )
    )

incdec_tab <- 
  with(join_dat, table(sed_ot, ach_ot, useNA = "ifany")) 

incdec_tab %>%
  kable(.)

# sum(incdec_tab)
# sum(diag(incdec_tab)) / sum(incdec_tab)

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
        med == "ach_score" ~ "Anticholinergic",
        med == "sed_score" ~ "Sedative",
        TRUE               ~ "ERROR!!!"
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
    cols = starts_with("ilr"), 
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
  dplyr::select(-starts_with("tu_")) %>% # keep only ilrs not time-use vars
  pivot_longer(
    ., 
    cols = starts_with("ilr"), 
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




# ---- pred_const_over_time_setup ----



get_mod_pred <- function(mod, dat) {
  pred_val <- predict(mod, newdata = dat, re.form = NA)
  pred_newd <- cbind.data.frame(pred_val, dat)
  pred_newd_w <- 
    pivot_wider(pred_newd, names_from = "ilr.no", values_from = "pred_val")
  ilr_cols <- grepl("ilr", colnames(pred_newd_w))
  time_use <- 1440 * unclass(ilrInv(pred_newd_w[, ilr_cols]))
  colnames(time_use) <- c("sl", "sed", "lp", "mv")
  return(bind_cols(pred_newd_w, time_use))
}

new_sed <- expand.grid(
  ilr.no = c("ilr1", "ilr2", "ilr3"),
  TrialStage = c("Baseline", "12 Months"),
  sed_score = seq(0, 9, 1)
)

new_ach <- expand_grid(
  ilr.no = c("ilr1", "ilr2", "ilr3"),
  TrialStage = c("Baseline", "12 Months"),
  ach_score = seq(0, 11, 1)
)

preds <- predict(mod_sed, newdata = new_sed, re.form = NA)
pred_df <- cbind.data.frame(preds, new_sed)

sed_preds <-
  get_mod_pred(mod_sed, new_sed) %>% 
  rename(medload = sed_score) %>%
  mutate(med = "Sedative load")

sed_preds %>%
  arrange(desc(TrialStage), medload) %>%
  select(TrialStage, med, medload, everything(), -starts_with("ilr")) %>%
  kable(., digits = 0)


(sed_preds %>%
  dplyr::filter(TrialStage == "12 Months", medload == 4) %>%
  select(6:9)) -
(sed_preds %>%
  dplyr::filter(TrialStage == "Baseline", medload == 2) %>%
  select(6:9))

preds <- predict(mod_ach, newdata = new_ach, re.form = NA)
pred_df <- cbind.data.frame(preds, new_ach)


ach_preds <-
  get_mod_pred(mod_ach, new_ach) %>% 
  rename(medload = ach_score) %>%
  mutate(med = "Anticholinergic load")

ach_preds %>%
  arrange(desc(TrialStage), medload) %>%
  select(TrialStage, med, medload, everything(), -starts_with("ilr")) %>%
  kable(., digits = 0)


all_pred <- 
  bind_rows(sed_preds, ach_preds) %>%
  select(-starts_with("ilr"))


all_pred <-
  inner_join(
    all_pred %>% filter(TrialStage == "12 Months"),
    all_pred %>% filter(TrialStage == "Baseline"),
    c("med", "medload")
  )

all_pred <-
  all_pred %>%
  mutate(
    change_Sleep = sl.x - sl.y,
    change_Sedentary = sed.x - sed.y,
    change_LightPA = lp.x - lp.y,
    change_MVPA = mv.x - mv.y
  )

all_pred <-
  all_pred %>%
  select(-matches("\\.(x|y)", perl = TRUE))


# ---- specific_preds ----

(all_pred %>%
  dplyr::filter(med == "Sedative load", medload == 4) %>%
  select(3:6)) -
(all_pred %>%
  dplyr::filter(med == "Sedative load", medload == 2) %>%
  select(3:6))
###  change_Sleep change_Sedentary change_LightPA change_MVPA
###     -3.821428         14.51255      -4.866429   -5.824697

(all_pred %>%
    dplyr::filter(med != "Sedative load", medload == 4) %>%
    select(3:6))
### change_Sleep change_Sedentary change_LightPA change_MVPA
###        -1.92             28.3          -15.6       -10.8

# ---- pred_const_over_time ----

all_pred_lng <-
  all_pred %>%
  pivot_longer(
    ., cols = starts_with("change"), names_to = "timeuse", values_to = "time"
  ) %>%
  mutate(
    timeuse =  gsub("change_", "", timeuse),
    med = fct_inorder(med),
    timeuse = fct_inorder(timeuse)
  )

all_pred_lng %>%
  ggplot(., aes(x = medload, y = time, group = timeuse)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(colour = timeuse)) +
  geom_point(aes(colour = timeuse)) +
  facet_wrap( ~ med) +
  scale_colour_manual(values = timeuse_col) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  labs(
    x = "Medication Load at 12-months",
    y = "Change in Activity (min/d)",
    colour = "Activity"
  ) +
  theme(text = element_text(family = "serif"))


# ---- pred_df ----



get_pred_diff <- function(mod, dat) {
  time_use <- get_mod_pred(mod, dat)
  time_use <- time_use[, c("sl", "sed", "lp", "mv")]
  return(time_use[2, ] - time_use[1, ])
}

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


cat(
  "This is expected change in minutes to the time-use composition\n",
  "when going from sed load = 2 to sed load = 4 from baseline to 12 months.\n"
)
get_pred_diff(mod_sed, newd_sed) %>% 
  kable(., digits = 1)


cat(
  "This is expected change in minutes to the time-use composition\n",
  "when going from anticholinergic load = 2 to sed load = 4\n",
  "from baseline to 12 months.\n"
)
# somewhat of an extrapolation
get_pred_diff(mod_ach, newd_ach) %>% 
  kable(., digits = 1)



