#===============================================================================
# 2021-02-23 -- outliving Demographic Research reproducibility package
# Replicate all the figures and calculations
# Marie-Pier Bergeron-Boucher, mpbergeron@sdu.dk (examples)
# Ilya Kashnitsky, ilya.kashnitsky@gmail.com (dataviz + code review)
#===============================================================================


# prepare session ---------------------------------------------------------
library(tidyverse)
library(magrittr)
library(patchwork)
library(hrbrthemes)
library(msm)


# ggplot theme setup ------------------------------------------------------

theme_set(
    theme_minimal(base_family = font_rc, base_size = 14)+
        theme(
            legend.position = "none",
            panel.grid.minor = element_blank(),
            line = element_line(lineend = "round")
        )
)

# colors
pal <- c("#ef5350", "#002171")

# Fig 1 -- French males and females -----------------------------------

df_fr <- read.table("dat/LifeTable_Sex_HMD.txt", skip = 2)

# plot
df_fr %>%
    filter(Age > 64) %>%
    group_by(Sex) %>%
    mutate(value = dx %>% prop.table) %>%
    ungroup() %>%
    ggplot(aes(Age, value, color = Sex, fill = Sex))+
    geom_path()+
    geom_ribbon(
        aes(ymax = value, ymin = 0, xmin = 65, xmax = 110),
        alpha = .25
    )+
    geom_hline(yintercept = 0, size = 1)+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal)+
    scale_x_continuous(breaks = c(65, seq(70, 110, 10)))+
    labs(x = "Age", y = "Probability")+
    annotate(
        "text", x = c(110, 80), y = .045,
        size = 7, color = pal,
        family = font_rc, hjust = 1,
        label = c("Females", "Males")
    )

fig_1 <- last_plot()

ggsave(
    "fig/figure-1-france.pdf", fig_1,
    width = 5, height = 3.5, device = cairo_pdf
)

# calculate the overlap
df_fr %>%
    filter(Age > 64) %>%
    group_by(Sex) %>%
    mutate(value = dx %>% prop.table) %>%
    ungroup() %>%
    transmute(Age, Sex, value) %>%
    pivot_wider(names_from = Sex) %>%
    mutate(
        diff = Males - Females,
        sign = diff > 0
    ) %>%
    group_by(sign) %>%
    summarise(area = diff %>% sum)


# Figure 2 -- Graphical representation of Phi ------------------------

df_graph <- df_fr %>%
    transmute(Age, Sex, dx = dx/1e5, lx = lx/1e5) %>%
    filter(Age > 64) %>%
    pivot_wider(names_from = Sex, values_from = dx:lx) %>%
    mutate(
        dx_f = dx_Females %>% prop.table
    )

df_graph %>%
    ggplot(aes(Age))+
    geom_path(aes(y = lx_Males), color = "#002171")+
    geom_hline(yintercept = 0, size = .5)+
    scale_x_continuous(breaks = c(65, seq(70, 110, 10)))+
    labs(x = "Age", y = "Survival")+
    annotate(
        "text", x = c(95), y = c(.6),
        size = 7, color = "#002171",
        family = "Palatino", fontface = 2, hjust = 1,
        label = c(expression(italic(l[1](x))))
    )

a <- last_plot()


df_graph %>%
    ggplot(aes(Age))+
    geom_path(aes(y = dx_f), color = "#ef5350")+
    geom_ribbon(
        aes(ymax = dx_f*lx_Males, ymin = 0, xmin = 65, xmax = 110),
        alpha = .5, color = NA, fill = "#26a69a"
    )+
    geom_ribbon(
        aes(ymax = dx_f, ymin = dx_f*lx_Males, xmin = 65, xmax = 110),
        alpha = .5, color = NA, fill = "#cdcdcd"
    )+
    geom_hline(yintercept = 0, size = .5)+
    scale_y_continuous(breaks = seq(0, .05, .01))+
    scale_x_continuous(breaks = c(65, seq(70, 110, 10)))+
    theme_minimal(base_family = font_rc, base_size = 14)+
    labs(x = "Age", y = "Probability")+
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank()
    )+
    annotate(
        "text", x = c(87.5, 92.5), y = c(.05, .004),
        size = 7, color = c("#ef5350", "#26a69a"),
        family = "Palatino", hjust = 1, fontface = "italic",
        label = c(
            expression(italic(d[2](x))),
            expression(italic(l[1](x)~d[2](x)))
        )
    )

b <- last_plot()

fig_2 <-  a + b + plot_annotation(tag_levels = "A")

ggsave(
    "fig/figure-2.pdf", fig_2,
    width = 8, height = 4, device = cairo_pdf
)



# Example 1: Lifespan Sex France --------------------------------------------

# This example uses French life tables for 2018, the same dataset as in Fig 1
# Life tables extracted from the HMD: MD. 2020. Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany), Available at www.mortality.org (Accessed on September 9)

head(df_fr)

# Check available population
unique(df_fr$Sex)

# Age range: from age 65 only
age <- 65:110
n <- length(age)

#### Function to simulate lifespan ####
sim <- function(nsim, Mx, s.age){
    l <- length(Mx) - 1
    lt <- as.data.frame(msm::rpexp(n = nsim, rate = Mx, 0:l))
    names(lt) <- "lt"
    lt$lt <- lt$lt + s.age
    return(lt)}

# number of simulations
nsim <- 1e6 # a million simulations


#### The data by sex ####

## Females

# Life table
LT_F <- df_fr[df_fr$Sex == "Females", ]

# Simulations
simF <- sim(nsim, LT_F$mx[age + 1], s.age = 65)


## Males

# Life table
LT_M <- df_fr[df_fr$Sex == "Males", ]

# Simulations
simM <- sim(nsim, LT_M$mx[age + 1], s.age = 65)


#### The dataset ####
FreqTableF <- prop.table(table(simF$lt))
dxF <- data.frame(
    Age = as.numeric(rownames(FreqTableF)),
    dx_females = as.numeric(FreqTableF)
)

FreqTableM <- prop.table(table(simM$lt))
dxM <- data.frame(
    Age = as.numeric(rownames(FreqTableM)),
    dx_males = as.numeric(FreqTableM)
)

dta_dx<-merge(dxF, dxM, by="Age", all.x = TRUE, all.y = TRUE)
dta_dx[is.na(dta_dx)] <- 0

dta_dx$Dx_females <- cumsum(dta_dx$dx_females)
dta_dx$Dx_males <- cumsum(dta_dx$dx_males)


#### phi based on random pairing ####
RP<-data.frame(
    Females=sample(simF$lt),
    Males=sample(simM$lt)
)

RP$MoutliveF <- 0
RP$MoutliveF[RP$Males > RP$Females] <- 1

phi_RP <- length(RP$MoutliveF[RP$MoutliveF == 1]) / nrow(RP)
phi_RP


#### phi based on eq. 6 ####
phi_eq6 <- sum(dta_dx$dx_males * dta_dx$Dx_females)
phi_eq6

#### phi based on eq. 1 ####
phi_eq1 <- sum(dta_dx$dx_females * (1 - dta_dx$Dx_males))
phi_eq1

#### phi based on discrete approximation (Bergeron-Boucher et al. 2020) ####
dxF <- LT_F$dx[age + 1] / sum(LT_F$dx[age + 1])
dxM <- LT_M$dx[age + 1] / sum(LT_M$dx[age + 1])
DxM <- cumsum(dxM)
lxM <- c(1, 1 - DxM[1:(n - 1)])

phi_discrete <- sum(dxF[1:(n - 1)] * lxM[2:n]) + sum(dxF * dxM) / 2
phi_discrete

#### Comparing approaches ####

# Random pairing
round(phi_RP, 4)
# Equation 1
round(phi_eq1, 4)
# Equation 6
round(phi_eq6, 4)
# Discrete
round(phi_discrete, 4)




# Example 2: Lifespan Race US -----------------------------------------------

#### R code to calculate phi by race in the United States, both sexes combined
#### Based on data from the CDC: CDC. 2019. "United States Life Tables, 2017."  National Vital Statistics Report 68 (7)


#### Load data ####
df_cdc <- read.table("dat/LifeTable_Race_CDC.txt", header = T, skip = 2)

head(df_cdc)

# Check available population
unique(df_cdc$population)


#### The data by race ####

## Non-hispanic white population

# Life table
LT_w <- df_cdc[df_cdc$population == " the non-Hispanic white population", ]

# dx
dxW <- LT_w$dx / LT_w$lx[1]

# lx
lxW <- LT_w$lx / LT_w$lx[1]


## Non-hispanic black population

# Life table
LT_b <- df_cdc[df_cdc$population == " the non-Hispanic black population", ]

# dx
dxB <- LT_b$dx / LT_b$lx[1]

# lx
lxB <- LT_b$lx / LT_b$lx[1]

#### phi: probability that non-hispanic black outlive non-hispanic white ####
n_race <-length(dxW)

# Phi
phi_race <- sum(dxW[1:(n_race - 1)] * lxB[2:n_race]) + sum(dxB * dxW) / 2
phi_race

# Complement of phi: probability that non-hispanic white outlive non-hispanic black

phi_race_comp <- sum(dxB[1:(n_race - 1)] * lxW[2:n_race]) + sum(dxB * dxW) / 2
phi_race_comp

# Sum check
phi_race + phi_race_comp




# Example 3: Lifespan Education Italy ---------------------------------------

#### R code to calculate phi by education degree for males in Italy
#### Based on data from ISTAT: Life tables by educational attainment for the year 2012. Available at http://dati.istat.it/Index.aspx?DataSetCode=DCIS_MORTALITA1&Lang=en (Accessed on September 8).


#### Load data ####
df_ed <- read.table("dat/LifeTable_Education_ISTAT.txt", header = T, skip = 2)

head(df_ed)

# Check available population
unique(df_ed$Degree)

# Age range: from age 30 only
age <- 30:90


#### The data by education level ####

## University degree

# Life table
LT_U <- df_ed[df_ed$Degree == "University", ]

# dx
dxU <- LT_U$dx[age + 1] / sum(LT_U$dx[age + 1])

# define n
n_ed <- length(dxU)

# lx
lxU <- c(1, 1 - cumsum(dxU)[1:(n_ed - 1)])


## Grade 8

# Life table
LT_G8 <- df_ed[df_ed$Degree == "Grade8", ]

# dx
dxG8 <- LT_G8$dx[age + 1] / sum(LT_G8$dx[age + 1])

# lx
lxG8 <- c(1, 1 - cumsum(dxG8)[1:(n_ed - 1)])


#### phi: probability that low educated outlive high educated ####

# Phi
phi_ed <- sum(dxU[1:(n_ed - 1)] * lxG8[2:n_ed]) + sum(dxU * dxG8) / 2
phi_ed


# Complement of phi: probability that high educated outlive low educated
phi_ed_comp <- sum(dxG8[1:(n_ed - 1)] * lxU[2:n_ed]) + sum(dxU * dxG8) / 2
phi_ed_comp


# Sum check
phi_ed + phi_ed_comp




# Example 4: Income and Sex US --------------------------------------------

#### R code to calculate phi for income distribution by sex in the US
#### Based on data from the US Census Bureau: U.S. Census Bureau. 2020. PINC-01. Selected Characteristics of People 15 Years and Over, by Total Money Income, Work Experience, Race, Hispanic Origin, and Sex. 2019 Total Work Experience.


#### Load data ####
df_inc<- read.table(
    "dat/Population_byIncomeSex_USCB.txt", header = T, skip = 2, sep = ""
)

head(df_inc)


#### The data by sex ####

## Females
females <- df_inc$Females

# distribution
dxF_inc <- females / sum(females)

# cumulative distribution
FxF_inc <- cumsum(dxF_inc)


## Males
males <- df_inc$Males

# distribution
dxM_inc <- males / sum(males)

# cumulative distribution
FxM_inc <- cumsum(dxM_inc)


#### phi: probability that females have a higher income than males ####
n_inc <- length(dxF_inc)

# Phi
phi_inc <-
    sum(dxF_inc[2:n_inc] * FxM_inc[1:(n_inc - 1)]) +
    sum(dxF_inc * dxM_inc) /
    2
phi_inc


# Complement of phi: probability that males have a higher income than females
phi_inc_comp <-
    sum(dxM_inc[2:n_inc] * FxF_inc[1:(n_inc - 1)]) +
    sum(dxF_inc * dxM_inc) /
    2
phi_inc_comp


# Sum check
phi_inc + phi_inc_comp




# Example 5: Populations Italy and Ireland --------------------------------

#### R code to calculate phi for population distribution differences between Italy and Ireland
#### Based on data from the EUROSTAT

#### Load data ####
df_eu <- read.csv("dat/Population_Italy_Ireland.txt", header = T, sep = ",")

head(df_eu)
dim(df_eu)


#### The data by country ####

## Irish population
pop_Ireland <- as.numeric(as.character(df_eu$Value[df_eu$GEO == "Ireland"]))

# dx
dxIR <- pop_Ireland / sum(pop_Ireland)

# Dx
DxIR <- cumsum(dxIR)


## Italian population
pop_Italy <- as.numeric(as.character(df_eu$Value[df_eu$GEO == "Italy"]))

# dx
dxIT <- pop_Italy / sum(pop_Italy)

# Dx
DxIT <- cumsum(dxIT)



#### phi: probability that an Irish be older than an Italian ####
n_eu <- length(dxIT)

# Phi
phi_eu <- sum(dxIR[2:n_eu] * DxIT[1:(n_eu - 1)]) + sum(dxIT * dxIR) / 2
phi_eu

# Complement of phi: probability that males have a higher income than females
phi_eu_comp <-
    sum(dxIT[2:n_eu] * DxIR[1:(n_eu - 1)]) +
    sum(dxIT * dxIR) /
    2
phi_eu_comp

# Sum check
phi_eu + phi_eu_comp




