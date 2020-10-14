library(ggplot2)
library(tidyverse)

num_vector <- rnorm(1000, mean = 10, sd = 5)
summary(num_vector)
mean_vector <- mean(num_vector)
sd_vector <- sd(num_vector)
quantile(num_vector, probs = seq(0, 1, 0.1))
# ggplot
hist(num_vector)
vector_df <- data.frame(NAME = num_vector)
head(vector_df)

ggplot(vector_df, aes(x = NAME)) +
    geom_histogram() +
    geom_vline(xintercept = mean_vector)

large_vector <- rnorm(1e7, mean = 10, sd = 5)
hist(large_vector, breaks = 10)

# For the standard normal distribution
# Plot the cumulative distribution function in the range .
# Plot the inverse cumulative distribution function for quantiles in 0.01 increment.
# Plot the density function in the range .
# What is the probability of observing a value greater than 2?
# What is the probability of observing a value between -2 and 2?
# What is the probability of observing a value more extreme than -2 or 2?


x_index = seq(-5,5,0.01)

dst <- tibble(
    x_coordinate = x_index,
    cdf_normal_distribution = pnorm(q = x_index, mean  = 0, sd = 1))

dst %>%
    ggplot(aes(x=x_coordinate, y=cdf_normal_distribution))+
    geom_point()

x_prob = seq(0,1,0.01)
dst <- tibble(
    x_prob = x_prob,
    inverse_cdf_normal = qnorm(p = x_prob, mean = 0, sd = 1))

dst %>%
    ggplot(aes(x = x_prob, y = inverse_cdf_normal)) +
    geom_point()

# density normal distribution
x_index = seq(-5,5,0.01)

dst <- tibble(
    x_coordinate = x_index,
    density_normal_distribution = dnorm(x = x_index, mean  = 0, sd = 1))

dst %>%
    ggplot()+
    geom_point(aes(x=x_coordinate, y=density_normal_distribution))

# What is the probability of observing a value greater than 2?

1-pnorm(2)
pnorm(2) - pnorm(-2)
1-(pnorm(2) - pnorm(-2))

# Statistical tests
# Use the summary() function to view some information about each column.
# Visualise the distribution of Sepal.Length , stratified by species.
# Is Sepal.Length length normally distributed? Overall? Within each species?
#  Is there a signi􀃘cant variation of Sepal.Length between the various species?

summary(iris)

# Visualise the distribution of Sepal.Length , stratified by species.

ggplot(iris, aes(Sepal.Length))+
    geom_histogram(color = "black")+
    facet_wrap(~Species, ncol = 1)

shapiro.test (iris$Sepal.Length)

species_setosa <- subset(iris, Species == "setosa",
                         select = Sepal.Length)
species_setosa

shapiro.test(species_setosa$Sepal.Length)

anova <- aov(formula = Sepal.Length ~ Species, data = iris)
summary(anova)

levels(iris$Species)

t.test(species_setosa$Sepal.Length, )

# Exercise Linear regression
summary(ChickWeight)
head(ChickWeight)
# Fit a linear mode to measure the effect of Time and Diet in the ChickWeight data set.
# Which diet leads to the fastest increase in body weight
# How much does weight increase per unit of Time for the top diet?
# Does the top diet drive an increase in body weight that is significantly faster than the next best diet?

formula <- formula(weight ~ Diet + Time)

lm(formula = formula, data = ChickWeight)

ggplot(ChickWeight, aes(x = Time, y = weight))+
    geom_point(aes(color = Diet))+
    geom_smooth(method = "lm")

# multiple test correction

log_counts <- read.csv("data/logcounts.csv", row.names = 1)
cell_md <- read.csv("data/cell_metadata.csv", row.names =1)

t.test(log_counts[1,],log_counts[2,])

gene1 <- data.frame("log_count" = as.numeric(log_counts[1,]), infection = cell_md$Infection)
t.test(formula = log_count~infection, data = gene1)

diff_exp <- function(gene_index, matrix, groups){
    gene_row <- data.frame(
        "log_count" = as.numeric(matrix[gene_index,]),
        infection = groups)
    test_result <- t.test(log_count~infection, gene_row)
    return(test_result[['p.value']])
}

diff_exp(3,log_counts,cell_md$Infection)

vapply(seq(1,nrow(log_counts)), diff_exp, numeric(1), matrix = log_counts, groups = cell_md$Infection)
names(p_values) <- rownames(log_counts)
p_values

# For each, build table with Fisher's test, how many belong in the signal pathway


