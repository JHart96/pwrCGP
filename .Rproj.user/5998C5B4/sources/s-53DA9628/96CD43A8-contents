# pwrIRGP - Power Analysis for Interaction Rate Social Networks using a Gamma-Poisson model

This package can be used to estimate the accuracy of observed social networks built from interaction rate data. It uses a Gamma-Poisson model of interaction data to estimate the correlation between a sampled network and the true underlying network. The package also includes methods to conduct power analysis on nodal regression and to estimate the point at which increases in sampling effort lead to diminishing returns in accuracy.

If using this package, please cite the paper, which includes more detailed information about the methods:

Hart, J. D. A., Franks, D. W., Brent, L. J. N., & Weiss, M. N. (2021). Accuracy and Power Analysis of Social Interaction Networks. BioRxiv, 2021.05.07.443094. https://doi.org/10.1101/2021.05.07.443094

## How to use

### Installation

To install this package you'll need the `devtools` library, then you can install `pwrIRGP` using the `install_github` function.

```{r}
devtools::install_github("JHart96/pwrIRGP")
```

### Example

First import the package.

```{r}
library(pwrIRGP)
```

Simulate undirected network data with 8 nodes, 10 mean units of observation time per dyad, a social differentiation of 2, and a mean interaction rate of 0.5 interactions per unit time. Extract the symmetric square matrices for use in the `net_cor` function. This won't be necessary if you're using your own data.

```{r}
set.seed(1)
sim_data <- simulate_data_gp(8, 10, 2, 0.5)

X <- sim_data$X # 8 x 8 symmetric matrix of integer observation counts.
D <- sim_data$D # 8 x 8 matrix of positive real-valued sampling times.
```

Let's have a look at what X looks like:
```{r}
X
```
```
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    0    0    2    0    2    2    0    0
[2,]    0    0    0    7    0    0    6   13
[3,]    2    0    0   12    5   14    1    0
[4,]    0    7   12    0    0    0    0    4
[5,]    2    0    5    0    0    1    0    1
[6,]    2    0   14    0    1    0    0    0
[7,]    0    6    1    0    0    0    0   11
[8,]    0   13    0    4    1    0   11    0
```

Use the two matrices `X` and `D` to estimate the correlation between the sampled network and the true, underlying network. This will provide both a summary table of several properties of the data as well as a QQ diagnostic plot to qualitatively verify that the data fit the Gamma-Poisson model.

```{r}
summary_obj <- net_cor(X, D)
summary_obj
```
```
                                Estimate     SE Lower CI Upper CI
Observed Social Differentiation    1.540     NA       NA       NA
Mean Interactions Per Dyad         0.274     NA       NA       NA
Sampling Effort                    2.700     NA       NA       NA
Interaction Rate                   0.274 0.1070    0.137     0.55
Social Differentiation             1.780 0.3460    1.230     2.58
Correlation                        0.946 0.0288    0.869     0.98
```

Now we have an estimate for underlying interaction rate and social differentiation. Also we have an estimate of how well the observed network correlates with the true network. Each of these estimates is accompanied by a confidence interval estimate, which defaults to the 95% confidence interval.

If you want to conduct power analysis for nodal regression, you will need to extract social differentiation, interaction rate, and sampling times from the summary object and the data matrices. You will also need to provide an `effect` value. This is the effect size (correlation coefficient) and reflects the true relationship between the response and predictors in the regression (the effect size we'd see with perfect, infinite sampling). 

```{r}
social_differentiation <- summary_obj[5, 1]
interaction_rate <- summary_obj[4, 1]
sampling_times <- D # Matrix of sampling times OR Single value of mean sampling times

# Calculate power of nodal regression for effect size r = 0.5
effect <- 0.5
pwr_nodereg(8, effect, social_differentiation, interaction_rate, sampling_times)
```
```
$nodes
[1] 8

$effect
[1] 0.5

$soc_diff
[1] 1.78

$int_rate
[1] 0.274

$power
[1] 0.238
```
This shows us that we could expect a power of 24% given the properties of the data and the true effect size = 0.5.

If a different type of analysis is being conducted such as network subsetting, the diminishing returns/elbow estimator could be used to determine if sufficient data are available:
```{r}
pwr_elbow(social_differentiation, rho_max=0.99) # Use rho_max=0.99 as in the paper.
```
```
$sampling_effort
[1] 1.34

$correlation
[1] 0.8996477
```

The elbow method says we need a sampling effort of 1.34 to reach the optimal level of correlation, which is roughly 90%. From our run of `net_cor` we know that sampling effort is 2.7, more than enough to get the optimal level of correlation. This indicates that the network could probably be subset to create two networks if required.

The elbow method should be used a general guideline, and careful consideration should be taken when using it on data with a wide confidence interval. In this case it is advisible to run the power analysis at both ends of the confidence interval to fully quantify the uncertainty.
