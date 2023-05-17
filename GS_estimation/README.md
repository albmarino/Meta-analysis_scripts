## Estimating expected C-values
table3 being Supplementary Table 3 and table1 being Supplementary Table 1, Weighted Least Square was calculated from C-values and assembly sizes in Supplementary Table 3 as follows:
```
library(dplyr)
m1 <- lm(Cvalue~Assembly_size, data=table3)
wt_m1 <- 1 / lm(abs(m1$residuals) ~ m1$fitted.values)$fitted.values^2
wls_m1 <- lm(Cvalue ~ Assembly_size, data = table3, weights=wt_m1)

#Predict values from assembly sizes with Quast_ContigN50 >= 50kb in table1
table1 <- filter(table1, Quast_ContigN50 >= 50kb)
as <- data.frame(Assembly_size = table1$Assembly_size)
corrected_as <- data.frame(Expected_Cvalue = predict(wls_m1,newdata = as))
```

## Testing reads effect
Method 1: z-test of the difference between the regression coefficients of WLS based on only SR and only LR
```
df_SR <- filter(table3, Reads_assembly %in% "SR")
df_LR <- filter(table3, Reads_assembly %in% c("LR", "LR-SR"))
m1_SR <- lm(Cvalue~Assembly_size, data=df_SR)
m1_LR <- lm(Cvalue~Assembly_size, data=df_LR)
wt_m1_SR <- 1 / lm(abs(m1_SR$residuals) ~ m1_SR$fitted.values)$fitted.values^2
wt_m1_LR <- 1 / lm(abs(m1_LR$residuals) ~ m1_LR$fitted.values)$fitted.values^2

wls_m1_SR <- lm(Cvalue ~ Assembly_size, data = df_SR, weights=wt_m1_SR)
wls_m1_LR <- lm(Cvalue ~ Assembly_size, data = df_LR, weights=wt_m1_LR)

compare.coeff <- function(coef_SR,se_SR,coef_LR,se_LR){
 return((coef_SR-coef_LR)/sqrt(se_SR^2+se_LR^2))
 }

coef_SR <- summary(wls_m1_SR)$coefficients[2,1]
se_SR <- summary(wls_m1_SR)$coefficients[2,2]
coef_LR <- summary(wls_m1_LR)$coefficients[2,1]
se_LR <- summary(wls_m1_LR)$coefficients[2,2]

pval <- 2*pnorm(-abs(compare.coeff(coef_SR,se_SR,coef_LR,se_LR)))
pval # 0.8127543
```

Method 2: test the effect of reads category (LR and SR) on the WLS slope
```
newheader <- "Reads_category"
newdf <- data.frame(matrix(nrow=0, ncol=length(newheader)))
for (i in table3$Reads_assembly) {
 if (i=="LR" || i == "LR-SR") {newdf <- rbind(newdf, "LR")}
 else if (i=="SR") {newdf <- rbind(newdf, "SR")}
 }
colnames(newdf) <- newheader
table3 <- cbind(table3, newdf)
m1 <- lm(Cvalue ~ Assembly_size + Reads_category + Assembly_size:Reads_category, data=table3))
wt_m1 <- 1 / lm(abs(m1$residuals) ~ m1$fitted.values)$fitted.values^2
wls_m1 <- lm(Cvalue ~ Assembly_size + Reads_category + Assembly_size:Reads_category, data = table3, weights=wt_m1)
summary(wls_m1)

######################
Coefficients:
    Estimate   Std. Error t value  Pr(>|t|)
(Intercept)                    -1.655e+07  1.281e+07  -1.292   0.1971
Assembly_size                   1.214e+00  2.159e-02  56.222   <2e-16 ***
Reads_categorySR                5.468e+07  2.064e+07   2.649   0.0084 **
Assembly_size:Reads_categorySR -5.914e-03  3.246e-02  -0.182   0.8555
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.441 on 384 degrees of freedom
Multiple R-squared:  0.9368, Adjusted R-squared:  0.9363
F-statistic:  1898 on 3 and 384 DF,  p-value: < 2.2e-16
####################
```
