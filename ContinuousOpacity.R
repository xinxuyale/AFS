# Continuous Opacity Correction Example

# Load in the example dataset. It contains five neighboring orders.
load("~/Desktop/AFS/examples/ContinuousOpacity.RData")

# Set tilde_q to be 0.9.
tilde_q <- .9
n1 <- dim(x1)[1]
n2 <- dim(x2)[1]
n3 <- dim(x3)[1]
n4 <- dim(x4)[1]
n5 <- dim(x5)[1]

# Save the indices of points whose intensities are in the top tilde_q quantile.
ind1 <- which(x1$intens>quantile(x1$intens, tilde_q))
ind2 <- which(x2$intens>quantile(x2$intens, tilde_q))
ind3 <- which(x3$intens>quantile(x3$intens, tilde_q))
ind4 <- which(x4$intens>quantile(x4$intens, tilde_q))
ind5 <- which(x5$intens>quantile(x5$intens, tilde_q))

# Record the start and end indices of each overlap region.
start1 <- which(x1$wv>=x2$wv[1])[1]
end1 <- which(x2$wv>x1$wv[n1])[1]-1
start2 <- which(x2$wv>=x3$wv[1])[1]
end2 <- which(x3$wv>x2$wv[n2])[1]-1
start3 <- which(x3$wv>=x4$wv[1])[1]
end3 <- which(x4$wv>x3$wv[n3])[1]-1
start4 <- which(x4$wv>=x5$wv[1])[1]
end4 <- which(x5$wv>x4$wv[n4])[1]-1

# Linearly adjust x2 to align with x1.
# Use the points whose intensities are in the top tilde_q quantile in the 
# overlap region to estimate the difference in slope and intercept.
temp1 <- x1$intens[start1:n1]
temp2 <- x2$intens[1:end1]
index1 <- which(temp1>=quantile(temp1, tilde_q))
index2 <- which(temp2>=quantile(temp2, tilde_q))
lm1 <- lm(temp1[index1] ~ x1$wv[start1:n1][index1])
lm2 <- lm(temp2[index2] ~ x2$wv[1:end1][index2])
offset <- lm2$coefficients[1] - lm1$coefficients[1]
slope <- lm2$coefficients[2] - lm1$coefficients[2]
x2$intens2 <- x2$intens - x2$wv*slope - offset  

# Linearly adjust x3 to align with x2.
# Use the points whose intensities are in the top tilde_q quantile in the 
# overlap region to estimate the difference in slope and intercept.
temp1 <- x2$intens2[start2:n2]
temp2 <- x3$intens[1:end2]
index1 <- which(temp1>=quantile(temp1, tilde_q))
index2 <- which(temp2>=quantile(temp2, tilde_q))
lm1 <- lm(temp1[index1] ~ x2$wv[start2:n2][index1])
lm2 <- lm(temp2[index2] ~ x3$wv[1:end2][index2])
offset <- lm2$coefficients[1] - lm1$coefficients[1]
slope <- lm2$coefficients[2] - lm1$coefficients[2]
x3$intens2 <- x3$intens-x3$wv*slope - offset  

# Linearly adjust x4 to align with x3.
# Use the points whose intensities are in the top tilde_q quantile in the 
# overlap region to estimate the difference in slope and intercept.
temp1 <- x3$intens2[start3:n3]
temp2 <- x4$intens[1:end3]
index1 <- which(temp1>=quantile(temp1, tilde_q))
index2 <- which(temp2>=quantile(temp2, tilde_q))
lm1 <- lm(temp1[index1] ~ x3$wv[start3:n3][index1])
lm2 <- lm(temp2[index2] ~ x4$wv[1:end3][index2])
offset <- lm2$coefficients[1] - lm1$coefficients[1]
slope <- lm2$coefficients[2] - lm1$coefficients[2]
x4$intens2 <- x4$intens - x4$wv*slope - offset  

# Linearly adjust x5 to align with x4.
# Use the points whose intensities are in the top tilde_q quantile in the 
# overlap region to estimate the difference in slope and intercept.
temp1 <- x4$intens2[start4:n4]
temp2 <- x5$intens[1:end4]
index1 <- which(temp1>=quantile(temp1, tilde_q))
index2 <- which(temp2>=quantile(temp2, tilde_q))
lm1 <- lm(temp1[index1] ~ x4$wv[start4:n4][index1])
lm2 <- lm(temp2[index2] ~ x5$wv[1:end4][index2])
offset <- lm2$coefficients[1] - lm1$coefficients[1]
slope <- lm2$coefficients[2] - lm1$coefficients[2]
x5$intens2 <- x5$intens - x5$wv*slope - offset 

# Fit another linear regression to get rid of the extra slope.
# The variable whole records the combined spectrum.
# Use the points saved in ind1 to ind5 to fit a linear model.
whole <- c(x1$intens, x2$intens2, x3$intens2, x4$intens2, x5$intens2)
index <- c(ind1, ind2+n1, ind3+n1+n2, ind4+n1+n2+n3, ind5+n1+n2+n3+n4)
wave <- c(x1$wv, x2$wv, x3$wv, x4$wv, x5$wv)
lm_whole <- lm(whole[index] ~ wave[index])
whole <- whole - wave*lm_whole$coefficients[2] 
whole <- whole - whole[1] + x1$intens[1]
# Save the corrected intensities to each order.
x1$intens2 <- whole[(1):(n1)]
x2$intens2 <- whole[(n1+1):(n1+n2)]
x3$intens2 <- whole[(n1+n2+1):(n1+n2+n3)]
x4$intens2 <- whole[(n1+n2+n3+1):(n1+n2+n3+n4)]
x5$intens2 <- whole[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]

# Show the comparison.
plot(x1$wv, x1$intens, type="l", xlim=c(4540, 4740))
lines(x2$wv, x2$intens)
lines(x3$wv, x3$intens)
lines(x4$wv, x4$intens)
lines(x5$wv, x5$intens)
abline(h=1, col="red")
lines(x1$wv, x1$intens2, col=alpha("green", alpha=0.5))
lines(x2$wv, x2$intens2, col=alpha("green", alpha=0.5))
lines(x3$wv, x3$intens2, col=alpha("green", alpha=0.5))
lines(x4$wv, x4$intens2, col=alpha("green", alpha=0.5))
lines(x5$wv, x5$intens2, col=alpha("green", alpha=0.5))




