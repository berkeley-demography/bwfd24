# Pop projection lecture. May 2024 

# Read in a bunch of Lx's, Fx's and Kx's from a few countries
# They're 5-year age groups so let's say our projection step is 5 years
# (i.e., "wide" rather than "long" format
# I'm assuming the radix is 1. If not, I have to adjust everything by the radix. 
# The data I'm reading in have 10 age groups, so 10 rows. I'm too lazy to generalize that.
# I should've looked at the length of each column to determine how many
# age groups I have and then created a nxn matrix for L but, you know. 

# This is in base R, not tidyverse, but I'm too lazy to re-do everything

options(digits=4)
dat =read.csv("projection-rates.csv")
head(dat)

# Define some functions to use 
# I wasn't smart enough to generalize the five-year groups to n-year groups
# I also wasn't smart enough to allow for a radix != 1.0

TFR <- function(f) return(5 * sum(f))       # This assumes 5-year rates
GRR <- function(f, ffab=0.4886) return(5 * sum(f) * ffab)
NRR <- function(l,f,ffab=0.4886) return( ffab*sum(l*f))
mu <- function(l,f, n) {
   x=seq_along(l)*n - (n/2)
   return( sum(x*l*f)/sum(l*f)) }


leslie <- function(l,f, ffab=0.4886) {
  n = length(l)-1
  L = matrix(0,n, n)
  L[1,] = ffab * l[1]*(f[1:n]+f[2:(n+1)]*l[2:(n+1)]/l[1:n])/2  # top row
  diag(L[2:n,1:(n-1)]) = l[2:n] / l[1:(n-1)]  # subdiagonal
return(L)
}

# Let's remind ourselves of what the input data look like

dat

# If you look at the columns, you'll see that we have the 1934 birth 
# cohort of US women's Lx and Fx values. We also have:
# Canada: low mortality
# Kenya: high mortality
# Brazil: low fertility
# Nigeria: high fertility
# Russia: "old" age structure
# Senegal: "young" age structure


# Examples of functions   =========

attach(dat)

TFR(Brazil.Fx) 
TFR(Niger.Fx)
TFR(USA1934.Fx)
plot(x,Niger.Fx,type="l",ylab="5Fx")
lines(x,USA1934.Fx,col="red")
lines(x,Brazil.Fx,col="blue")

NRR(USA1934.Lx,USA1934.Fx, 0.4877) # in this case, we know the Ffab = 0.4877
mu(USA1934.Lx,USA1934.Fx, 5)

NRR(Canada.Lx,Brazil.Fx)
mu(Kenya.Lx,Niger.Fx, 5)
NRR(Canada.Lx, Niger.Fx)


# Let's create some Leslie matrices from
# "hybrid" populations with low mortality and low fertility
# and with high mortality and high fertility


CB = leslie(Canada.Lx,Brazil.Fx)  # We can create Leslie matrices using a hybrid
CB   # You can see all the structural zeroes, the top row elements, and the subdiag elements

UU= leslie(USA1934.Lx,USA1934.Fx, 0.4877) #These are the US rates with US Ffab
UU

# ================ Projection code here, using UU and initial Russian Kx
A = UU                       # For this projection, use the US 1934 rates

K = matrix(0,nr=10,nc=51)    # Create a big empty matrix to hold our projection results

K[1:10,1] = Russia.Kx[1:10]  # Here's where we initialize the population in year 0.
# K[1:10,1] = c(0,0,0,100,100, 0, 0, 0, 0, 0)
K[,2] = A %*% K[,1]  # Project the initial pop 5 years into the future
K[,3] = A %*% K[,2]  # And another 5 years

K[,4] = A %*% K[,3]  # And again

K[,1:10]

# I'm tired of doing that one step at a time. Let's do all 50 steps in a loop

for(i in 1:50) K[,i+1] = A %*% K[,i]  # Here's where we do the entire projection
K                            # That was quick

totpop = colSums(K)  # Let's sum up the age groups to get a total pop for each 5-year step
totpop


year=seq(0,250,by=5)
plot(year,totpop)
plot(year, totpop, log="y")

# That looks pretty linear, especially at the right hand side
# Let's calculate the growth rate for the total pop between each 5-year step
# Take the log(P(t+5)/P(t))/5

log(totpop[2]/totpop[1]) /5
log(totpop[3]/totpop[2]) /5

# that's boring. Let's calculate all of them all at once. 

# notice: the log of a quotient is the difference of the logs, so 
# a faster way to calculate all of the growth rates is to take the logs,
# subtract each one from the next (that is, get the difference) and then
# divide that by n. In this case, n=5

r = diff(log(totpop))/5
r

plot(r)
plot(r, type="l") # sometimes it's easier to see if we connect the dots


# We can plot the pop pyramid for any given projection step
barplot(K[,1], horiz=T)  # This is the original population
barplot(K[,2], horiz=T)  # This is the pop after 5 years (1 projection step)

# That's boring. Let's plot them all, one by one
devAskNewPage(T); for(i in 1:51) barplot(K[,i], horiz=T, main=i)

# So, it looks like the population is growing but the pop pyramid is staying
# the same. 
# We can use either sweep() or scale() to get the percentage distribution
# That gives us the population pyramid

pyr <- sweep(K, 2, totpop, "/")  # or: scale(K, center=F, scale=totpop)

# Let's save those projection results. The projection used US Lx, US Fx, 
# and the Russian age structure, so let's call that UUR
# If we do a projection with Canadian Lx, Brazilian Fx, and Senegalese
# age structure, we'll call that CBS.


K.UUR <- K
totpop.UUR <- totpop
r.UUR <- r
pyr.UUR <- pyr


# We can do the same thing for the same A = UU but Senegal as the initial age structure
A = UU
K[,1] = Senegal.Kx[1:10]              # Initialize
for(i in 1:50) K[,i+1] = A %*% K[,i]  # Project
totpop = colSums(K)                   # Get total pop for each 5-year step
r = diff(log(totpop))/5               # calc 5-year growth rates
pyr <- sweep(K, 2, totpop, "/")  # or: scale(K, center=F, scale=totpop)

K.UUS <- K
totpop.UUS <- totpop
r.UUS <- r
pyr.UUS <- pyr

# Hmmmm. Compare the ending r and pyr with UUK. 

# Let's try an odd initial pop
K[,1] = c(0,1,2,3,4,5,6,7,8,9)
for(i in 1:50) K[,i+1] = A %*% K[,i]  
devAskNewPage(T); for(i in 1:51) barplot(K[,i], horiz=T, main=i)

# Two more odd initial pops
K[,1] = round(100*runif(10))
for(i in 1:50) K[,i+1] = A %*% K[,i]  
for(i in 1:51) barplot(K[,i], horiz=T, main=i)

K[,1] = c(0, 0, 0, 0, 0, 25, 25, 0, 0, 0)
for(i in 1:50) K[,i+1] = A %*% K[,i]  
for(i in 1:51) barplot(K[,i], horiz=T, main=i)


# Let's look at another Leslie matrix
A=CB
K[,1] = Senegal.Kx[1:10]
for(i in 1:50) K[,i+1] = A %*% K[,i]  
for(i in 1:51) barplot(K[,i], horiz=T, main=i)

totpop = colSums(K)                   # Get total pop for each 5-year step
r = diff(log(totpop))/5               # calc 5-year growth rates
pyr <- sweep(K, 2, totpop, "/")  # or: scale(K, center=F, scale=totpop)

K.CBS <- K
totpop.CBS <- totpop
r.CBS <- r
pyr.CBS <- pyr

# What do you think it looks like if we started with a different 
# initial pop?


# ================ Ignore the part below this for now

devAskNewPage(F)

##### Reproductive  value
A=UU
K[,1]=rep(1,10)
for(i in 1:50) K[,i+1] = A %*% K[,i]
z0 <- sum(K[,51])
z0

x=matrix(0,10,10); 
diag(x) = 1
y=x+1
y
z=rep(0,10)
for(j in 1:10) {K[,1]=y[,j]; for(i in 1:50) K[,i+1] = A %*% K[,i]; z[j]=sum(K[,51]) }
z
plot(z-z0, main="Additional pop at step 50 of 1 additional person at different ages")

##### Momentum & Keyfitz scenario
NRR.USA= NRR(USA1934.Lx, USA1934.Fx, 0.4877)
UU2 = leslie(USA1934.Lx, USA1934.Fx/NRR.USA, 0.4877)
A = UU; for(i in 1:15) K[,i+1] = A %*% K[,i]
A = UU2; for(i in 16:50) K[,i+1] = A %*% K[,i]
plot(colSums(K), type="l")
abline(v=16, lty=3)
 
 
####### Eigendecomposition
eigen.CB = eigen(CB)
r.CB = log(Re(eigen.CB$value[1]))/5
v = Re(eigen.CB$vectors[,1])
pyr.CB = v/sum(v)

 