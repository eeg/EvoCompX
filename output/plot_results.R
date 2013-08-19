myyellow <- rgb(0.750, 0.750, 0.188) #BFBF30
myorange <- rgb(0.832, 0.367, 0)     #D45E00
mygreen <- rgb(0.161, 0.718, 0.008)  #29B702
myblue <- rgb(0, 0.445, 0.695)       #0071B1
mypurple <- rgb(0.443, 0.035, 0.667) #7109AA
# mycol <- c(myyellow, myorange, mygreen, myblue, mypurple)
mycol <- c(myorange, mygreen, myblue, mypurple)

# mv ../input/*.dat .
# or make a simlink
# setwd("../test/03-dispersal_first/")
# setwd("/home/emma/plastic-cd/EvoCompX/input/")

# slopes of the environmental and optimum phenotype functions
env.slope <- 1
opt.slope <- 0.5
bbar <- 0.1
abar.slope <- env.slope * (opt.slope - bbar)

# number of species
nsp <- 2

# (could get above from params file instead)

# open an output file
pdf("results.pdf", width=7, height=14)
par(mfrow=c(3,1))

#--------------------------------------------------
# Changes over time
#-------------------------------------------------- 

num <- list()
zbar <- list()
abar <- list()
for (i in seq(nsp))
{
    num[[i]] <- t(read.table(paste("num", i, ".dat", sep="")))
    zbar[[i]] <- t(read.table(paste("zbar", i, ".dat", sep="")))
    abar[[i]] <- t(read.table(paste("abar", i, ".dat", sep="")))

    zbar[[i]][which(zbar[[i]] < -999)] <- NA
    abar[[i]][which(abar[[i]] < -999)] <- NA
}
# swap rows and columns while reading in data
# so each time slice is a column and each row is a cell

tm <- ncol(num[[1]])
x <- seq(nrow(num[[1]])) - 1
lty <- c(2, rep(1, tm-1))

# plot abundances
plot(NA, ylim=c(0, max(sapply(num, max))), xlim=c(0, max(x)),
     xlab="spatial position", ylab="number of individuals")
title("changes over time")
for (i in seq(nsp))
{
    lines(x, num[[i]][,tm], lwd=3, col=mycol[i])
    matlines(x, num[[i]][], col=mycol[i], lty=lty)
}

# plot mean phenotypes
plot(x, x*opt.slope, type="l", xlab="spatial position", ylab="mean phenotype")
for (i in seq(nsp))
{
    lines(x, zbar[[i]][,tm], lwd=3, col=mycol[i])
    matlines(x, zbar[[i]][], col=mycol[i], lty=lty)
}

# plot breeding values
plot(x, x*abar.slope, type="l", xlab="spatial position", ylab="mean breeding value")
for (i in seq(nsp))
{
    lines(x, abar[[i]][,tm], lwd=3, col=mycol[i])
    matlines(x, abar[[i]][], col=mycol[i], lty=lty)
}

#--------------------------------------------------
# The final state of the system
#-------------------------------------------------- 

dat.num <- read.table("num_final.dat")
dat.zbar <- read.table("zbar_final.dat")
dat.abar <- read.table("abar_final.dat")

dat <- cbind(dat.num, dat.zbar, dat.abar)
names(dat) <- paste(c(rep("num", nsp), rep("zbar", nsp), rep("abar", nsp)),
                    rep(as.character(seq(nsp)), 2), sep="")
for (i in (nsp+1):ncol(dat))
{
    dat[which(dat[,i] < -999), i] <- NA
}

x <- seq(nrow(dat)) - 1

# number of individuals
plot(x, dat[["num1"]], type="l", lwd=3, col=mycol[1], xlab="spatial position",
        ylab="number of individuals")
if (nsp > 1)
{
    for (i in 2:nsp)
    {
        lines(x, dat[[paste("num", i, sep="")]], lwd=3, col=mycol[i])
    }
}
title("final state")

# mean phenotype
plot(x, x*env.slope*opt.slope, type="l", xlab="spatial position", ylab="mean phenotype")
for (i in 1:nsp)
{
    lines(x, dat[[paste("zbar", i, sep="")]], lwd=3, col=mycol[i])
}

# mean breeding value
plot(x, x*abar.slope, type="l", xlab="spatial position", ylab="mean breeding value")
for (i in 1:nsp)
{
    lines(x, dat[[paste("abar", i, sep="")]], lwd=3, col=mycol[i])
}

#--------------------------------------------------
# Asymmetry of displacement
#--------------------------------------------------

# These are a bit confusing.  There must be a clearer way to look at asymmetry.
# Also, there must be a way that doesn't assume the optimum is known.  
#
# Probably need to compare sympatry with allopatry.  But comparing slopes
# assumes that the gradient is linear.
#
# Or maybe something with abundances?  Work out expected amount of displacement
# for given abundance difference, and look for deviations from that?
# Old standby instead is dominance-divergence plot.

# Dominance-divergence plot
plot.dom.div <- function(n1, n2, z1, z2, ylab.type)
{
    dom <- log(n1/n2)
    xlab <- "ln(num1/num2)"

    if (sum(z1) > sum(z2))
    {
        div <- z1 - z2
        ylab <- sprintf("%s%d - %s%d", ylab.type, 1, ylab.type, 2)
    } else {
        div <- z2 - z1
        ylab <- sprintf("%s%d - %s%d", ylab.type, 2, ylab.type, 1)
    }

    x0 <- which.min(div)
    x1 <- min(abs(c(1,100) - x0))
    i <- seq(x0-x1, x0+x1)

    plot(dom[i], div[i], xlab=xlab, ylab=ylab, type="l", lwd=3, main="dominance-divergence")
    abline(v=0, lty=2)
}

if (nsp == 2)
{
    x <- seq(nrow(dat)) - 1
    x0 <- which.min(abs(dat.num[,1] - dat.num[,2]))

    # mean phenotype

    opt.z <- x * env.slope * opt.slope
    z1 <- dat.zbar[,1]
    z2 <- dat.zbar[,2]

    plot(x, abs((opt.z - z1) - (z2 - opt.z)), type="l", lwd=3, main="mean phenotype diff")

    plot(x, opt.z - z1, type="l", lwd=3, col=mycol[1], main="mean phenotype clines")
    lines(x+2*x0-100, rev(z2 - opt.z), lwd=3, col=mycol[2])
    abline(h=0)

    plot.dom.div(dat.num[,1], dat.num[,2], z1, z2, "zbar")

    # mean breeding value

    opt.a <- x * abar.slope
    a1 <- dat.abar[,1]
    a2 <- dat.abar[,2]

    plot(x, abs((opt.a - a1) - (a2 - opt.a)), type="l", lwd=3, main="mean breeding value diff")

    plot(x, opt.a - a1, type="l", lwd=3, col=mycol[1], main="mean breeding value clines")
    lines(x+2*x0-100, rev(a2 - opt.a), lwd=3, col=mycol[2])
    abline(h=0)

    plot.dom.div(dat.num[,1], dat.num[,2], a1, a2, "abar")
}


junk <- dev.off()  # suppress useless message
