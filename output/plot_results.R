myyellow <- rgb(0.750, 0.750, 0.188) #BFBF30
myorange <- rgb(0.832, 0.367, 0)     #D45E00
mygreen <- rgb(0.161, 0.718, 0.008)  #29B702
myblue <- rgb(0, 0.445, 0.695)       #0071B1
mypurple <- rgb(0.443, 0.035, 0.667) #7109AA
mycol <- c(myorange, mygreen, myblue, mypurple)

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
