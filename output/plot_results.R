myorange = rgb(0.832, 0.367, 0)    # D45E00
myblue = rgb(0, 0.445, 0.695)      # 0071B1
mygreen = rgb(0.161, 0.718, 0.008) # 29B702
mycol = c(myorange, myblue, mygreen)

#--------------------------------------------------
# mv ../input/*.dat .
# or make a simlink
#-------------------------------------------------- 

# number of species
nsp = 3

# slope of the optimum phenotype function
opt.slope = 0.5


#--------------------------------------------------
# The final state of the system
#-------------------------------------------------- 

dat = read.table("final.dat")
names(dat) = paste(append(rep("num", nsp), rep("zbar", nsp)), 
                   rep(as.character(seq(nsp)), 2), sep="")

x = seq(length(dat[,1]))

# number of individuals
plot(x, dat$num1, type="l", lwd=3, col=mycol[1], ylab="number of individuals")
lines(x, dat$num2, lwd=3, col=mycol[2])
lines(x, dat$num3, lwd=3, col=mycol[3])

# mean phenotype
plot(x, x*opt.slope, type="l", ylab="mean phenotype")
lines(x, dat$zbar1, lwd=3, col=mycol[1])
lines(x, dat$zbar2, lwd=3, col=mycol[2])
lines(x, dat$zbar3, lwd=3, col=mycol[3])


#--------------------------------------------------
# Changes over time
#-------------------------------------------------- 

plotme = function(i, dat, col)
{
    lines(x, dat[i,], col=col)
}

num1 = read.table("num1.dat")
num2 = read.table("num2.dat")
num3 = read.table("num3.dat")
zbar1 = read.table("zbar1.dat")
zbar2 = read.table("zbar2.dat")
zbar3 = read.table("zbar3.dat")

t = length(num1[,1])
x = seq(length(num1[1,]))

plot(x, num1[t,], type="l", lwd=3, col=mycol[1], xlab="number of individuals")
lapply(seq(t-1), plotme, num1, mycol[1])
lines(x, num2[t,], lwd=3, col=mycol[2])
lapply(seq(t-1), plotme, num2, mycol[2])
lines(x, num3[t,], lwd=3, col=mycol[3])
lapply(seq(t-1), plotme, num3, mycol[3])

plot(x, x*opt.slope, type="l", ylab="mean phenotype")
lines(x, zbar1[t,], lwd=3, col=mycol[1])
lapply(seq(t-1), plotme, zbar1, mycol[1])
lines(x, zbar2[t,], lwd=3, col=mycol[2])
lapply(seq(t-1), plotme, zbar2, mycol[2])
lines(x, zbar3[t,], lwd=3, col=mycol[3])
lapply(seq(t-1), plotme, zbar3, mycol[3])
