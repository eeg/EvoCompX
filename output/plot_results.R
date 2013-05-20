myyellow = rgb(0.750, 0.750, 0.188) #BFBF30
myorange = rgb(0.832, 0.367, 0)     #D45E00
mygreen = rgb(0.161, 0.718, 0.008)  #29B702
myblue = rgb(0, 0.445, 0.695)       #0071B1
mypurple = rgb(0.443, 0.035, 0.667) #7109AA
mycol = c(myyellow, myorange, mygreen, myblue, mypurple)

#--------------------------------------------------
# mv ../input/*.dat .
# or make a simlink
#-------------------------------------------------- 

# number of species
nsp = 5

# slope of the optimum phenotype function
opt.slope = 0.5


# open an output file
pdf("results.pdf", width=8.5, height=11)
par(mfrow=c(2,1))


#--------------------------------------------------
# The final state of the system
#-------------------------------------------------- 

dat.num = read.table("num_final.dat")
dat.zbar = read.table("zbar_final.dat")
dat = cbind(dat.num, dat.zbar)
names(dat) = paste(append(rep("num", nsp), rep("zbar", nsp)), 
                   rep(as.character(seq(nsp)), 2), sep="")

x = seq(length(dat[,1]))

# number of individuals
plot(x, dat$num1, type="l", lwd=3, col=mycol[1], xlab="spatial position",
        ylab="number of individuals")
lines(x, dat$num2, lwd=3, col=mycol[2])
lines(x, dat$num3, lwd=3, col=mycol[3])
lines(x, dat$num4, lwd=3, col=mycol[4])
lines(x, dat$num5, lwd=3, col=mycol[5])

title("final state")

# mean phenotype
plot(x, x*opt.slope, type="l", xlab="spatial position", ylab="mean phenotype")
lines(x, dat$zbar1, lwd=3, col=mycol[1])
lines(x, dat$zbar2, lwd=3, col=mycol[2])
lines(x, dat$zbar3, lwd=3, col=mycol[3])
lines(x, dat$zbar4, lwd=3, col=mycol[4])
lines(x, dat$zbar5, lwd=3, col=mycol[5])


#--------------------------------------------------
# Changes over time
#-------------------------------------------------- 

# swap rows and columns while reading in data
# so each time slice is a column and each row is a cell
num1 = t(read.table("num1.dat"))
num2 = t(read.table("num2.dat"))
num3 = t(read.table("num3.dat"))
num4 = t(read.table("num4.dat"))
num5 = t(read.table("num5.dat"))
zbar1 = t(read.table("zbar1.dat"))
zbar2 = t(read.table("zbar2.dat"))
zbar3 = t(read.table("zbar3.dat"))
zbar4 = t(read.table("zbar4.dat"))
zbar5 = t(read.table("zbar5.dat"))

t = length(num1[1,])
x = seq(length(num1[,1]))

# plot abundances
plot(x, num1[,t], type="l", lwd=3, col=mycol[1], xlab="spatial position",
        ylab="number of individuals")
matlines(x, num1[], col=mycol[1], lty=1)
lines(x, num2[,t], lwd=3, col=mycol[2])
matlines(x, num2[], col=mycol[2], lty=1)
lines(x, num3[,t], lwd=3, col=mycol[3])
matlines(x, num3[], col=mycol[3], lty=1)
lines(x, num4[,t], lwd=3, col=mycol[4])
matlines(x, num4[], col=mycol[4], lty=1)
lines(x, num5[,t], lwd=3, col=mycol[5])
matlines(x, num5[], col=mycol[5], lty=1)

title("changes over time")

# plot mean phenotypes
plot(x, x*opt.slope, type="l", xlab="spatial position", ylab="mean phenotype")
lines(x, zbar1[,t], lwd=3, col=mycol[1])
matlines(x, zbar1[], col=mycol[1], lty=1)
lines(x, zbar2[,t], lwd=3, col=mycol[2])
matlines(x, zbar2[], col=mycol[2], lty=1)
lines(x, zbar3[,t], lwd=3, col=mycol[3])
matlines(x, zbar3[], col=mycol[3], lty=1)
lines(x, zbar4[,t], lwd=3, col=mycol[4])
matlines(x, zbar4[], col=mycol[4], lty=1)
lines(x, zbar5[,t], lwd=3, col=mycol[5])
matlines(x, zbar5[], col=mycol[5], lty=1)


dev.off()




#--------------------------------------------------
# TEMP
#-------------------------------------------------- 

svg("turnover.svg", width=9, height=4)
par(mfrow=c(1,2))

dat.num = read.table("num_final.dat")
dat.zbar = read.table("zbar_final.dat")
dat = cbind(dat.num, dat.zbar)
names(dat) = paste(append(rep("num", nsp), rep("zbar", nsp)), 
                   rep(as.character(seq(nsp)), 2), sep="")

x = seq(length(dat[,1]))

# number of individuals
plot(x, dat$num1, type="l", lwd=3, col=mycol[1], xlab="space",
        ylab="abundance", ylim=c(0,9))
lines(x, dat$num2, lwd=3, col=mycol[2])
lines(x, dat$num3, lwd=3, col=mycol[3])
lines(x, dat$num4, lwd=3, col=mycol[4])
lines(x, dat$num5, lwd=3, col=mycol[5])

# mean phenotype
plot(x, x*opt.slope, type="l", xlab="space", ylab="mean phenotype")
lines(x, dat$zbar1, lwd=3, col=mycol[1])
lines(x, dat$zbar2, lwd=3, col=mycol[2])
lines(x, dat$zbar3, lwd=3, col=mycol[3])
lines(x, dat$zbar4, lwd=3, col=mycol[4])
lines(x, dat$zbar5, lwd=3, col=mycol[5])

dev.off()
