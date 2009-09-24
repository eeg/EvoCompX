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

title("final state")

# mean phenotype
plot(x, x*opt.slope, type="l", xlab="spatial position", ylab="mean phenotype")
lines(x, dat$zbar1, lwd=3, col=mycol[1])
lines(x, dat$zbar2, lwd=3, col=mycol[2])
lines(x, dat$zbar3, lwd=3, col=mycol[3])


#--------------------------------------------------
# Changes over time
#-------------------------------------------------- 

# swap rows and columns while reading in data
# so each time slice is a column and each row is a cell
num1 = t(read.table("num1.dat"))
num2 = t(read.table("num2.dat"))
num3 = t(read.table("num3.dat"))
zbar1 = t(read.table("zbar1.dat"))
zbar2 = t(read.table("zbar2.dat"))
zbar3 = t(read.table("zbar3.dat"))

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

title("changes over time")

# plot mean phenotypes
plot(x, x*opt.slope, type="l", xlab="spatial position", ylab="mean phenotype")
lines(x, zbar1[,t], lwd=3, col=mycol[1])
matlines(x, zbar1[], col=mycol[1], lty=1)
lines(x, zbar2[,t], lwd=3, col=mycol[2])
matlines(x, zbar2[], col=mycol[2], lty=1)
lines(x, zbar3[,t], lwd=3, col=mycol[3])
matlines(x, zbar3[], col=mycol[3], lty=1)


dev.off()
