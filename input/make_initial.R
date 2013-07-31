#--------------------------------------------------
# Modeled on make_initial.py
#--------------------------------------------------

library(stringr)

get.best.abar <- function(x, sp)
{
    env.slope <- 1  # slope of the environment function, C
    (params$opt_slope - params$bbar[sp]) * env.slope * (x - 1)
}

undef <- -9999 # this must match the value in input.h

parfile <- (commandArgs(TRUE))
tmp <- read.table(parfile, sep="=", as.is=T)
params <- str_split(str_trim(tmp$V2), "\\s+")
names(params) <- str_trim(tmp$V1)
i <- which(!names(params) %in% c("initial_num", "initial_abar"))
params[i] <- lapply(params[i], as.numeric)

for (varname in c("V_u", "delta", "bbar"))
{
    if (length(params[varname] == 1))
        params[[varname]] <- rep(params[[varname]], params$num_species)
}

sp.init <- matrix(NA, nrow=params$num_species, ncol=4)
colnames(sp.init) <- c("x.start", "x.stop", "abun", "abar.off")
#                              sp1  sp2
sp.init[, "x.start"]     <- c(   1,  91)
sp.init[, "x.stop"]      <- c(  10, 100)
sp.init[, "abun"]        <- c(   4,   6)
sp.init[, "abar.off"]    <- c(-0.5, 0.3)

num.in  <- matrix(0, nrow=params$space_size, ncol=params$num_species)
abar.in <- matrix(undef, nrow=params$space_size, ncol=params$num_species)
for (i in seq(nrow(sp.init)))
{
    idx <- sp.init[i, "x.start"] : sp.init[i, "x.stop"]
    num.in[idx, i]  <- sp.init[i, "abun"]
    abar.in[idx, i] <- get.best.abar(idx, i) + sp.init[i, "abar.off"]
}

write.table( format(num.in, digits=8), file= "num.in", row.names=F,
            col.names=F, sep="\t", quote=F)
write.table(format(abar.in, digits=8), file="abar.in", row.names=F,
            col.names=F, sep="\t", quote=F)
