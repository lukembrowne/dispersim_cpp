
library(ggplot2)
library(dplyr)
### Basing equation off Harms et al. 2000

## Assuming intercept (alpha = 0) - is this a safe assumption??

## Output is predicted number of recruits of a specific species ID and genotype, given number of seeds and cndd and gndd values

## Bgi = strength of GNDD - slope of log log relationship: 1 = no NDD
## Bci = strength of CNDD

## Sgi = density of seeds of given species and genotype
## Sci = density of seeds of given species
harms <- function(bgi, bci, sgi, sci, verbose = FALSE){
  
  ## Effectively halves the strength of the slopes from Harms et al. 2000
 #R = (bgi * log(sgi) + bci * log(sci) + log(sgi/sci)) / 2
 
  # Original harms equation, but modified to return genotype densities
 R = bci * log(sci) + log(sgi/sci)
 
 ### If slope of a  line is still 1, then something is wrong because that's not NDD - rare species doesn't have a proportional advantage of common species, slope needs to be less than 1!! Need to have higher proportional loss at higher densities

 if(verbose == TRUE){
   cat("\n ----------------------- \n")
  cat("bgi =", bgi, ", bci = ", bci, " \n")
  cat("Original number of seeds per genotype:")
  print(sgi)
  cat("\n Remaining number of Recruits per genotype:")
  print(exp(R))
  cat("\n Percent loss of seeds:")
  print((sgi - exp(R))/sgi)
  cat("\n ----------------------- \n")
 }

 loss <- (sgi - exp(R))/sgi
 
 # return(R)
 return(loss)
 
}

## Generate data of number of genotypes for species XX
#genotypes <- as.numeric(table(sample.int(20, size = 500, replace = TRUE)))
genotypes <- 0:1000

harms(bgi = 1, bci = .8,
      sgi = 1:10, sci = sum(genotypes), verbose = TRUE)


## Checking percent loss of seeds based on gndd and cndd values
loss_test <- data.frame(mean_cndd = NA, mean_gndd = NA, percent_loss = NA)
row = 1
for(i in seq(-.2, .9, length.out = 100)){
for(j in 1){
  # Log sequence
  #for(i in exp(seq(log(0.001), log(1), length.out = 20))){
    #for(j in exp(seq(log(0.001), log(1), length.out = 20))){
  
 
    
    loss_test[row, "mean_cndd"] <- i
    loss_test[row, "mean_gndd"] <- j
    loss_test[row, "percent_loss"] <- harms(bgi = j, bci = i,
                                       sgi = 1:20, 
                                       sci = sum(genotypes))[1]
    row = row + 1
  }
}

# Plots of percent loss 
ggplot(loss_test, aes(x = mean_cndd, y = percent_loss)) +
  geom_line() + geom_point() + theme_bw(15)



####### Example plot of all genotypes within one plot

plot(0:7, 0:7, type = "n", ylab = "log(recruits)", xlab = "log(seeds)", las = 1)
abline(a = 0, b = 1)

points(log(genotypes), (harms(bgi = 1, bci = 1,
                                 sgi = genotypes, sci = sum(genotypes))), pch = 19)


points(log(genotypes), (harms(bgi = .5, bci = 1,
                                 sgi = genotypes, sci = sum(genotypes))), pch = 19,
       col = "red")

points(log(genotypes), (harms(bgi = 1, bci = 0.5,
                                 sgi = genotypes, sci = sum(genotypes))), pch = 19,
       col = "skyblue2")

points(log(genotypes), (harms(bgi = 0.5, bci = 0.5,
                                 sgi = genotypes, sci = sum(genotypes))), pch = 19,
       col = "orange")

legend("topleft", legend = c("bgi = 1.0, bci = 1.0", 
                             "bgi = 0.5, bci = 1.0", 
                             "bgi = 1.0, bci =  0.5",
                             "bgi = 0.5, bci =  0.5"),
       fill = c("black", "red", "skyblue2", "orange"), bty = "n")



### Write a function that generates data in plots for a single genotype and randomize this and plot the lines

sample_harms <- function(bgi, bci, n, max, verbose = FALSE){
  
  gen_seeds <- NA
  gen_recruits <- NA
  cons_seeds <- NA
  cons_recruits <- NA
  
  for(x in 1:n){
    gen_dens = sample(1:max, size = 1) # Sample density of seeds for one genotype
    if(gen_dens == max){
      cons_dens = max
      cons_seeds[x] <- log(cons_dens)
      cons_recruits[x] <- log(sum(exp(harms(bgi = bgi, bci = bci, 
                                            sgi = c(gen_dens), 
                                            sci = cons_dens,
                                            verbose = verbose))))
      
      next
      
    } else {
      cons_dens = sample(gen_dens:max, size = 1)
    }
    
     #cat("Gen_dens:", gen_dens, "| cons_dens: ", cons_dens, "\n")
    
    gen_seeds[x] <- log(gen_dens)
    gen_recruits[x] <- harms(bgi = bgi, bci = bci, 
                      sgi = gen_dens, 
                      sci = cons_dens,
                      verbose = verbose)
    cons_seeds[x] <- log(cons_dens)
    cons_recruits[x] <- log(sum(exp(harms(bgi = bgi, bci = bci, 
                           sgi = c(gen_dens, cons_dens - gen_dens), 
                           sci = cons_dens,
                           verbose = verbose))))
    
  }
    
    return(list(gen_seeds = gen_seeds, 
                gen_recruits = gen_recruits,
                cons_seeds = cons_seeds,
                cons_recruits = cons_recruits))
}


sum(exp(harms(bgi = bgi, bci = bci, 
                  sgi = c(gen_dens, 5), 
                  sci = cons_dens,
                  verbose = verbose)))



sample_harms(bgi = 0.5, bci = 0.5, n = 10, max = 1000 )

sample_harms(bgi = .98, bci = .98, n = 10, max = 1000, verbose = T )


plot(c(.0, .5, .75, .9, 1), c(0.9989236, 0.96375, 0.778174, .486988, 0 ),
     xlab = "cndd", ylab = "percent loss of seeds", las  = 1, 
     pch = 21, bg = "skyblue2")


####### Example plot
plot(0:6, 0:6, type = "n", ylab = "log(recruits)", xlab = "log(seeds)", las = 1)
abline(a = 0, b = 1)

n = 1000 ## Number of points
max = 100 # Max number of individuals

## One to one
scen_11 <- sample_harms(bgi = 1, bci = 1, n = n, max = max)
points(scen_11$gen_seeds, scen_11$gen_recruits, pch = 19)
(lm1 <- lm(gen_recruits ~ 0 + gen_seeds, data = scen_11))
abline(lm1, lwd = 2)

## 0.5 1.0
scen_051 <- sample_harms(bgi = 0.5, bci = 1, n = n, max = max)
points(scen_051$gen_seeds, scen_051$gen_recruits, pch = 19, col = "red")
(lm2 <- lm(gen_recruits ~ 0 + gen_seeds, data = scen_051))
abline(lm2, lwd = 2, col = "red")

## 1.0 0.5
scen_105 <- sample_harms(bgi = 1, bci = 0.5, n = n, max = max)
points(scen_105$gen_seeds, scen_105$gen_recruits, pch = 19, col = "skyblue2")
(lm3 <- lm(gen_recruits ~ 0 + gen_seeds, data = scen_105))
abline(lm3, lwd = 2, col = "skyblue2")


## .5 .5
scen_0505 <- sample_harms(bgi = 0.5, bci = 0.5, n = n, max = max)
points(scen_0505$gen_seeds, scen_0505$gen_recruits, pch = 19, col = "orange")
(lm4 <- lm(gen_recruits ~ 0 + gen_seeds, data = scen_0505))
abline(lm4, lwd = 2, col = "orange")

scen_0505$gen_recruits <- log(exp(scen_0505$gen_recruits) + 1)
scen_0505$gen_seeds <- log(exp(scen_0505$gen_seeds) + 1)

(lm4 <- lm(gen_recruits ~ 0 + gen_seeds, data = scen_0505))
(lm4 <- lm(gen_recruits ~ gen_seeds, data = scen_0505))
points(scen_0505$gen_seeds, scen_0505$gen_recruits, pch = 19, col = "orange")
abline(lm4, lwd = 2, col = "orange")



### Conspecific perspective
scen_0505 <- sample_harms(bgi = 1, bci = 1.0, n = n, max = max)
#points(scen_0505$cons_seeds, scen_0505$cons_recruits, pch = 19, col = "grey80")
(lm4 <- lm(cons_recruits ~ 0 + cons_seeds, data = scen_0505))
abline(lm4, lwd = 2, col = "grey80")



scen_0505 <- sample_harms(bgi = 1, bci = 0.5, n = n, max = max)
#points(scen_0505$cons_seeds, scen_0505$cons_recruits, pch = 19, col = "purple")
(lm4 <- lm(cons_recruits ~ 0 + cons_seeds, data = scen_0505))
abline(lm4, lwd = 2, col = "purple")

### Conspecific
scen_0505 <- sample_harms(bgi = 0.5, bci = 0.5, n = n, max = max)
#points(scen_0505$cons_seeds, scen_0505$cons_recruits, pch = 19, col = "brown")
(lm4 <- lm(cons_recruits ~ 0 + cons_seeds, data = scen_0505))
abline(lm4, lwd = 2, col = "brown")

scen_0505 <- sample_harms(bgi = 0.5, bci = 1.0, n = n, max = max)
#points(scen_0505$cons_seeds, scen_0505$cons_recruits, pch = 19, col = "darkblue")
(lm4 <- lm(cons_recruits ~ 0 + cons_seeds, data = scen_0505))
abline(lm4, lwd = 2, col = "darkblue")


legend("topleft", legend = c("bgi = 1.0, bci = 1.0", 
                             "bgi = 0.5, bci = 1.0", 
                             "bgi = 1.0, bci =  0.5",
                             "bgi = 0.5, bci =  0.5",
                             "bgi = 1.0, bci =  1.0 - species level",
                             "bgi = 1.0, bci =  0.5 - species level",
                             "bgi = 0.5, bci =  0.5 - species level",
                             "bgi = 0.5, bci =  1.0 - species level"),
       fill = c("black", "red", "skyblue2", "orange", "grey80", "purple", "brown",
                "darkblue"), bty = "n")

### Can force intercepts to go through 0 with ~ 0 + seeds
## But sometimes lines will have an intercept less than 0 because they are losing fractions of a seedling when there is 1 seedling for a genotype and there is conspecific density dependense

### Thoughts on graph..

# is cool and makes sense that bgi and bci have equivalent effects at genotype and species level

# All slopes are less than 1, which is good

# 0.5, 0.5 is strongest for both species and genotypes, which is good


## Not sure if i should present just straight lines for example graph instead of points..
# can just use abline(a = 0, b = ##) to draw the lines using sample_harms and lm() to make sure the slopes are right. Seems like slope is usually (bci + bgi) / 2

## Make two separate panels --
# 1) Plot of genotypes within one species
# 2) Plot of slope of species within a community, with varying levels of bci and bgi

# With ranges of bci and bgi from 1, 0.5, 0

## With Bgi = 1, Bci = -1 = overall slope of 0 (lowest range of harms and lamanna studies)
# but bgi -1, bci -1 would be overall -1 slope



#############################
### CNDD PLOT
  # CNDD Plot ---------------------------------------------------------------
  
  pdf("./figures/Figure 1 - NDD schematic.pdf", width = 12, height = 4)
  
  par(mfrow = c(1, 2), mar = c(5.1, 4.1, 2.1, 8.1))
  
  axis_vals <- c(.01, 1, 5, 10, 25, 50, 100)
  
  x_pos = 350 ## For labels
  max_size = 150
  n = 1000
  
  genotypes <- 1:100
  
  text_cex = .75

### AT SPECIES LEVEL  
  plot(seq(0, log(max_size), by = .1), seq(0, log(max_size), by = .1), type = "n", 
       ylab = expression('Recruit density (R'['i']*')'), xlab = expression('Seed density (S'['i']*')'), las = 1,
       bty = "l", xaxt = "n", yaxt = "n", yaxt = "n", xaxs="i", yaxs = "i")
  axis(side = 1, at = log(axis_vals), labels = axis_vals)
  axis(side = 2, at = log(axis_vals), labels = axis_vals, las = 1)
 # abline(a = 0, b = 1)
  
  mtext("a", side = 3, adj = -.1)
  
  #### GNDD = 1, CNDD = 1
  harms_out = sample_harms(bgi = 1, bci = 1, n = n, max = max_size)
 # points(harms_out$cons_seeds, harms_out$cons_recruits, pch = ".")
  lm1 <- lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1)
  abline(lm1, lty = 2)
  b_overall <- as.character(sprintf("%.2f", round(lm1$coefficients, 2)))
  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits),
       labels = bquote('b'['con']*' = 1.00, b'['gen']*' = 1.00'),
       xpd = NA , cex = text_cex)
  
  
  #### GNDD = .5, CNDD = 1
  harms_out = sample_harms(bgi = .5, bci = 1, n = n, max = max_size)
 # points(harms_out$cons_seeds, harms_out$cons_recruits, pch = ".")
  lm1 <- lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1)
  abline(lm1)
  b_overall <- as.character(sprintf("%.2f", round(lm1$coefficients, 2)))
  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits), 
       labels = bquote('b'['con']*' = 1.00, b'['gen']*' = 0.50'),
       xpd = NA, cex = text_cex)

  #### GNDD = 1, CNDD = .5
  harms_out = sample_harms(bgi = 1, bci = .5, n = n, max = max_size)
 # points(harms_out$cons_seeds, harms_out$cons_recruits, pch = ".")
  lm1 <- lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1)
  abline(lm1)
  b_overall <- as.character(sprintf("%.2f", round(lm1$coefficients, 2)))
  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits), 
       labels = bquote('b'['con']*' = 0.50, b'['gen']*' = 1.00'),
       xpd = NA, cex = text_cex)
  
  
  # #### GNDD = .5, CNDD = .5
 #  harms_out = sample_harms(bgi = .5, bci = .5, n = n, max = max_size)
 # # points(harms_out$cons_seeds, harms_out$cons_recruits, pch = ".")
 #  abline(lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1))
 #  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits), labels = expression('b'['con']*' = 0.5, b'['gen']*' = 0.5'), xpd = NA)

  
  #### GNDD = .1, CNDD = 1
  harms_out = sample_harms(bgi = .1, bci = 1.0, n = n, max = max_size)
  #points(harms_out$cons_seeds, harms_out$cons_recruits, pch = ".")
  lm1 <- lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1)
  abline(lm1)
  b_overall <- as.character(sprintf("%.2f", round(lm1$coefficients, 2)))
  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits), 
       labels = bquote('b'['con']*' = 1.00, b'['gen']*' = 0.10'),
       xpd = NA, cex = text_cex)

  
  #### GNDD = 1, CNDD = 0.1
  harms_out = sample_harms(bgi = 1, bci = .1, n = n, max = max_size)
 # points(harms_out$cons_seeds, harms_out$cons_recruits, pch = ".")
  lm1 <- lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1)
  abline(lm1)
  b_overall <- as.character(sprintf("%.2f", round(lm1$coefficients, 2)))
  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits), 
       labels = bquote('b'['con']*' = 0.10, b'['gen']*' = 1.00'),
       xpd = NA, cex = text_cex)
  
  # #### GNDD = .1, CNDD = .1
  harms_out = sample_harms(bgi = .1, bci = .1, n = n, max = max_size)
 # points(harms_out$cons_seeds, harms_out$cons_recruits, pch = 19)
  lm1 <- lm(harms_out$cons_recruits ~ harms_out$cons_seeds - 1)
  abline(lm1)
  b_overall <- as.character(sprintf("%.2f", round(lm1$coefficients, 2)))
  text(x = c(log(x_pos)), y = max(harms_out$cons_recruits), 
     #  labels = bquote('b'['con']*' = 0.10, b'['gen']*' = 0.10'*', b'['total']*' = '*.(b_overall)),
     labels = bquote('b'['con']*' = 0.10, b'['gen']*' = 0.10'),
       xpd = NA, cex = text_cex)
  
  
  
##### AT GENOTYPE LEVEL
  
 
  plot(seq(0, log(max_size), by = .1), seq(0, log(max_size), by = .1), type = "n", 
       ylab = expression('Recruit density (R'['gi']*')'), xlab = expression('Seed density (S'['gi']*')'), las = 1,
       bty = "l", xaxt = "n", yaxt = "n", xaxs="i", yaxs = "i")
  axis(side = 1, at = log(axis_vals), labels = axis_vals)
  axis(side = 2, at = log(axis_vals), labels = axis_vals, las = 1)
  # abline(a = 0, b = 1)
  mtext("b", side = 3, adj = -.1)
  
  
  #### GNDD = 1, CNDD = 1
  harms_out <- (harms(bgi = 1, bci = 1,
                      sgi = genotypes, sci = max(genotypes)))
  lines(log(genotypes), harms_out, lty = 2)
  text(x = c(log(x_pos)), y = max(harms_out),
       labels = bquote('b'['con']*' = 1.00, b'['gen']*' = 1.00'),
       xpd = NA , cex = text_cex)
  
  #### GNDD = .5, CNDD = 1
  harms_out <- (harms(bgi = .5, bci = 1,
                      sgi = genotypes, sci = max(genotypes)))
  lines(log(genotypes), harms_out, lty = 1)
  text(x = c(log(x_pos)), y = max(harms_out),
       labels = bquote('b'['con']*' = 1.00, b'['gen']*' = 0.50'),
       xpd = NA , cex = text_cex)
  
  #### GNDD = 1, CNDD = .5
  harms_out <- (harms(bgi = 1, bci = 0.5,
                      sgi = genotypes, sci = max(genotypes)))
  lines(log(genotypes), harms_out, lty = 1)
  text(x = c(log(x_pos)), y = max(harms_out),
       labels = bquote('b'['con']*' = 0.50, b'['gen']*' = 1.00'),
       xpd = NA , cex = text_cex)
  
  
  #### GNDD = .1, CNDD = 1
  harms_out <- (harms(bgi = .1, bci = 1,
                      sgi = genotypes, sci = max(genotypes)))
  lines(log(genotypes), harms_out, lty = 1)
  text(x = c(log(x_pos)), y = max(harms_out),
       labels = bquote('b'['con']*' = 1.00, b'['gen']*' = 0.10'),
       xpd = NA , cex = text_cex)
  
  
  #### GNDD = 1, CNDD = 0.1
  harms_out <- (harms(bgi = 1, bci = 0.1,
                      sgi = genotypes, sci = max(genotypes), verbose = T))
  lines(log(genotypes), harms_out, lty = 1)
  text(x = c(log(x_pos)), y = max(harms_out),
       labels = bquote('b'['con']*' = 0.10, b'['gen']*' = 1.00'),
       xpd = NA , cex = text_cex)
  
  
  # #### GNDD = .1, CNDD = .1
  harms_out <- (harms(bgi = .1, bci = 0.1,
                      sgi = genotypes, sci = max(genotypes), verbose = TRUE))
  lines(log(genotypes), harms_out, lty = 1)
  text(x = c(log(x_pos)), y = max(harms_out),
       labels = bquote('b'['con']*' = 0.10, b'['gen']*' = 0.10'),
       xpd = NA , cex = text_cex)
  
  

dev.off()
  

