################
# Gender Style #
################

library(dplyr)
library(VIM)
library(lme4)
library(lmerTest)
library(effects)

rm(list = ls())

set.seed(1979)
options(scipen = 10)
rm(list = ls())

G <- read.csv("gender_dataset_18801999.csv", header=TRUE)

G <- droplevels(filter(G, Inputfile != "input/NVT_1979_94_Hirs.txt")) #all NAs

G$Namen_p.z <- scale(G$Namen_p)
summary(mm.prop <- lmer(Namen_p.z ~ GENDER + (1|AUTHORS), data=G))
plot(allEffects(mm.prop))

G <- dplyr::select(G, Word_per_doc, Namen_p, Let_per_wrd_nsam, Morf_per_wrd_zn, Wrd_freq_zn_log, Alg_nw_p, Alg_ww_p, Abstr_ww_p, MTLD_lem, Pv_Frog_per_zin, D_level, Wrd_per_dz, Mv_inbed_per_zin, Bijzin_per_zin, Inputfile, YEAR, AUTHORS, GENDER)
mice_plot <- aggr(G, col=c("forestgreen","red"), numbers=TRUE, sortVars=TRUE)

G <- dplyr::select(G, !c(Mv_inbed_per_zin))
mice_plot <- aggr(G, col=c("forestgreen","red"), numbers=TRUE, sortVars=TRUE)

G$Namen_p.z <- as.vector(scale(G$Namen_p))

G$Let_per_wrd_nsam.z <- as.vector(scale(G$Let_per_wrd_nsam))
G$Morf_per_wrd_zn.z <- as.vector(scale(G$Morf_per_wrd_zn))
G$Wrd_freq_zn_log.z <- as.vector(scale(G$Wrd_freq_zn_log))
G$Alg_nw_p.z <- as.vector(scale(G$Alg_nw_p))
G$Alg_ww_p.z <- as.vector(scale(G$Alg_ww_p))
G$Abstr_ww_p.z <- as.vector(scale(G$Abstr_ww_p))
G$MTLD_lem.z <- as.vector(scale(G$MTLD_lem))

G$Pv_Frog_per_zin.z <- as.vector(scale(G$Pv_Frog_per_zin))
G$D_level.z <- as.vector(scale(G$D_level))
G$Wrd_per_dz.z <- as.vector(scale(G$Wrd_per_dz))
G$Multiple_Sub.z <- as.vector(scale(G$Bijzin_per_zin))

summary(manova.model.overall <- manova(cbind(Let_per_wrd_nsam.z, Morf_per_wrd_zn.z, Wrd_freq_zn_log.z, Alg_nw_p.z, Alg_ww_p.z, Abstr_ww_p, MTLD_lem.z, Pv_Frog_per_zin.z, D_level.z, Wrd_per_dz.z) ~ GENDER, data=G))

summary(mm.size.lex <- lmer(Let_per_wrd_nsam.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.size.morph <- lmer(Morf_per_wrd_zn.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.size.freq <- lmer(Wrd_freq_zn_log.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.div.lex <- lmer(MTLD_lem.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.alg.n <- lmer(Alg_nw_p.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.alg.v <- lmer(Alg_ww_p.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.abs.v <- lmer(Abstr_ww_p.z ~ GENDER + (1|AUTHORS), data=G))

summary(mm.syn.pv <- lmer(Pv_Frog_per_zin.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.syn.Dlevel <- lmer(D_level.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.size.syn <- lmer(Wrd_per_dz.z ~ GENDER + (1|AUTHORS), data=G))
summary(mm.mult.sub <- lmer(Multiple_Sub.z ~ GENDER + (1|AUTHORS), data=G))

G$year.c <- as.vector(scale(G$YEAR, center=TRUE, scale=FALSE))

summary(mm.size.lex.y <- lmer(Let_per_wrd_nsam.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.size.lex.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Size, average number of characters per word", main="")

summary(mm.size.morph.y <- lmer(Morf_per_wrd_zn.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.size.morph.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Size, average number of morphemes per word", main="")

summary(mm.size.freq.y <- lmer(Wrd_freq_zn_log.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.size.freq.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Word frequency", main="")

summary(mm.div.lex.y <- lmer(MTLD_lem.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.div.lex.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Lexical diversity", main="")

summary(mm.alg.n.y <- lmer(Alg_nw_p.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.alg.n.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="General nouns", main="")

summary(mm.alg.v.y <- lmer(Alg_ww_p.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.alg.v.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="General verbs", main="")

summary(mm.abs.v.y <- lmer(Abstr_ww_p.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.abs.v.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Abstract verbs", main="")

summary(mm.syn.pv.y <- lmer(Pv_Frog_per_zin.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.syn.pv.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Number of finite verbs", main="")

summary(mm.syn.Dlevel.y <- lmer(D_level.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.syn.Dlevel.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="D-level", main="")

summary(mm.size.syn.y <- lmer(Wrd_per_dz.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.size.syn.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Clause length", main="")

summary(mm.mult.sub.y <- lmer(Multiple_Sub.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.mult.sub.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Multiple subordination", main="")

range(G$YEAR)

summary(mm.names.y <- lmer(Namen_p.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.names.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Proportion of proper names", main="")

TTR <- read.csv("TTR_gender.csv", header=TRUE)
TTR <- dplyr::select(TTR, Inputfile, TTR_inhwrd, TTR_inhwrd_zonder_abw)
G <- left_join(G, TTR, by="Inputfile")

G$TTR.G <- G$TTR_inhwrd *  G$Word_per_doc / sqrt(G$Word_per_doc)
G$TTR.G.z <- as.vector(scale(G$TTR.G))

cor.test(G$TTR.G.z, G$MTLD_lem.z)

summary(mm.TTRG.y <- lmer(TTR.G.z ~ GENDER * year.c + (1|AUTHORS), data=G))
plot(allEffects(mm.TTRG.y), lines=list(multiline=TRUE, col="black", lty=c(1,2)), confint=list(style="bands"), xlab="Year (centered, 0 = 1939)", ylab="Guiraud TTR", main="")

#Early years (1880-1939) vs. Later years (1940-1999)#
GE <- droplevels(filter(G, YEAR < 1940))
GL <- droplevels(filter(G, YEAR >= 1940))

#Lexical measures#
summary(mm.size.lex.e <- lmer(Let_per_wrd_nsam.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.size.morph.e <- lmer(Morf_per_wrd_zn.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.size.freq.e <- lmer(Wrd_freq_zn_log.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.div.lex.e <- lmer(MTLD_lem.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.alg.n.e <- lmer(Alg_nw_p.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.alg.v.e <- lmer(Alg_ww_p.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.abs.v.e <- lmer(Abstr_ww_p.z ~ GENDER + (1|AUTHORS), data=GE))

summary(mm.size.lex.l <- lmer(Let_per_wrd_nsam.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.size.morph.l <- lmer(Morf_per_wrd_zn.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.size.freq.l <- lmer(Wrd_freq_zn_log.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.div.lex.l <- lmer(MTLD_lem.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.alg.n.l <- lmer(Alg_nw_p.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.alg.v.l <- lmer(Alg_ww_p.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.abs.v.l <- lmer(Abstr_ww_p.z ~ GENDER + (1|AUTHORS), data=GL))

#Syntactic measures#

summary(mm.syn.pv.e <- lmer(Pv_Frog_per_zin.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.syn.Dlevel.e <- lmer(D_level.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.size.syn.e <- lmer(Wrd_per_dz.z ~ GENDER + (1|AUTHORS), data=GE))
summary(mm.mult.sub.e <- lmer(Multiple_Sub.z ~ GENDER + (1|AUTHORS), data=GE))

summary(mm.syn.pv.l <- lmer(Pv_Frog_per_zin.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.syn.Dlevel.l <- lmer(D_level.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.size.syn.l <- lmer(Wrd_per_dz.z ~ GENDER + (1|AUTHORS), data=GL))
summary(mm.mult.sub.l <- lmer(Multiple_Sub.z ~ GENDER + (1|AUTHORS), data=GL))

dev.off()
