getwd()
install.packages("stats")
require(stats)
my_data <- read.delim("coastalBMIsPhylum.txt")
ZetaAnalysis <- read.delim("coastalBMIsPhylum.txt")
my_data
plot(ZetaAnalysis)
ZetaAnalysis
plot(ConditionScore, zeta_1)
data('coastalBMIsPhylum')
data('ConditionScore')


plot(x=ConditionScore, y=zeta_1) 

zetaAnalysis <- read.delim('coastalBMIsPhylum.txt', sep = '\t', header = T)
plot(zetaAnalysis$ConditionScore, zetaAnalysis$zeta_1)

plot(zetaAnalysis$ConditionScore, zetaAnalysis$zeta_1,
     xlab="Condition Score",
     ylab="zeta_1")

aov(y ~ A + B + A:B, data=mydataframe)

a = aov(ConditionScore~zeta_1, data=zetaAnalysis)
summary(a)


a = aov(zeta_1~ConditionScore, data=zetaAnalysis)
summary(a)

zetaAnalysis <- read.delim("coastalBMIsPhylum.txt")
!is.na(zetaAnalysis$SQO_BRI)
zetaAnalysis<- zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]
nrow(zetaAnalysis)


a = aov(zeta_2_DD~stratumType*year*SQO_BRI+stratumType:year:SQO_BRI, data=zetaAnalysis)
summary(a)

citation("zetadiv")
install.packages("zetadiv")


summary(aov(zeta_2_DD~stratumType*year*SQO_BRI,data=zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]))
