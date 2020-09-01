#testing character dis hypothesis
#in Xylotrupes beetles

library(ape);
library(phytools);
my.tree <- read.nexus('XgPopTre.nex');

#plot tree to map
states <- read.table('coor.txt', header = TRUE);
row.names(states) <- states[,1];
lat <- NULL;
lon <- NULL;
for(i in 1:length(my.tree$tip.label)){
	index.sample <- which(states$taxa == my.tree$tip.label[i]);
	lat <- c(lat, states[i, 2]);
	lon <- c(lon, states[i, 3]);
}
coords <- states[,2:3];


phylo.to.map(my.tree, coords, ylim = c(-20,30), xlim = c(70,170), fsize = 1.5, cex.points=c(2,2));


#remove tips that we do not know the horn character state
my.treeED <- drop.tip(my.tree,c('LAMG','MALX','SIBB','PAGB','KLUB','EJAX','GUMF','WETF','BUTL','MUNA',"AMBG",'XSP'));
char.state <- read.csv('state_BAM.csv');
states <- char.state[,2];
names(states) <- char.state[,1];

ERAnc <- ace(states, my.treeED, type = 'discrete', model = 'ER');
ARDAnc <- ace(states, my.treeED, type = 'discrete', model = 'ARD');

#the resulted lik is not sig diff between models; P = 0.4906723
1-pchisq(2*abs(ERAnc$loglik - ARDAnc$loglik),1);

#garbage test see Luke Harmon 2019's book chapter 8.5
#the log likelihood value of the garbage model is
#n0*log(p) + (n-n0)*log(1-p)
#where n0 is 37, the total number n is 69
37*log(37/69) + 32*log(1-(37/69));
#the resulted lnL is -47.64548, which is much worse that those based on ER (-35.1983) and ARD (-34.96077) models

#plot empirical tree
states <- char.state[,2];
names(states) <- char.state[,1];
fitER<-rerootingMethod(my.treeED,states,model='ER');

xx <- matrix(c(1,1,1,1,2,2), nrow = 2, ncol = 3);
quartz.options(height = 6, width = 6, dpi = 72);
layout(xx);
par(oma=c(1,1,1,1));

plotTree(my.treeED,setEnv=TRUE,offset=0.5,fsize=0.5, show.tip.label = FALSE)
nodelabels(node=as.numeric(rownames(fitER$marginal.anc)),pie=fitER$marginal.anc,piecol=c('white', 'black'),cex=0.5)
habQ <- ERAnc$index.matrix;
habQ[2,1] <- 0.26;
habQ[1,2] <- 0.26;
library(diagram);
plotmat(habQ, main = 'ML estimated transition rates', name = c('absent','present'), box.col = c('white','black'), txt.col = c('black','white'));



library(phytools);
#use simmap to retrieve Q
Q <- matrix(c(-0.26, 0.26, 0.26, -0.26),2,2);#2 states

num.dis.char <- NULL;
#perform 1000 simulations
	BEXBcount <- 0;
	XWBOcount <- 0;
	FFcount <- 0;
	SScount <- 0;
	GIcount <- 0;
	TTcount <- 0;
	WWcount <- 0;
	LEPULcount <- 0;
	PHASPcount <- 0;
	char.dis.count <- NULL;

for(i in 1:1000){
	
	temp.tree <- sim.history(my.tree,Q);
	count.temp <- 0;
	
	if(temp.tree$states['BE']!=temp.tree$states['XB']){BEXBcount <- BEXBcount + 1};
	if(temp.tree$states['XWIK']!=temp.tree$states['BO']){XWBOcount <- XWBOcount + 1};
	if(temp.tree$states['G']!=temp.tree$states['XIN']){GIcount <- GIcount + 1};	
	if(temp.tree$states['SOF']!=temp.tree$states['XF']){FFcount <- FFcount + 1};
	if(temp.tree$states['XFS']!=temp.tree$states['SOS']){SScount <- SScount + 1};
	if(temp.tree$states['XFT']!=temp.tree$states['XGI']){TTcount <- TTcount + 1};
	if(temp.tree$states['WETF']!=temp.tree$states['GW']){WWcount <- WWcount + 1};	
	if(temp.tree$states['PLE']!=temp.tree$states['XPUL']){LEPULcount <- LEPULcount + 1};
	if(temp.tree$states['XPHA']!=temp.tree$states['XSP']){PHASPcount <- PHASPcount + 1};	


	if(temp.tree$states['BE']!=temp.tree$states['XB']){count.temp <- count.temp + 1};
	if(temp.tree$states['BO']!=temp.tree$states['XWIK']){count.temp <- count.temp + 1};
	if(temp.tree$states['XIN']!=temp.tree$states['G']){count.temp <- count.temp + 1};
	if(temp.tree$states['SOS']!=temp.tree$states['XFS']){count.temp <- count.temp + 1};
	if(temp.tree$states['SOF']!=temp.tree$states['XF']){count.temp <- count.temp + 1};
	if(temp.tree$states['XGI']!=temp.tree$states['XFT']){count.temp <- count.temp + 1};
	if(temp.tree$states['GW']!=temp.tree$states['WETF']){count.temp <- count.temp + 1};
	if(temp.tree$states['PLE']!=temp.tree$states['XPUL']){count.temp <- count.temp + 1};
	if(temp.tree$states['XPHA']!=temp.tree$states['XSP']){count.temp <- count.temp + 1};


	char.dis.count <- c(char.dis.count,count.temp);
}
#to should the number of certain sympatric taxa evolving distinct phenotypic states
cat(c('the number of BE and XB evolving distinct phenotypes is', BEXBcount));
cat(c('the number of XWIK and BO evolving distinct phenotypes is', XWBOcount));
cat(c('the number of G and XIN evolving distinct phenotypes is', GIcount));
cat(c('the number of SOF and SOF evolving distinct phenotypes is', FFcount));
cat(c('the number of XFS and SOS evolving distinct phenotypes is', SScount));
cat(c('the number of XFT and XGI evolving distinct phenotypes is', TTcount));
cat(c('the number of WETF and GW evolving distinct phenotypes is', WWcount));
cat(c('the number of PLE and XPUL evolving distinct phenotypes is', LEPULcount));
cat(c('the number of XPHA and XSP evolving distinct phenotypes is', PHASPcount));

hist(char.dis.count, xlim = c(-1,10), xlab = 'Expected number of sympatric taxa \n exhibiting different male phenotypes', main = NULL);
abline(v = 6, col = 'red', lwd = 2, lty = 2);
abline(v = 4, col = 'red', lty = 2, lwd = 2);