# LTBI prevalence compare
model_ltbi = list('India'=c(7.01,11.23,16.26), 'Pakistan'=c(6.08,9.94,14.60), 'Philippines'=c(26.53, 34.05, 41.44),
             'China'=c(5.39, 10.33, 16.35), 'Indonesia'=c(18.61,26.03,34.39))

setwd('C:/Users/rrag0004/Models/SNAP_TB_BMC_MED/')

ltbi_houben = read.csv('data/LTBI_Houben.csv', header=TRUE)
colnames(ltbi_houben)[1]='country'

png(filename = 'ltbi_compare.png', res=150,width = 10, height=10,units = 'in')

plot(0,0,type='n',xlab='',ylab='LTBI prevalence (%)',cex.lab=1.5, cex.axis=1.3, xlim=c(0.5,5.5), ylim=c(0,100), xaxt='n')

delta = 0.1
x = 0
col_houb = rgb(red = 57,green = 106,blue = 177,maxColorValue = 255)
col_snap = rgb(red = 217,green = 124,blue = 48,maxColorValue = 255)
for (country in c('India', 'Indonesia', 'China', 'Philippines', 'Pakistan')){
  x = x+1
  
  # Houben's LTBI estimates
  houben_data = ltbi_houben[ltbi_houben$country == country,]
  houben_mean = 100*houben_data$n_ltbi / houben_data$population
  houben_low = 100*houben_data$n_ltbi_low / houben_data$population
  houben_high = 100*houben_data$n_ltbi_high / houben_data$population
  
  segments(x0=x-delta,x1=x-delta,y0=houben_low,y1=houben_high,lwd=2, col=col_houb)
  points(x=x-delta, y=houben_mean, cex=1.5,pch=16, col=col_houb)
  mtext(text = country, side = 1, cex=1.3, at=x, las=1, line=1.5)
  
  # SNAP-TB estimates
  segments(x0=x+delta,x1=x+delta,y0=model_ltbi[[country]][1],y1=model_ltbi[[country]][3],lwd=2, col=col_snap)
  points(x=x+delta, y=model_ltbi[[country]][2], cex=1.5,pch=16, col=col_snap)

  if (country == 'China'){
    segments(x0=x,x1=x,y0=13,y1=20,lwd=2, col='black')
  }
    
  legend(x=3.5, y=90,legend = c('Houben & Dodd 2016', 'SNAP-TB model'),col = c(col_houb,col_snap),lty = 0,pch = 16,cex=1.3)  
}


dev.off()

