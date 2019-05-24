
plot_calibration_posteriors <- function(country){
  directory = paste('C:/Users/rrag0004/Models/SNAP_TB_BMC_MED/calibrated_models/lhs_calibration_',
                    country, sep='')
  setwd(directory)
  # setwd('C:/Users/rrag0004/Models/SNAP_TB_BMC_MED/calibrated_models/lhs_calibration_Indonesia')
  
  table = read.csv('all_lhs_values.csv',header = TRUE)
  
  par_names = c('TB mortality multiplier', 'g (>= 15 y.o.)', 'infectiousness switching age',
                'g (5 to 14 y.o.)', 'g (0 to 4 y.o.)', 'transmission probability')
  
  filename = paste('calibration_graph_', country, '.png',sep='')
  png(filename = filename, res=150,width = 15, height=10,units = 'in')
  par(mfrow=c(2,3))
  
  fileConn<-file("calib.txt") 
  str = ''
  for (j in c(6, 3, 1, 5, 4, 2)){
    par_vals=table[,j]
    n = table[,8]
    
    sample_par = c()
    for (i in 1:length(n)){
      sample_par = c(sample_par, rep(par_vals[i],n[i]))
    }
    hist(sample_par, main=par_names[j], breaks=9, xlab='', cex.lab=1.3, cex.main=1.5, cex.axis =1.3)
    qs = quantile(sample_par,probs = c(0.025,0.5,0.975),names = FALSE)
    str = paste(str, "\n",par_names[j], ': ', qs[2],' (',qs[1],'-',qs[3],')',sep='')
  }
  writeLines(str,fileConn)
  close(fileConn)
  dev.off()
}

plot_calibration_ranges <- function(countries){
  directory = 'C:/Users/rrag0004/Models/SNAP_TB_BMC_MED/calibrated_models'
  setwd(directory)
  tables = list()
  for (country in countries){
    filepath = paste('lhs_calibration_',country,'/all_lhs_values.csv', sep='')
    tables[[country]] = read.csv(filepath,header = TRUE)
  }
  par_names = c('TB mortality multiplier', 'g (>= 15 y.o.)', 'infectiousness switching age',
                'g (5 to 14 y.o.)', 'g (0 to 4 y.o.)', 'transmission probability')
  
  x_lims = list('1'=c(0.8, 1.2), '2'=c(0.5, 1), '3'=c(10, 20), 
            '4'=c(0.5, 1), '5'=c(0.5, 1), '6'=c(0, 0.005))  # mortality_mult, g_adult, infectiouss age, g_teen, g_child, p_trans
  
  filename = 'calibration_ranges.png'
  png(filename = filename, res=150,width = 15, height=10,units = 'in')
  par(mfrow=c(2,3))
  counter = 0
  for (j in c(6, 3, 1, 5, 4, 2)){  # for each parameter
    counter = counter + 1
    plot(0,0,type='n',xlab='', ylab='',xlim=x_lims[[as.character(j)]], ylim=c(0.5, 5.5), main=par_names[j],
         cex.main=1.5, cex.axis =1.3, yaxt='n', bty='n')
    
    h = 0
    for (country in countries){
      h = h+1
      table = tables[[country]]
      par_vals=table[,j]
      n = table[,8]
      sample_par = c()
      for (i in 1:length(n)){
        sample_par = c(sample_par, rep(par_vals[i],n[i]))
      }
      qs = quantile(sample_par, probs = c(0.025,0.5,0.975),names = FALSE)
      segments(x0=qs[1],x1 = qs[3],y0=h,y1 = h,lwd=2)
      points(x=qs[2],y=h,cex=2,pch=16)
      if (counter %in% c(1,4)){
        mtext(text = country,side = 2,at = h,adj=0, las=1,line = 3, cex=1.1)
      }
    }
  }
  
  dev.off()
}


countries = c('India', 'Indonesia', 'China', 'Philippines', 'Pakistan')
for (country in countries){
  plot_calibration_posteriors(country)
}

plot_calibration_ranges(countries)

