# Create a trumpet curve

#get dimensions of temp
num_col<-ncol(temp)
num_row<-nrow(temp)

#create t1 matrix with 3 columns(min,mean and max) for each depth (number of columns in temp)
t1=temp[1:num_col,1:3]

for (column in 1:num_col){
  # for each column in t1 substitute the min, mean and max of time series for each column in temp
  
  t1[column,1]=min(temp[1:num_row,column])
  print(column)
  t1[column,2] = mean(temp[1:num_row,column])
  t1[column,3] = max(temp[1:num_row,column])
}

plot (t1[,1],type="l",ylim=c(min(temp),max(temp)))
lines(t1[,2],col=2)
lines(t1[,3],col=3)


  
