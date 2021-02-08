

```{r}
# github.monte
which.github = "local";  # local ... remote 

github.local = "C:/_git_/github/RachelJudd/WSU_STATS419_Spring2021/";

```


# Handshake

library(plotrix);
library(pracma);

is.wholenumber=function(x,tol=.Machine$double.eps)
{
  abs(x-round(x))<tol;
}
handShake = function(n=1, plotMe=TRUE)
{
  if (n<1) {stop("n must be greater than or equal to 1");}
  if (!is.wholenumber(n)) {stop("n must be a whole number");}
  h=n*(n-1)/2;
  if (plotMe)
  {
    plot(0,0,type="n",main=paste0(h," handShakes for ",n," People"),xlim=c(-6,6),ylim=c(-6,6),asp=1,xlab="",ylab="",xaxt="n",yaxt="n",bty="n");
    r=5;
    xs=ys=c();
    angle=360/n;
    draw.circle(0,0,r);
    x=0;y=5;xn=yn=NULL;
    points(x,y);
    xs=c(xs,x);
    ys=c(ys,y);
    for (p in 2:n)
    {
      # xn=x+r*sin(deg2rad(angle));
      # yn=sqrt(r^2-xn^2);
      xn=x*cos(deg2rad(angle))+y*sin(deg2rad(angle));
      yn=y*cos(deg2rad(angle))-x*sin(deg2rad(angle));
      points(xn,yn,cex=p);
      x=xn;
      y=yn;
      xs=c(xs,x);
      ys=c(ys,y);
      
    }
    print(xs);
    print(ys);
    for (source in 1:n)
    {
      for (target in 1:n)
      { 
        print("pi")
        segments(xs[source],ys[source],xs[target],ys[target])
      }
    }
    
  }
  h;
}


# Alphabet Matrix

countLetterInString=function(str,letter)
{
  nchar(as.character(str)) -nchar( gsub(letter, "", str,fixed=TRUE))
}
AlphabetCounter = function(str)
{
  str=gsub("[[:space:]]", "", str)
  str=tolower(str)
  df=data.frame(matrix(0, nrow=1, ncol=27, byrow=TRUE))
  colnames(df)=c(letters,"OTHER")
  for (letter in letters)
  {
    idx = which(letters == letter)
    df[1,idx] = countLetterInString(str,letter)
    str=gsub(letter,"",str,fixed=TRUE)
  }
  df[1,27]=nchar(str)
  df;
}
# readChar()
```

```{r}
path.declaration = paste0(github.local,"datasets/declaration/");
final = readChar(paste0(path.declaration,"final.txt"),99999)
draft = readChar(paste0(path.declaration,"draft.txt"),99999)
df=NULL
df=rbind(df,AlphabetCounter(draft))
df=rbind(df,AlphabetCounter(final))
rownames(df)=c("draft","final")


```


```{r}
prop=df #proportion
prop[1,]=prop[1,]/sum(prop[1,])
prop[2,]=prop[2,]/sum(prop[2,])
rowSums(prop)
```


```{r}
myMatrix1 = matrix(c(0.06675487,0.01407265,0.02706279,0.03127255,0.1311042,0.02958865,0.01972576,0.05304306,0.06916045,0.002285302,0.002405581,0.03175367,0.02357469,0.06567236,0.07914361,0.02165023,0.0008419533,0.06687515,0.06879962,0.09514073,0.02934809,0.01166707,0.01563628,0.001443349,0.01299014,0.0003608371,0.02862641,
                    0.07113095,0.01413690,0.02708333,0.03764881,0.1281250,0.02678571,0.01934524,0.05193452,0.06681548,0.002380952,0.002083333,0.03392857,0.02157738,0.07202381,0.07604167,0.02053571,0.0008928571,0.06309524,0.07098214,0.09464286,0.03065476,0.01101190,0.01443452,0.001636905,0.01205357,0.0005952381,0.02842262), nrow=2, byrow=TRUE);
prop2 = prop
dp = data.frame( draft = c(0.06675487,0.01407265,0.02706279,0.03127255,0.1311042,0.02958865,0.01972576,0.05304306,0.06916045,0.002285302,0.002405581,0.03175367,0.02357469,0.06567236,0.07914361,0.02165023,0.0008419533,0.06687515,0.06879962,0.09514073,0.02934809,0.01166707,0.01563628,0.001443349,0.01299014,0.0003608371,0.02862641), 
                 final = c(0.07113095,0.01413690,0.02708333,0.03764881,0.1281250,0.02678571,0.01934524,0.05193452,0.06681548,0.002380952,0.002083333,0.03392857,0.02157738,0.07202381,0.07604167,0.02053571,0.0008928571,0.06309524,0.07098214,0.09464286,0.03065476,0.01101190,0.01443452,0.001636905,0.01205357,0.0005952381,0.02842262));
rownames(dp)=c(letters,"OTHER");

prop2 = do.call(rbind, dp);

barplot(myMatrix1, names.arg=(rownames(dp)), cex.names=0.78, cex.axis=0.9, beside = TRUE, ylim = c(0,0.14),axisnames=TRUE, axis.lty=1, las=2, legend.text = colnames(dp), xlab="Alphabet Count", ylab="Proportion")
```



# Compute 3x3 Determinant
```{r}

computeDeterminant2=function(myMatrix)
{
  nrow=nrow(myMatrix);
  ncol=ncol(myMatrix);
  if (nrow != ncol) {stop("Matrix must be square");}
  if (nrow != 2) {stop("Matrix must be 2 by 2");}
  a=myMatrix[1,1]
  b=myMatrix[1,2]
  c=myMatrix[2,1]
  d=myMatrix[2,2]
  
  a*d-b*c;
}
computeDeterminant3=function(myMatrix)
{
  nrow=nrow(myMatrix);
  ncol=ncol(myMatrix);
  if (nrow != ncol) {stop("Matrix must be square");}
  if (nrow != 3) {stop("Matrix must be 3 by 3");}
  a=myMatrix[1,1]
  b=myMatrix[1,2]
  c=myMatrix[1,3]
  d=myMatrix[2,1]
  e=myMatrix[2,2]
  f=myMatrix[2,3]
  g=myMatrix[3,1]
  h=myMatrix[3,2]
  i=myMatrix[3,3]
  #https://www.chilimath.com/lessons/advanced-algebra/determinant-3x3-matrix/
  a*computeDeterminant2(matrix(c(e,f,h,i), nrow=2, byrow=TRUE)) - b*computeDeterminant2(matrix(c(d,f,g,i), nrow=2, byrow=TRUE)) + c*computeDeterminant2(matrix(c(d,e,g,h), nrow=2, byrow=TRUE))
}
# compute33Determinant=function()
```

```{r}
myMatrix2 = matrix(c(1,2,3,
                    4,5,6,
                    7,8,9), nrow=3, byrow=TRUE);
zeroIsh(det(myMatrix2))
zeroIsh(computeDeterminant3(myMatrix2))
```

```{r}
myMatrix3 = matrix(c(4,8,25,
                    4,56,6,
                    7,45,9), nrow=3, byrow=TRUE);
zeroIsh(det(myMatrix3))
zeroIsh(computeDeterminant3(myMatrix3))
```


```{r}
github.remote = "https://raw..../";
# 
# if(which.github == "remote")
#   {
#   include.me = paste0( github.remote, "functions/functions-intro.R");
#   library(devtools);
#   source_url(include.me);
#   } else {
#           include.me = paste0( github.local, "functions/functions-intro.R");
#           source(include.me);
#           }
```

