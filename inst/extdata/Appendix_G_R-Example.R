library(mvtnorm)
library(MASS)
library(abind)

B0=matrix(c(670.25868,1480.24821,1938.53579,2466.25469,2837.84888,3003.52391,
            3055.38674,3132.93838,3141.18638,3159.72524,
            767.98833,1592.50266,2463.79447,3019.71976,3374.72689,3553.61387,3602.27898,
            3627.28386,3645.5656,NA,
            740.57952,1615.79681,2345.85028,2910.52511,3201.5226,3417.71335,3506.58672,
            3529.00243,NA,NA,
            862.11956,1754.90405,2534.77727,3270.85361,3739.88962,4003.00219,4125.30694,
            NA,NA,NA,
            840.94172,1859.02531,2804.54535,3445.34665,3950.47098,4185.95298,NA,NA,NA,NA,
            848.00496,2052.922,3076.13789,3861.03111,4351.57694,NA,NA,NA,NA,NA,
            901.77403,1927.88718,3003.58919,3881.41744,NA,NA,NA,NA,NA,NA,
            935.19866,2103.97736,3181.75054,NA,NA,NA,NA,NA,NA,NA,
            759.32467,1584.91057,NA,NA,NA,NA,NA,NA,NA,NA,
            723.30282,NA,NA,NA,NA,NA,NA,NA,NA,NA),10,10,byrow=TRUE)
dnom=c(39161.,38672.4628,41801.048,42263.2794,41480.8768,40214.3872,43598.5056,
       42118.324,43479.4248,49492.4106)

# Identify model to be used
#   Berquist for the Berquist-Sherman Incremental Severity
#   CapeCod for the Cape Cod
#   Hoerl for the Generalized Hoerl Curve Model with trend
#   Wright for the Generalized Hoerl Curve with individual accident year levels
#   Chain for the Chain Ladder model
model="Berquist"

# Toggle graphs off if desired
graphs=TRUE

# Toggle simulations off if desired
simulation=TRUE

# Input (B0) is a development array of cumulative averages with a the
# exposures (claims) used in the denominator appended as the last column.
# Assumption is for the same development increments as exposure increments
# and that all development lags with no development have # been removed.
# Data elements that are not available are indicated as such.  This should
# work (but not tested for) just about any subset of an upper triangular
# data matrix.  Another requirement of this code is that the matrix contain
# no columns that are all zero.

# Set tau to have columns with entries 1 through 10
tau=t(array((1:10),c(10,10)))

# Calculate incremental average matrix
A0=cbind(B0[,1],(B0[,(2:10)]+0*B0[,(1:9)])-
           (B0[,(1:9)]+0*B0[,(2:10)]))

# Generate a matrix to reflect exposure count in the variance structure
logd=log(matrix(dnom,10,10))

# Set up matrix of rows and columns, makes later calculations simpler
r=row(A0)
c=col(A0)

# msk is a mask matrix of allowable data, upper triangular assuming same
# development increments as exposure increments, msn picks off the first
# forecast diagonal, msd picks off the to date diagonal
msk=(10-r)>=c-1
msn=(10-r)==c-2
msd=(10-r)==c-1

# Amount paid to date
ptd=rowSums(B0*msd,na.rm=TRUE)

# *********************************************************************
# *   START OF MODEL SPECIFIC CODE                                    *
# *********************************************************************

if(model=="Berquist") {
  #
  # g - Assumed loss emergence model, a function of the parameters a.
  # Note g must be matrix-valued with 10 rows and 10 columns

  # g itself
  # Basic design is for g to be a function of a single parameter vector, however
  # in the simulations it is necessary to work on a matrix of parameters, one
  # row for each simulated parameter, so g.obj must be flexible enough to handle
  # both.
  # Here g.obj is the Berquist-Sherman incremental severity model
  g.obj=function(theta) {
    if(is.vector(theta))
    {outer(exp(theta[11]*(1:10)),theta[1:10])}
    else
    {array(exp(outer(theta[,11],(1:10))),c(nrow(theta),10,10))*
        aperm(array(theta[,(1:10)],c(nrow(theta),10,10)),c(1,3,2))}
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 11 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad=function(theta) {
    abind(
      aperm(array(rep(exp(theta[11]*(1:10)),10*10*10),c(10,10,10)),c(2,1,3))*
        outer((1:10),outer(rep(1,10),(1:10)),"=="),
      outer((1:10)*exp(theta[11]*(1:10)),theta[1:10]),
      along=1)
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 11 (=length(theta)), 11, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess=function(theta)  {
    aa=aperm(outer(diag(rep(1,10)),
                   array((1:10)*exp((1:10)*theta[11]),c(10,1))),c(1,4,3,2))
    abind(abind(array(0,c(10,10,10,10)),aa,along=2),
          abind(aperm(aa,c(2,1,3,4)),
                array(outer((1:10)^2*exp((1:10)*theta[11]),theta[(1:10)]),c(1,1,10,10)),
                along=2),
          along=1)
  }

  # Set up starting values.  Essentially start with classical chain ladder
  # ultimate estimates and estimate trend from that and incremental average
  # start values based on incrementals from classic chain ladder
  ptd=((!ptd==0)*ptd)+(ptd==0)*mean(ptd)
  tmp=c(
    (colSums(B0[,2:10]+0*B0[,1:9],na.rm=TRUE)/
       colSums(B0[,1:9]+0*B0[,2:10],na.rm=TRUE)),
    1)
  yy=1/(cumprod(tmp[11-(1:10)]))[11-(1:10)]
  xx=yy-c(0,yy[1:9])
  ww=t(array(xx,c(10,10)))
  uv=ptd/((10==rowSums(msk))+(10>rowSums(msk))*rowSums(msk*ww))
  tmp=na.omit(data.frame(x=1:10,y=log(uv)))
  trd=0.01
  trd=array(coef(lm(tmp$y~tmp$x))[2])[1]
  a0=c((xx*mean(uv/(exp(trd)^(1:10)))),trd)
}

if(model=="CapeCod") {
  #
  # g - Assumed loss emergence model, a function of the parameters a.
  # Note g must be matrix-valued with 10 rows and 10 columns

  # g itself
  # Basic design is for g to be a function of a single parameter vector, however
  # in the simulations it is necessary to work on a matrix of parameters, one
  # row for each simulated parameter, so g.obj must be flexible enough to handle
  # both.
  # Here g.obj is nonlinear and based on the Kramer Chain Ladder parmaterization
  g.obj=function(theta) {
    if(is.vector(theta))
    {theta[1]*outer((c(1,theta[2:10])),c(1,theta[11:19]))}
    else
    {array(theta[,1],c(nrow(theta),10,10))*
        array(cbind(1,theta[,2:10]),c(nrow(theta),10,10))*
        aperm(array(cbind(1,theta[,11:19]),c(nrow(theta),10,10)),c(1,3,2))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 19 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad=function(theta) {
    abind(
      outer(c(1,theta[2:10]),c(1,theta[11:19])),
      theta[1]*
        abind(aperm(array(t(outer((2:10),(1:10),"==")),c(10,9,10)),c(2,1,3))*
                aperm(array(c(1,theta[11:19]),c(10,10,9)),c(3,2,1)),
              aperm(array(t(outer((2:10),(1:10),"==")),c(10,9,10)),c(2,3,1))*
                aperm(array(c(1,theta[2:10]),c(10,10,9)),c(3,1,2)),
              along=1),
      along=1)
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 19 (=length(theta)), 19, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess=function(theta)  {
    a1=abind(aperm(array(t(outer((2:10),(1:10),"==")),c(10,9,10)),c(2,1,3))*
               aperm(array(c(1,theta[11:19]),c(10,10,9)),c(3,2,1)),
             aperm(array(t(outer((2:10),(1:10),"==")),c(10,9,10)),c(2,3,1))*
               aperm(array(c(1,theta[2:10]),c(10,10,9)),c(3,1,2)),
             along=1)
    a2=theta[1]*(aperm(array(1:10,c(10,10,9,9)),c(3,4,1,2))==
                   array(2:10,c(9,9,10,10)))*
      (aperm(array(1:10,c(10,10,9,9)),c(4,3,2,1))==
         aperm(array(2:10,c(9,9,10,10)),c(2,1,4,3)))
    abind(abind(array(0,c(10,10)),a1,along=1),
          abind(a1,
                abind(abind(array(0,c(9,9,10,10)),a2,along=2),
                      abind(aperm(a2,c(2,1,3,4)),array(0,c(9,9,10,10)),along=2),
                      along=1),along=1),along=2)
  }

  # Set up starting values based on development factors for columns and
  # relative sizes for rows
  ptd=((!ptd==0)*ptd)+(ptd==0)*mean(ptd)
  tmp=c(
    (colSums(B0[,2:10]+0*B0[,1:9],na.rm=TRUE)/
       colSums(B0[,1:9]+0*B0[,2:10],na.rm=TRUE)),
    1)
  yy=1/(cumprod(tmp[11-(1:10)]))[11-(1:10)]
  xx=yy-c(0,yy[1:9])
  ww=t(array(xx,c(10,10)))
  uv=ptd/((10==rowSums(msk))+(10>rowSums(msk))*rowSums(msk*ww))
  a0=c((uv[1]*xx[1]),(uv[2:10]/uv[1]),(xx[2:10]/xx[1]))
}

if(model=="Hoerl") {
  #
  # g itself
  # Basic design is for g to be a function of a single parameter vector, however
  # in the simulations it is necessary to work on a matrix of parameters, one
  # row for each simulated parameter, so g.obj must be flexible enough to handle
  # both.
  # Here g.obj is Wright's operational time model with trend added

  g.obj=function(theta) {
    if(is.vector(theta))
    {exp(theta[1]+
           colSums(abind(tau,abind(tau^2,log(tau),along=0.5),along=1)*
                     array(theta[c(2,3,4)],c(3,10,10)))+
           theta[5]*array((1:10),c(10,10)))}
    else
    {exp(array(theta[,1],c(nrow(theta),10,10))+
           colSums(abind(aperm(array(tau,c(10,10,nrow(theta))),c(3,1,2)),
                         abind(aperm(array(tau^2,c(10,10,nrow(theta))),c(3,1,2)),
                               aperm(array(log(tau),c(10,10,nrow(theta))),c(3,1,2)),along=0.5),along=1)*
                     aperm(array(theta[,c(2,3,4)],c(nrow(theta),3,10,10)),c(2,1,3,4)))+
           array(theta[,5],c(nrow(theta),10,10))*
           aperm(array((1:10),c(10,nrow(theta),10)),c(2,1,3)))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 5 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad=function(theta) {
    abind(array(1,c(10,10)),
          abind(tau,abind(tau^2,
                          abind(log(tau),
                                array((1:10),c(10,10)),along=0.5),
                          along=1),
                along=1),
          along=1)*aperm(array(g.obj(theta),c(10,10,5)),c(3,1,2))
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 5 (=length(theta)), 5, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess=function(theta)  {
    aa=aperm(
      array(abind(array(1,c(10,10)),
                  abind(tau,abind(tau^2,
                                  abind(log(tau),
                                        array((1:10),c(10,10)),along=0.5),
                                  along=1),
                        along=1),
                  along=1),c(5,10,10,5)),
      c(4,1,2,3))
    aa*aperm(aa,c(2,1,3,4))*aperm(array(g.obj(theta),c(10,10,5,5)),c(3,4,1,2))
  }

  # Base starting values on classic chain ladder forecasts and inherent trend
  ptd=((!ptd==0)*ptd)+(ptd==0)*mean(ptd)
  tmp=c(
    (colSums(B0[,2:10]+0*B0[,1:9],na.rm=TRUE)/
       colSums(B0[,1:9]+0*B0[,2:10],na.rm=TRUE)),
    1)
  yy=1/(cumprod(tmp[11-(1:10)]))[11-(1:10)]
  xx=yy-c(0,yy[1:9])
  ww=t(array(xx,c(10,10)))
  uv=ptd/((10==rowSums(msk))+(10>rowSums(msk))*rowSums(msk*ww))
  tmp=na.omit(data.frame(x=1:10,y=log(uv)))
  trd=0.01
  trd=array(coef(lm(tmp$y~tmp$x))[2])[1]
  tmp=na.omit(data.frame(x1=c(tau),x2=c(tau^2),x3=log(c(tau)),
                         y=c(log(outer(uv,xx)))))
  ccs=array(coef(lm(tmp$y~tmp$x1+tmp$x2+tmp$x3)))[1:4]
  a0=c(c(ccs),trd)
}

if(model=="Wright") {
  #
  # g itself
  # Basic design is for g to be a function of a single parameter vector, however
  # in the simulations it is necessary to work on a matrix of parameters, one
  # row for each simulated parameter, so g.obj must be flexible enough to handle
  # both.
  # Here g.obj is Wright's operational time model with separate level by
  # accident year

  g.obj=function(theta) {
    if(is.vector(theta))
    {exp(array(theta[1:10],c(10,10))+theta[11]*tau+theta[12]*tau^2+
           theta[13]*log(tau))}
    else
    {exp(array(theta[,1:10],c(nrow(theta),10,10))+
           array(theta[,11],c(nrow(theta),10,10))*
           aperm(array(tau,c(10,10,nrow(theta))),c(3,1,2))+
           array(theta[,12],c(nrow(theta),10,10))*
           aperm(array(tau^2,c(10,10,nrow(theta))),c(3,1,2))+
           array(theta[,13],c(nrow(theta),10,10))*
           aperm(array(log(tau),c(10,10,nrow(theta))),c(3,1,2)))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 13, 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad=function(theta) {
    abind(array(outer(1:10,1:10,"=="),c(10,10,10)),
          abind(tau,abind(tau^2,log(tau),
                          along=0.5),
                along=1),
          along=1)*aperm(array(g.obj(theta),c(10,10,13)),c(3,1,2))
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 13, 13, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess=function(theta)  {
    aa=array(
      abind(array(outer(1:10,1:10,"=="),c(10,10,10)),
            abind(tau,abind(tau^2,log(tau),
                            along=0.5),
                  along=1),
            along=1),
      c(13,10,10,13))
    aperm(aa,c(4,1,2,3))*aperm(aa,c(1,4,2,3))*
      aperm(array(g.obj(theta),c(10,10,13,13)),c(3,4,1,2))
  }

  # Base starting values on classic chain ladder forecasts
  ptd=((!ptd==0)*ptd)+(ptd==0)*mean(ptd)
  tmp=c(
    (colSums(B0[,2:10]+0*B0[,1:9],na.rm=TRUE)/
       colSums(B0[,1:9]+0*B0[,2:10],na.rm=TRUE)),
    1)
  yy=1/(cumprod(tmp[11-(1:10)]))[11-(1:10)]
  xx=yy-c(0,yy[1:9])
  ww=t(array(xx,c(10,10)))
  uv=ptd/((10==rowSums(msk))+(10>rowSums(msk))*rowSums(msk*ww))
  tmp=na.omit(data.frame(x1=c(tau),x2=c(tau^2),x3=log(c(tau)),
                         y=c(log(outer(uv,xx)))))
  ccs=array(coef(lm(tmp$y~tmp$x1+tmp$x2+tmp$x3)))[1:4]
  a0=c(log(uv/rowSums(exp(ccs[2]*tau+ccs[3]*tau^2+ccs[4]*log(tau)))),ccs[2:4])
}

if(model=="Chain") {
  #
  # g itself
  # Basic design is for g to be a function of a single parameter vector, however
  # in the simulations it is necessary to work on a matrix of parameters, one
  # row for each simulated parameter, so g.obj must be flexible enough to handle
  # both.
  # Here g.obj is a version of the Cape Cod model but with the restriction
  # that the expected cumulative averge payments to date equal the actual
  # average paid to date.  Because of this restriction the incremental averages
  # are expressed as a percentage times the expected ultimate by row.
  # Formulae all assume a full, square development triangle.

  g.obj=function(theta) {
    if(is.vector(theta))
    {th=t(array(c(theta[1:9],(1-sum(theta[1:9]))),c(10,10)))
    uv=ptd/((10==rowSums(msk))+(10>rowSums(msk))*rowSums(msk*th))
    th*array(uv,c(10,10))}
    else
    {th=aperm(array(cbind(theta[,1:9],(1-rowSums(theta[,1:9]))),
                    c(nrow(theta),10,10)),c(1,3,2))
    mska=aperm(array(msk,c(10,10,nrow(theta))),c(3,1,2))
    ptda=t(array(ptd,c(10,nrow(theta))))
    uva=ptda/((10==rowSums(mska,dims=2))
              +(10>rowSums(mska,dims=2))*rowSums(mska*th,dims=2))
    th*array(uva,c(nrow(theta),10,10))
    }
  }

  v1=aperm(
    array((1:10),c(10,10,9)),
    c(3,2,1))==
    array((1:9),c(9,10,10))
  v2=aperm(
    array(msk[,1:9],c(10,9,10)),
    c(2,1,3))&
    aperm(
      array(msk,c(10,10,9)),
      c(3,1,2))
  v2[,1,]=FALSE
  rsm=rowSums(msk)

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 9 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad=function(theta) {
    th=t(array(c(theta,(1-sum(theta))),c(10,10)))
    psm=rowSums(msk*th)
    psc=rowSums(th[,1:9]*(1-msk[,1:9]))
    uv=ptd/((10==rsm)+(10>rsm)*psm)
    uva=aperm(array(uv,c(10,10,9)),c(3,1,2))
    thj=aperm(
      array(
        outer(
          (1/psm),c(theta,1-sum(theta))),
        c(10,10,9)),
      c(3,1,2))
    xx=uva*(v1-v2*thj)
    xx[,,10]=-uva[,,1]*((1-v2[,,1])+v2[,,1]*t(array((1-psc)/psm,c(10,9))))
    xx[,1,10]=-uv[1]
    xx
  }

  d1=array((1:9),c(9,9,10,10))
  d2=aperm(
    array((1:9),c(9,9,10,10)),
    c(2,1,3,4))
  d3=aperm(
    array((1:10),c(10,9,9,10)),
    c(2,3,1,4))
  d4=aperm(
    array((1:10),c(10,9,9,10)),
    c(2,3,4,1))
  rsma=aperm(
    array(rsm,c(10,9,9,10)),
    c(2,3,1,4))
  mm1=(d1<=rsma)&(d2<=rsma)
  mm2=((d4==d1)&(d2<=rsma))|((d4==d2)&(d1<=rsma))
  mm3=(d1==d4)&(d2==d4)&(d1<=rsma)&(d2<=rsma)
  mm5=(((d1>rsma)&(d2<=rsma))|((d2>rsma)&(d1<=rsma)))&!((d1>rsma)&(d2>rsma))

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 9 (=length(theta)), 9, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess=function(theta)  {
    th=t(array(c(theta,(1-sum(theta))),c(10,10)))
    psm=rowSums(msk*th)
    psc=rowSums(th[,1:9]*(1-msk[,1:9]))
    uv=ptd/(((10==rowSums(msk))+(10>rowSums(msk))*psm)^2)
    psc=rowSums(th[,1:9]*(1-msk[,1:9]))
    uva=aperm(array(uv,c(10,10,9,9)),c(3,4,1,2))
    thj=aperm(
      array(
        outer(
          (1/psm),2*c(theta,1-sum(theta))),
        c(10,10,9,9)),
      c(3,4,1,2))
    xx=uva*(thj*mm1-mm3-mm2)
    xx[,,,10]=uva[,,,1]*(mm5[,,,1]+mm1[,,,1]*
                           aperm(
                             array(2*(1-psc)/psm,c(10,9,9)),
                             c(2,3,1)))
    xx[,,1,]=0
    xx
  }

  # Starting values essentially classical chain ladder

  tmp=c(
    (colSums(B0[,2:10]+0*B0[,1:9],na.rm=TRUE)/
       colSums(B0[,1:9]+0*B0[,2:10],na.rm=TRUE)),
    1)
  yy=1/(cumprod(tmp[11-(1:10)]))[11-(1:10)]
  a0=(yy-c(0,yy[1:9]))[1:9]

}

# *********************************************************************
# *   END OF MODEL SPECIFIC CODE                                      *
# *********************************************************************

# Negative loglikelihood function, to be minimized.  Note that the general
# form of the model has parameters in addition to those in the loss model,
# namely the power for the variance and the constant of proprtionality that
# varies by column.  So if the original model has k parameters with 10
# columns of data, the total objective function has k+11 parameters
l.obj=function(a,A) {
  npar=length(a)-2
  e=g.obj(a[1:npar])
  v=exp(-outer(logd[,1],rep(a[npar+1],10),"-"))*(e^2)^a[npar+2]
  t1=log(2*pi*v)/2
  t2=(A-e)^2/(2*v)
  sum(t1+t2,na.rm=TRUE)}

# Gradient of the objective function
l.grad=function(a,A) {
  npar=length(a)-2
  p=a[npar+2]
  Av=aperm(array(A,c(10,10,npar)),c(3,1,2))
  e=g.obj(a[1:npar])
  ev=aperm(array(e,c(10,10,npar)),c(3,1,2))
  v=exp(-outer(logd[,1],rep(a[npar+1],10),"-"))*(e^2)^p
  vv=aperm(array(v,c(10,10,npar)),c(3,1,2))
  dt=rowSums(g.grad(a[1:npar])*((p/ev)+(ev-Av)/vv-p*(Av-ev)^2/(vv*ev)),
             na.rm=TRUE,dims=1)
  yy=1-(A-e)^2/v
  dk=sum(yy/2,na.rm=TRUE)
  dp=sum(yy*log(e^2)/2,na.rm=TRUE)
  c(dt,dk,dp)
}

# Hessian of the objective function
# e is the expectated value matrix
# v is the matrix of variances
# A, e, v all have shape c(10,10)
# The variables _v are copies of the originals to shape c(npar,10,10), paralleling
# the gradient of g.
# The variables _m are copies of the originals to shape c(npar,npar,10,10),
# paralleling the hessian of g
l.hess=function(a,A) {
  npar=length(a)-2
  p=a[npar+2]
  Av=aperm(array(A,c(10,10,npar)),c(3,1,2))
  Am=aperm(array(A,c(10,10,npar,npar)),c(3,4,1,2))
  e=g.obj(a[1:npar])
  ev=aperm(array(e,c(10,10,npar)),c(3,1,2))
  em=aperm(array(e,c(10,10,npar,npar)),c(3,4,1,2))
  v=exp(-outer(logd[,1],rep(a[npar+1],10),"-"))*(e^2)^p
  vv=aperm(array(v,c(10,10,npar)),c(3,1,2))
  vm=aperm(array(v,c(10,10,npar,npar)),c(3,4,1,2))
  g1=g.grad(a[1:npar])
  gg=aperm(array(g1,c(npar,10,10,npar)),c(4,1,2,3))
  gg=gg*aperm(gg,c(2,1,3,4))
  gh=g.hess(a[1:npar])
  dtt=rowSums(
    gh*(p/em+(em-Am)/vm-p*(Am-em)^2/(vm*em))+
      gg*(1/vm+4*p*(Am-em)/(vm*em)+p*(2*p+1)*(Am-em)^2/(vm*em^2)-p/em^2),
    dims=2,na.rm=TRUE)
  dkt=rowSums((g1*(Av-ev)+p*g1*(Av-ev)^2/ev)/vv,na.rm=TRUE)
  dtp=rowSums(g1*(1/ev+(log(ev^2)*(Av-ev)+(p*log(ev^2)-1)*(Av-ev)^2/ev)/vv),
              na.rm=TRUE)
  dkk=sum((A-e)^2/(2*v),na.rm=TRUE)
  dpk=sum(log(e^2)*(A-e)^2/(2*v),na.rm=TRUE)
  dpp=sum(log(e^2)^2*(A-e)^2/(2*v),na.rm=TRUE)
  m1=rbind(array(dkt),c(dtp))
  rbind(cbind(dtt,t(m1)),cbind(m1,rbind(cbind(dkk,c(dpk)),c(dpk,dpp))))
}

# End funciton specificaitons now on to the minimization

# Get starting values for kappa and p parameters, default 10 and 1
ttt=c(10,1)

# For starting values use fitted objective function and assume variance for a
# cell is estimated by the square of the difference between actual and expected
# averages.  Note since log(0) is -Inf we need to go through some machinations
# to prep the y values for the fit
E=g.obj(a0)
yyy=(A0-E)^2
yyy=logd+log(((yyy!=0)*yyy)-(yyy==0))
sss=na.omit(data.frame(x=c(log(E^2)),y=c(yyy)))
ttt=array(coef(lm(sss$y~sss$x)))[1:2]
a0=c(a0,ttt)

set.seed(1) # to check reproducibility with new code
max=list(iter.max=10000,eval.max=10000)

# Actual minimization
mle= nlminb(a0,l.obj,gradient=l.grad,hessian=l.hess,
            scale=abs(1/(2*((a0*(a0!=0))+(1*(a0==0))))),A=A0,control=max)

# mean and var are model fitted values, stres standardized residuals
npar=length(a0)-2
p=mle$par[npar+2]
mean=g.obj(mle$par[1:npar])
var=exp(-outer(logd[,1],rep(mle$par[npar+1],10),"-"))*(mean^2)^p
stres=(A0-mean)/sqrt(var)
g1=g.grad(mle$par[1:npar])
gg=aperm(array(g1,c(npar,10,10,npar)),c(4,1,2,3))
gg=gg*aperm(gg,c(2,1,3,4))
meanv=aperm(array(mean,c(10,10,npar)),c(3,1,2))
meanm=aperm(array(mean,c(10,10,npar,npar)),c(3,4,1,2))
varm=aperm(array(var,c(10,10,npar,npar)),c(3,4,1,2))

# Masks to screen out NA entries in original input matrix
s=0*A0
sv=aperm(array(s,c(10,10,npar)),c(3,1,2))
sm=aperm(array(s,c(10,10,npar,npar)),c(3,4,1,2))

# Calculate the information matrix using second derivatives of the
# log likelihood function

# Second with respect to theta parameters
tt=rowSums(sm+gg*(1/varm+2*p^2/(meanm^2)),dims=2,na.rm=TRUE)

# Second with respect to theta and kappa
kt=p*rowSums(sv+g1/meanv,na.rm=TRUE)

# Second with respect to p and theta
tp=p*rowSums(sv+g1*log(meanv^2)/meanv,na.rm=TRUE)

# Second with respect to kappa
kk=(1/2)*sum(1+s,na.rm=TRUE)

# Second with respect to p and kappa
pk=(1/2)*sum(s+log(mean^2),na.rm=TRUE)

# Second with respect to p
pp=(1/2)*sum(s+log(mean^2)^2,na.rm=TRUE)

# Create information matrix in blocks
m1=rbind(array(kt),c(tp))
inf=rbind(cbind(tt,t(m1)),cbind(m1,rbind(c(kk,pk),c(pk,pp))))

# Variance-covariance matrix for parameters, inverse of information matrix
vcov=solve(inf)

# Initialize simulation array to keep simulation results
sim=matrix(0,0,11)
smn=matrix(0,0,11)
spm=matrix(0,0,npar+2)

# Simulation for distribution of future amounts
# Want 10,000 simulations, but exceeds R capacity, so do
# in batches of 5,000
nsim=5000
smsk=aperm(array(c(msk),c(10,10,nsim)),c(3,1,2))
smsn=aperm(array(c(msn),c(10,10,nsim)),c(3,1,2))

if(simulation) {for (i in 1:5) {
  # Randomly generate parameters from multivariate normal
  spar=rmvnorm(nsim,mle$par,vcov)

  # Arrays to calculate simulated means
  esim=g.obj(spar)

  # Arrays to calculate simulated variances
  ksim=exp(aperm(outer(array(spar[,c(npar+1)],c(nsim,10)),log(dnom),"-"),c(1,3,2)))
  psim=array(spar[,npar+2],c(nsim,10,10))
  vsim=ksim*(esim^2)^psim

  # Randomly simulate future averages
  temp=array(rnorm(nsim*10*10,c(esim),sqrt(c(vsim))),c(nsim,10,10))

  # Combine to total by exposure period and in aggregate
  # notice separate array with name ending in "n" to capture
  # forecast for next accounting period
  sdnm=t(matrix(dnom,10,nsim))
  fore=sdnm*rowSums(temp*!smsk,dims=2)
  forn=sdnm*rowSums(temp*smsn,dims=2)

  # Cumulate and return for another 5,000
  sim=rbind(sim,cbind(fore,rowSums(fore)))
  smn=rbind(smn,cbind(forn,rowSums(forn)))
  spm=rbind(spm,spar)
}}

model
summary(sim)
summary(smn)

# Scatter plots of residuals & Distribution of Forecasts
if(graphs) {x11(title=model)

  # Prep data for lines for averages in scatter plots of standardized residuals
  ttt=array(cbind(c(r+c-1),c(stres)),c(length(c(stres)),2,19))
  sss=t(array((1:19),c(19,length(c(stres)))))

  # Plotting
  par(mfrow=c(2,2))
  plot(na.omit(cbind(c(r+c-1),c(stres))),
       main="Standardized Residuals by CY",xlab="CY",
       ylab="Standardized Residual",pch=18)
  lines(na.omit(list(x=(1:19),y=colSums(ttt[,2,]*
                                          (ttt[,1,]==sss),na.rm=TRUE)/colSums((ttt[,1,]==sss)+
                                                                                0*ttt[,2,],na.rm=TRUE))),col="red")
  plot(na.omit(cbind(c(c),c(stres))),
       main="Standardized Residuals by Lag",xlab="Lag",
       ylab="Standardized Residual",pch=18)
  lines(na.omit(list(x=c[1,],y=colSums(stres,na.rm=TRUE)/
                       colSums(1+0*stres,na.rm=TRUE))),col="red")
  qqnorm(c(stres))
  qqline(c(stres))
  if(simulation) {proc=list(x=(density(sim[,11]))$x,
                            y=dnorm((density(sim[,11]))$x,
                                    sum(matrix(c(dnom),10,10)*mean*!msk),
                                    sqrt(sum(matrix(c(dnom),10,10)^2*var*!msk))))
  truehist(sim[,11],ymax=max(proc$y),
           main="All Years Combined Future Amounts",xlab="Aggregate")
  lines(proc)}}

# Summary of mean, standard deviation, and 90% confidence interval from
# simulation, similar for one-period forecast
sumr=matrix(0,0,4)
sumn=matrix(0,0,4)

for (i in 1:11) {
  sumr=rbind(sumr,c(mean(sim[,i]),sd(sim[,i]),quantile(sim[,i],c(.05,.95))))
  sumn=rbind(sumn,c(mean(smn[,i]),sd(smn[,i]),quantile(smn[,i],c(.05,.95))))
}


