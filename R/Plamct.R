# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Title cal_B
#'
#' @param q order of polynomial function
#' @param p number of spline basis functions
#' @param kjn knots number of spline
#'
#' @return list
#' @export
#'
#' @examples
#' cj1 = cal_B(q = 2,p = 8,kjn = 5)[[1]]
cal_B<-function(q,p)
{
  cj = list()
  for(k in 1:q)
  {
    cjk_row = p-k
    cjk_col = p-k+1
    c1 = cbind(rep(0,cjk_row),diag(cjk_row))
    c2 = cbind(diag(cjk_row),rep(0,cjk_row))
    cj[[k]] = c1-c2
  }
  return(cj)
}


#' Title cal_H
#'
#' @param p number of spline basis functions
#' @param q order of polynomial function
#' @param J0 number of additive components that you want to test
#'
#' @return list
#' @export
#'
#' @examples
#' H = cal_H(p = 8,q = 2,J0 = 1)
cal_H<-function(p,q,J0)
{
  Hj = list()
  projection_matrix = diag(p-q)-(rep(1,p-q)%*%t(rep(1,p-q)))/(p-q)
  for(i in 1:J0)
  {
    Hj[[i]] = projection_matrix[1:(p-q-1),]
  }
  if(J0==1)
  {
    H = Hj[[1]]
  }else {
    H = cbind(Hj[[1]],matrix(0,nrow = nrow(Hj[[1]]),ncol = (J0-1)*ncol(projection_matrix)))
    for(i in 2:J0)
    {
      h1 = cbind(matrix(0,nrow = nrow(Hj[[i]]),ncol = (i-1)*ncol(projection_matrix)),Hj[[i]],matrix(0,nrow(Hj[[i]]),ncol = (J0-i)*ncol(projection_matrix)))
      H = rbind(H,h1)
    }
  }
  return(H)
}


#' Title cal_C
#'
#' @param p number of spline basis functions
#' @param D number of linear parameters
#' @param J number of additive functions
#' @param ind which additive component you want to test
#'
#' @return matrix
#' @export
#'
#' @examples
#' C = cal_C(p = 8,D = 3,J = 3,ind = 1)
cal_C<-function(p,D,J,ind)
{
  # for (i in 1:J0)
  # {
  C = cbind(matrix(0,nrow = p,ncol = D+(ind-1)*(p-1)),diag(1,p,p),matrix(0,nrow = p,ncol = (J-ind)*(p-1)))
  # }
  return(C)
}

#' Title cal_A
#'
#' @param Bj matrix bj
#' @param C  matrix c
#'
#' @return matrix
#' @export
#'
#' @examples
#' A = cal_A(Bj,C)
cal_A<-function(Bj,C)
{
  return(Bj%*%C)
}

#' Title myfit
#'
#' @param y response
#' @param design_x covariate x and b-spline basis
#' @param tau quantile level
#'
#' @return
#' @export list
#'
#' @examples
#'result = myfit(y,design_x,tau = 0.5)
myfit<-function(y,design_x,tau)
{
  result = rq(y~design_x-1,tau)
  return(result)
}



#' Title knots:calcu uu1
#'
#' @param p number of spline basis functions
#' @param degree spline order
#' @param boundary.knots bourdary.knots of spline
#' @param Knots_z1 knots_z1
#'
#' @return vector
#' @export
#'
#' @examples
#'uu1 = knots(p,degree,bound,z.knots.z1)
uu1<-function(p,degree,boundary.knots,Knots_z1)
{
  m = p+degree
  uu1=vector(length = p+degree+1)
  uu1[1:(degree+1)]=boundary.knots[1]
  uu1[(m-degree+1):(m+1)]=boundary.knots[2]
  uu1[(degree+2):(m-degree)]=Knots_z1
  return(uu1)
}


#' Title coef.order
#'
#' @param uu1 knots of spline
#' @param p number of spline basis functions
#' @param order derivative order now calcu :1-q
#' @param degree spline order
#'
#' @return vector
#' @export
#'
#' @examples
#' coef.order1 = coef.order(uu1,p,q,order=1,degree)
coef.order<-function(uu1,p,order,degree)
{
  index = 1:(p-order)
  coef.order = 1/(uu1[index+degree+1]-uu1[index+order])
  return(coef.order)
}

#' Title knots
#'
#' @param z additive covariate
#' @param u sequence u
#'
#' @return vector
#' @export
#'
#' @examples
#'knots_z = knots(z,u)
knots<-function(z,u)
{
  col = ncol(z)
  knots = list()
  for (i in 1:col) {
    knots[[i]] = as.numeric(quantile(z[,i], u))
  }
  return(knots)
}

#' Title z.basis
#'
#' @param z additive covariate
#' @param u vector u
#' @param ind index
#' @param p number of spline basis functions
#'
#' @return matrix
#' @export
#'
#' @examples
#'basis_z = z.basis(z,u,ind,p)
z.basis<-function(z,u,ind,p)
{
  col = ncol(z)
  basis = list()
  knots_z = knots(z,u)
  for (i in 1:col) {
    index = ifelse(ind==i,0,1)
    basis[[i]] = bs(z[,i], knots=knots_z[[i]], intercept=T, degree=degree)[,1:(p-index)]
  }
  basis.z = basis[[1]]
  for (j in 2:col) {
    basis.z = cbind(basis.z,basis[[j]])
  }
  return(basis.z)
}

#' Title kn.sic
#'
#' @param y response
#' @param x linear covariate
#' @param z additive covariate
#' @param n sample size
#' @param tau quantile level
#' @param degree order of spline
#' @param ind index indicate which additve function is testing now
#'
#' @return a number
#' @export
#'
#' @examples
#'kjn = kn.sic(y,x,z,n,tau,degree,ind = 1)
kn.sic<-function(y,x,z,n,tau,degree,ind)

{
  SIC = vector(length = maximum)
  for (kn in 0:maximum)
  {
    p = kn+degree+1
    u = seq(0, 1, length=kn+2)[-c(1,kn+2)]
    basis_z = z.basis(z,u,ind,p)
    design_total = cbind(x,basis_z)
    pkn = ncol(design_total)
    result = rq(y~design_total-1,tau)
    res = result$resid
    sic = log(sum(res*(tau-1*(res<0)))) + 0.5*log(n)*pkn/n
    SIC[kn] = sic
  }
  kn = max(which(SIC == min(SIC)))
  return(kn)
}


#' Title z.spline
#'
#' @param y response
#' @param x linear covariate
#' @param z additive covariate
#' @param n sample size
#' @param p number of spline basis functions
#' @param tau quantile level
#' @param kjn number of spline knots
#' @param degree order of spline
#' @param ind index
#'
#' @return a list
#' @export
#'
#' @examples
#'z.sp = z.spline(y,x,z,n,p,tau,kjn,degree,ind = 1)
z.spline<-function(y,x,z,n,p,tau,kjn,degree,ind)
{
  l = list()
  u = seq(0, 1, length=kjn+2)[-c(1,kjn+2)]
  knots = knots(z,u)
  basis.z = z.basis(z,u,ind,p)
  for (i in 1:ncol(z)) {
    if(ind == i)
    {
      knots.z = knots[[i]]
      boundary.knots = attr(bs(z[,i], knots=knots.z, intercept=T, degree=degree),'Boundary.knots')
    }
  }
  l$z.basis = basis.z
  l$knots.z = knots.z
  l$bound = boundary.knots
  return(l)
}


#' Title Tn
#'
#' @param H matrix H
#' @param A matrix A
#' @param Bj matrix Bj
#' @param gamma matrix gamma
#' @param coef_z coefficient of spline
#'
#' @return a number
#' @export
#'
#' @examples
#'stat.norm = Tn(H1,A1,B_z1,gamma1,coef1)
Tn<-function(H,A,Bj,gamma,coef_z)
{
  tmp = eigen(H%*%A%*%gamma%*%t(A)%*%t(H))
  tn = tmp$vectors%*%diag(1/sqrt(abs(tmp$values)))%*%t(tmp$vectors)%*%H%*%Bj%*%coef_z
  len = length(tn)
  stat.norm=((t(tn)%*%tn)-len)/sqrt(2*len)
  return(stat.norm)
}


#' Title calcu_Bz
#'
#' @param q polynomial order
#' @param p number of spline basis functions
#' @param degree order of spline
#' @param z.sp list
#'
#' @return a matrix
#' @export
#'
#' @examples
#' B_z = calcu_Bz(q,p,degree,z.sp)
calcu_Bz<-function(q,p,degree,z.sp)
{
  B = list()
  coef = list()
  for (i in 1:q) {
    B[[i]] = cal_B(q,p)[[i]]
  }
  uu1 = uu1(p,degree,z.sp$bound,z.sp$knots.z)
  for(j in 1:q){
    coef[[j]] = coef.order(uu1,p,order=j,degree)
  }
  B_z = B[[q]]*coef[[q]]
  if(q>=2)
  {
    for (t in (q-1):1) {
      B_z = B_z%*%(B[[t]]*coef[[t]])
    }
  }
  return(B_z)
}


#' Title poly.test
#'
#' @param result quantile regression result
#' @param design_x x and z_spline
#' @param tau quantile level
#' @param n sample size
#' @param J number of additive functions
#' @param D number of linear part
#' @param p number of spline basis functions
#' @param q polynomial order
#' @param degree order of spline
#' @param B.num wildbootstrap number
#' @param weight wildbootstrap weight
#' @param z.sp list
#' @param ind index
#'
#' @return a list
#' @export
#'
#' @examples
#'r = poly.test(result,design_x,tau,n,J,D,p,q,degree,B.num,weight=NULL,z.sp,ind = 1)
poly.test<-function(result,design_x,tau,n,J,D,p,q,degree,B.num,weight=NULL,z.sp,ind)
{
  r = list()
  coef_z = result$coefficients[(D+(ind-1)*(p-1)+1):(D+(ind-1)*(p-1)+p)]
  H = cal_H(p,q,1)
  C = cal_C(p,D,J,ind = ind)
  B_z = calcu_Bz(q,p,degree,z.sp)
  # B_z%*%coef_z
  # H%*%B_z%*%coef_z
  A = cal_A(B_z,C)
  cov.boot = boot.wild(tau,n,B.num,weight = NULL,result,design_x)
  r = list(H=H,A=A,B_z =B_z,cov.boot=cov.boot,coef_z=coef_z)
  return(r)
}


#' Title cal_pvalue
#'
#' @param r a list
#' @param J0 number of additive functions you want to test
#'
#' @return a number
#' @export
#'
#' @examples
#' pvalue = cal_pvalue(r,J0)
cal_pvalue<-function(r,J0)
{
  cov.boot = list()
  H = list()
  A = list()
  B_z = list()
  coef = list()
  gamma = list()
  for (i in 1:J0) {
    cov.boot[[i]] = r[[i]]$cov.boot
    H[[i]] = r[[i]]$H
    A[[i]] = r[[i]]$A
    B_z[[i]] = r[[i]]$B_z
    coef[[i]] = r[[i]]$coef_z
    gamma[[i]] = cov(t(cov.boot[[i]]))
  }
  gamma1 = as.matrix(bdiag(gamma))
  H1 = as.matrix(bdiag(H))
  A1 = as.matrix(bdiag(A))
  B_z1 = as.matrix(bdiag(B_z))
  coef1 = as.vector(unlist(coef))
  stat.norm = Tn(H1,A1,B_z1,gamma1,coef1)
  pvalue = (1-pnorm(abs(stat.norm)))*2
  return(pvalue)
}


#' Title score
#'
#' @param TAU  quantile level
#' @param e_hat estimated error
#'
#' @return a number
#' @export
#'
#' @examples
#'s = score(tau,result$residuals)
score<-function(TAU,e_hat)
{
  score = TAU - ifelse(e_hat<0,1,0)
  return(score)
}


#' Title resample
#'
#' @param weight weight
#' @param res residuals
#' @param fit fitted value
#'
#' @return vector ynew
#' @export
#'
#' @examples
#'y.boot=apply(w.boot,2,resample,res=residuals.b,fit=result$fitted)
resample<-function(weight,res,fit)
{
  re = fit + weight*abs(res)
  return(re)
}


#' Title regre
#'
#' @param y response
#' @param design_x design matrix
#'
#' @return list
#' @export
#'
#' @examples
#'result.boot = apply(y.boot,2,regre,design_x=design_x)
regre<-function(y,design_x)
{
  return(rq(y~design_x-1,tau))
}

#' Title coeff
#'
#' @param x a list
#'
#' @return
#' @export
#'
#' @examples
#' coef.boot = sapply(result.boot,coeff)
coeff<-function(x)
{
  return(x$coefficients)
}


#' Title boot.wild
#'
#' @param tau quantile level
#' @param n sample size
#' @param B.num wildbootstrap number
#' @param weight weight
#' @param result quantile regression result
#' @param design_x design matrix
#'
#' @return a matrix
#' @export
#'
#' @examples
#' cov.boot = boot.wild(tau,n,B.num,weight = NULL,result,design_x)
boot.wild<-function(tau,n,B.num,weight = NULL,result,design_x)
{
  boot.num = B.num
  if(is.null(weight))
  {
    # default weights
    w.boot=matrix(sample(c(-2*tau,2*(1-tau)),size=n*boot.num,prob=c(tau,1-tau),replace=T),n,boot.num)
  }else{
    # specified weights
    boot.n=dim(Weight)[2]
    w.boot=Weight
  }
  cc = akj(result$residuals,z=0)$dens
  w=diag(design_x%*%solve(t(design_x)%*%design_x)%*%t(design_x))
  residuals.b=result$residuals +2*score(tau,result$residuals)*w/cc
  y.boot=apply(w.boot,2,resample,res=residuals.b,fit=result$fitted)
  result.boot = apply(y.boot,2,regre,design_x=design_x)
  coef.boot = sapply(result.boot,coeff)
  return(coef.boot)
  # return(cov(t(coef.boot)))
}

#' Title Plamct.test
#'
#' @param y response y
#' @param x linear covariate
#' @param z additve covariate
#' @param n sample size
#' @param tau quantile level
#' @param degree order of spline
#' @param ind index
#' @param J number of additive functions
#' @param D number of linear part
#' @param q polynomial order
#' @param B.num wildbootstrap number
#'
#' @return a number
#' @export
#'
#' @examples
#' a= Plamct.test(y,x,z,n,tau,degree,J,D,q,ind,B.num)
Plamct.test<-function(y,x,z,n,tau,degree,J,D,q,ind,B.num)
{
  r = list()
  # y = x%*%beta+m1_null(z1)-mean(m1_null(z1))+m2(z2)-mean(m2(z2))+m3(z3)-mean(m3(z3))+rnorm(n,0,1)-qnorm(tau)
  for (j in 1:J0) {
    kjn = kn.sic(y,x,z,n,tau,degree,ind = ind[j])
    p = kjn+degree+1
    z.sp = z.spline(y,x,z,n,p,tau,kjn,degree,ind = ind[j])
    design_x = cbind(x,z.sp$z.basis)
    result = myfit(y,design_x,tau)
    r[[j]] = poly.test(result,design_x,tau,n,J,D,p,q,degree,B.num,weight=NULL,z.sp,ind = ind[j])
  }
  pvalue = cal_pvalue(r,J0)
  return(pvalue)
}

