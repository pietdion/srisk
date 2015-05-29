NB. =================================
NB.  All below here is old stuff
NB. ================================

dpl=:('xticpos 03 04 05 06 07 08 09 10 11 12 13 14 15'&plot)@(yrmo&;)

stop


ld=: 3 : 0
 'u v'=.P"1 y
 x=.(>:@i.% >:)#v
 ((covar"2 v,:"1 (u>"1 0 x))%"1-:x*1-x);x
)

pl=: 3 : 0
  pd 'new;reset'
  tics=.'-20 -15 -10 -5 0 5 10 15'
  pd 'pensize 2;type line;color black;xticpos ',tics,';yticpos ',tics
  for_i. i.#y do.
    pd i{y
  end.
  pd 'show'
return.
   pd 'reset;new'
   tics=.'-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30'
   pd 'pensize 2;type dot;color red;xticpos ',tics,';yticpos ',tics
   pd rm;ri
   pd 'type line;color green'
   pd LRMES
   pd 'show'
)



pdraw=: 4 : '(<.x*#y){/:~y'        NB. percentile draws from y
zdraw=: 4 : '(Phi x) pdraw y'   NB. zeta draws from y
   
covzeta=:cov@:|:@:(zeta"1)@|:

fitarch=: 3 : 0"1
  'y' Rset y 
  R 'r=garch(y,c(1,1),include.mean=TRUE)'
  sig=.}.{."1 >{:{.Rget 'r$fitted.values'
  R 'r=summary(r)'
  wrs 'a0 a1 b1'=.coef['coef stdev tval pval'=.|:>{:{.Rget'r$coef'
  'yn sign'=.{:y,.0,sig
  uf=.yf=.sigf=.''
  for_t. i.250 do.
    sigf=.sigf,sign=.%:a0+(a1**:yn)+b1**:sign
    yf=.yf,yn=.sign*rnorm 1
  end.
  (}.y)%sig
)

tarch=: 3 : 0   NB. threshold model;library(rugarch);library(tseries) 
   'y' Rset y 
   R 'm=list(armaOrder=c(0,0),include.mean=TRUE)'
   R 'v=list(garchOrder=c(1,1),model ="gjrGARCH")'
   R 'u=ugarchspec(mean.model=m,variance.model=v)'
   R 'fit=ugarchfit(u,y,solver.control=list(trace=0));print(fit)'
   R 'res=fit@fit'
   (n)=.{."1>{:{.matcoef[n=a.>{.>{:{:wrs matcoef=.Rget 'res$matcoef'
   sigma=:Rget 'res$sigma'
)


pctarch=: 3 : 0
  sig=.|: tarch"1 y
  'mu s'=.(mean,:sd) sig 
  'U D V'=.svd cor sig
   -U tp |:|(sig-"1 mu)
NB. 'U D V'=.svd cov y
NB. tarch"1 |:U tp"2 1 (y-"1 mu)
)
 

plot pctarch zeta"1 cba,nab,anz,:mqg


getfin=: 3 : 0
  'y' Rset y
  NB. Rcmd 'library(tseries)'
  R 'y=get.hist.quote(c(y),quote="Adj",provider=c("yahoo"),method = NULL,origin ="1899-12-30",compression = "d",retclass = c("zoo"),quiet = FALSE, drop = FALSE)'
  R 'print(summary(y))'
  ,>{:1{ Rget 'y'
)

lag=: 4 : 'x(}.,:-@[}.])y' 


stop

sgarch=: 0 : 0   NB. Variance model applied to y
   var.mod <- list(model ='sGARCH', garchOrder=c(1,1))
   mean.mod <- list(armaOrder = c(1,1))
   spec <- ugarchspec(var.mod,mean.mod)
   garch <- ugarchfit(spec, y, solver.control=list(trace=0))
   print(garch@fit)
)

ap=: 0 : 0   NB. threshold model applied to y
   var.mod <- list(model ='apARCH', garchOrder=c(1,1))
   mean.mod <- list(armaOrder = c(1,1))
   f.pars <- list(delta = 1)
   spec <- ugarchspec(var.mod,mean.mod,fixed.pars=f.pars)
   garch <- ugarchfit(spec, y, solver.control=list(trace=0))
   print(garch@fit)
)
 
mug=: 0 : 0
 # GARCH-M -- Specify y and x
   mean.mod.garch.m <- list(armaOrder = c(1,1), external.regressors = matrix(x), archm=TRUE)
   spec <- ugarchspec(variance.model <- var.mod, mean.model <- mean.mod.garch.m)
   garch.m <- ugarchfit(spec=spec, data=y, solver.control=list(trace=0))
   print(garch.m)
 # Adding external regressors
   mean.mod.xreg <- list(armaOrder = c(1,1), external.regressors = matrix(x))
   spec <- ugarchspec(variance.model <- var.mod, mean.model <- mean.mod.garch.m)
   garch.xreg <- ugarchfit(spec=spec, data=y)
   garch.xreg
   print(garch.xreg)
)

rmgarch=: 0 : 0 
  data(dji30ret) 
  data <- dji30ret[, 1:3, drop = FALSE] 
  uspec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "std") 
  cspec <- cgarchspec(uspec = multispec( replicate(3, uspec)))
  cfit <- cgarchfit(cspec, data =data, spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", kernel = "epanech"),fit.control = list(eval.se = TRUE, trace = TRUE, stationarity = TRUE),solver = "solnp", solver.control = list(), out.sample = 0, parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), fit = NULL, VAR.fit = NULL) 
)

mkr=:".@(10j6&":)   NB. make real

pc=: 3 : 0
  'U D V'=.svd cor y
  getd D
)

pdcomp=: 3 : 0 
  'U D V'=. svd ln tarch"1 zeta"1 banks=.cba,wbc,nab,anz,:mqg
   (2&{."1 U) mp tarch"1 [-2{. D mp |:V
)

 pdcomp''


stop

  
simarch=: 3 : 0
  'n a b'=.y
  eps=.rnorm n+100
  s2=.1[z=.{.eps
  for_t. i.n+10 do.
    z=.z,(t{eps)*%:s2=.1+(a*<:s2)+b*<:*:{:z
  end.
  <.n*pnorm (-n){.z
)


fitrugarch=: 3 : 0"1
  'y' Rset y=.y-mu=.mean y  NB. percentages
  Rcmd 'r=ugarchfit(y,spec=ugarch(spec))'
)


 



lc=: 3 : 0
 series=.;:'anz cba mqg nab wbc banks areit asx audusd comm dspread tspread bills'
 'U D V'=. svd z=.zsc"1 r=.>".L:0 series
 K=. D mp |: V
 sigt=.fitarch"1 }:K
)

  
accum=:*/\@:>:@:%&100
 


stop

rdfire=: 3 : 0
  dat=:readcsv '~/documents/research/chi/fire.csv'
  ({.dat)=:|:}.dat
  date=:(todayno@:getdate)"1 date
  state=:s:state
  fire=:".fire
  loss=:,(". :: 0:)"1 house_loss
  temp=:,".tmax
  ffdi=:".ffdi
  wind=:".wind
  rain=:,".rainc6
)

cum=:(*/\)@:>:

pstock=: 3 : 0
  n=.#y=.cum r=.r%100['r s h'=.y
  Z=.,:1,h$0 [ G=.,:s,0
  T=.1(<h,h)}(h(<0,h)}(_1|.id h),.0),0
  H=.((3-%3),1) (0 1;h,1)} ((h+1),2)$0
  start=.(('';(h+1)$_));T,.}."1 H
  D=:start,y;"0 2 (Z,.G),T,.H
  'a cova i covi'=: |:xx=:PR KF D
  r=.{:"1 a[covr=.(<1 1){"2 sqr cova
  plot u=.r%covr+*:r
  return.
  s=.0[w=.c=.100
  for_t. i.#y do.
    pt=.t{y
    ut=.t{u 
    wt=.({:c)+({:s)*pt
    st=.(wt*ut)%pt
    ct=.(1-ut)*wt
    w=.w,wt[s=.s,st[c=.c,ct
    wrs t,ut,st,ct,wt,100*pt
  end.
  plot u
) 

sh=: 3 : 0
  n=.#r=.y%100
  c=.s=.p=.1
  for_t. i.n do.
    p=.p,({:p)*1+(t-1){r
    gain=.({:s)*({:p)*(t-1){r
    c=.c,({:c)+0.5*gain
    s=.s,({:s)-0.5*gain%{:p
  end.
  plot (+:p),:c+s*p
)

test=: 3 : 0
  x=. rnorm 200 3
 'mu C'=.(mean;cov) x
 (C+*/~ mu)%.mu
NB. +/mu%.(C+mu*"0 1 mu)-lambda*id#C
)


test ''

NB.  pstock (nab);(^10);1

stop

sh %&100 rnorm 100

stop
NB. rddat''

phi=:(%&(%:&o.2))@^@-@-:@*:     NB.  standard normal density
rhot=:%@%:@>:@%@*:              NB.  

Zsc=:zsc@:(Phi^:_1)@:P

cusph=:3 : 0  NB. cubic spline model with state starting conditions x
   h=: 1"0 P q[c=.3-%:3['z q s'=.y
   start=.(_;_ _);(0 1 1,c),:0 0 1 1
   D_t=.(1 1 0 0 0),(0 1 1 0,c),:0 0 1 0 1
   D=:start,z;"0 2 (q,.s*h) (0 0;0 3)}"1 2 D_t
)

figx=: 3 : 0
  pd 'reset;sub ',":2##y
  pd 'show'[res=:panel"1/~>y 
  header=:<'     b       r       l0',:'  t-stat     R       lf'
  header,res
)

panel=: 4 : 0
  pd 'new;pensize 1;yticpos -3 -2 -1 0 1 2 3 ;xticpos -3 -2 -1 0 1 2 3'
  if. x-:y do. <'' return. end.
  l0=.(#q)*1+^.1-*:r=.mean*/'q z'=.Zsc"1 'x y'=.(x,:y)/:"1 x[n=.#x 
  'nds b C ell'=.res_KFR['a d'=.|:ad['ad covad int covint'=.|:SMI cusph z;q;s=.^5
  tstat=.b%sig*c['nh det sig'=.nds[b=.{.b[c=.{.|getd C['l0f lif'=.ell
NB.  r_t=.rhot (psi=.b+(#d)*(phi q)*d)%s*sig*h
  pd L:_1 (<q;z);~'type dot;color black;'
  pd L:_1 (<q;r*q);~'type line;color red'
  pd L:_1 (<q;a+b*q);~'type line;color green'
NB. pd L:_1 (<(r*q);~a+b*q);~'type line;color purple'
NB.  pd L:_1 (<q;a);~'type line;color yellow'
NB.    pd L:_1 (<q;((a+b*q)%s*sig)[(r*q%%:1-*:r));~'type line;color red'
NB.  pd L:_1 (<q;r*q);~'type line;color black'
NB.  pd L:_1 (<q;r_t);~'type line;color purple'
  <8j3":L:0 (b,tstat),.(r,R),.l0,(#q)*1+^.1-*:R=.%:1-*:s*sig
)

NB. rdfire ''

fcop=: 3 : 0
  chisq=.+/*:rnorm 'df n'=.y
  chisq,:eps=.chisq+0.31*rnorm n
)


comp=: 3 : 0
  'U D V'=.svd Z=.>Zsc L:0 y
  wrs %:1-% getd %.cov |: Z
  c=.-{.D mp|:V
  figx c;|.y
)


NB.figx  (fcop 1 1000)](,<@comp)anz;cba;nab;wbc;mqg
comp anz;cba;nab  NB.;wbc;mqg;areit;asx;dspread;tspread
stop

simnorm=: 3 : 0
  'n rho'=.y
  (chol 2 2 $ 1,rho,rho,1) mp rnorm 2,n
)

fcop=: 3 : 0
  chisq=.+/*:rnorm 'df n'=.y
  chisq,:eps=.chisq+1.5*rnorm n
)

dplot=: 'dot;pensize 3'&plot
pcop=:dplot@(P L:0)  NB. Plot copula
pncop=:dplot@(qnorm@:P L:0)

cuspf=:3 : 0  NB. cubic spline model with forcing x
   h=.3-%:3['y x s'=.y 
   start=.(_;_ _);(0 1 1,h),:0 0 1 1
   D_t=.(1 1 0,s,0),(0 1 1 0,h),:0 0 1 0 1
   D=:start,y;"0 2 x (<0 0)}"0 2 D_t
)

cusp=:3 : 0  NB. cubic spline model
   h=.3-%:3['y s'=.y
   start=.('';_ _);(1 1,h),:0 1 1
   D_t=.(1 0,s,0),(1 1 0,h),:0 1 0 1
   D=:start,y;"0 2 D_t
)

sreg=: 3 : 0
  rho=.{:{.cor|:'q z'=.|: y
  'muq muz sdq sdz'=.(mean,sd) y
  CoV=.sdq%muq
  b=.rho*sdz%sdq
  a=.muz-b*muq
  a,b,rho,CoV
)


treg=: 3 : 0
 'x y'=.P"1 y
 'q z'=.(-@-:@# {. ])"1 Phi^:_1 'u v'=.(x,:y)/:"1 x
  'a b rho CoV'=.|:sreg\.q,.z
   plot   _40}."1 rho,a,:b*q
)

simnorm=: 3 : 0
  'n rho'=.y
  (chol 2 2 $ 1,rho,rho,1) mp rnorm 2,n
)




treg fcop 1 1000
  
NB. treg fcop 1 1000
  
stop


R=:Phi^:_1@P

fig1=: 3 : 0
  'a b'=.y
  q=.Phi^:_1 u=./:~ P runif n=.1000
  ssig=.%:1-(3%~*:a)+*:b
  z=.(Ez=.(a*1-+:u)+(b*q))+ssig*rnorm#q
  wrs b,{:{.cor z,.q
  pd 'reset;sub 3 4;pensize 1'
  pd L:_2 ('type dot';(<q;z));('type line';(<q;Ez))
  pd L:_2 'new';('type line';(<u;Ez));('type dot';(<u;z))
  pd L:_1 'new;type dot';(<u;v=.Phi z)
  pd L:_1 'new;type dot';(<(P z); v)
  pd 'show'
)

NB.fig1 0 0.8
NB. stop


sim=: 3 : 0   NB.  draw randomly from empirical copula
  n=.1000
  u=.pnorm ({., {.+ 20*{:)"1 rnorm n,2
  chi=. (_1++:?n#2)*  %:+/"1 *:rnorm n,3$u
  chi=.chi,.(_1++:?n#2)*  %:+/"1 *:rnorm n,3
  dplot ;/|:pnorm chi*<:+:u
)





finds=: 4 : 0
  'q z'=.Phi^:_1 'u v'=.(P"1 x,:y)/:"1 x
  {:@>{:@KFR@KF@cuspm"1 z;q;s=.^8
)
   
NB. cba finds anz 



fitpair=: 4 : 0
  'q z'=.qnorm 'u v'=.|:(|:/:{.) P"1 x,:y
  sm=.smooth@((z;q)&(,<))
NB.  for_i. 7+-:-:-:i.10 do.
NB.    zhat=.sm wrs s=.^i
NB.    wrs 10j3": ell_KFR,(%>:*:s*sig_KFR%beta),beta=.{.beta_KFR end.
  pd 'reset;new;type dot;pensize 2'
  pd u;z
  pd 'type line'
  pd  u;mu=:sm 1600
  pd 'show'
)

anz fitpair *:cba

stop



  

stop
  



C=: 3 : 0
  'x y p q'=.y
  'u v'=. P L:0 x;y
  mean (u<:p)*.(v<:q)
)

br=: 4 : 0
  h=:x
  x=:/:~rnorm 1000
  Fx=:+/\P x
  fx=: ((h}.Fx) - (-h)}.Fx)%(h}.x)-(-h)}.x
)
   

dydx=:}:@{.,:((%~/)@:(dif"1))  
F=:(,:P)@:(/:~)
f=:dydx@F
If=:({.,:{.*{:)@:f
G=:{.@f ,:(}:@{:@F - {:@If)


cop=: 3 : 0
  'u v'=.P L:0  [ 5&{. L:0 cba;anz
NB.  suv=.u (,"0/)&(/:~) v
NB. c=.(|.+/(u,.v) (*/@:<:)"1"1 3 suv)%#u
  'su sv'=./:~ L:0 u;v
  'm1 m2'=.(u;v) (<:"0 1)L:0 (su;sv) 
  
) 

update=: 4 : 0
   'q x y '=.x [ qij=.y
   (C x;y;(q+qij);q)-*:q
)

s=: 4 : 0
   q=.x
   qij=.((q;y)&update)^:_[0
   1<.qij%q*1-q
)

S=: 4 : 0        NB. Sensitivity matrix:   x=q , y=boxed list of series
   (x&s)@,"0/~y   
)
  

NB. F=:(;P)@/:~                NB.  x;F(x) from data series 
NB. f=:%~/@:>@(dif L:0)@F      NB.  f(x) from data series 
u=:>:@i. % >:              NB.  uniform grid or quantiles
ma=: 100&(mean\)           NB.  moving average

rf=: 3  : 0                NB. estimate geom mean of ratio
   u=.ma"1 u <:#>{. y  
   u; %/ ^ ma"1 ^. >  f L:0 y
)
   
panel1=: 3 : 0
   pd 'new;pensize 1'
   if. ({.-:{:)y do.
      pd 'type line' 
      pd ({.;ln@(1:-{:))>F>{.y
   else. 
      pd 'type dot'
      pd u=.P L:0 y
      u=./:~>{.u
      pd 'type line;color red'
      pd u;u s"0 2 y end.   
)


fig1=: 3 : 0
   wrs p=.":2#$y
   pd 'reset;sub ',p
   (panel1@:,)"0/~y
   pd 'show'
   pd 'pdf 320 300 ',fdir,'fig1.pdf' 
)

NB. fig1 cba;anz;mqg;wbc;banks

panel2=: 3 : 0
   pd 'new;pensize 1'
   if. ({.-:{:)y do.
      pd 'type line' 
      pd ({.;ln@(1:-{:))>rf>{.y
   else. 
      pd 'type dot'
      pd u=.P L:0 y
      u=./:~>{.u
      pd 'type line;color red'
      pd u;u s"0 2 y end.   
)


fig2=: 3 : 0
   p=.":2#$y
   pd 'reset;sub ',p
   (panel2@:,)"0/~y
   pd 'show'
   pd 'pdf 320 300 ',fdir,'fig2.pdf' 
)

NB.  fig2 wbc;anz

sens=: 3 : 0
  q=.0.8+(0.95-0.8)*u 4 
  q S"0 1 y
)

NB. sens anz;wbc;mqg;cba