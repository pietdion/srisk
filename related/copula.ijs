NB.  Copula stuff stuff

load '~user\kfs.ijs'
load '~user\minimize.ijs'

dir=:'c:\Documents and Settings\pdejong\my documents\research\apra2\data\'

NB.  GM & GE data

gmge=:>{:L:0 iread dir,'GMGE'
NB. gmge=:simf 2 2 11370

P=:(/:@/:%#)"1      NB. Prob transform -- percentile ranks
srt=:/:&/:{[        NB. sort  x. using y.

rcor=:cor@|:@:P     NB. rank correlation matrix

simf=: 3 : 0
  'lambda sigma N'=.y
  f=.(% runif N)^ lambda
  f,:f+sigma*rnorm N
) 

fig1=:3 : 0
   lambda=. 0.1 2 3 
   sigma=. 3 2 1
   pd 'reset;sub ', 2j0": ($lambda), $sigma
   opt=:'titlefont Symbol 15;pensize 2;ticminor 1;type dot;'
   ticpos=.": 0 0.2 0.4 0.6 0.8 1
   opt=:opt,'yticpos ',ticpos,';xticpos ',ticpos  
   panel"1 (,/lambda ,"0"0 1 sigma),.10000
   pd 'show'
NB.   pd 'eps 315 315 c:\fig1.eps'
)

panel=:3 : 0
  pd 'new;',opt,';title l=',(":{.}:y),', s=',(":{:}:y)
  pd ;/u=.P f=.simf y
  n=.1
  wrs {.{:cor |:f*n*u^<:n
  
)

simd=:3 : 0
NB.  'Y R'=.y
  n=.y['p N'=.$Y=.10000*gmge
  R=.2 2$ 1 _0.3 _0.3 1
  phi=.n&([*(^<:)~)[n&(([*^@*)%(-&1)@^@[)
  U=.pnorm (sqrt R) mp rnorm p,N
  Y=.Y ((rank@] { /:~@[))"1 U
  cov=.(mean@:*)&(-mean)
  c=.Y cov"1 phiU=:phi U
  d=.Y cov"1 phiup=:phi P +/Y 
  rho=. c % (sigy=.sd"1 Y)*sigphi=.sd"1 phiU 
  rhop=.d % sigy*sigphiup=.sd phiup
  w=.(%+/) sigy [ v=.(%+/) c
  out=.;:'meany sigy c d v rho rhop ''1-rhop/rho'' ''v*(1-rhop/rho)'' w ''rho-rhop'' ''w*(rho-rhop)'''
  out=.out,:,.L:0 (meany=.mean"1 Y);sigy;c;d;v;rho;rhop;(1-rhop%rho);(v*1-rhop%rho);w;(rho-rhop);w*rho-rhop
  out,(+/meany);(+/sigy);(+/L:0 c;d),(+/v);('meanphi';mean phiup);('sigphi';sigphiup);'';(+/v*1-rhop%rho);(+/w);'';+/w*rho-rhop
)


  

   


NB. Choose from x according to y (uniform)
M=:(<.@(* #)~ { /:~@[)"1 0


copcor=: 4 : 0
  'omega tau'=.x[s=.y
  {.{: cor |: (^ - omega * (| tau - P s))*s   
)


ell=: 4 : 0
  'a b'=.x
  'rho n'=.y
  s=.simf rho,n,0,1
  100 qt (+/b%"0 1 (1-s)^"1 0 %a)-+/b
)

snorm=: 3 : 0"1
  a=.0.10
  <:%(%a)%:1-y
)


rat=: 3 : 0
   rho=.(>:@i.%])y
   cor=.%@%:@>:@*:@^. rho
   (cor,:rho);~rho
)

fig=:3 : 0
   rho=.0.1 0.5 0.75 0.95 [ N=.1000 [ t=.0 [ df=.1
   pd 'reset;sub 2 2'
   opt=:'pensize 2;ticminor 1;type dot;'
   ticpos=.":0 1  NB. 0.2 0.4 0.6 0.8 1
   opt=:opt,'yticpos ',ticpos,';xticpos ',ticpos  
   panel"1 rho,"0 1 N,t,df
   pd 'show'
NB.   pd 'eps 315 315 c:\briefcase\apra2\fig3.eps '  NB. 480 360 default 
)


fig3=: 3 : 0
   rho=.0.1 0.5 0.75 0.9 [ N=.10000 [ t=.0 [ df=.1
   pd 'reset;sub 2 2'
   opt=:'pensize 2;ticminor 1;type dot;'
   ticpos=.": 0 0.2 0.4 0.6 0.8 1
   opt=:opt,'yticpos ',ticpos,';xticpos ',ticpos
   panel 0.5,N,0,1 
   panel 0.5,N,1,-1  
   panel 0.5,N,0,-2
   panel 0.5,N,2,2
   pd 'show'
   pd 'eps 420 315 c:\briefcase\apra2\fig3.eps '  NB. 480 360 default 
   pd 'clip'
) 

tau=: 3 : 0    NB.  Kendall's tau.
  ((*:-]){:$y) %~ (+/^:2@:*@:*&>)/~ <"2 -/~"1 y
) 


gibbs=: 3 : 0 
  'df N'=.y.
  srcor=.{.@{:@rcor 
  r=.srcor gmge [ n=.{:$gmge [ s=.-4.24
  for. i.N do. 
    s=.s,({:s)+r-srcor simf ({:s),n,1 2 
  end.
  fe=.,:rnorm 2,n
  for. i.N do.
  fe=.fe,({.,:{:-{.) (simf ({:s),n,1 2) srt"1 gmge
  end.
  (mean,:sd) {."2 }.fe
)

gof=: 3 : 0         NB. freq count
  'dim n'=.$data['data cuts'=.y.
  cp=.<.cuts*P"1 data 
  'freq index'=.(#/.~;~.)/:~|:cp
  obs=.|.freq (<"1 index)} (dim#cuts)$0
)

test=: 3 : 0
  'data cuts t df'=.y.
  'dim n'=.$data
  s=._4.2[N=.100000[t=.1 
  e=.n*(%(+/^:2))gof (simf s,N,t,df);cuts
  o=.gof gmge;cuts
  chisq=.(+/^:2 e%~*:o-e)
  nu=.*:<:cuts
  (chisq-nu)%%:+:nu
  (+/+/e%~*:o-e)
)
  



NB.  read data

'optime claim injcode accmonth repmonth finmonth legal'=:iread dir,'suncorp'
AFG=:iread dir,'AFG'
Y=:>.&0 <.&1 {: >_2{. iread dir,'ozmortality'  NB. age x time

quantile=:0.5&+@:rank % #

sim=: 3 : 0
  n=.y.
  z1=.invz u1=.runif n
  lambda=. %:-<:*:u1
  z2=.(u1*z1) + lambda * rnorm n
  u2=.pnorm z2
  'dot' plot u1;u2
)



TCE=: mean\.@{:     NB.  Tail conditional expectation
TCSD=: sd\.@{:      NB.  Tail conditional std dev

tcp=:3 : 0
  u=.{.{.y.
  mean u<:{:"1 y.
)




NB.  Routine to pick u quantiles from vector x
NB.  Use linear interpolation

upick=: 4 : 0"0 1
u=.(i.%])x.
'n s'=:(#;/:~)y.
r=.<.u*n
f=.(u*n)-r
({.+f&*@(-/))((>:,:])r){(,{:)s
)

invz=: 3 : 0  
    p=.y.
NB.   Coefficients in rational approximations.
   a_1 =. -3.969683028665376e1
   a_2 =.  2.209460984245205e2
   a_3 =. -2.759285104469687e2
   a_4 =.  1.383577518672690e2
   a_5 =. -3.066479806614716e1
   a_6 =.  2.506628277459239e0

   b_1 =. -5.447609879822406e1
   b_2 =.  1.615858368580409e2
   b_3 =. -1.556989798598866e2
   b_4 =.  6.680131188771972e1
   b_5 =. -1.328068155288572e1

   c_1 =. -7.784894002430293e_3
   c_2 =. -3.223964580411365e_1
   c_3 =. -2.400758277161838e0
   c_4 =. -2.549732539343734e0
   c_5 =.  4.374664141464968e0
   c_6 =.  2.938163982698783e0

   d_1 =.  7.784695709041462e_3
   d_2 =.  3.224671290700398e_1
   d_3 =.  2.445134137142996e0
   d_4 =.  3.754408661907416e0

NB.   Define break-points.

   p_low  =. 0.02425
   p_high =. 1 - p_low

   if.   (p <: p_low) +. (p>:p_high) do.
      if. p>p_high do. pd =. 1-p else. pd=.p end.
      q =. %:(-2*ln(pd))
      x =. (c_6+q*c_5+q*c_4+q*c_3+q*c_2+q*c_1) % 1+q*d_4+q*d_3+q*d_2+q*d_1
      if. p>:p_high do. x=.-x end.
   else.
      q =. p - 0.5
      r =. q*q
      x =. (q*a_6+r*a_5+r*a_4+r*a_3+r*a_2+r*a_1)%1+r*b_5+r*b_4+r*b_3+r*b_2+r*b_1
   end.
)

normit=: invz@quantile

NB.  --- Now the real stuff

F_star=:quantile"1       NB. quantile transformation of rows 
N_star=:invz@F_star      NB. normit transform of rows

qcor=: 3 : 0   NB. compute correlations in each quantile quadrant
  m=.y.
  q=.|:>quantile L:0  (setd;claim) /:L:0 setd
  locsd=.<.m*q
  r=.locsd ({.@{:@cor@:((invz"1)&.|:))/. q
  qc=.r (<"1 ~. locsd)} (m,m) $ 0
  'surface' plot qc
  qc
)

pwr=: 3 : 0
  r=.0.9999
  data=.setd,:claim
  n=.{:$data
  z=. (sqrt 2 2 $ 1,r,r,1)mp nrand 2,n
  u=. <.n*1&chiprob *:z
  u{"1 data
NB. 'dot' plot ;/u
)

ct=: 3 : 0
   d=./:~|:quantile"1 [10{."1 setd,:claim
   xley=:*/@:<:    NB. x. less than or equal to y. in all components  
   nl=.+/@(xley"1 {:)  NB. number less than or equal to last
   p=.(#d)%~ nl\d
  ((,<./,:*/)@|:d),p
)

utd=: 3 : 0   NB.  Upper tail dependence
   d=:/:~|:normit"1 y.
   xgey=.*/@:>:    NB.  All x. greater than y. 
   phi=.mean@:(xgey"1 {.@{.)
   rho=.{.@{:@cor
   psi=.{:@mean
   ({."1 d),|:((psi,rho,phi)\.d)
)

bt=: 3 : 0
  a=.y.
  u=.urand 10000
  alpha=.1-a*ln 1-u [ beta=.1-a*ln u
  v=.,betarand"1  alpha,.beta,.1
  'dot' plot ;/x=:u,:v
  return.
  label_name.
  u=. 0.7
  plot hist betarand (-a*ln 1-u),(-b*ln u),10000 
)

clay=: 3 : 0  NB. Simulate from bivariate clayton
  'delta n'=.y.
  gamma=.gammarand (%delta),n
  v=. urand n,d=.2
  u=.|:(1-(ln v) % gamma)^-%delta
)

aclay=: 3 : 0
   'delta alpha beta n'=.y.
   v=.clay delta,n
   u=.urand 2,n
   >./"2  (v^%alpha,beta),:"1 u^%1-alpha,beta
)
   

frac=: i.@>: % ]
span=: (<./,>./)@,
support=:frac@[ ({.@] + [ * ({:-{.)@])  span@]  
Ex=:mean"_1@:(] <:"0 1 support)
R=:mean@:(*/~"1) @:|:@:(normit"1)
rho=:>@(mean L:0)@</.@|.@R 

NB. Age-at-death modelling

NB. U=: 1- */\"1 [ 1- |: Y   NB.  F_*(i) at i=i.101
NB. Z=: invz U


optime=:optime%100
logcl=:ln claim             NB. log of claim size
setd=:finmonth-accmonth     NB.  settlement delay
repd=:repmonth-accmonth     NB.  reporting delay


logit=: ln@(% -@<:)

intact=.|:@(~.)@:(,/)@:(*"1"1 _)&|: NB.  Interact the columns of two matrices.

powp=: (ln@[)`(<:@^%])@.(0&~:@])



gres=: 3 : 0   NB.  Resample from the given marginals
    normits=.normit"1 y.  NB. y. is 2xn  data.
    rho=.{:{.cor|:normits
    distns=./:~"1 y.
    A=.sqrt 2 2 $ 1,rho,rho,1
    index=. <.(#|:distns) * nprob A mp nrand 2,n
    index {" 1 distns
)  

spl=: 3 : 0
   m=.y.
   sc=.normit setd,.claim
   'stick' plot ;/|:(~.,. (#/.~))@:<.@(m&*) sc
)  

   

Psi=: 4 : 0  NB.  Not sure if still correct
  rho=.x.
  'x y'=. y.
  (1-rho^2)%~^-((rho^2*(x^2+y^2))-2*rho*x*y)%2*1-rho^2
)

NB.  Rejection sampling with Gaussian copula

gcop=: 3 : 0   NB.  This has mistake since 
   rho=.0.8[n=.100000
   lnu=.ln urand n
   'normitx normity'=.nrand 2,n
   rhs=. (-*:normitx-rho*normity)%+:1-*:rho
   (lnu<rhs)# normitx,.normity
)

gcop2=: 3 : 0
   rho=.0.8[n=.100
   A=.sqrt 2 2 $ 1,rho,rho,1
   'normitd normitc'=.A mp nrand 2,n
   c=.claim invnormit"1 0 normitc
   d=.setd invnormit"1 0 normitd
)
   



plp=: 3 : 0
  rho=.y.
  unigrid=.(>:@i.%>:)100
  grid=. (<@,)"0 / ~ unigrid
  vals=. > rho Psi L:0 grid
  pd 'type surface'
  pd 'edgecolor 64 64 64'
  pd 'backcolor lightgray'
  pd vals
  pd 'show'
)

cop=: 3 : 0
    's c'=.normit"1 y.
    'dot' plot s;c
return.

    'dot' plot 500 sd\ s /: c

return.  
   'dot' plot ;/normit"1 y.
return.
    a=. 1 100 1000 mean\"0 1 c/: s
   'dot' plot a, 1 100 1000 sd\"0 1 c/: s
)
       


    

dplot=: 'point'&plot        NB.  dot plot

optime2=:quantile setd          NB.  Note: differs from optime because of days

                   NB.  seven variables in the following order
                   NB.  operational time, payout, injury code, acc2, rep2,  legalrep
                   NB.  injury codes 1 and 2 are low, 3 4 5 6 are higher.  9 is unknown
                   NB.  accmonth:  accident date converted to month number  (nominal)
                   NB.  repmonth:  report date converted to month number  (1989 jan = 1)
                   NB.  finmonth:  finalization month
                   NB.  legal:  legal representation 1=yes, 0=no

opt=:'labelfont Arial 40';'captionfont Arial 40';'border 1';'backcolor white'
opt=:opt,'border 1';'titlefont Ariel 55';'ticminor 0'
opt=:opt,<'symbolfont Arial 40'

monthplot=: 3 : 0
   toyear=.  (1991&+)@(%&12) 
   pd 'reset'
   pd 'new'
   pd 'backcolor white'

   pd 'new 20 520 460 460';'type point';opt
   pd 'title accident month'
   pd 'xcaption  year'
   pd 'xticpos 1992 1996 2000'
   pd 'ycaption log claim'
   pd 'yticpos 0 4 8 12 16'
   pd (toyear accmonth);ln claim

   pd 'new 20 20 460 460';'type point';opt   
   pd 'title finalization month'
   pd 'xcaption year'
   pd 'xticpos 1992 1996 2000'
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd (toyear finmonth);ln claim

   pd 'new 520 520 460 460';'type point';opt
   pd opt 
   pd 'title report month'
   pd 'xcaption year'
   pd 'xticpos 1992 1996 2000'
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd (toyear repmonth);ln claim

   pd 'new 520 20 460 460';'type point';opt
   pd 'title injury code'
   pd 'xcaption code'
NB.   pd 'xticpos '
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd (0,injcode);0,ln claim

    pd 'show'
    pd 'clip'
)
   
delplot=: 3 : 0 
   pd 'reset'
   pd 'new'
   pd 'backcolor white'

   pd 'new 20 520 460 460';opt
   pd 'title reporting delay'
   pd 'xcaption  months'
NB.   pd 'xticpos 0 50 100'
   pd 'ycaption log claim'
   pd 'yticpos 0 4 8 12 16'

    pd (repmonth-accmonth);ln claim

   pd 'new 20 20 460 460';opt   
   pd 'title reporting delay'
   pd 'xcaption quantentile'
NB.   pd 'xticpos 0 0.50 1'
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd (quantile repmonth-accmonth);ln claim

   pd 'new 520 520 460 460'
   pd opt 
   pd 'title settlement delay'
   pd 'xcaption months'
NB.   pd 'xticpos 0 50 100'
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd (finmonth-repmonth);ln claim

   pd 'new 520 20 460 460';opt
   pd 'title settlement delay'
   pd 'xcaption quantentile'
NB.   pd 'xticpos '
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd optime;ln claim

    pd 'show'
    pd 'clip'
)
   
pp=: 3 : 0    NB.  pp plot for normal,gamma
  data=.y.
  'mu sigma'=.(mean,sd) data
  nu=.*: mu%sigma  NB.  Shape of gamma
  pdata=. quantile data
  pnorm=. nprob (data-mu)%sigma
  pgamma=.(+:nu) chiprob (+:nu)*data % mu
  pdata;pnorm,pgamma,:pdata
)
   
regress=: 3 : 0
  'y X'=.y.
  'n p'=.$X
  X=.1,.X
  b=.,y%.X
  e=.y-yhat=.X mp b
  s=.%:(tp e)%n-p+1
  covb=: (*:s)*%. tp X
  seb=.%:getd covb
  halfcovpred=.(sqrt covb) sqr X
  sepred=.%:+/*:halfcovpred
  pred=. X mp b
  R2=.1- e %&(*:@sd)   y
  ,.(e,:sepred);(b,:seb);R2
)

dummy=:|:@}.@=  NB. generate dummy variables from categorical: first level is base case

   


 fr=:|:@:".@('m'&fread)@(dir&,)         NB.  Read data into 2 x something array

figa=: 3 : 0
   pd 'reset'
   pd 'new'
opt=:'labelfont Arial 40';'captionfont Arial 40';'border 1';'backcolor white'
opt=:opt,'border 1';'titlefont Ariel 55';'ticminor 0'

ny=.#y.
   mu=.mean y.
   sig=.sd y.
   sy=./:~y.
   max=.{:sy
   min=.{.sy
   range=.max-min
   nbins=.>.%:ny
   binwidth=.range%nbins
   cuts=. min+binwidth*>:i.nbins
   freq=.(({. , }. - }:)+/ sy<:/cuts)%(ny*binwidth)   NB.  Note integral is 1.
   midpoints=.cuts-0.5*binwidth
   normal=.(%sig*%:2p1)*(^--:*:(midpoints-mu)%sig)                               
   shape=.*:mu%sig
   scale=.mu%*:sig
   gamma=.((scale^shape)*(midpoints^<:shape)*^-scale*midpoints)%!<:shape

   pd opt
   pd 'xcaption log claim size'
NB.  pd 'xticpos 20 40 60 80 100'
   
   pd 'ycaption frequency'
NB.  pd 'yticpos 0 0.1 0.2 0.30 0.4'
    pd  freq,normal,:gamma
  
    pd 'show'
    pd 'clip'
*:mu%sig
)

logvsop=: 3 : 0
   pd 'reset'
   pd 'new'
opt=:'labelfont Arial 40';'captionfont Arial 40';'border 1';'backcolor white'
opt=:opt,'border 1';'titlefont Ariel 55';'ticminor 0'
   pd opt
   pd 'type point'
   pd 'xcaption operational time'
NB.  pd 'xticpos 20 40 60 80 100'
   
   pd 'ycaption log claim'
NB.  pd 'yticpos 0 0.1 0.2 0.30 0.4'
    pd  optime;ln claim
  
    pd 'show'
    pd 'clip'
)

boxplot=: 3 : 0
   'inj cl' =.sort injcode,:claim
   'n code mu x stdev'=.|: inj (#,mean,sd)/. inj,.ln cl
   low=. mu-+:stdev%%:n
   hi=.mu++:stdev%%:n

      pd 'reset'
   pd 'new'
opt=:'labelfont Arial 40';'captionfont Arial 40';'border 1';'backcolor white'
opt=:opt,'border 1';'titlefont Ariel 55';'ticminor 0'
   pd opt
   pd 'type point'
   pd 'xcaption injury code'
NB.  pd 'xticpos 20 40 60 80 100'
   
   pd 'ycaption log claim'
NB.  pd 'yticpos 0 0.1 0.2 0.30 0.4'
NB.    pd  low,mu,:hi
    pd (0,inj);(0,ln cl)
  
    pd 'show'
    pd 'clip'
)

   

NB.  Sort on the basis of first row

sort=: |:@(|: /: {.)

tr=: 2 : '({. y.) ,: x. @ {: y.'   NB.  Adverb to define any transformation 


NB.  plot data

pl=: ('dot'&plot)@:(<"1)                   NB.  dot plot of 2 x ... array with x sorted.

NB.  residual from regression together with x variable  

res=: {.,:(({: - (1&,"0)@{. mp {: regn {.) % sd@{:)          NB.  x and residuals from linear regression

pow=: <:@^ % ]                                               NB.  power transform

normq=:(1&-)@:(norm"0)@:zsc@:(/:~)                           NB.  Ordered normal quantiles
normqq=:plot@:({.;])@:((>:@:i.%])@#@{.,normq"1)              NB.  Normal qq plot of y.


chisc=: +:@:(* mean % *:@sd)                                NB.  Chi scores for gamma data
chidf=: +:@*:@(mean % sd)                                   NB.  Chi df for gamma data
chipr=: chiprob@,                                            NB.  lower chi probs with x.=df, y.=argument
gammaq=:(chidf chipr"0 ])@:chisc@:(/:~)                      NB.  Ordered gamma quantiles
gammaqq=:plot@:({.;])@:((>:@:i.%])@#@{.,gammaq"1)
 
plotqq=: plot@:({.;])@:((>:@:i.%])@#,gammaq,:(normq@:^.))    NB.  plot log normal and gamma qq, where y.=data

hist=: 3 : 0                                                 NB.  histogram plot
   ny=.#y.
   mu=.mean y.
   sig=.sd y.
   sy=./:~y.
   max=.{:sy
   min=.{.sy
   range=.max-min
   nbins=.>.%:ny
   binwidth=.range%nbins
   cuts=. min+binwidth*>:i.nbins
   freq=.(({. , }. - }:)+/ sy<:/cuts)%(ny*binwidth)   NB.  Note integral is 1.
   midpoints=.cuts-0.5*binwidth
   normal=.(%sig*%:2p1)*(^--:*:(midpoints-mu)%sig)                               
   shape=.*:mu%sig
   scale=.mu%*:sig
   gamma=.((scale^shape)*(midpoints^<:shape)*^-scale*midpoints)%!<:shape
   (freq);~midpoints 
)

claimhist=: 3 : 0
   pd 'reset'
   pd 'new'
   pd 'backcolor white'

   pd 'new 20 520 460 460';opt
   pd 'title claim frequency'
   pd 'xcaption  log claim size'
NB.   pd 'xticpos 0 50 100'
   pd 'ycaption frequency'
NB.   pd 'yticpos 0 4 8 12 16'
   pd hist claim powp 0

   pd 'new 20 20 460 460';opt   
   pd 'title pp-plot'
   pd 'xcaption quantentile'
NB.   pd 'xticpos 0 0.50 1'
   pd 'ycaption quantentile'
NB.   pd 'yticpos 0 4  8  12 16'
   pd pp claim powp 0

   pd 'new 520 520 460 460'
   pd opt 
   pd 'title settlement delay'
   pd 'xcaption months'
NB.   pd 'xticpos 0 50 100'
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd (finmonth-repmonth);ln claim

   pd 'new 520 20 460 460';opt
   pd 'title settlement delay'
   pd 'xcaption quantentile'
NB.   pd 'xticpos '
   pd 'ycaption log claim'
   pd 'yticpos 0 4  8  12 16'
   pd optime;ln claim

    pd 'show'
    pd 'clip'
)


initplot=: 3 : 0
   pd 'reset'
   pd 'new'
   pd 'backcolor white'

   pd 'new 20 520 460 460';opt
   pd 'title claim size distribution'
   pd 'xcaption  log claim size'
NB.   pd 'xticpos 0 50 100'
   pd 'ycaption frequency'
NB.   pd 'yticpos 0 4 8 12 16'
   pd 'type point'
   pd hist claim powp 0

   pd 'new 520 520 460 460'
   pd opt 
   pd 'title legal representation'
   pd 'xcaption legal representation'
NB.   pd 'xticpos 0 50 100'
   pd 'ycaption frequency/1000'
   pd 'xlabel 0 1'
   pd 'type sbar'
   pd %&1000 {."1 >{:"1 >freq legal

   pd 'new 20 20 460 460';opt   
   pd 'title injury level'
   pd 'xcaption injury code'
   pd 'xticpos 1 2 3 4 5 6 9'
   pd 'ycaption frequency/1000'
NB.   pd 'yticpos 0 4  8  12 16'
   pd 'type sbar'
   pd %&1000 {."1 >{:"1 >freq injcode


   pd 'new 520 20 460 460';opt
   pd 'title settlement delay'
   pd 'xcaption months'
   pd 'ycaption proportion'
NB.   pd 'yticpos 0 4  8  12 16'
NB.   pd 'yticpos 0 4  8  12 16'
   pd 'type point'
   pd ((i.%])#setd);~/:~setd

    pd 'show';'clip'
)

confband=:{. ([,+,:-)+:@{:

finplot=: 3 : 0
   xo=./:~finmonth
   yo=.(ln claim)/:finmonth
   'yho syho'=:>{.regress yo;,.xo
   xm=.~.xo
   ym=.|:xo (mean,sd,#)/.yo
   yhm=.~.yho
   syhm=.~.syho

   pd 'reset';'new';'backcolor white';opt

   pd 'new 20 520 460 460';opt
   pd 'title log claim'
   pd 'xcaption finalization month'
   pd 'ycaption log claim'
   pd 'yticpos 8 9 10'
   pd  xm;(confband yhm,:syhm),confband ({.,:1&{%%:@{:)ym

   pd 'new 520 520 460 460';opt 
   pd 'title claim'
   pd 'xcaption finalization month'
   pd 'ycaption claim/1000'
   pd  xm;%&1000^( confband yhm,:syhm),confband ({.,:1&{%%:@{:)ym

   pd 'show';'clip'
)  


compare=: 3 : 0
NB. s=.|:(|:y.)/:0{y.
   'optime total legrep'=:   <"1 (0 1 6{ y.)
   optime=.optime%100
   prec=.%1.5167e_4+(_3.07638e_4*optime)+(1.67989e_4**:optime)+(_2.81629e_6*legrep)
   plog=.^8.2446+(4.2016*optime)+(_0.7290**:optime)+0.2397*legrep
   pnorm=.^(7.6329+(4.0054*optime)+(_0.8308**:optime)+0.4826*legrep) +-:*:1.1370
   pid=.2479.98+(10010.65*optime)+(95906.8**:optime)+3270.77*legrep
   'dot' plot optime;prec,pnorm,:plog
) 

cov=.6.336e_5 _0.0002 0.0002 _1.592e_5
cov=.cov,_0.0002  0.0013  _0.0012 _1.797e_6 
cov=.cov, 0.0002 _0.0012 0.0013 _7.618e_7
cov=.4 4 $cov,_1.592e_5 _1.797e_6 7.618e_7 2.597e_5

beta=:2.6337 1.2946 _0.2118 0.1217
xt=:1 0.6 0.36 1


NB.   Neural Network stuff

nn=:;`(%.@>:@^@-@(>@{: (+/ . *) (1&,)@nn@}:))@.((1&<)@#)   NB.  Compute output from NN inputs

ann=: 3 : 0                       
  Y=.0 1 1 0                      NB.  training set outputs
  X=.0 0,0 1,1 0,:1 1             NB.  training set inputs
  act=. %.@>:@^@-                 NB.  logistic activation
  top=.1 2 ,: 2 2                 NB.  network topology
  topb=. top +"1 (0 1)            NB.  network topology including biases
  Xt=.< |: X                      NB.  boxed transposed X 
  nopar=.+/  */"1 topb            NB.  number of pars including biases
  lsq=: +/@:*:@(Y&-)@,@nn@(,&Xt)@(topb&tpar)    NB.  SSE for given pars
  minres=: lsq mins init    NB.  SSE simplex minimization iterations 
  topb tpar > 2{ {: minres  NB.  network weights
)
                               
tpar=: 4 : 0         NB.  x. is topb , y. is pars
    z=.  (+/\) (*/"1) x.
    cuts=. +/@(0&,@}: =/ i.@{:)z
    boxed=. cuts <;.1 y.
    x. <@($ >)"1 0 boxed
)
  
nopt=: 23.70 _17.5 _17.5 7.9 _14 16 1.2 19 _9

init=: 10 _10 _10 10 _10 10 3 20 _10

opt2=: 23.70 _17.5 _17.5 7.9 _14 _9 1.2 19 6


ta=: 3 : 0
  Ex=. 2 0.79
  p2sd=. 4.83 2.02
  sd=.-:p2sd-Ex
  stdi=.(%&sd)@:(-&Ex)
  v=.Ex,.sd
  v=.v,.sd%Ex
  v=.v,. stdi 4.74 1.86
  v=.v,. stdi 5.09 2.13 
  v=.v,. stdi 4.62 1.92
  v=.v,. stdi 5.79 2.58
)

taylor=: 3 : 0
  u=. 0.6 0.9 0.99 1
  ph=. 0.75 1.12 1.87 4.67
  X =. u^"0 1[ 1 2 3 4
  ph%.X
)