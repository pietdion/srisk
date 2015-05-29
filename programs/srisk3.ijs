load '~user/piet.ijs'
load 'stats/r/rserve csv dates tables/tara'

figdir=:'/users/pietdejong/documents/research/srisk/figures/'

datadir=:'~/documents/research/srisk/data/'
time=: 6!:2

dplot=:'dot;pensize 2'&plot
R=:Rcmd
tocor=:((<0 1)&{)@:(diag@:%@:%:@:getd (mp mp [) ])
E=:mean
covar=:mean@(*/)@:((-mean)"1)
range=:>./-<./
stats=:mean,sd,:range

cte=:(%@] * <:)& 0.05
min12=:(12&*)@(^&11)@(1&-)
phi=:min12
Ex=:+/@:(* %"1 +/@:])                 NB.   (mxT) Ex mxT where weighting is on right
put=:*>&0                            

lo=:_15"0
ln=: (lo`^.@.(>&0))"0                     	NB. log or zero

rddat=: 3 : 0
  dat=:readcsv datadir,'cifrdatdaily.csv'
  ({.dat)=:|: (".)L:0  '-_'&charsub L:0 }.dat
  k=:0.08['year month day'=:|:ymd=:todate 1&todayno date
  dyear=:ymd tsdiff 2000 1 1   
  prices=:cba,anz,nab,wbc,mqg,boq,ben,:aba
  debt=:cbad,anzd,nabd,wbcd,mqgd,boqd,bend,:abad
  newmonth=:1,0~:dif month   NB. assumes first is new month
  ind=:(year>:2003)*.newmonth
  'dyr yr mo da'=.|:ind#dyear,.ymd
  banks=:'cba','anz','nab','wbc','mqg','boq','ben',:'aba'
NB.  'ldwit pn'=:|:(1&ct)"1 banks      NB. (mx3xT);mx2xNxT  -- m=#y
  'ldwit pn'=:iread datadir,'data.j'
  'l_it d_it w_it'=:1 0 2|:ldwit        NB. 3xmxT
  Ex_d=:Ex&d_it                         NB. debt weighted averaging
  Ex_w=:Ex&w_it                         NB. equity weighted averaging
  'nu_oit temp'=:1 2 0 3|:pn            NB. NxmxT;NxmxT 
  nu_omt=:{."2 temp                     NB. NxT
  {.dat
)

Init=: 0 : 0
   rddat''
NB.  run  R
NB.  Assume .First=function(){library(Rserve)}
NB.  Rserve(args="--no-save")
   Ropen''
   Rlib 'rmgarch'   NB.  'quantmod' 
   sr''
)

dcccore=: 4 : 0        NB. library(rmgarch)
  'S h'=.6000 22['y' Rset r=.y   
  R 'm=list(armaOrder=c(0,0),include.mean=TRUE)' 
  R 'v=list(garchOrder=c(1,1),model="gjrGARCH")'
  R 'u=ugarchspec(mean.model=m,variance.model=v,distribution.model="norm")'
  R 'c=dccspec(uspec=multispec(replicate(ncol(y),u)),dccOrder=c(1,1),distribution="mvnorm")'
  R 'fit=dccfit(c,y);model=fit@model;mfit=fit@mfit;#print(names(model));#print(names(mfit))'  
  NB.  R 'fcst=dccforecast(fit,n.ahead=22);print(fcst)'
  'mu sig rho Qbar'=.Rget L:0 'model$mu';'model$sig';'unlist(mfit$R,use.names=FALSE)';'mfit$Qbar'
  'epsi epsm'=:|:eps=.(r-mu)%sig[rho=.1{|:((#r),4)$rho
  veps=.((epsi-rho*epsm)%%:1-*:rho),.epsm  NB. Tx2
  a=.ai,am[o=.oi,om[b=.bi,bm[g=.gi,gm[mu=.mui,mum['ai bi gi mui oi am bm gm mum om aj bj'=.{:"1 Rget 'mfit$coef'
  vepsf=.(S,h,2)$,(?(S*h)##veps){veps[9!:1[3^10  NB. choose past veps given random seeds
  rhot=.tocor"2 Q['rf sigt Q'=.(S&#)@:,: L:0 ('',.'');({:sig);2 2$1,rhot,1,~rhot=.{:rho    
  for_t. i.h do.
    'veit emt'=.|:t{"2 vepsf
    et=.((rhot*emt)+(%:1-rhot^2)*veit),.emt
    rf=.rf,"2 1 mu+"1 sigt*et
    sigt=.%:o+"1(*:sigt)*b+"1(a+"1 g*"1(et<0))*et^2              NB. Sx2
    rhot=.tocor"2 Q=.(Qbar*1-aj+bj)+"2 (aj*(*/~)"1 et)+"2 bj*Q   NB. S, Sx2x2
  end.
  |:(/:{:"1) +/"2 rf  NB. Sx2 pairs (nu_i,nu_m) sorted  according nu_m
)

dcc=:0&dcccore
dccold=:1&dcccore

ct=: 4 : 0   NB. actual (x=0) or predicted (x~:0) capital shortfall for a single firm
  'd w r'=:"."1 y,"1 0 'dw '
  if. x=0 do. d,w,:r return. end.
  l=.(-/"1 ln ind#d,.w)+logit k              NB. T
  nu=.(ind#i.#r) (dcc@:{.)"0 2 r,.asx        NB. Tx2xS predictive distributions
  (l,|:ind#d,.w);1 2 0|:nu%100               NB. 2xT;2xSxT
)

dwk=:0&ct   NB. T x 3 matrix d,w,k     
pcs=:1&ct   NB. T x S predicted capital shortfalls with rows ordered accoring market percentil

NB. %~/time"1  'dcc cba,.asx',:'Rget ''unlist(mfit$R,use.names=FALSE)''[dcc cba,.asx'
NB. xxx=:pcs L:0  ;:'cba wbc nab anz mqg'

NB. xxx=:pcs 'cba'



NB. xxx=:stats"2 [ 0 1 dcc"0 2 cba,.asx 
NB. time '0 dcc cba,.asx'
NB. time '1 dcc cba,.asx'

NB. pcs 'cba'

NB. xxx=:1 ct 'cba'



sr=: 3 : 0    
  phiu_omt=:phi (P"1)&.|:nu_omt         NB. NxT <- TxN <- NxT 
  sig_phi=:E sd phiu_omt                NB. ''  
  Es=:E@:(*"1&phiu_omt)                 NB. Stressed expectation
  r_oit=:nu_oit-"2 l_it                 NB. NxmxT 
  p_oit=:put 1-^r_oit                   NB. NxmxT
  p_ot=:Ex_d"2 p_oit                    NB. NxT
  mu_it=:E p_oit                        NB. mxT
  s_it=:(Es p_oit)-mu_it                NB. mxT
  q_it=:E nu_oit <:"2 l_it              NB. mxT
  e_it=:put 1- Es^r_oit                 NB. mxT  -- engle put
  pi_it=:(%"1+/)d_it                      NB. mxt
  'qb_t mub_t sb_t'=:Ex_d"2 q_it,mu_it,:s_it NB. T;T
  nu_ot=:^. Ex_w"2 ^ nu_oit             NB. NxT
  l_t=:^. Ex_w ^ l_it                   NB. T
  p_ot=:put 1-^nu_ot-"1 l_t             NB. NxT 
  mu_t=:E p_ot                           NB. T
  s_t=:(Es p_ot)-mu_t                    NB. T
  q_t=:E nu_ot <:"1 l_t                 NB. T
  yrmo=:(yr-2000)+(mo-1)%12['yr mo da'=:|:ind#ymd
)


figprices=: 3 : 0
   pd 'reset'
   opt=. 'yticpos 0 3 6 9;xticpos 1 3 5 7 9 11 13 15'
   opt=. opt,';xcaption year-2000'
   pd 'sub 1 2'
   pd 'new'
   pd opt
   pd 'key cba anz nab wbc'
   pd 'keyfont Arial 6'
   pd 'keypos top left'
   pd 'keystyle thin' 
   pd dyear;*/\"1 ^%&100 [ 4{.prices
   pd 'new'
   pd opt
   pd 'key mqg boq ben aba'
   pd 'keyfont Arial 6'
   pd 'keypos top left'
   pd 'keystyle thin'
   pd dyear;*/\"1 ^%&100 [_4{.prices
   pd 'pdf 350 150 ',figdir,'prices.pdf'
)

NB.  figprices''


figbloglev=: 3 : 0
   pd 'reset'
   opt=. 'yticpos -1 -0.5 0 0.5 1;xticpos 3 4 5 6 7 8 9 10 11 12 13 14 15'
   opt=. opt,';xcaption year-2000'
   pd 'sub 1 2'
   pd 'new'
   pd opt
   pd 'key cba anz nab wbc'
   pd 'keyfont Arial 6'
   pd 'keypos top left'
   pd 'keystyle thin'    
   pd dyr;4{.l_it
   pd 'new'
   pd opt
   pd 'key mqg boq ben aba'
   pd 'keyfont Arial 6'
   pd 'keypos top left'
   pd 'keystyle thin'
   pd dyr;_4{.l_it
   pd 'pdf 350 150 ',figdir,'bloglev.pdf'
)

NB. figbloglev''

figsimulation=: 3 : 0
   months=.'first secnd'=.72 142   NB. jan09 nov14
   'nu_sm1 nu_sm2'=.|:months {"1 nu_smt
   'nu_s11 nu_s12'=.|:{."2 months {"1 nu_sit  NB. choose CBA
   'l_11 l_12'=.{.months {"1 l_it             NB. choose CBA
   'cbaprice asxprice'=.((*/\)@:>:@:(%&100))"1 cba,:asx
   'nu_1t nu_mt'=.|:dif ln ind#cbaprice,.asxprice 
   pd 'reset'
   opt=.'xticpos -0.8 -0.6 -0.4 -0.2 0 0.2'
   opt=.opt,';yticpos -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6'
   opt=.opt,';xcaption market log return'
   opt=.opt,';ycaption cba log return'
   pd 'sub 1 2'
   pd 'new;type dot'
   pd opt
   pd nu_sm1;nu_s11
   pd 'type line;color red' 
   pd nu_sm1;(#nu_sm1)#l_11
   pd 'type dot;color black;pensize 2'
   pd 2&# L:0 (first{nu_mt);first{nu_1t
   pd 'type dot;color yellow;pensize 2'
   pd 2&# L:0 [ 0;0
   pd 'new;type dot;pensize 1;color blue'
   pd opt
   pd nu_sm2;nu_s12
   pd 'type line;color red'
   pd nu_s12;(#nu_s12)#l_12
   pd 'color black;type dot;pensize 2'
   pd  2&# L:0 (secnd{nu_mt);secnd{nu_1t
   pd 'type dot;color yellow;pensize 2'
   pd 2&# L:0 [ 0;0
   pd 'pdf 350 150 ',figdir,'simulation.pdf'
)

NB. figsimulation''

figdefault=: 3 : 0
   pd 'reset'
   opt=.'xticpos 3 5 7 9 11 13 15'
   pd 'sub 3 2'
   pd 'new'
   pd opt
   pd 'yticpos 0 0.25 0.5 0.75 1'
   pd yrmo;4{.q_it
   pd 'new'
   pd opt
   pd 'yticpos 0 0.25 0.5 0.75 1'
   pd yrmo;_4{.q_it
   pd 'new'
   pd opt
   pd 'yticpos 0   0.2 0.4 0.6'
   pd yrmo;4{.mu_it
   pd 'new'
   pd opt
   pd 'yticpos 0   0.2 0.4 0.6'
   pd yrmo;_4{.mu_it
   pd 'new'
   pd 'key cba anz nab wbc'
   pd 'keyfont Arial 6'
   pd 'keypos top left'
   pd 'keystyle thin'
   pd opt
   pd 'yticpos 0  0.1  0.2 0.3 '
   pd yrmo;0 1 2 3{s_it
   pd 'new'
   pd 'key mqg boq ben aba'
   pd 'keyfont Arial 6'
   pd 'keypos top left'
   pd 'keystyle thin'
   pd opt
   pd 'yticpos 0 0.1   0.2 0.3 '
   pd yrmo;4 5 6 7{s_it
   pd 'pdf 350 250 ',figdir,'default.pdf'
)

NB. figdefault ''

figmuqs=: 3 : 0
   pd 'reset'
   pd 'sub 2 3'
   pd 'new;type line'
   pd 'key ave tot'
   pd 'keyfont Arial 6'
   pd 'keypos tl'
   pd 'keystyle thin'
   pd 'xticpos 3 5 7 9 11 13 15'
   pd 'yticpos 0 0.1 0.2 0.3' 
   pd  yrmo;mub_t,:mu_t
   pd 'new;type line'
   pd 'xticpos 3 5 7 9 11 13 15'
   pd 'yticpos 0 0.2 .4 .6 .8 1 '
   pd yrmo;qb_t,:q_t
   pd 'new;type line'
   pd 'xticpos 3 5 7 9 11 13 15'
   pd 'yticpos 0 0.05 0.1 .15 0.2'
   pd yrmo;sb_t,:s_t
   pd 'new;type dot;pensize 1'
   pd 'xticpos 0 .1 .2 .3; yticpos 0 0.1 0.2 0.3'
   pd  mub_t;mu_t
   pd 'new;type dot;pensize 1'
   pd 'yticpos 0 0.2 .4 .6 .8 1 ;xticpos 0 .2 .4 .6 .8 1 '
   pd  qb_t;q_t
   pd 'new;type dot;pensize 1'
   pd 'xticpos 0 .05 .1 .15 .2; yticpos 0 0.05 0.1 .15 0.2'
   pd  sb_t;s_t
   pd 'pdf 350 150 ',figdir,'muqs.pdf'
)

NB.  figmuqs''


figsysstress=: 3 : 0
   pd 'reset'
   pd 'sub 1 2' 
   ticpos=.''
   opt=.'yticpos 0 0.05 0.1 0.15 0.2 0.25'
   opt=.opt,';xticpos 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7'
   opt=.opt,';xcaption background stress;ycaption systemic stress'
   opt=.opt,';pensize 1;type dot'
   pd 'new'
   pd opt
  pd 'key cba anz nab wbc'
   pd 'keyfont Arial 6'
   pd 'keypos top right'
   pd 'keystyle thin'
   pd"1  ;/ 4&{."2 mu_it,:s_it
   pd 'new'
   pd opt
  pd 'key mqg boq ben aba'
   pd 'keyfont Arial 6'
   pd 'keypos top right'
   pd 'keystyle thin'
   pd"1  ;/ _4&{."2 mu_it,:s_it
   pd 'pdf 350 150 ',figdir,'sysstress.pdf'
)
  
figsysstress''

stop





  