load '~user/piet.ijs'
load 'stats/r/rserve csv dates tables/tara'

figdir=:'/users/pietdejong/documents/research/srisk/figures/'
datadir=:'~/documents/research/srisk/data/'

initR=:0 : 0  
NB.  run  R and assuming .First=function(){library(Rserve)}
NB.  in R start up Rserve
NB.  Rserve(args="--no-save")
NB.  Open R in J
   Ropen''
NB.  Load relevant libraries
   Rlib L:0 'rmgarch';'xtable'   NB.  xtable required for latex output of Tables
)

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
Ex=:+/@:(* %"1 +/@:])                 NB.   (mxT) Ex mxT where y is weighting
put=:*>&0                            

lo=:_15"0
ln=: (lo`^.@.(>&0))"0                        NB. log or zero

ct=: 3 : 0   NB. predicted capital shortfall distribution for a single firm
  'd w r'=:"."1 y,"1 0 'dw '
  nu=.(ind#i.#r) (dcc@:{.)"0 2 r,.asx        NB. Tx2xN predictive distributions
  l=.(-/"1 ln ind#d,.w)+logit k              NB. T  adjusted log--leverage
  (l,|:ind#d,.w);1 2 0|:nu%100               NB. 3xT;2xNxT
)

dcc=: 3 : 0        NB. library(rmgarch)
  'N h'=.6000 22['y' Rset r=.y   
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
  vepsf=.(N,h,2)$,(?(N*h)##veps){veps[9!:1[3^10  NB. choose past veps given random seeds
  rhot=.tocor"2 Q['rf sigt Q'=.(N&#)@:,: L:0 ('',.'');({:sig);2 2$1,rhot,1,~rhot=.{:rho    
  for_t. i.h do.
    'veit emt'=.|:t{"2 vepsf
    et=.((rhot*emt)+(%:1-rhot^2)*veit),.emt
    rf=.rf,"2 1 mu+"1 sigt*et
    sigt=.%:o+"1(*:sigt)*b+"1(a+"1 g*"1(et<0))*et^2              NB. Nx2
    rhot=.tocor"2 Q=.(Qbar*1-aj+bj)+"2 (aj*(*/~)"1 et)+"2 bj*Q   NB. N, Nx2x2
  end.
  |:(/:{:"1) +/"2 rf  NB. Nx2 pairs (nu_i,nu_m) sorted  according nu_m
)


init=: 3 : 0
  dat=:readcsv datadir,'cifrdatdaily.csv'
  ({.dat)=:|: (".)L:0  '-_'&charsub L:0 }.dat
NB. 
  banks=:'cba','anz','nab','wbc','mqg','boq','ben',:'aba'
NB.
  k=:0.08['year month day'=:|:ymd=:todate dayno=.1&todayno date
  dyear=:(dayno - 1&todayno 20000101)%365.12
  prices=:"."1 banks[debt=:"."1 banks ,"1 0 'd'
  newmonth=:1,0~:dif month              NB. assumes first obs is new month
  ind=:(year>:2003)*.newmonth           NB. indicators of new months starting 2003
  yrmo=:(yr-2000)+(mo-1)%12['dyr yr mo da'=:|:ind#(ymd tsdiff 2000 1 1),.ymd
NB. estimation and simulation
  if. -.fexist datadir,'data.j' do. (|:ct"1 banks) iwrite datadir,'data.j' end. 
  'ldw_it pn'=:iread datadir,'data.j'   NB. (mx3xT);mx2xNxT  -- m=#banks
NB.
  'l_it d_it w_it'=:1 0 2|:ldw_it       NB. 3xmxT
  Ex_d=:Ex&d_it                         NB. debt weighted averaging
  Ex_w=:Ex&w_it                         NB. equity weighted averaging
  'nu_oit temp'=:1 2 0 3|:pn            NB. NxmxT;NxmxT 
  nu_omt=:{."2 temp                     NB. NxT 
  phiu_omt=:phi (P"1)&.|:nu_omt         NB. NxT <- TxN <- NxT 
  sig_phi=:E sd phiu_omt                NB. ''  
  Es=:E@:(*"1&phiu_omt)                 NB. Stressed expectation
  r_oit=:nu_oit-"2 l_it                 NB. NxmxT -- adjusted returns
  p_oit=:put 1-^r_oit                   NB. NxmxT -- S^+/(kd) -- put
  p_ot=:Ex_d"2 p_oit                    NB. NxT -- debt weighted put
  mu_it=:E p_oit                        NB. mxT -- BASRISK/(kd)
  s_it=:(Es p_oit)-mu_it                NB. mxT -- PSIRISK/(kd)
  q_it=:E nu_oit <:"2 l_it              NB. mxT -- P(S^+>0)
  e_it=:put 1- Es^r_oit                 NB. mxT -- Engle put
  pi_it=:(%"1+/)d_it                    NB. mxt -- debt proportions and sum
  'qb_t mub_t sb_t'=:Ex_d"2 q_it,mu_it,:s_it NB. T;T -- debt weighted
  nu_ot=:^. Ex_w"2 ^ nu_oit             NB. NxT -- pooled return
  l_t=:^. Ex_w ^ l_it                   NB. T  -- pooled adj log lev
  p_ot=:put 1-^nu_ot-"1 l_t             NB. NxT -- pooled put
  mu_t=:E p_ot                          NB. T -- pooled BASRISK
  s_t=:(Es p_ot)-mu_t                   NB. T -- pooled PSIRISK
  q_t=:E nu_ot <:"1 l_t                 NB. T -- pooled default prob
  S=:k*d_it*"2 (1-^nu_oit-"2 l_it)      NB. NxmxT -- S
  {.dat
)

major=:'CBA ANZ NAB WBC'[minor=:'MQG BOQ BEN ABA'

key=:3 : 0
   pd 'key ',y
   pd 'keyfont Arial 6;keypos top left;keystyle thin'
)
   
figprices=: 3 : 0
   pd 'reset'
   opt=. 'yticpos 0 3 6 9;xticpos 1 3 5 7 9 11 13 15'
   opt=. opt,';xcaption years since 2000'
   pd 'sub 2 2'
   pd 'new'
   pd opt
   key major
   pd dyear;*/\"1 ^%&100 [ 4{.prices
   pd 'new'
   pd opt
   key minor
   pd dyear;*/\"1 ^%&100 [_4{.prices
   opt=. 'yticpos -1 -0.5 0 0.5 1;xticpos 3 4 5 6 7 8 9 10 11 12 13 14 15'
   opt=. opt,';xcaption years since 2000'
   pd 'new'
   pd opt
   key major   
   pd dyr;4{.l_it
   pd 'new'
   pd opt
   key minor
   pd dyr;_4{.l_it
   pd 'pdf 350 300 ',figdir,'prices2.pdf'
)


NB. figprices''

figbloglev=: 3 : 0
   pd 'reset'
   opt=. 'yticpos -1 -0.5 0 0.5 1;xticpos 3 4 5 6 7 8 9 10 11 12 13 14 15'
   opt=. opt,';xcaption years since 2000;axes 1'
   pd 'sub 1 2'
   pd 'new'
   pd opt
   key 'CBA ANZ NAB WBC'    
   pd dyr;4{.l_it
   pd 'new'
   pd opt
   key 'MQG BOQ BEN ABA'
   pd dyr;_4{.l_it
   pd 'pdf 350 150 ',figdir,'bloglev.pdf'
)

NB. figbloglev'' 

figCBA=: 3 : 0
   months=.'first secnd'=.72 142             NB. jan09 nov14 
   'l_it d_it w_it'=.{.ldw_it                NB. assumes CBA is first bank
   'nu_om1 nu_om2'=.|:months {"1 nu_omt      NB. market simulations for relevant months
   'nu_o11 nu_o12'=:|:months {"1 {."2 nu_oit NB. assumes CBA is first bank
   'l_11 l_12'=.months{l_it
   'S_o11 S_o12'=.1-^(nu_o11,:nu_o12)-"1 0 l_11,l_12
   'cbaprice asxprice'=.((*/\)@:>:@:(%&100))"1 cba,:asx
   'nu_it nu_mt'=.|:dif ln ind#cbaprice,.asxprice 
   pd 'reset'
   opt=.'xticpos -0.8 -0.6 -0.4 -0.2 0 0.2'
   opt=.opt,';yticpos -1.5 -1 -0.5 0 0.5 1'
   opt=.opt,';xcaption market return'
   opt=.opt,';ycaption S/(kd);axes 1'
   pd 'sub 1 2'
   pd 'new'
   pd  opt
   pd 'title January 2009'
   pd 'type dot;color red'
   pd nu_om1;S_o11
   pd 'type dot;color black;pensize 2'
   pd 2&# L:0 (first{nu_mt);1-^first{nu_it-l_11
   pd 'new;type dot;pensize 1;color red'
   pd opt
   pd 'title December 2014'
   pd nu_om2;S_o12
   pd 'color black;type dot;pensize 2'
   pd  2&# L:0 (secnd{nu_mt);xx=:1-^secnd{nu_it-l_12
   pd 'pdf 350 150 ',figdir,'figCBA.pdf'
)

NB. figCBA''

figdefault=: 3 : 0
   pd 'reset'
   opt=.'xticpos 3 5 7 9 11 13 15'
   pd 'sub 2 2'
   pd 'new'
   pd opt
   pd 'yticpos 0   0.2 0.4 0.6 0.8'
   pd yrmo;4{.mu_it
   pd 'new'
   pd opt
   pd 'yticpos 0 0.2 0.4 0.6 0.8'
   pd yrmo;_4{.mu_it
   pd 'new'
   key major
   pd opt
   pd 'yticpos 0  0.1  0.2 0.3 '
   pd yrmo;0 1 2 3{s_it
   pd 'new'
   key minor
   pd opt
   pd 'yticpos 0 0.1   0.2 0.3 '
   pd yrmo;4 5 6 7{s_it
   pd 'pdf 350 250 ',figdir,'default.pdf'
)

NB. figdefault ''

figPut=: 3 : 0
   'EPdkd ESg0'=.E L:0 (put S%"2 k*d_it);S>0
   pd 'reset'
   opt=. 'xticpos 4 6 8 10  12 14;xcaption years since 2000'
   pd 'sub 2 2'
   pd 'new;type line;yticpos  0  0.2 0.4  0.6 0.8'
   pd 'ycaption E(S+)/kd ;xticpos 4  6  8  10  12  14'
   key major 
   pd dyr;4{.EPdkd
   pd 'new;type line;ycaption ;yticpos  0  0.2 0.4  0.6 0.8;xticpos 4  6  8  10  12  14'
   key minor
   pd dyr;_4{.EPdkd
   pd 'new;type line'
   pd 'yticpos 0 0.2 0.4 0.6 0.8 1'
   pd opt
   pd 'ycaption P(S>0)'
   pd dyr;4{.ESg0
   pd 'new;type line'
   pd 'ycaption ;yticpos 0 0.2 0.4 0.6 0.8 1'
   pd opt
   pd dyr;_4{.ESg0
   pd 'pdf 350 250 ',figdir,'figPut.pdf'
)

NB. figPut''

figStress=: 3 : 0
   'Spdkd Sg0'=.(put S%"2 k*d_it);S>0
   Estressed=. |:>E L:0 (nu_omt<_0.1) (<@:#"1)&:|: Spdkd
   inc=.Estressed - E Spdkd
   pd 'reset'
   opt=. 'type line;xticpos 4 6 8 10  12 14;xcaption years since 2000;yticpos 0 0.100 0.200 0.300 0.400 0.5 0.6' 
   pd 'sub 2 2'
   pd 'new'
   pd opt
   pd 'ycaption baseline risk'
   pd dyr;4{.E Spdkd
   key major
   pd 'new'
   key minor
   pd opt
   pd dyr;_4{.E Spdkd
   pd 'new'
   pd opt
   pd 'ycaption systemic risk'
   key major 
   pd dyr;4{.inc
   pd 'new'
   pd opt
   key minor
   pd dyr;_4{.inc
   pd 'pdf 350 250 ',figdir,'figStress.pdf'
)

NB. figStress''


figmuqs=: 3 : 0
   pd 'reset'
   pd 'sub 2 2'
   pd 'new;type line'
   pd 'key BRISK BRISK* BASABS'
   pd 'keyfont Arial 6'
   pd 'keypos tl'
   pd 'keystyle thin'
   pd 'xticpos 3 5 7 9 11 13 15'
   pd 'yticpos -0.1 0 0.1 0.2 0.3'
   pd  yrmo;(,-/)mub_t,:mu_t
   pd 'new;type line'
   pd 'key RISK RISK* PSIABS'
   pd 'keyfont Arial 6'
   pd 'keypos tl'
   pd 'keystyle thin'
   pd 'xticpos 3 5 7 9 11 13 15'
   pd 'yticpos -0.1 0 0.1 0.2 0.3'
   pd yrmo;(,-/)sb_t,:s_t
   pd 'new;type dot;pensize 1'
   pd 'xticpos 0 .1 .2 .3; yticpos 0 .1 .2 .3'
   pd  mub_t;mu_t
   pd 'new;type dot;pensize 1'
   pd 'xticpos 0 .1 .2 .3; yticpos 0 .1 .2 .3'
   pd  sb_t;s_t
   pd 'pdf 350 300 ',figdir,'muqs.pdf'
)

NB.  figmuqs''


figsysstress=: 3 : 0
   pd 'reset'
   pd 'sub 1 2' 
   ticpos=.''
   opt=.'yticpos 0 0.05 0.1 0.15 0.2 0.25'
   opt=.opt,';xticpos 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7'
   opt=.opt,';xcaption BRISK;ycaption PRISK'
   opt=.opt,';pensize 1;type dot'
   pd 'new'
   pd opt
   key 'CBA ANZ NAB WBC'
   pd"1  ;/ 4&{."2 mu_it,:s_it
   pd 'new'
   pd opt
   key 'MQG BOQ BEN ABA'
   pd 'keyfont Arial 6'
   pd 'keypos top right'
   pd 'keystyle thin'
   pd"1  ;/ _4&{."2 mu_it,:s_it
   pd 'pdf 350 150 ',figdir,'sysstress.pdf'
)
  
NB figsysstress''


table=: 3 : 0
   months=.'first secnd'=.72 142                      NB. jan09 nov14
   mus_it=:pi_it*mu_it %"1 mub_t                      NB. mxT -- %BASRISK
   ss_it =:pi_it*s_it %"1 sb_t                        NB. mxT -- %PSIRISK
   mat=:l_it,pi_it,mus_it,:ss_it                      NB. 4xmxT
   row=:(Ex_d l_it),(%&1e11 +/d_it),mub_t,:sb_t       NB. mxT -- Total
   mat=:mat ,"2 1 row                                 NB. 4x(m+1)xT
   row=:l_t,(%&1e11 +/d_it),mu_t,:s_t                 NB. mxT -- Pooled
   mat=:mat ,"2 1 row                                 NB. 4x(m+2)xT
   mat=:mat ,"2 1 
   mat=:100*,./months{|:mat                           NB. 2x(m+2)x4
  'mat' Rset wrs mat
   R 'rnames=c("CBA","ANZ","NAB","WBC","MQG","BOQ","BEN","ABA","wAve","Pooled")'
   R 'cnames=c("Aloglev","dperc","RISK","RISK","Aloglev","dperc","RISK","PRISK")'
   R 'dimnames(mat)=list(rnames,cnames)'
   R 'print(xtable(mat))'
)

NB. table''






  