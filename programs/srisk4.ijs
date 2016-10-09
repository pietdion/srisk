load '~user/piet.ijs'
load 'stats/r/rserve csv dates tables/tara'

figdir=:'/users/pietdejong/documents/research/srisk/figures/'
datadir=:'~/documents/research/srisk/data/'

initR=:0 : 0  
NB.  run  R and assuming .First=function(){library(Rserve)}
NB.  in R start up Rserve
NB.  Rserve(args="--no-save")
   Ropen''  NB.  Open R in J
   Rlib L:0 'rmgarch';'xtable'   NB.  xtable required for latex output of Tables
   init''   NB.  
)

dplot=:'dot;pensize 2'&plot
R=:Rcmd
tocor=:((<0 1)&{)@:(diag@:%@:%:@:getd (mp mp [) ])
E=:mean
covar=:mean@(*/)@:((-mean)"1)
range=:>./-<./
stats=:mean,sd,:range
logit2=:(_10"_&>.)@ logit
plus=:*>&0

cte=:(%@] * <:)& 0.05
min12=:(12&*)@(^&11)@(1&-)
phi=:min12
Wa=:+/@:(* %"1 +/@:])              NB.  Weighted ave: (mxT) Wa mxT where y is weighting

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
  'ldw_it pn'=:iread datadir,'data.j'   NB. (mx3xT);mx2xNxT  -- m=#banks  -- pn see below                                     
  'l_it d_it w_it'=:1 0 2|:ldw_it       NB. 3xmxT -- leverage, debt, equity
  'r_oit r_ot'=:({.;({."2)@:{:)1 2 0 3|:pn  NB. NxmxT;NxT -- firm and market return
  S_oit=:1-^r_oit-"2 l_it               NB. NxmxT -- shortfall
  psi_ot=:phi (P"1)&.|:r_ot             NB. NxT -- stress function
  Es=:E@:(*"1&psi_ot)                   NB. stressed expectation
  brisk=:E@plus
  srisk=:Es@plus
  strain=:srisk-brisk
  D=:(Wa&d_it)"2                        NB. debt weighted averaging
  SRISK=:D@plus@Es                      NB. Engle SRISK
  Sum=:D@:plus,:plus@:D                 NB. compute S^* and S^+
  rho=:1 :'-/ln u"2 Sum y'              NB. resilience  eg (Es,:E) rho 
  DPTR=:1 :'(plus (%,:])&u >&0) y'      NB. Default probability and Tail risk
  ABSORB=:1 :'%~/u"2 y'                 
  all=:brisk,strain,srisk, E DPTR , Es DPTR
  {.dat
)

major=:'CBA ANZ NAB WBC'[minor=:'MQG BOQ BEN ABA'

key=:3 : 0
   pd 'key ',y
   pd 'keyfont Arial 6;keypos top left;keystyle thin'
)
   

fig1=: 3 : 0
   pd 'reset'
   pd 'sub 1 2'
   pd 'new'
   opt=. 'yticpos 0 3 6 9;xticpos 1 3 5 7 9 11 13 15'
   opt=. opt,';xcaption years since 2000'
   opt=. opt,';ycaption stock price'
   pd opt
   pd 'title major banks'
   key major
   pd dyear;*/\"1 ^%&100 [ 4{.prices
   pd 'new'
   pd opt
   pd 'title minor banks'
   key minor
   pd dyear;*/\"1 ^%&100 [ _4{.prices
   pd 'pdf 350 150 ',figdir,'fig1.pdf'
)

NB. fig1''


fig2=: 3 : 0
   pd 'reset'
   pd 'sub 1 2'
   pd 'new'
   opt=. 'yticpos -1 -0.5 0 0.5 1;xticpos 3 5 7  9  11 13  15'
   opt=. opt,';xcaption years since 2000;axes 1;ycaption default index'
   pd opt
   pd 'title major banks'
   key major    
   pd dyr;4{.l_it
   pd 'new'
   pd opt
   pd 'title minor banks'
   key minor    
   pd dyr;_4{.l_it
   pd 'pdf 350 150 ',figdir,'fig2.pdf'
)

NB. fig2''


figCBA=: 3 : 0
   months=.'first secnd'=.72 142             NB. jan09 nov14 
   'l_it d_it w_it'=.{.ldw_it                NB. assumes CBA is first bank
   'r_om1 r_om2'=.|:months {"1 r_ot      NB. market simulations for relevant months
   'r_o11 r_o12'=:|:months {"1 {."2 r_oit NB. assumes CBA is first bank
   'l_11 l_12'=.months{l_it
   'S_o11 S_o12'=.|:months {"1 {."2 S_oit
   'cbaprice asxprice'=.((*/\)@:>:@:(%&100))"1 cba,:asx
   'r_it r_mt'=.|:dif ln ind#cbaprice,.asxprice 
   pd 'reset'
   opt=.'xticpos -0.8 -0.6 -0.4 -0.2 0 0.2'
   opt=.opt,';yticpos -1.5 -1 -0.5 0 0.5 1'
   opt=.opt,';xcaption market return'
   opt=.opt,';axes 1'
   pd 'sub 1 2'
   pd 'new'
   pd  opt
   pd 'title January 2009'
   pd 'type dot;color red'
   pd r_om1;S_o11
   pd 'type dot;color black;pensize 2;ycaption S'
   pd 2&# L:0 (first{r_mt);1-^first{r_it-l_11
   pd 'new;type dot;pensize 1;color red'
   pd opt
   pd 'title December 2014'
   pd r_om2;S_o12
   pd 'color black;type dot;pensize 2;ycaption S'
   pd  2&# L:0 (secnd{r_mt);1-^secnd{r_it-l_12
   pd 'pdf 350 150 ',figdir,'figCBA.pdf'
)

NB. figCBA''


fig4=: 3 : 0
   plotdata=. (brisk;strain) S_oit
   pd 'reset'
   opt=.'xticpos 3 5 7 9 11 13 15'
   pd 'sub 2 2'
   pd 'new'
   pd opt
   key major
   pd 'yticpos 0   0.2 0.4 0.6 0.8;ycaption base risk;title major banks'
   pd yrmo; 4{. brisk S_oit
   pd 'new'
   pd opt
   key minor
   pd 'yticpos 0 0.2 0.4 0.6 0.8;ycaption base risk;title minor banks'
   pd yrmo;_4{. >{.plotdata
   pd 'new'
   pd opt
   pd 'yticpos 0 0.05 0.1 0.15 0.2 0.25;ycaption strain;xcaption years since 2000'
   pd yrmo;4{. strain S_oit
   pd 'new'
   pd opt
   pd 'yticpos 0 0.05 0.1 0.15 0.2 0.25;ycaption  strain;xcaption years since 2000'
   pd yrmo;_4{. >{:plotdata
   pd 'pdf 350 250 ',figdir,'fig4.pdf'
)

NB. fig4 ''


fig5=: 3 : 0
   plotdata=.(brisk ;"1 strain) S_oit
   pd 'reset'
   pd 'sub 1 2' 
   ticpos=.''
   opt=.'yticpos 0 0.05 0.1 0.15 0.2 0.25'
   opt=.opt,';xticpos 0  0.2  0.4 0.6 0.8'
   opt=.opt,';pensize 1;type dot'
   pd 'new'
   pd opt
   key major
   pd 'title major banks'
   pd 'keyfont Arial 6'
   pd 'keypos top right'
   pd 'xcaption base risk;ycaption strain' 
   pd"1  (4&{. plotdata) 
   pd 'new'
   pd opt
   key minor
   pd 'title minor banks'
   pd 'keyfont Arial 6'
   pd 'keypos top right'
   pd 'keystyle thin'
   pd 'xcaption base risk;ycaption strain' 
   pd"1  (_4&{. plotdata)
   pd 'pdf 350 150 ',figdir,'fig5.pdf'
)
  
NB. fig5''


fig6=: 3 : 0
   agg=.Sum S_oit
   pd 'reset'
   opt=.'xticpos 2008 2009 2010 2011'
   tly=.((36&{.)@(_84&{.))"1 L:0
   y=.2000+yrmo
   pd 'sub 2 2'
   pd 'new'
   pd opt
   key 'non-pooled pooled'
   pd 'keypos top right'
   pd 'yticpos 0 .10 .20 .30 ;ycaption base risk'
   pd  tly y;brisk"2 agg
   pd 'new'
   pd opt
   key 'non-pooled pooled'
   pd 'keypos top right'
   pd 'yticpos 0 .10 .2 .3 ;ycaption base risk '
   pd tly y;strain"2 agg
   pd 'new'
   pd opt
   key 'non-pooled SRISK pooled'
   pd 'keypos top right'
   pd 'yticpos 0 .1 .2 .3 .4;ycaption SRISK'
   pd tly y;(srisk {.agg), (SRISK S_oit) ,: srisk {: agg
   pd 'new'
   pd opt
   pd 'yticpos 0 .2 .4 .6 .8 1 ;ycaption absorbability'
   key 'srisk base-risk'
   pd 'keypos top right'
   pd tly y; ((srisk ABSORB),:(brisk ABSORB)) agg
   pd 'pdf 350 250 ',figdir,'fig6.pdf'
)

NB. fig6''

fig=: 3 : 0
 ". 'fig',(1j0":y),'''',''''
)


table=: 3 : 0
   pi_it=.(%"1 +/)d_it                         NB. 
   prop=.(pi_it&*)@:(%"1 D)                    NB. proportion contributions
   mat=.|:prop"2 (brisk,strain,:srisk) S_oit   NB. Txmx3
   mat=.(|:l_it,:pi_it),."2 mat                NB. Txmx5
   resil=.-/@:ln                                NB. resilience
   m2=.(],resil)"2[2 0 1 |:(brisk,strain,:srisk)"2 Sum S_oit   NB. Tx3x3
   mat=.mat,"2 (0,."2 [ 0,."2 m2)  NB.Tx(m+2)x5
  'mat' Rset wrs y{100*mat
   R 'rnames=c("CBA","ANZ","NAB","WBC","MQG","BOQ","BEN","ABA","non--pooled","pooled","resilience")'
   R 'cnames=c("default","debt","basis","strain","srisk")'
   R 'dimnames(mat)=list(rnames,cnames)'
   R 'print(xtable(mat))'
)

table 72









  