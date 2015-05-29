NB. coclass 'kfs'

o=:>@{                           NB. indexed open box
bigk=:1e11                       NB. big k -- cannot be less than 1e11
eps=:1000 * % bigk               NB. epsilon
                                               

ptnz=:(>&eps)@|                  NB. point to numbers with magnitude bigger than eps
gnz=:ptnz # ]                    NB. get non zero numbers with zero as defined in ptnz
iom=:$#:(i. <./)@,               NB. array index of minimum

desmat=:(^/ i.@>:)L:0            NB. boxed design matrices from boxed lists in x  -- y is power
kpb=:>`(kp/@:>)@.(*@<:@#)        NB. kronecker product of boxed list of matrices 
reo=:,@((</.)"2^:(<:@$@$))@i.@cd NB. reorder list to put lowest powers first in kp
   cd=:;@(({:@$)L:0)             NB. last dimensions in boxed arrays
dra=:;^:L.                       NB. deep raze -- raze out of lowest level -- monadic !
poldn=:] (dra@(>:@[{.reo@:]){"1 kpb@:]) desmat  NB. polynomial design matrix in "usual" order
qpol=:,@]%.(poldn&2)@[           NB. regn vec from quad - (x1;x2;x3;..) qpol y - ($>i{x) = (i{$y)
                                 NB. !!! x1 slowest, then x2 etc.
nbv =: -@-:@({.}.)               NB. -0.5b in fit f(t)=c+b't+t'At with fit using qpol
Amat=:sym@(>:@[ }. ])            NB. A matrix in quad fit -- t=(z,y,x) rather than (x,y,z)
nth=:#@[ (nbv %. Amat) qpol      NB. optimal t from quadratic fit ie t=-0.5A^{-1}b
                                 NB. use (x1;x2;...) nth y where y as in qpol
  

pow=: 4 : 'x&mp^:(<:y) x'     NB. multiply x by itself integer y times
ric=:$@[$(,@]%.((kp~+kp)id@#@[)) NB. solve xv+vx'=y  for v
ri2=:$@[$(,@]%.(id@*:@# - kp~)@[)NB. solve xvx'+y=v  for v
fp=:((#@]{."1[)mp]),"1#@]}."1[   NB. (A,"1 B) fp X is ((A mp X),"1 B)
ioi=:_&= # i.@$                  NB. indexes of infinity
zinf=:0&(ioi@]})                 NB. zero the infinite entries


NB.  Regression from lower triangular square root of sum of sq matrix

NB.  Use :  regnlt y where (tp y) = (ssq matrix) with y square upper triangular
NB.  Out :  (rss;b;S) where  (sqr S) = cov(b) (excluding sig sq)
NB.          b=reg coefs, (tp rss) = residual ssq 

regnlt=: 3 : 0
    C=. |: linv }: }:"1 y
    beta=.-(}:{:y) sqr C
    q=.{:{:y
    q;beta;C
)


pinv =: 3 : 0
c =.   eps&< | getd y
linv&.(c&#@:(c&#"1)) y
)

onclean =: 1 : 0
'r c' =. (+./"1 ; +./) * y
u&.(r&#@:(c&#"1)) y
)



                            NB.  Kalman Filter iteration

KF =: 3 : 0
   'gamma alpha0 D0'=.(>@{.,{:)@{. y               
   BsA0  =.|:@((ioi{id@#),zinf) (gamma,alpha0)  
   ngamma=. # gamma   
   na0   =. # alpha0                                 
   Bmat  =. ngamma {. BsA0                   
   r=.,:(((ngamma+na0){."1 D0) mp BsA0);(ngamma+na0)}."1 D0
   t=.1 
   while. t< #y do.
     'At Pt'=.2{.{:r                                
     'yt Dt'=.t{y                                  
      nyt=. # yt
      nat=. # At                
      VtsAtp =. ((ngamma+nat) {."1 Dt) mp  Bmat,At  NB.  (V_t,A_t+) before dealing with y_t 
      Vt  =.yt (-@}:@],(-{:))"0 1 (nyt{.VtsAtp)     NB.  Now "subtract" Vt from yt
      QRmat =. QRt (ngamma }."1 Dt) fp Pt
      Ft=.  nyt ([ {."1 {.) QRmat
      if. (0 = */ * getd Ft) do. 'KF warning: cov(y) singular' 1!:2 (2) end. 
      Fit =. linv Ft
      Ptp1 =. nyt ([ }."1 }.) QRmat
      Kt =. (nyt ([ {."1 }.) QRmat)  mp Fit
      Atp1=.(nyt }. VtsAtp) + Kt mp Vt
      t=.>:t
      r=.(}:r),(At;Pt;Vt;Fit;Kt),:Atp1;Ptp1
   end.
)


KFR =: 3 : 0                                          NB.  Uses the output of KF
   V =.> 2{"1 y
   Fi=.> 3{"1 y
   a=.  gnz ;  getd"_1   Fi                           NB.  Get nonzero
   nh=.# a                                            NB.  Number nonzero
   d=.-  +: +/ ^. | a                                 NB.  _2 * Sum ln|inv(F_t)|
   'q beta C'=. regnlt QRt |: (,/) Fi mp"_1  V
   sig=.|q%%:nh                                       NB.  sigma hat
   ell0=.d+nh*>:+:^.sig
   if. (p=.#beta)=0 do. elli=.ell0  
   else.  detp=. +: ^. | -/ . * C
          elli =.(d-detp)+(nh-p)*>:^.(nh**:sig)%(nh-p) end.
   res_KFR=:(nh,d,sig);beta;C;(ell0,elli)
) 


gr=:{@:((<@:({.@>@] +({:-{.)@:>@] * (i.@>: % ])@[))"0)  NB. x grid from boxed intervals in y

scovb=.std@(*:@{:@>@{. * sqr@>@(2&{))

ell0=:{.@>@{:@KFR
elli=:{:@>@{:@KFR
                           

KFS =:  KF ([ ,"1 SF) ]                  NB. Kalman Filter Smoother -- note allignment of output

SF=: 4 : 0                               NB.   x is KF output
    r=.,:@(((0&*)&.>)@(2&{.)@{:) x      NB.   y is cell array 
    ngamma=.#>{.>{.{.y  
    t=.<:# y                          
    whilst. t>0 do.
       'Rt Nt'=.{.r
       'At Pt Vt Fit Kt'=.(t-1){x
       'yt Dt'=.t{y
       nat =.# At
       nyt =.# yt
       natp1=.(#Dt)-nyt
       Zt    =. ((0,ngamma),:(nyt,nat)) ];.0 Dt
       Tt    =. ((nyt,ngamma),:(natp1,nat)) ];.0 Dt
       Lt    =. Tt - Kt mp Zt
       Rtm1  =.  (Lt tp Rt) + (Zt tp (tp Fit) mp Vt)
       Ntm1  =.  nat {."1 QRt (Zt tp (|: Fit)) ,"1  Lt tp Nt         NB.  Check this !!
       t=.<:t
       r=.(Rtm1;Ntm1),r
    end.
)

dA =: 4 : 0     NB.  Diffuse adjustment
  'A P' =. x
  's2 b covb'=. y
  if. (#b)> 0 do.
     P=. P + (}:"1 A) qp covb
     A=. A mp b,1
  end. 
  A; s2*P
)

dAsqrt =: 4 : 0  NB. Diffuse adjustment in square root form
  'A P' =. x
  'sig b C'=. y
  if. ((#b)>0) do.                NB. *.((#A)>0) do.
     P=. P qr (}:" 1 A) mp C
     A=. A mp b,1
  end. 
  A; sig*P
)

NB.  Prediction factoring in the diffuse effects

PR=: 3 : 0                               NB.  Uses output of KF
    'sc beta C ell'=.KFR r=. y
    'nh d sig'=. sc
    sc=.sig;beta;C
    'At Pt'=.0 1{ {: r
    out=.,:(At;Pt) dAsqrt sc
    t=.<:#y
    whilst. t>0 do.
      t=.<:t
      'At Pt Vt Fit Kt'=.t{r
      aP=.(At;Pt) dAsqrt sc
      Ft=. linv Fit
      vF=.(Vt;Ft) dAsqrt sc
      out=.(aP,vF),out
   end.
) 

NB.  Smoothing and interpolation factoring in the diffuse effects

SMI=: 3 :  0
    'sc b C ell'=.KFR r=.KFS y 
    'nh d sig'=.sc
    sc=.(*:sig);b;sqr C
    'At Pt Rt Nt'=. 0 1 5 6{ {: r
    Pt=.sqr Pt
    out=.,:(At;Pt) dA sc
    t=.<:#y 
    whilst. t>0 do.
      t=.<:t
      'At Pt Vt Fit Kt Rtm1 Ntm1'=.t{r
      Pt=. sqr Pt
      Fit=. tp Fit
      Ntm1 =. sqr Ntm1
      At=.At + Pt mp Rtm1  
      Pt=.Pt - Pt qp Ntm1
      aP=.(At;Pt) dA sc
      Mit=.%. Fit + Kt qtp Nt
      Ut=.Mit mp (Fit mp Vt) -(Kt tp Rt) 
      uM=.(Ut;Mit) dA sc
      Rt=.Rtm1
      Nt=.Ntm1
      out=.(aP,uM),out
    end.
    out=.}:out
)


GCV=: 3 : 0                             NB.  CV and GCV
   ell=.>{:KFR y
   interp=. _2&{"1 r=. SMI y
   cov=.{:"1 r
   CV=.+/ ; (sqr L:0) interp
   num=. +/ ,  interp ((sqr@%.)&>)"0 cov
   den=. *: +/ ; ((trace@%.) L:0) cov
   CV , (num % den) , ell
)
 
NB.  Error smoothing, factoring in the diffuse effects. 

ESM=:4 : 0                     
    'sc b C ell'=.KFR y          NB.  x is ssm -- y is KF output
    ngamma=.#>{.>{.{.x 
    'nh d sig'=.sc
    sc=.(*:sig);b;sqr C
    s=.y SF x             NB. (R_0;N_0),...,:(R_n;N_n)
    out=.1 1 $<''
    t=.#y 
    while. t>1 do.
      t=.<:t
      'Vt Fit Kt'=.2 3 4{(t-1){y
      'yt Dt'=.t{x
      'Rt Nt'=.t{s
      Fit=.tp Fit 
      Nt=. sqr Nt
      nyt=.#yt
      nat=.#>{.(t-1){s
      GtsHt =.(ngamma+nat)}."1 Dt    
      Gt=.nyt{.GtsHt
      Ht=.nyt}.GtsHt
      Jt=.Ht - Kt mp Gt
      Et=.(Gt tp Fit mp Vt) + Jt tp Rt
      Ct=.(id@# - ]) (Gt qtp Fit) +  Jt qtp Nt
      'et Ct'=.(Et;Ct) dA sc
      out=.(et;Ct;(GtsHt mp et);(GtsHt qp Ct);yt-Gt mp et),out
   end.
   out=.}:out
)

DIAG=: 3 : 0   NB.  Computes chi^2  and component z-stats  - uses output of KF
     'sc b C ell'=.KFR r=.y
     'nh d sig'=.sc
     sc=.sig;b;C
     s=.}.r SF y
     out=.1 1$<''
     t=.0
     while. t<<:#y do.
       'Vt Fit'=.2 3{t{r
       'Rt Nt'=.t{s
       Ft=.linv Fit
       'vt Ft'=.(Vt;Ft) dAsqrt sc
       'rt Nt'=.(Rt;Nt) dAsqrt sc
       Fit=.tp linv Ft
       Nit=.tp linv Nt
       chit=.vt , Nit mp rt
       sdchit=. %:(getd sqr Ft), getd Nit
       chisqt=. (vt qtp Fit) + rt qtp Nit 
       out=.out,chisqt;chit%sdchit
       t=.>:t
     end.
     out=.}.out
)

DIAGE=: 3 : 0            NB. standardized E(\eps|y)
    r=.ESM y
    sig=.{: > {. KFR y
    e=.> {."1 r
    c=.> getd L:0 (1&{"1) r
    c=. %: (*: sig) - c
    e=. e%"0 c 
)





      

                           NB.  ARMA(x,y) - the max(p,q) form

arma=:|:@({.,}:@}.@id@{:@$,{:)@pa      NB.  Pearlman's ARMA(x,y) representation Dt matrix
  pa=:1&,"1@({.,:+/)@,:&vec            NB.  Preliminary set up for arma

ARMA=: 4 : 0                           NB.  Use phi ARMA theta to do KF on ARMA model
                                       NB.  y data is in bj
  n=.12                                NB.  Sample size -- vary to suit
  DT=.x arma y                       NB.  Set up D_t matrix as per KFS writeup 
  y=.1.1 4.3 3.0 _0.7 1.6 3.2 0.3 _1.9 _0.3 _0.3 0.8 2.0   NB.  Box Jenkings sec 7.1.5 data 
  yarma=.(n $ y) (;"0 2)  DT           NB.  Right argument for KF
  T0=.}. }:"1 DT                       NB.  The T_0 matrix is taken as T
  H=.,. }. {:"1 DT                     NB.  This is H in Pearlman specification
  H0=.QRt  T0 ri2 (sqr H)              NB.  Square root of P_1=TP_1T'+HH'
  D0=.((#T0)#0);(T0,"1 H0 )            NB.  Left argument for KF -- note missing gamma
  D0 KF yarma                          NB.  Apply KF
)

uneq=: 4 : 0                           NB.  Work out T(k) and H(k) for unequal spacing
    p=.#y                             NB.  y is the (T,H) array
    T=.p {."1 y
    H=.p }."1 y                       NB.  H is definitely a matrix
    q=.{:$H                            NB.  column dimension of H
    k=.1
    r=.,:((id p) ,"1 (0*H))  
    while. k<:x do.
      Tk=.T mp p{."1 {:r
      Hk=.q{."1 QRt (T mp p}."1 {: r),"1 H
      Hk=.  ]`(_1&*)  @. ((_1&=)@:<./@:*@:getd) Hk   NB.  Change sign if negative
      k=.k+1
      r=.r,(Tk,"1 Hk) 
    end.   
) 


('z';'piet';'kfs') copath 'base'

   
