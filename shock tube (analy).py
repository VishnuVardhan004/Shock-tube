from numpy import*
import matplotlib.pyplot as plt


l=1                              #length of shock tube
n=1000                            #total no of cells 
delx=1/1000
Y=1.4                             #gamma
T=1
er=10                          


x=zeros(n)                      
for i in range(0,n):            
    x[i]=(i+0.5)*delx



#Initial conditions
rl=1;rr=0.125                                                   #density
pl=1;pr=0.1                                                     #pressure
ul=0;ur=0                                                       #velocity
el=(pl/(rl*(Y-1)))+(ul*ul/2);er=(pr/(rr*(Y-1)))+(ur*ur/2)
hl=el+(pl/rl);hr=er+(pr/rr)
al=sqrt((Y-1)*(hl-((ul*ul)/2)));ar=sqrt((Y-1)*(hr-((ur*ur)/2)))
UL=zeros(4);UR=zeros(4)
UL[1]=rl;UR[1]=rr                         
UL[2]=rl*ul;UR[2]=rr*ur           
UL[3]=rl*el;UR[3]=rr*er
U=zeros((4,n+2));F=zeros((4,n+2))

#analytical
A=2/(Y+1);B=(Y-1)/(Y+1)
ps=0.5*(pl+pr)
if al>abs(ul) and ar>abs(ur):                                         #subsonic case,waves travelling in both directions
    while er>pow(10,-6):
        SR=RR=SL=RL=0                                          
        if ps>pr:
            FR=(ps-pr)*sqrt((A/rr)/(ps+(B*pr)))                             ;SR=1                         #right moving shock
            FRP=sqrt((A/rr)/(ps+(B*pr)))*(1-((ps-pr)/(2*((B*pr)+ps))))
        else:
            FR=((pow((ps/pr),((Y-1)/(2*Y))))-1)*2*ar/(Y-1)                  ;RR=1                    #right moving rarefaction
            FRP=pow((ps/pr),-(Y-1)/(2*Y))/(rr*ar)
        if ps>pl:
            FL=(ps-pl)*sqrt((A/rl)/(ps+(B*pl)))                              ;SL=1                 #left moving shock
            FLP=sqrt((A/rl)/(ps+(B*pl)))*(1-((ps-pl)/(2*((B*pl)+ps))))   
        else:
            FL=((pow((ps/pl),((Y-1)/(2*Y))))-1)*2*al/(Y-1)                    ;RL=1            #left moving rarefaction
            FLP=pow((ps/pr),-(Y-1)/(2*Y))/(rr*ar)
        F=FL+FR+ur-ul
        psn=ps-(F/(FLP+FRP))
        er=abs(psn-ps)/abs(0.5*(psn+ps))
        ps=psn
        
        
    if SR==1:
        print("right moving shock")
        us=ur+FR
        rrs=rr*((B+(ps/pr))/(((ps*B)/pr)+1))
    if RR==1:
        print("right moving rarefaction")
        us=ur+FR
        rrs=rr*pow((ps/pr),(1/Y))
    if SL==1:
        rls=rl*((B+(ps/pl))/(((ps*B)/pl)+1))
        print("left moving shock")
    if RL==1:
        rls=rl*pow((ps/pl),(1/Y))
        print("left moving rarefaction")   
    
    #shock speed
    ssr=ssl=srrt=srrh=srlt=srlh=0
    ars=sqrt(Y*ps/rrs)                                                  #speed of sound at right of contact discontinuity
    als=sqrt(Y*ps/rls)
    if SR==1:
        qsr=-sqrt((ps-pr)/((1/rr)-(1/rrs)))                           #qsr=rho*u relative  for right moving shock
        ssr=ur-(qsr/rr)                                               #ssr=speed of shock moving 
        RMW=[ssr,ssr]                                                       #RMV=right moving wave speed
    if SL==1:
        qsl=sqrt((ps-pl)/((1/rl)-(1/rls)))
        ssl=ul-(qsl/rl)
        LMW=[ssl,ssl]
    if RR==1:
        srrt=us+ars                                                           #srrt=speed of  rarefaction's tail moving right ;tail is considered the one close to time axis
        srrh=ur+ar
        RMW=[srrt,srrh]
    if RL==1:
        srlt=us-als
        srlh=ul-al
        LMW=[srlt,srlh]

    #cal values at a given location
    t=float(input("enter the time of computation"))
    r=zeros(n);u=zeros(n);p=zeros(n)
    for i in range(0,n):        
        if (x[i]-(l/2))>0:                                                             #(x[i]-(l/2)) can never be zero
            if (x[i]-(l/2))>(RMW[0]*t) and  (x[i]-(l/2))>(RMW[1]*t):
                r[i]=rr
                p[i]=pr
                u[i]=ur 
            if ((x[i]-(l/2)))<(RMW[0]*t) and ((x[i]-(l/2)))<(RMW[1]*t):
                if (x[i]-(l/2))>(us*t):
                    r[i]=rrs
                else:
                    r[i]=rls     
                p[i]=ps
                u[i]=us
            if  (x[i]-(l/2))>(RMW[0]*t) and ((x[i]-(l/2)))<(RMW[1]*t):
                r[i]=rr*pow((A-((B/ar)*(ur-((x[i]-(l/2))/t)))),(2/(Y-1)))
                u[i]=A*(-ar+(0.5*(Y-1)*ur)+((x[i]-(l/2))/t))
                p[i]=pr*pow((A-((B/ar)*(ur-((x[i]-(l/2))/t)))),((2*Y)/(Y-1)))
                j=j+1

        if (x[i]-(l/2))<0:    
            if (x[i]-(l/2))<(LMW[0]*t) and  (x[i]-(l/2))<(LMW[1]*t):
                r[i]=rl
                p[i]=pl
                u[i]=ul 
            if ((x[i]-(l/2)))>(LMW[0]*t) and ((x[i]-(l/2)))>(LMW[1]*t):
                if (x[i]-(l/2))<(us*t):
                    r[i]=rls
                else:
                    r[i]=rrs     
                p[i]=ps
                u[i]=us
            if  ((x[i]-(l/2)))<(LMW[0]*t) and ((x[i]-(l/2)))>(LMW[1]*t):
                r[i]=rl*(A+((B/al)*(ul-((x[i]-(l/2))/t))))**(2/(Y-1))
                u[i]=A*(ar+(0.5*(Y-1)*ul)+((x[i]-(l/2))/t))
                p[i]=pl*((A+((B/al)*(ul-(((x[i]-(l/2))/t)))))**((2*Y)/(Y-1)))
                
        print(i,p[i])

print("ps=",ps,"us=",us,"rls=",rls,"rrs=",rrs)
print("RMW=",RMW,"LMW=",LMW)
print("************************************************************************")
plt.plot(x,r,color='black',linestyle='-',label='density')
#print(r)