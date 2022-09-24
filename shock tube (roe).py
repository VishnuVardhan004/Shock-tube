from numpy import*
from matplotlib.pyplot import*

l=1                               #length of shock tube
n=1000                            #total no of cells 
delx=l/1000
Y=1.4                             #gamma
T=1
eps = 10**(-6)                          


x=zeros(n+2)                      # 0, n+1 are ghost cells
for i in range(1,n+1):            # includes 1,n
    x[i]=(i-0.5)*delx
x[0]=-0.5*delx
x[n+1]=1+(0.5*delx)

#Initial conditions
rl=1;rr=0.125                                                   #density
pl=1;pr=0.1                                                     #pressure
ul=0;ur=0                                                       #velocity
el=(pl/(rl*(Y-1)))+(ul*ul/2);er=(pr/(rr*(Y-1)))+(ur*ur/2)
hl=el+(pl/rl);hr=er+(pr/rr)
UL=zeros(4);UR=zeros(4)
UL[1]=rl;UR[1]=rr                         
UL[2]=rl*ul;UR[2]=rr*ur           
UL[3]=rl*el;UR[3]=rr*er
U=zeros((4,n+2));F=zeros((4,n+2))
for i in range(1,int((n+2)/2)):
    U[1][i]=rl
    U[2][i]=rl*ul
    U[3][i]=rl*el
    F[1][i]=rl*ul
    F[2][i]=(rl*ul*ul)+pl
    F[3][i]=((rl*el)+pl)*ul

U[1][0]=UL[1]; U[2][0]=-UL[2];U[3][0]=UL[3]                                        #ghost cells
F[1][0]=rl*ul;F[2][0]=(rl*ul*ul)+pl;F[3][0]=ul*((rl*el)+pl)   

for i in range(int((n+2)/2),n+1):        
    U[1][i]=rr
    U[2][i]=rr*ur
    U[3][i]=rr*er
    F[1][i]=rr*ur
    F[2][i]=(rr*ur*ur)+pr
    F[3][i]=((rr*er)+pr)*ur

U[1][n+1]=UR[1];U[2][n+1]=-UR[2];U[3][n+1]=UR[3];                                   #ghost cell
F[1][0]=rr*ur;F[2][0]=(rr*ur*ur)+pr;F[3][0]=ur*((rr*er)+pr)




#iterations
ht=zeros(n+1); ut=zeros(n+1); at=zeros(n+1)
lambt=zeros((4,n+1));K=zeros((4,4,n+1))
del2=zeros(n+1);del1=zeros(n+1);del3=zeros(n+1);DEL=zeros((4,n+1))
Ff=zeros((4,n+3))


t=0
while t<T:
    for i in range(0,n+1):
        

        ht[i]=((sqrt(U[1][i])*((Y*U[3][i]/U[1][i])-((Y-1)*(U[2][i]/U[1][i])*(U[2][i]/U[1][i])/2)))+(sqrt(U[1][i+1])*((Y*U[3][i+1]/U[1][i+1])-((Y-1)*(U[2][i+1]/U[1][i+1])*(U[2][i+1]/U[1][i+1])/2))))/(sqrt(U[1][i])+sqrt(U[1][i+1]))
        ut[i]=((sqrt(U[1][i])*(U[2][i]/U[1][i]))+(sqrt(U[1][i+1])*(U[2][i+1])/U[1][i+1]))/(sqrt(U[1][i])+sqrt(U[1][i+1]))
        at[i]=sqrt((Y-1)*(ht[i]-((ut[i]*ut[i])/2)))

        
        lambt[1][i]=ut[i]-at[i];lambt[2][i]=ut[i];lambt[3][i]=ut[i]+at[i]
        for j in range(1,4):                                                           #entropy fix
            if (abs(lambt[j][i])<eps):
                if(lambt[j][i] !=0):
                    lambt[j][i] = 0.5*((lambt[j][i]/eps) + eps)



        
        K[1][1][i]=1;                         K[1][2][i]=1;                             K[1][3][i]=1
        K[2][1][i]=ut[i]-at[i];               K[2][2][i]=ut[i];                         K[2][3][i]=ut[i]+at[i]
        K[3][1][i]=ht[i]-(ut[i]*at[i]);       K[3][2][i]=ut[i]*ut[i]/2;                 K[3][3][i]=ht[i]+(ut[i]*at[i])  
        

        del2[i]=((Y-1)/(at[i]*at[i]))*(((U[1][i+1]-U[1][i])*(ht[i]-(ut[i]*ut[i])))+(ut[i]*(U[2][i+1]-U[2][i]))-(U[3][i+1]-U[3][i]))
        del1[i]=(1/(2*at[i]))*(((U[1][i+1]-U[1][i])*(at[i]+ut[i]))-(U[2][i+1]-U[2][i])-(at[i]*del2[i]))
        del3[i]=U[1][i+1]-U[1][i]-del1[i]-del2[i]
        
        DEL[1][i]=del1[i];DEL[2][i]=del2[i];DEL[3][i]=del3[i]

        
        Ff[1][i]=0
        Ff[2][i]=0
        Ff[3][i]=0
        for j in range(1,4):
            if lambt[j][i]<0:
                Ff[1][i]+=(lambt[j][i]*DEL[j][i]*K[1][j][i])
                Ff[2][i]+=(lambt[j][i]*DEL[j][i]*K[2][j][i])
                Ff[3][i]+=(lambt[j][i]*DEL[j][i]*K[3][j][i])
            
        Ff[1][i]=F[1][i]+Ff[1][i]
        Ff[2][i]=F[2][i]+Ff[2][i]
        Ff[3][i]=F[3][i]+Ff[3][i]

    
       


    delt=(0.8*delx)/amax(abs(lambt)) 



    for i in range(1,n+1):    
        U[1][i]=U[1][i]-((delt/delx)*(Ff[1][i]-Ff[1][i-1]))
        U[2][i]=U[2][i]-((delt/delx)*(Ff[2][i]-Ff[2][i-1]))
        U[3][i]=U[3][i]-((delt/delx)*(Ff[3][i]-Ff[3][i-1]))
    
    #boundary conditions updation
    U[1][n+1]=U[1][n];             U[1][0]=U[1][1]
    U[2][n+1]=-U[2][n];            U[2][0]=-U[2][1]
    U[3][n+1]=U[3][n];             U[3][0]=U[3][1]
    
    
    for i in range(0,n+1):
        F[1][i]=U[2][i]
        F[2][i]=(U[3][i]*(Y-1))+(((3-Y)/2)*U[2][i]*U[2][i]/U[1][i])
        F[3][i]=(U[2][i]/U[1][i])*((Y*U[3][i])-((Y-1)*0.5*U[2][i]*U[2][i]/U[1][i]))
    
    t=t+delt

    




plot(U[1])
for i in range(1,n+1):
    print(i,U[1][i])
print(delt)             
