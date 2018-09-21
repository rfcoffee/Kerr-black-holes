(**********************  Calculate the transmission coefficient for s=1 mode  ************************)

Needs["DifferentialEquations`NDSolveProblems`"];
Needs["DifferentialEquations`NDSolveUtilities`"];
Needs["FunctionApproximations`"];
ClearAll[s,w,rp,rm,M,a,astar];

(****************** Common prameters *************************)
rp=1;
astar=0.99999;
M=rp*(1-Sqrt[1-astar^2])/astar^2;
rm=(M*astar)^2/rp;
a=Sqrt[rp*rm];
(** precision parameters **)
r0=rp+10^(-4);(*10^(-7);*)
startstep=10^(-6);(*10^(-11);*)
rmax=1000;(*2000;*)
plotrmax=1000;(*2000;*)
numsteps=Infinity;
stepmax=0.1;

(*****************angular eigenvalues Elm***********************)
s=-1;
m=1.0+10^-22;
L=1.0+10^-15;
Nw=49;
apb=Max[Abs[m],Abs[s]];
amb=m*s/apb;
H[L_]:=(L^2-apb^2)*(L^2-s^2)*(L^2-amb^2)/2/(L-1/2)/L^3/(L+1/2);
f0=L*(L+1);
f1=-2*s^2*m/L/(L+1);
f2=H[L+1]-H[L]-1;
f3=2*s^2*m*(H[L]/(L-1)/L^2/(L+1)-H[L+1]/(L+2)/(L+1)^2/L);
f4=4*s^4*m^2*(H[L+1]/(L+2)^2/(L+1)^4/L^2-H[L]/(L-1)^2/L^4/(L+1)^2)+1/2*(H[L+1]^2/(L+1)+H[L+1]*H[L]/(L+1)/L-H[L]^2/L)+1/4*((L-1)*H[L]*H[L-1]/L/(L-1/2)-(L+2)*H[L+1]*H[L+2]/(L+1)/(L+3/2));
f5=8*s^6*m^3*(H[L]/(L-1)^3/L^6/(L+1)-H[L+1]/(L+2)^3/(L+1)^6/L^3)+s^2*m*(3*H[L]^2/(L-1)/L^3/(L+1)-(7*L^2+7*L+4)*H[L]*H[L+1]/(L-1)/L^3/(L+1)^3/(L+2)-3*H[L+1]^2/(L+2)/(L+1)^3/L+1/2*((3*L+7)*H[L+1]*H[L+2]/(L+3)/(L+3/2)/(L+1)^3/L-(3*L-4)*H[L]*H[L-1]/(L-2)/(L-1/2)/L^3/(L+1)));

term1=16*s^8*m^4*(H[L+1]/(L+2)^4/(L+1)^8/L^4-H[L]/(L-1)^4/L^8/(L+1)^4);
term2=4*s^4*m^2*(3*H[L+1]^2/(L+2)^2/(L+1)^5/L^2+(11*L^4+22*L^3+31*L^2+20*L+6)*H[L]*H[L+1]/(L-1)^2/L^5/(L+1)^5/(L+1)^2-3*H[L]^2/(L-1)^2/L^5/(L+1)^2+1/2*((3*L^2-8*L+6)*H[L]*H[L-1]/(L-2)^2/(L-1)/(L-1/2)/L^5/(L+1)^2-(3*L^2+14*L+17)*H[L+1]*H[L+2]/(L+3)^2/(L+2)/(L+3/2)/(L+1)^5/L^2));
term3=1/4*(2*H[L+1]^3/(L+1)^2+(2*L^2+4*L+3)*H[L]^2*H[L+1]/L^2/(L+1)^2-(2*L^2+1)*H[L+1]^2*H[L]/(L+1)^2/L^2-2*H[L]^3/L^2+(L+2)*(3*L^2+2*L-3)*H[L]*H[L+1]*H[L+2]/4/(L+3/2)^2/(L+1)^2/L-(L-1)*(3*L^2+4*L-2)*H[L+1]*H[L]*H[L-1]/4/(L-1/2)^2/L^2/(L+1)+(L+2)*H[L+2]^2*H[L+1]/4/(L+3/2)^2/(L+1)^2-(L-1)^2*H[L-1]^2*H[L]/4/(L-1/2)^2/L^2+(L-1)*(7*L-3)*H[L]^2*H[L-1]/4/(L-1/2)^2/L^2-(L+2)*(7*L+10)*H[L+1]^2*H[L+2]/4/(L+3/2)^2/(L+1)^2+(L+3)*H[L+1]*H[L+2]*H[L+3]/12/(L+3/2)^2/(L+1)-(L-2)*H[L]*H[L-1]*H[L-2]/12/(L-1/2)^2/L);
f6=term1+term2+term3;

(****************** s=-1 **********************************)
Clear[k,λ,S,A,B,α,β,P,Q,Rhor,χhor,R,χ];
C1=(I λ)/(2 w);
C2=-((a m-a^2 w-2 I a m M w)/(2 w));
Tri[r_]:=r^2-2*M*r+a^2;
x[r_]:=r+2*M/(rp-rm)*(rp*Log[r-rp]-rm*Log[r-rm])
k=w-m*a/(1+a^2);
S[r_]:=r*Exp[I*w*x[r]+C1/r+C2/r^2];
A[r_]:=-(s+1)Tri'[r]/Tri[r];
B[r_]:=-((r^2+a^2)^2*w^2 - 4 a M r w m + a^2*m^2 + 2 I a (r-M) m s - 2 I M (r^2-a^2) w s + (2 I r w s-λ) Tri[r])/Tri[r]^2;
α[r_]:=A[r]-2*S'[r]/S[r];
β[r_]:=(A[r]*S'[r] + B[r]*S[r]-S''[r])/S[r];
P[r_]:=α[r]+β'[r]/β[r];
Q[r_]:=α'[r]+β[r]-α[r]*β'[r]/β[r];
Rhor[r_]:=Tri[r]^(-s)*Exp[-I*k*x[r]];
χhor[r_]:=(Rhor'[r]*S[r]-Rhor[r]*S'[r])/(S[r]^2);

(*parameters*)
Elm=N[f0]+N[f1]*a*w+N[f2]*a^2*w^2+N[f3]*(a*w)^3+N[f4]*(a*w)^4+N[f5]*(a*w)^5+N[f6]*(a*w)^6;
λ=-(-1)*(-1+1)+a^2*w^2+Elm;
BB=Sqrt[(Elm+a^2*w^2-2*a*m*w)^2+4*m*a*w-4*a^2*w^2];
dw=0.01;
list1=dw*Range[1,Nw];
list2=list1;

For[i=1,i<=Nw,i++,
  w=dw*(0+i);

  solR2=NDSolve[{Rm''[r]==A[r]*Rm'[r]+B[r]*Rm[r],Rm[r0]==Rhor[r0],Rm'[r0]==Rhor'[r0]},Rm,{r,r0,rmax},StartingStepSize->startstep,Method->{"StiffnessSwitching"},AccuracyGoal->10,MaxSteps->numsteps,MaxStepSize->stepmax];

  ϵ=Sqrt[M^2-a^2]/(4*M);
  list2[[i]]=128*w*k*M^3*(k^2+4*ϵ^2)/(BB^2)/((1/rmax*Abs[Evaluate[Rm[rmax]/.solR2]])^2+128*w*k*M^3*(k^2+4*ϵ^2)/(BB^2))];

Export["out_s1.dat",Transpose[{list1,Flatten[list2]}],"Table"];
Exit[];
