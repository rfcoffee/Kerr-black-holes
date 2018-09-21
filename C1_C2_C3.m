(*%%%%%%%%%%%%  Derive the C1, C2, C3 parameters for s = ∓2, ±1, ±0.5 for chosen S - function %%%%%%%%%%%%%%%%%%%%%%*)

Tri[r_]:=r^2-2*M*r+a^2
x[r_]:=Integrate[(r^2+a^2)/(Tri[r]),r]
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=-2           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3]
s=-2;
RR[r_]:=r^3*Exp[I*w*x[r]+C1/r+C2/r^2+C3/r^3]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,2}]
Solve[{-2 I C1 w-λ==0,-2 C1-4 I a m-6 M+6 I a^2 w-4 I C2 w+4 I C1 M w-4 a m M w+2 M λ==0,C1^2-2 C2+10 C1 M+4 I a m M+12 M^2-6 I C3 w+8 I C2 M w+a^2 (6+m^2-4 I C1 w-12 I M w-λ)==0},{C1,C2,C3}]
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=2           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3];
s=2;
RR[r_]:=r^(-1)*Exp[-I*w*x[r]+C1/r+C2/r^2+C3/r^3]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,2}]
Solve[{-4+2 I C1 w-8 I M w-λ==0,-2 C1+4 I a m+10 M+2 I a^2 w+4 I C2 w-4 I C1 M w-4 a m M w+2 M λ==0,C1^2-2 C2+2 C1 M-4 I a m M-4 M^2+6 I C3 w-8 I C2 M w+a^2 (-2+m^2+4 I C1 w+4 I M w-λ)==0},{C1,C2,C3}]

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=1           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3]
s=1;
RR[r_]:=r^(-1)*Exp[-I*w*x[r]+C1/r+C2/r^2]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,2}]
Solve[{-2+2 I C1 w-4 I M w-λ==0,2 I a m+4 M+2 I a^2 w+4 I C2 w-4 I C1 M w-4 a m M w+2 M λ==0},{C1,C2}]
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=-1           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3]
s=-1;
RR[r_]:=r*Exp[I*w*x[r]+C1/r+C2/r^2]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,2}]
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=0.5          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3]
s=0.5;
RR[r_]:=r^(-1)*Exp[-I*w*x[r]+C1/r]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,2}]
Simplify[Solve[{-1+(2 I) C1 w-(2 I) M w-1 λ==0},{C1}]]
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=-0.5           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3]
s=-0.5;
RR[r_]:=1/r*Exp[-I*w*x[r]+C1/r]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,2}]
Solve[{(-2 I) C1 w-λ==0},{C1}]
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        s=0          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
Clear[s,C1,C2,C3]
s=0;
RR[r_]:=r^(-1)*Exp[-I*w*x[r]+C1/r+C2/r^2+C3/r^3]
termone[r_]:=Tri[r]^2*RR''[r]
termtwo[r_]:=(s+1)*Tri[r]*Tri'[r]*RR'[r]
termthree[r_]:=((r^2+a^2)^2*w^2 -4*a*M*r*w*m +a^2*m^2 + 2*I*a*(r-M)*m*s -2*I*M*(r^2-a^2)*w*s+(2*I*r*w*s-λ)*Tri[r])*RR[r]
termall[r_]:=(termone[r]+termtwo[r]+termthree[r])/RR[r]
Series[Simplify[termall[r]],{r,∞,4}]
