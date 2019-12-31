points = 10;
GrmFukui = zeros(points, 1);
GrmRight = zeros(points, 1);
GrmYoshiSquared = zeros(points, 1);
for i=1:points
rng(i)
Fmin = rand();
Ref = rand();
Gs = rand();
Gof = rand();
Bs = rand();
Bof = rand();
Gamax = rand();
Reg = rand();
Gog = rand();
Bog = rand();

c = ((Gof^2)*Ref+(Gog^2)*M*Reg+(Bof^2)*Ref+(Bog^2)*M*Reg)/(Ref + M*Reg);
M = (Fmin - 1 + Ref/Gs * ((Gs-Gof)^2+(Bs-Bof)^2))/(1-1/Gamax - Reg/Gs*( (Gs-Gog)^2+ (Bs-Bog)^2) );
ShouldBeZero = Gs/(M*Reg+Ref) * (M*(1/Gamax - 1) + Fmin - 1) + Gs^2 + Bs^2 - 2*Gs*(Gof*Ref+Gog*M*Reg)/(M*Reg+Ref)...
    -2*Bs*(Bof*Ref+Bog * M *Reg)/(Ref+M*Reg) + c;


Gm = (  M*(2*Reg*Gog - 1/Gamax + 1) + 2*Ref*Gof - Fmin +1  )/(2*(M*Reg+Ref) );
Bm = (M*Reg*Bog+Bof*Ref)/((M*Reg+Ref));

ShouldAlsoBeZero= Gs^2+Bs^2-2*Gs*Gm-2*Bs*Bm+c;


GrmFukui(i) = 1/(2*(M*Reg+Ref) ) * ( (M*(1/Gamax - 1) + Fmin - 1)^2 ...
    -4* (M*Reg*Gog + Ref*Gof)*(M*(1/Gamax - 1) + Fmin - 1)...
    - 4*M*Reg*Ref*(Gog^2+Bog^2+Gof^2+Bof^2-4*Bof*Bog) );
GrmYoshiSquared(i) = 1/(4*(M*Reg+Ref)^2 ) * ( (M*(1/Gamax - 1) + Fmin - 1)^2 ...
    -4* (M*Reg*Gog + Ref*Gof)*(M*(1/Gamax - 1) + Fmin - 1)...
    - 4*M*Reg*Ref*(Gog^2+Bog^2+Gof^2+Bof^2-2*(Bof*Bog+Gof*Gog)) );
GrmRight(i) = (Gs-Gm)^2 + (Bs-Bm)^2;

end

indices= 1:points;
semilogy(indices, abs(GrmFukui-GrmRight), indices, abs(GrmYoshiSquared-GrmRight));
legend("Fukui", "Yoshi");
title("Error vs (Gs-Gm)^2 + (Bs-Bm)^2");

