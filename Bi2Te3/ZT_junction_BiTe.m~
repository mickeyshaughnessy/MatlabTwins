%This script computes ZT for the single twin junction
%It uses the expression
%R_net(E) = R_dif(E) + N*R_twin(E)
%and R = 1/G = 1/(2e^2/h * T(E))
%For T(E), I use either the output of PWCOND or the expression lambda(E)/L*M(E)
%The input to the ZT integrals in T_net(E) = 2e^2/h*1/R_net(E) =
%1/(l/lambda(E)+(1/T_twin(E)-1/T_pure(E)))
%The temperature is modified in several places: the T in ZT, the T in the
%FD functions, and the value of the bulk lattice conductivity
clear all
%load ZT.mat
load pure
load twin
load DOSpure
load DOStwin
DOSpure(:,1) = DOSpure(:,1)- 6.2833+0.01; %Ef pure, +0.1 to center Ef in midgap at zero doping
DOStwin(:,1) = DOStwin(:,1)- 6.272+0.01; %Ef twin, +0.1 to center Ef in midgap at zero doping
pure(:,1) = pure(:,1)-0.025; %these center the Fermi energy at midgap
twin(:,1) = twin(:,1)-0.025;

Emin = -3.025;
Emax = 2.125;
deltaE = 0.01;
Esteps = (Emax - Emin)/deltaE;
E = linspace(Emin, Emax, Esteps);
Ttwin = interp1(twin(:,1), twin(:,2), E, 'linear');
Tpure = interp1(pure(:,1), pure(:,2), E, 'linear');
DOSp =  interp1(DOSpure(:,1), DOSpure(:,2), E, 'linear');
DOSp = 2*DOSp %raw data comes from 30 atom VASP calculation, multiply by 2 to normalize to 60 atom unit cell
%DOSpure = interp1(E-EF25,DOS25(:), E, 'spline');
%lambda_0 = 30E-9; %[m]
%lambda_1 = 30E-9;
ncarrier_0 = 2.8E16; %[cm^-3]
%mu_min = 68.5;
%mu_max = 1414;
%Nr = 9.2E16;
%alpha = 0.711;
%lambda = 18E-6 %[m] from Datta and Jeong paper on phonons
lambda = 18E-9 %[m] from Datta and Jeong paper on electrons

kp_exp = 1 %[W/mK], thermal conductivity from the literature
sigma = 3.838e9 %[W/m^2K] %kapitza conductance
%fact1 = lambda_0/2.8e19
%fact2 = lambda_0/333; %this emprical factor converts between mobility and mean free path

N = 1; %number of twins, so L = twin-twin distance
%L = 100E-9 %m length of unit cell
%L = 7.52644E-9 %m length of unit cell
%L = 4*7.52644E-9 %m length of unit cell
%volcellj = 4*volcell %cm^3 since we are using a 48-atom cell which is 4 times as long as the 12-atom one
%Acellm = Acell*0.0001 %Cell area in meters
A = 16.8344e-20 %m^2 %60.1165 (a.u.)^2
L = 4*30.980205783e-10 %m,  4 for 60 atom cell
volcell = A*L %m^3

Tr = 1./((L/lambda - N)./(Tpure+0.00001) + N./(Ttwin+0.00001)); %This is the net transmission with a twin (inverse of net resistance)
Tp = (lambda/L)*Tpure;  %p stands for pure (no twin)

kp_units = kp_exp/(1+kp_exp/(sigma*L));%[W/mK]
kp = kp_units*A/L; %(W/K)
kp_exp_nodim = kp_exp*A/L %(W/K), for pure ZT

q = 1.602176565e-19; %Coulomb
h = 6.626068e-34; %m^2 kg / s
kb = 1.3806503e-23; %m^2 kg s-2 K-1
Doprange = 100;
Tp(200)
Tpure(200)
lambda
G = zeros(Doprange,1); Gp = zeros(Doprange,1);
S = G; 
ke = S; 
ZT=kp; 
Sp = Gp; 
kep = Sp; 
ZTp=kep; 
%lambda1 = zeros(Doprange,1);
%lambda2 = zeros(Doprange,1);
%mu = zeros(Doprange,1);
ncarrier = zeros(Doprange,1);
pcarrier = zeros(Doprange,1);
c = zeros(Esteps,Doprange);

dFDdE = zeros(Esteps,Doprange);
%dNdE = zeros(Esteps,Doprange);
FD = zeros(Esteps,Doprange);

DOP = linspace(-1,1,Doprange);  %eV

%lambda = 0.5e-6%1.0e-6; %cm, electron mean free path
%lambdap = 1e-5 %cm, phonon mean free path

for Dop = 1:Doprange 
temp = 300;
kbT = kb*temp*6.24150974e18; %in eVs
    
%This section computes the FD and BE distributions as functions of T and
%doping
c(:,Dop) = (E'-DOP(Dop))./kbT; %the c's are the arguments to the distribution functions, the 0.25 is a shift to center the distribution about Ef
dFDdE(:,Dop) = -(1/(kbT))*(1./(2+exp(c(:,Dop))+exp(-c(:,Dop))));
%dNdE(:,Dop) = -(1/(kbT))*(exp(Ep./(kbT)))./(exp(Ep./(kbT)).^2); %BE distribution
FD(:,Dop) = 1./(1+exp(c(:,Dop)));

for i=1:Esteps  %compute the transport integrals
    if (E(i)>0)
        ncarrier(Dop) = ncarrier(Dop)+DOSp(i)*FD(i,Dop)*deltaE; %ncarrier density
        %pcarrier(Dop) = pcarrier(Dop)+DOSpure(i)*(-FD(i,Dop))*deltaE; %pcarrier density
    end
    if(E(i)<0)
        pcarrier(Dop) = pcarrier(Dop)+DOSp(i)*(1-FD(i,Dop))*deltaE; %pcarrier density
        %ncarrier(Dop) = ncarrier(Dop)+DOSpure(i)*(FD(i,Dop)-1)*deltaE; %ncarrier density
    end
end

end
ncarrier = ncarrier./volcell*1E-6; %convert volcell from m^-3 to cm^-3
pcarrier = pcarrier./volcell*1E-6;

for Dop = 1:Doprange
%mu(Dop) = mu_min+(mu_max-mu_min)/(1+(ncarrier(Dop)/Nr)^alpha);

%lambda1(Dop) = fact1*ncarrier(Dop)*mu(Dop); %m
%lambda2(Dop) = fact2*mu(Dop); %m
%Tr = 1./((L/lambda2(Dop) - N)./Tpure + N./Ttwin); %This is the net transmission with a twin (inverse of net resistance)
%Tp = (lambda2(Dop)/L)*Tpure;  %p stands for pure (no twin)

%Si thermal conductivity = 1.4W/(cm-K)

sum1 = 0;
sumc=0; %These sums are for computing the transport integrals
sumS = 0;
sumK = 0;
sum1p = 0;
sumcp=0; %These sums are for computing the transport integrals
sumSp = 0;
sumKp = 0;
%phonSum = 0;E
for i=1:Esteps  %compute the transport integrals
    sum1 = sum1 + (-dFDdE(i,Dop)*Tr(i))*deltaE;%G
    sumS = sumS +(-dFDdE(i,Dop)*Tr(i)*(E(i)-DOP(Dop))/kbT)*deltaE; %S
    sumK = sumK +(-dFDdE(i,Dop)*Tr(i)*((E(i)-DOP(Dop))/kbT)^2)*deltaE;%K_e
    sum1p = sum1p + (-dFDdE(i,Dop)*Tp(i))*deltaE;%G
    sumSp = sumSp +(-dFDdE(i,Dop)*Tp(i)*(E(i)-DOP(Dop))/kbT)*deltaE; %S
    sumKp = sumKp +(-dFDdE(i,Dop)*Tp(i)*((E(i)-DOP(Dop))/kbT)^2)*deltaE;%K_e
    %phonSum = phonSum+(-dNdE25(i,Dop)*(Ep(i)/kbT)^2*DOMppure(i))*deltaEp; phonSum_111 = phonSum_111+(-dNdE111(i,Dop)*(Ep(i)/kbT)^2*DOMp111(i))*deltaEp; %K_p
end
%lambda = 5E-9 %m
%Tr = 1./((L/lambda - N)./Tpure + N./Ttwin); %This is the net transmission with a twin (inverse of net resistance)
%Tp = (lambda/L)*Tpure;  %p stands for pure (no twin)

G(Dop) = sum1*(2*q^2/h); %(1/Ohms)
S(Dop) = (sumS/sum1) *kb/q; %(V/K)
ke(Dop) = 2*temp*(sumK-(sumS*sumS/sum1))*kb*kb/h; % (W/K) 
ZT(Dop)= S(Dop).*S(Dop).*G(Dop)*temp./(ke(Dop)+kp); %dimensionless
Gp(Dop) = sum1p*(2*q^2/h); 
Sp(Dop) = (sumSp/sum1p) *kb/q;
kep(Dop) = 2*temp*(sumKp-(sumSp*sumSp/sum1p))*kb*kb/h;
ZTp(Dop)= Sp(Dop).*Sp(Dop).*Gp(Dop)*temp./(kep(Dop)+kp_exp_nodim);


end

carrier = ncarrier-pcarrier %change x-axis to effective carrier density
%ZT_exp = 24970*400*400*1e-12*300/148
%close all;
LorNum = ke(Doprange/2)/(G(Doprange/2)*temp);
figure(1);
semilogx(carrier,G*L/A,':o','linewidth',8), xlim([1e18,1e21]), ylabel('Conductivity ($\Omega^{-1}$m$^{-1}$)','Interpreter', 'latex', 'fontsize',32,'linewidth',2); %xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32);
hold all;
semilogx(carrier,Gp*L/A,'linewidth',8), xlim([1e18,1e21]), ylabel('Conductivity ($\Omega^{-1}$m$^{-1}$)','Interpreter', 'latex', 'fontsize',32,'linewidth',2); %xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32);
%plot(1.7E19,24970, '.', 'markersize', 70) %This is from "Transport Properties of Si" by Weber and Gmelin (1991)
set(gca, 'fontsize', 30);
h_legend= legend('Twin Bi$_2$Te$_3$','Pure Bi$_2$Te$_3$');
set(h_legend,'Interpreter', 'latex','fontsize',16);
% figure;
% semilogx(ncarrier,G*L/A*0.01,'linewidth',8), xlim([1e12,1e21]), ylabel('Conductivity ($\Omega^{-1}$cm$^{-1}$)','Interpreter', 'latex', 'fontsize',32,'linewidth',2); %xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32);
% hold all;
% semilogx(ncarrier,Gp*L/A*0.01,'linewidth',8), xlim([1e12,1e21]), ylabel('Conductivity ($\Omega^{-1}$cm$^{-1}$)','Interpreter', 'latex', 'fontsize',32,'linewidth',2); %xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32);
% %plot(1.7E19,24970, '.', 'markersize', 70) %This is from "Transport Properties of Si" by Weber and Gmelin (1991)
% set(gca, 'fontsize', 30);
%h_legend= legend('Single twin Si','Pure Si','Experiment');
%set(h_legend, 'fontsize',30);
figure;
semilogx(carrier,S*1e6,':o','linewidth',8), xlim([1e18,1e21]), ylabel('S ($\mu$ V/K)', 'Interpreter', 'latex','fontsize',32,'linewidth',2); %xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
hold all;
semilogx(carrier,Sp*1e6,'linewidth',8), xlim([1e18,1e21]), ylabel('S ($\mu$ V/K)', 'Interpreter', 'latex','fontsize',32,'linewidth',2); %xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
%plot(1.7E19, 400, '.', 'markersize', 70) %This is from "Transport Properties of Si" by Weber and Gmelin (1991)
set(gca, 'fontsize', 30);
h_legend= legend('Twin Bi$_2$Te$_3$','Pure Bi$_2$Te$_3$');
set(h_legend,'Interpreter', 'latex','fontsize',16);

figure;
semilogx(carrier,(ke+kp_exp_nodim)*(L/A),':o','linewidth',8),xlim([1e18,1e21]), ylabel('$\kappa_e$+$\kappa_{ph}$ (W/m-K)','Interpreter', 'latex', 'fontsize',32,'linewidth',2);% xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
hold all;
%semilogx(ncarrier,(ke)*(L/A),'linewidth',8),xlim([1e18,1e21]);% xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
semilogx(carrier,(kep+kp_exp_nodim)*(L/A),'linewidth',8),xlim([1e18,1e21]), ylim([-5,30]), ylabel('$\kappa_e$+$\kappa_{ph}$ (W/m-K)','Interpreter', 'latex', 'fontsize',32,'linewidth',2);% xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
%plot(1.7E19, 148, '.', 'markersize', 70) %This is from "Transport Properties of Si" by Weber and Gmelin (1991)
set(gca, 'fontsize', 30);
h_legend= legend('Twin Bi$_2$Te$_3$','Pure Bi$_2$Te$_3$');
set(h_legend,'Interpreter', 'latex','fontsize',16);
%figure;
%semilogx(ncarrier,ke*(L/Acellm)),xlim([1e18,3e20]), ylabel('K$_e$ (W/cm-K)','Interpreter', 'latex', 'fontsize',32,'linewidth',2), xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',2);
%hold all;
%semilogx(ncarrier,kep*(L/Acellm)),xlim([1e18,3e20]), ylabel('K$_e$ (W/cm-K)','Interpreter', 'latex', 'fontsize',32,'linewidth',2), xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',2);
%set(gca, 'fontsize', 20);
%h_legend= legend('Twin Bi$_2$Te$_3$','Pure Bi$_2$Te$_3$');
%set(h_legend,'Interpreter', 'latex','fontsize',16);

figure;
semilogx(carrier, ZT,':o','linewidth',8),xlim([1e18,1e21]),ylabel('ZT', 'Interpreter', 'latex','fontsize',32,'linewidth',2);%xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
hold all;
semilogx(carrier, ZTp,'linewidth',8),xlim([1e18,1e21]),ylabel('ZT', 'Interpreter', 'latex','fontsize',32,'linewidth',2);%xlabel('n$_e$ (cm$^{-3}$)','Interpreter', 'latex','fontsize',32,'linewidth',8);
%plot(1.7E19, ZT_exp, '.', 'markersize', 70) %This is computed from the above three
set(gca, 'fontsize', 30);
h_legend= legend('Twin Bi$_2$Te$_3$','Pure Bi$_2$Te$_3$');
set(h_legend,'Interpreter', 'latex','fontsize',16);

figure;
plot(E, Ttwin,'linewidth',4),xlim([Emin,Emax]),ylim([0,4]),ylabel('Transmission (2e$^2$/h)', 'Interpreter', 'latex','fontsize',32,'linewidth',2),xlabel('E-E$_F$ (eV)','Interpreter', 'latex','fontsize',32,'linewidth',8);
hold all;
plot(E, Tpure,'linewidth',4),xlim([Emin,Emax]),ylim([0,4]),ylabel('Transmission (2e$^2$/h)', 'Interpreter', 'latex','fontsize',32,'linewidth',2),xlabel('E-E$_F$ (eV)','Interpreter', 'latex','fontsize',32,'linewidth',8);
set(gca, 'fontsize', 30);
h_legend= legend('Twin Bi$_2$Te$_3$','Pure Bi$_2$Te$_3$');
set(h_legend,'Interpreter', 'latex','fontsize',16);
%figure;
%plot(Ej, Ttwin,'linewidth',4),xlim([-1.75,2.25]),ylim([0,4]),ylabel('Transmission (2e$^2$/h)', 'Interpreter', 'latex','fontsize',32,'linewidth',2),xlabel('E-E$_F$ (eV)','Interpreter', 'latex','fontsize',32,'linewidth',8);
%hold all;
%plot(Ej, Tpure,'linewidth',4),xlim([-1.75,2.25]),ylim([0,4]),ylabel('Transmission (2e$^2$/h)', 'Interpreter', 'latex','fontsize',32,'linewidth',2),xlabel('E-E$_F$ (eV)','Interpreter', 'latex','fontsize',32,'linewidth',8);
%plot(E-EF111, DOM111, 'b--','linewidth',2)
%plot(E-EF25, DOM25, 'g-','linewidth',2)
%set(gca, 'fontsize', 20);
%h_legend= legend('Single twin Si','Pure Si','Periodic Twins');
%set(h_legend, 'fontsize',20);
%save ZTjunction.mat