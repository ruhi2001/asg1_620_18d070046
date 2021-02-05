%%Parameters

%constants
q = 1.6e-19;             %charge
eps_0 = 8.85e-12; 
VkBT = 26e-3;            %KBT/q at 300k

%semiconductor
Nsub = -1e18*1e6;        %negative for NMOS and positive for PMOS
k_si = 12;               %Dielectric Constant
ni = 1.5e10*1e6;         %Intriensic carrier concentration
Eg = 1.1*q;              %Bandgap
L = 1e-6;                %substrate thickness
eps_si = k_si*eps_0;     %Permitivity

if Nsub < 0
    psub = abs(Nsub); nsub = ni^2/abs(Nsub);  
elseif Nsub > 0
    psub = ni^2/Nsub; nsub = Nsub;
else
    psub = ni; nsub = ni;
end

%Oxide
tox = 10e-9;            %Thickness
k_ox = 4;               %Dielectric constant
Eg_ox = 9*q;            %Bandgap
eps_ox = k_ox*eps_0;    %Permittivity
Ec_off = 3.1*q;         %Conduction band offset
C_ox = eps_ox/tox;

%Flatband/Metal related
tm = 10e-9;                                    %Thickness
phi_b = -sign(Nsub)*VkBT*log(abs(Nsub)/ni);   %(Ei - Ef)/q of semiconductor
Vfb = sign(Nsub)*Eg/(2*q) - phi_b;            %n+/nmos, p+/pmos

%solver
dx = 1e-9;    %meshing size
NRmax = 10;   %Maximum NR iterations
tol = 1e-6;    %Error tolerance in potential

%% Input

 Vg = -0.3;
 
 if Vg < Vfb
     F = @(psi_v) Vfb - Vg + psi_v -(1/C_ox)*sqrt(2*eps_si*VkBT*q)*sqrt(psub*(exp(-psi_v/VkBT)+(psi_v/VkBT)-1) + nsub*(exp(psi_v/VkBT)-(psi_v/VkBT)-1));
     psi_s = fsolve(F,-0.5);
 elseif Vg >Vfb
     F = @(psi_v) Vfb - Vg + psi_v +(1/C_ox)*sqrt(2*eps_si*VkBT*q)*sqrt(psub*(exp(-psi_v/VkBT)+(psi_v/VkBT)-1) + nsub*(exp(psi_v/VkBT)-(psi_v/VkBT)-1));
     psi_s = fsolve(F,0.5);
 else
     psi_s = phi_b;
 end
 
 

%% mesh
x = 0:dx:L;       %semiconductor mesh
N = length(x);    %number of points

 
%% Newton Raphson - Matrices
A = (-2*diag(ones(1,N),0) + 1*diag(ones(1,N-1),1) + 1*diag(ones(1,N-1),-1));
A(1,:) = 0; A(1,1) = 1;
A(N,:) = 0; A(N,N) = 1;

b = zeros(N,1);
b(1) = psi_s; b(N) = 0;

psi = zeros(N,1);


%% Newton-Raphson - Iterations

for i = 1:NRmax
    
    %Generate b
    p = psub*exp(-psi/VkBT);        %psi(x), p(psi)
    n = nsub*exp(psi/VkBT);
    rho = q*(Nsub +p - n);          %Nsub = Nd-Na
    b = -rho/eps_si*dx^2;
    b(1) = psi_s;  b(N) = 0;        %Dirichlet Boundary Conditions
    f = A*psi -b;
 
    %Jacobian calculations
    delp = -1/VkBT*p;
    deln = 1/VkBT*n;
    delrho = q*(delp - deln);
    delb = -delrho/eps_si*dx^2;
    delb(1) = 0;  delb(N) = 0;
    J = A - diag(delb);              %delb gets added yo diagonal entries
    dV = -J\f;                       %potential update across x
    
    if max(abs(dV)) < tol
        break;
    end
    
    psi = psi + dV;
end


%% Addition of oxide

Es = (psi(1) - psi(2))/dx;          %E = -dpsi/dx
Eox = k_si/k_ox*Es;                 %in absence of oxide charges
xox = -tox:dx:0;                    %oxide mesh
psiox = psi(1) - Eox*xox;           %Potential vector - Continuity
Vox = Eox*tox;                      %Drop across oxide


%% Band Diagram
Ei = -q*psi;                 %Ei in bulk is taken as zero reference
Ec = Ei + Eg/2;
Ev = Ei - Eg/2;
Ef = (Ei(end) - q*phi_b)*ones(size(Ei));  %Because phi_b = (Ei-Ef)/q
a = Ef(1,1)
b = Ei(1,1)
Ecox = -q*psiox + Eg/2 + Ec_off;
Evox = Ecox - Eg_ox;


Em = Ef(end) - q*Vg;                      %Because q*Vg = En -Ef
xm = -tox - tm:dx:-tox;                   %Metal mesh
Em = Em*ones(size(xm));

figure(1);
plot([xm xox x]*1e9, [Em'; Ecox'; Ec]/q,'r',[xm xox x]*1e9,[Em';Evox';Ev]/q,'b',x*1e9,Ei/q,'m--',x*1e9,Ef/q,'k');
xlim([-30,150]);%ylim([-7,4]);
axes = gca;
axes.LineWidth = 1;
axes.FontSize = 18;
axes.FontWeight = 'bold';
axes.Box = 'on';
xlabel('position(nm)');
ylabel('Energy(eV)');
lines=axes.Children;
set(lines,'LineWidth', 2);
leg = legend('E_C','E_V','E_i','E_F');
leg.Orientation = 'horizontal';

%% Carrier Density plot
figure(2);
p = psub*exp(-psi/VkBT);
n = nsub*exp(psi/VkBT);
semilogy(x*1e9,n*1e-6,'r',x*1e9,p*1e-6,'b')
axes = gca;
axes.LineWidth = 1;
axes.FontSize = 18;
axes.FontWeight = 'bold';
axes.Box = 'on';
axes.YTick = [1 1e5 1e10 1e15 1e20];
axes.YTickLabel = {'10^0' '10^{5}' '10^{10}' '10^{15}' '10^{20}'};
leg = legend('n','p');
leg.Orientation = 'horizontal';
xlim([0 150])
xlabel('position(nm)');
ylabel('carrier density(/cm^3)');
lines=axes.Children;
set(lines,'LineWidth', 2);

