clear
clc

Nb=15;
%% Problem Definition
global NFE;
NFE=0;
model=CreateModel(Nb);
CostFunction=@(p,model,Nb) MyNPW(p,model,Nb); % Cost Function
%% Self PSO

nVar=model.N; % Number of Decision Variables                                                                 %???????????????????
T=3;
VarSize=[1 nVar];
VarMin=model.ps(1,:);
VarMax=model.ps(2,:);
%% PSO Parameters
MaxIt=100;
nPop=5;
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;
% wdamp=1;
% c1=chi*phi1;   
% c2=chi*phi2;   
% VelMax=0.1*(VarMax-VarMin);
% VelMin=-VelMax;
% w=chi;         
% wdamp=1;       
% c1=chi*phi1;   
% c2=chi*phi2;   
w_min=9.5;
w_max=9.5;
w=w_max;
c1=5;
c2=5;
% VelMax=0.1*(VarMax-VarMin);
% VelMin=-VelMax;
VelMax=6*2000;
VelMin=0;

% P_PVmin=1;
% P_PVmax=2000;
% P_WTmin=1;
% P_WTmax=2000;
% P_Bamin=-0.3203;                                                       %???????????????????
% P_Bamax=0.3203;
% P_DG1min=1;
% P_DG1max=2000;
% P_DG2min=1;
% P_DG2max=2000;
% P_UGmin=1;
% P_UGmax=20000;

% P_PVmin=1;
% P_PVmax=200000;
% P_WTmin=1;
% P_WTmax=400000;
% P_Bamin=-500*0.3203;                                                       %???????????????????
% P_Bamax=500*0.3203;
% P_DG1min=1;
% P_DG1max=80000000;
% P_DG2min=1;
% P_DG2max=80000000;
% P_UGmin=1;
% P_UGmax=200000000;


% ps_min=[0.9*ones(1,Nb),(-180)*ones(1,Nb),P_PVmin,P_WTmin,P_Bamin,P_DG1min,P_DG2min,P_UGmin];
% ps_max=[1.1*ones(1,Nb),180*ones(1,Nb),P_PVmax,P_WTmax,P_Bamax,P_DG1max,P_DG2max,P_UGmax];

% V_magmin=0.9;
% V_magmax=1.1;
% V_angmin=-180;
% V_angmax=180;
%% Initialization
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Sol=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
% empty_particle.Valuefunction=[];
empty_particle.Best.Sol=[];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=-inf;

% qbest_Fitaux=zeros(nPop,1,T);
% sbest_Fitaux=zeros(1,1,T);

ind_T=1;

for i=1:nPop
%     particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
%     particle(i).Position(1:Nb)=ps_min+(ps_max-ps_min).*rand(1,nVar);
%     particle(i).Position(1:Nb)=V_magmin+(V_magmax-V_magmin)*rand(1,Nb,1);
%     particle(i).Position((Nb+1):(2*Nb))=V_angmin+(V_angmax-V_angmin)*rand(1,Nb,1);
%     particle(i).Position((2*Nb+1))=P_PVmin+(P_PVmax-P_PVmin)*rand(1,1,1);
%     particle(i).Position((2*Nb+2))=P_WTmin+(P_WTmax-P_WTmin)*rand(1,1,1);
%     particle(i).Position((2*Nb+3))=P_Bamin+(P_Bamax-P_Bamin)*rand(1,1,1);
%     particle(i).Position((2*Nb+4))=P_DG1min+(P_DG1max-P_DG1min)*rand(1,1,1);
%     particle(i).Position((2*Nb+5))=P_DG2min+(P_DG2max-P_DG2min)*rand(1,1,1);
%     particle(i).Position((2*Nb+6))=P_UGmin+(P_UGmax-P_UGmin)*rand(1,1,1);
%     particle(i).Position((2*Nb+7))=P_PVmin+(P_PVmax-P_PVmin)*rand(1,1,1);
%     particle(i).Position((2*Nb+8))=P_WTmin+(P_WTmax-P_WTmin)*rand(1,1,1);
%     particle(i).Position((2*Nb+9))=P_Bamin+(P_Bamax-P_Bamin)*rand(1,1,1);
%     particle(i).Position((2*Nb+10))=P_DG1min+(P_DG1max-P_DG1min)*rand(1,1,1);
%     particle(i).Position((2*Nb+11))=P_DG2min+(P_DG2max-P_DG2min)*rand(1,1,1);
%     particle(i).Position((2*Nb+12))=P_UGmin+(P_UGmax-P_UGmin)*rand(1,1,1);
    particle(i).Position=CreateRandomSolution_1(model,Nb);
    particle(i).Velocity=zeros(VarSize);
%     [particle(i).Cost,particle(i).Valuefunction]=Fitness_Val(particle(i).Position);
    [particle(i).Cost, particle(i).Sol]=CostFunction(particle(i).Position,model,Nb);
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.Sol=particle(i).Sol;
    if particle(i).Best.Cost>GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
end
BestCost=zeros(MaxIt,1);
nfe=zeros(MaxIt,1);
%% PSO Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        particle(i).Position = particle(i).Position + particle(i).Velocity; 
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
%         [particle(i).Cost,particle(i).Valuefunction]=CostFunction(particle(i).Position);
        [particle(i).Cost, particle(i).Sol]=CostFunction(particle(i).Position,model,Nb);
        if particle(i).Cost>particle(i).Best.Cost
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Sol=particle(i).Sol;
            if particle(i).Best.Cost>GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
        end
    end
    BestCost(it)=GlobalBest.Cost;
%     BestCost(it)= particle(i).Best.Sol;
%     nfe(it)=nPop*(it+1);
    nfe(it)=NFE;
    
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost = ' num2str(BestCost(it))]);
%     w=w*wdamp;
    w=w_max-((w_max-w_min)/MaxIt)*it;
%     Result=particle.Valuefunction;
end

% disp([' Time ' num2str(toc)])

figure(1)
plot(BestCost)
xlabel('Iteration');
ylabel('Best Net Present Worth');
% legend('');
title('Output');

Final_a=[(GlobalBest.Position(1,31))/2;(GlobalBest.Position(1,32))/4;(GlobalBest.Position(1,33))/65.85;(GlobalBest.Position(1,34));(GlobalBest.Position(1,35));(GlobalBest.Position(1,43))/51.35;(GlobalBest.Position(1,44))/599;(GlobalBest.Position(1,45))/20.5497]