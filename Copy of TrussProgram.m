clear all
clc
close all

%Problem Data
EE=71e9;
AA=0.01*0.01;
NN=6;
NE=9;

% X,Y,Fx,Fy, u, v
Nodes=[0       ,0          ,0   ,0   ,1  ,2 ;
       2.0     ,0          ,0   ,0   ,3  ,4 ;
       2.5     ,0.5*sqrt(3),-100,0   ,5  ,6 ;
       1.5     ,0.5*sqrt(3),0   ,0   ,7  ,8 ;
       0.5     ,0.5*sqrt(3),0   ,0   ,9  ,10;
       1.0     ,0          ,0   ,0   ,11 ,12];

%      Node#1, Node#2, Modulus of Elasticity , Area
Elements=[1 ,5 ,EE ,AA;
          1 ,6 ,EE ,AA;
          5 ,4 ,EE ,AA;
          4 ,3 ,EE ,AA;
          6 ,2 ,EE ,AA;
          2 ,3 ,EE ,AA;
          2 ,4 ,EE ,AA;
          6 ,4 ,EE ,AA;
          6 ,5 ,EE ,AA];

BCs=[1,2,4]; %Fixed degrees of freedom
%Creating Empty Global matrix and force vector 
KGlobal=zeros(2*NN);
FGlobal=zeros(2*NN,1);

%Looping on all elements
for ii=1:NE
    %Identifying the current element nodes and data
    Node1=Elements(ii,1);
    Node2=Elements(ii,2);    
    Stiff=Elements(ii,3);
     Area=Elements(ii,4);
    %Getting the nodes' coordinates 
       x1=Nodes(Node1,1);
       x2=Nodes(Node2,1);
       y1=Nodes(Node1,2);
       y2=Nodes(Node2,2);
    %Getting the nodes' degrees of freedom
       u1=Nodes(Node1,5);
       v1=Nodes(Node1,6);
       u2=Nodes(Node2,5);
       v2=Nodes(Node2,6);
       UV=[u1,v1,u2,v2];
    %Creating the stiffness matrix of each element
    LL=sqrt((x2-x1)^2+(y2-y1)^2); %Element length 
    CC=(x2-x1)/LL; %Cos(Theta)
    SS=(y2-y1)/LL; %Sin(Theta)
    %Rotation MAtrix
    TT(:,:,ii)=[CC,-SS,0 ,0  ;
        SS,CC ,0 ,0  ;
        0 ,0  ,CC,-SS;
        0 ,0  ,SS,CC];
    %Element stiffness matrix in local coordinates
    Kl=[1,0,-1,0; 0,0,0,0; -1,0,1,0; 0,0,0,0]*Stiff*Area/LL;
    %Element stiffness matrix in GLOBAL coordinates
    KK(:,:,ii)=TT(:,:,ii)'*Kl*TT(:,:,ii);
    %Assembling the global stiffness matrix 
    KGlobal(UV,UV)=KGlobal(UV,UV)+KK(:,:,ii);
end
%Filling the global force vector from nodes' data
for ii=1:NN
    FGlobal(Nodes(ii,5))=Nodes(ii,3);
    FGlobal(Nodes(ii,6))=Nodes(ii,4);
end

%Separating Auxiliary equations
KAux=KGlobal(BCs,:);
KAux(:,BCs)=[];
%Getting reduced stiffness and force vector 
%   by applying boundary conditions
KReduced=KGlobal;
KReduced(BCs,:)=[];
KReduced(:,BCs)=[];
FReduced=FGlobal;
FReduced(BCs)=[];

%Solving for the displacements
Displacements=KReduced\FReduced
%Evaluating the reactions from the auxiliary equations
Reactions=KAux*Displacements

%To obtain the element forces ...
BCsC=[1:2*NN]';  %Complementary boundary conditions
BCsC(BCs)=[];
Displ=zeros(2*NN,1);
%Storing the resulting displacements in their respective location
Displ(BCsC)=Displacements; 

%Looping on the elements
for ii=1:NE
    Node1=Elements(ii,1);
    Node2=Elements(ii,2);    
       u1=Nodes(Node1,5);
       v1=Nodes(Node1,6);
       u2=Nodes(Node2,5);
       v2=Nodes(Node2,6);
       UV=[u1,v1,u2,v2];
    %Evaluating the loval force vector 
    ff(:,ii)=TT(:,:,ii)*KK(:,:,ii)*Displ(UV);
end
%Element foces are tension when f1 is negative
ElementForces=-ff(1,:)'
