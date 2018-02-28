clear all
clc
close all

%Problem Data
EE=71e9;
AA=0.01*0.01;
II=AA*0.01*0.01/12;
NN=6;
NE=9;

% X,Y,Fx,Fy,M, u, v, v'
Nodes=[0       ,0          ,0   ,0   ,0   ,1  ,2 ,3 ;
       2.0     ,0          ,0   ,0   ,0   ,4  ,5 ,6 ;
       2.5     ,0.5*sqrt(3),-100,0   ,0   ,7  ,8 ,9 ;
       1.5     ,0.5*sqrt(3),0   ,0   ,0   ,10 ,11,12;
       0.5     ,0.5*sqrt(3),0   ,0   ,0   ,13 ,14,15;
       1.0     ,0          ,0   ,0   ,0   ,16 ,17,18];

% Node#1, Node#2, Modulus of Elasticity , Area, Second moment of area 
Elements=[1 ,5 ,EE ,AA,II;
          1 ,6 ,EE ,AA,II;
          5 ,4 ,EE ,AA,II;
          4 ,3 ,EE ,AA,II;
          6 ,2 ,EE ,AA,II;
          2 ,3 ,EE ,AA,II;
          2 ,4 ,EE ,AA,II;
          6 ,4 ,EE ,AA,II;
          6 ,5 ,EE ,AA,II];

BCs=[1,2,5]; %Fixed degrees of freedom
%Creating Empty Global matrix and force vector 
KGlobal=zeros(3*NN);
FGlobal=zeros(3*NN,1);

%Looping on all elements
for ii=1:NE
    %Identifying the current element nodes and data
    Node1=Elements(ii,1);
    Node2=Elements(ii,2);    
    Stiff=Elements(ii,3);
     Area=Elements(ii,4);
      IIe=Elements(ii,5);
    %Getting the nodes' coordinates 
       x1=Nodes(Node1,1);
       x2=Nodes(Node2,1);
       y1=Nodes(Node1,2);
       y2=Nodes(Node2,2);
    %Getting the nodes' degrees of freedom
       u1=Nodes(Node1,6);
       v1=Nodes(Node1,7);
       t1=Nodes(Node1,8);
       u2=Nodes(Node2,6);
       v2=Nodes(Node2,7);
       t2=Nodes(Node2,8);
       UV=[u1,v1,t1,u2,v2,t2];
    %Creating the stiffness matrix of each element
    LL=sqrt((x2-x1)^2+(y2-y1)^2); %Element length 
    CC=(x2-x1)/LL; %Cos(Theta)
    SS=(y2-y1)/LL; %Sin(Theta)
    %Rotation MAtrix
    TT(:,:,ii)=[CC,-SS,0 , 0, 0 ,0  ;
                SS, CC,0 , 0, 0 ,0  ;
                0 , 0 ,1 , 0, 0 ,0  ;
                0 , 0 ,0 ,CC,-SS,0  ;
                0 , 0 ,0 ,SS, CC,0
                0 , 0 ,0 , 0, 0 ,1  ];
    %Element stiffness matrix in local coordinates
    Kl=[Area, 0           , 0       ,-Area, 0            , 0        ;
         0  , 12*IIe/LL/LL, 6*IIe/LL, 0   , -12*IIe/LL/LL, 6*IIe/LL ;
         0  , 6*IIe/LL    , 4*IIe   , 0   , -6*IIe/LL    , 2*IIe    ;
       -Area, 0           , 0       , Area, 0            , 0        ;
         0  ,-12*IIe/LL/LL,-6*IIe/LL, 0   ,  12*IIe/LL/LL,-6*IIe/LL ;
         0  , 6*IIe/LL    , 2*IIe   , 0   , -6*IIe/LL    , 4*IIe    ]*Stiff/LL;
    %Element stiffness matrix in GLOBAL coordinates
    KK(:,:,ii)=TT(:,:,ii)'*Kl*TT(:,:,ii);
    %Assembling the global stiffness matrix 
    KGlobal(UV,UV)=KGlobal(UV,UV)+KK(:,:,ii);
end
%Filling the global force vector from nodes' data
for ii=1:NN
    FGlobal(Nodes(ii,6))=Nodes(ii,3);
    FGlobal(Nodes(ii,7))=Nodes(ii,4);
    FGlobal(Nodes(ii,8))=Nodes(ii,5);
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
BCsC=[1:3*NN]';  %Complementary boundary conditions
BCsC(BCs)=[];
Displ=zeros(3*NN,1);
%Storing the resulting displacements in their respective location
Displ(BCsC)=Displacements; 

%Looping on the elements
for ii=1:NE
    Node1=Elements(ii,1);
    Node2=Elements(ii,2);    
       u1=Nodes(Node1,6);
       v1=Nodes(Node1,7);
       t1=Nodes(Node1,8);
       u2=Nodes(Node2,6);
       v2=Nodes(Node2,7);
       t2=Nodes(Node2,8);
       UV=[u1,v1,t1,u2,v2,t2];
    %Evaluating the loval force vector 
    ff(:,ii)=TT(:,:,ii)*KK(:,:,ii)*Displ(UV);
end
%Element foces are tension when f1 is negative
ElementForces=-ff(1,:)'
