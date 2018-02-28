clear all
close all
clc
%Problem Data
LengthX=1;
LengthY=1;
Nx=10;
Ny=10;
Ne=Nx*Ny;
Nn=(Nx+1)*(Ny+1);
%Length of the elements in x and y-directions
Lx=LengthX/Nx;
Ly=LengthY/Ny;
%Evaluating the stiffness matrix 
T1=CalcTRect(Lx,Ly);
%KK=CalcRect(T1,Lx,Ly);
KK=[2*(Lx*Lx+Ly*Ly) ,   Lx*Lx-2*Ly*Ly,-  Lx*Lx-Ly*Ly  ,-2*Lx*Lx+Ly*Ly  ; ...
       Lx*Lx-2*Ly*Ly, 2*(Lx*Lx+Ly*Ly),-2*Lx*Lx+Ly*Ly  ,-  Lx*Lx-Ly*Ly  ; ...
    -  Lx*Lx-Ly*Ly  ,-2*Lx*Lx+Ly*Ly  , 2*(Lx*Lx+Ly*Ly),   Lx*Lx-2*Ly*Ly; ...
    -2*Lx*Lx+Ly*Ly  ,-  Lx*Lx-Ly*Ly  ,   Lx*Lx-2*Ly*Ly, 2*(Lx*Lx+Ly*Ly)]/Lx/Ly/6;
[Nodes,Elements,BCs]=GenerateMesh(Nx,Ny,Lx,Ly);
%Creating Empty Global matrix
KGlobal=zeros(Nn,Nn);
FGlobal=zeros(Nn,1);
%Looping on all elements
for ii=1:Ne
    FF=zeros(4,1);
    if ii==1
        FF=0*Lx/2*[0;0;1;1];
    elseif ii==Nx
        FF=Lx/2*[1;1;0;0];
    endif
    %Identifying the current element nodes and data
    Node1=Elements(ii,1);
    Node2=Elements(ii,2);    
    Node3=Elements(ii,3);
    Node4=Elements(ii,4);    
    UV=[Node1,Node2,Node3,Node4];
    %Assembling the global stiffness matrix 
    KGlobal(UV,UV)=KGlobal(UV,UV)+KK;
    FGlobal(UV)=FGlobal(UV)+FF;
end
BCsC=[1:Nn]';  %Complementary boundary conditions
BCsC(BCs)=[];
%The matrix that will be multiplied by the boundary values
KAux=KGlobal(:,BCs);
KAux(BCs,:)=[];
FAux=FGlobal(BCs);

%Getting reduced stiffness and force vector 
%   by applying boundary conditions
KReduced=KGlobal;
KReduced(BCs,:)=[];
KReduced(:,BCs)=[];
FReduced=FGlobal;
FReduced(BCs)=[];
%The right-hand-side vector
 
FReduced=FReduced-KAux*Nodes(BCs,3);
%FReduced((Ny+1)*(Nx+1)-1)=1;
%Solving for the function values
Displacements=KReduced\FReduced;
%Saving the solution values in the Nodes' registry
Nodes(BCsC,3)=Displacements;

for ii=1:Nx+1
    Xx=(ii-1)*Lx; %x-coordinate
    for jj=1:Ny+1
        Yy=(jj-1)*Ly; %y-coordinate
        NN=(jj-1)*(Nx+1) + ii; %node number
        aa(ii,jj)=Nodes(NN,3); %saving the data
        
    endfor
endfor
contour(aa)
grid