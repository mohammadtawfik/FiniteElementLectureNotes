clear all
close all
clc
%Problem Data
LengthX=1;
LengthY=1;
Nx=2;
Ny=2;
Ne=Nx*Ny;
Nn=(Nx+1)*(Ny+1);
%Length of the elements in x and y-directions
Lx=LengthX/Nx;
Ly=LengthY/Ny;
%Generate the mesh 
[Nodes,Elements,BCs]=GenerateMesh(Nx,Ny,Lx,Ly);
%Creating Empty Global matrix
KGlobal=zeros(Nn,Nn);
%Looping on all elements
for ii=1:Ne
    %Identifying the current element nodes and data
    Node1=Elements(ii,1);
    Node2=Elements(ii,2);    
    Node3=Elements(ii,3);
    Node4=Elements(ii,4);    
    UV=[Node1,Node2,Node3,Node4];
    Xnodes=[Nodes(Node1,1);Nodes(Node2,1);Nodes(Node3,1);Nodes(Node4,1)];
    Ynodes=[Nodes(Node1,2);Nodes(Node2,2);Nodes(Node3,2);Nodes(Node4,2)];
    T1=CalcTxy(Xnodes,Ynodes);
    KK=CalcQuad(T1,Xnodes,Ynodes);
    %Assembling the global stiffness matrix 
    KGlobal(UV,UV)=KGlobal(UV,UV)+KK;
end
BCsC=[1:Nn]';  %Complementary boundary conditions
BCsC(BCs)=[];
%The matrix that will be multiplied by the boundary values
KAux=KGlobal(:,BCs);
KAux(BCs,:)=[];
%Getting reduced stiffness and force vector 
%   by applying boundary conditions
KReduced=KGlobal;
KReduced(BCs,:)=[];
KReduced(:,BCs)=[];
%The right-hand-side vector 
FReduced=-KAux*Nodes(BCs,3);
%Solving for the function values
Displacements=KReduced\FReduced;
%Saving the solution values in the Nodes' registry
Nodes(BCsC,3)=Displacements