%Clearing the memory and display
clear all
clc
%Problem Data
NE=2; %number of elements
Length=2.0; %beam length
Width=0.02; %beam width
Thickness=0.01; %beam thickness
Modulus=71e9; %Modulus of Elasticity Aluminum (GPa)
Rho=2700; %Density (Kg/m^3)
Alpha=22.5e-6; %Thermal Expansion coefficientt
%Cross-section area
Area=Width*Thickness;
%Second moment of area
Imoment=Width*Thickness*Thickness*Thickness/12;
Le=Length/NE; %Element Length
%Element stiffness matrix
Ke=Modulus*Imoment*[12,6*Le,-12,6*Le; ...
          6*Le,4*Le*Le,-6*Le,2*Le*Le; ...
          -12,-6*Le,12,-6*Le; ...
          6*Le,2*Le*Le,-6*Le,4*Le*Le]/Le/Le/Le;
%Global stiffness and mass matrix assembly
%Initializing an empty matrix
KGlobal=zeros(2*(NE+1),2*(NE+1));
%Assembling the global matrix
for ii=1:NE
    KGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
                  KGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Ke;
end
%For a cantilever beam the first and second degree of freedom are fixed
%Obtaining the auxiliary equations
KAux=KGlobal(1:2, 3:2*(NE+1));
%Applying the boundary conditions
KGlobal(1:2,:)=[];
KGlobal(:,1:2)=[];
%force Vector
FGlobal=zeros(2*NE,1); %This is the empty force fector
FGlobal(2*NE-1)=1; %Adding a single point load at the tip
%Obtainning the solution displacement field
WW=inv(KGlobal)*FGlobal
%Obtaining the support reaction
Reactions= KAux*WW