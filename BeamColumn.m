%Clearing the memory and display
clear all
clc
%Problem Data
NE=5; %number of elements
Length=1.0; %beam length
Width=0.02; %beam width
Thickness=0.001; %beam thickness
Modulus=71e9; %Modulus of Elasticity Aluminum (GPa)
Rho=2700; %Density (Kg/m^3)
Alpha=22.5e-6; %Thermal Expansion coefficientt
%Cross-section area
Area=Width*Thickness;
%Second moment of area
Imoment=Width*Thickness*Thickness*Thickness/12;
Le=Length/NE; %Element Length
%Element stiffness matrix
Ke=Modulus*Imoment*[12     ,6*Le    ,-12  ,6*Le; ...
                    6*Le   ,4*Le*Le ,-6*Le,2*Le*Le; ...
                    -12    ,-6*Le   ,12   ,-6*Le; ...
                    6*Le   ,2*Le*Le ,-6*Le,4*Le*Le]/Le/Le/Le;
%Element Geometricstiffness matrix
Kg=[36  ,3*Le   ,-36  ,3*Le; ...
    3*Le,4*Le*Le,-3*Le,-Le*Le; ...
    -36 ,-3*Le  ,36   ,-3*Le; ...
    3*Le,-Le*Le ,-3*Le,4*Le*Le]/30/Le;
%Global matrices assembly
%Initializing empty matrices
KGlobal=zeros(2*(NE+1),2*(NE+1));
KGGlobal=zeros(2*(NE+1),2*(NE+1));
%Assembling the global matrix
for ii=1:NE
    KGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
                  KGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Ke;
    KGGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
                  KGGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Kg;
end
%Applying the boundary conditions
BCS=[1,2*(NE+1)-1]; %Simply Supported BCS
KGlobal(BCS,:)=[];
KGlobal(:,BCS)=[];
KGGlobal(BCS,:)=[];
KGGlobal(:,BCS)=[];
%Evaluating the eigenvalues
sort(eig(inv(KGGlobal)*KGlobal))