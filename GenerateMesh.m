function [Nodes,Elements,BCs]=GenerateMesh(Nx,Ny,Lx,Ly)

Ne=Nx*Ny;
Nn=(Nx+1)*(Ny+1);


%Creating the Nodes' registry
Nodes=zeros(Nn,3); % x,y,function value
BCsCounter=0;
for ii=1:Nx+1
    Xx=(ii-1)*Lx; %x-coordinate
    for jj=1:Ny+1
        Yy=(jj-1)*Ly; %y-coordinate
        NN=(jj-1)*(Nx+1) + ii; %node number
        Nodes(NN,:)=[Xx,Yy,0]; %saving the data
        if and(ii==1,jj==1) %setting the boundary values
            BCsCounter=BCsCounter+1;
            BCs(BCsCounter,1)=NN;
            Nodes(NN,3)=0;
%        elseif or(or(jj==1,jj==Ny+1),ii==Nx+1)
%            BCsCounter=BCsCounter+1;
%            BCs(BCsCounter,1)=NN;
%            Nodes(NN,3)=0;
        endif
    endfor
endfor
%Creating the Elements' registry
Elements=zeros(Ne,4);
for ii=1:Nx
    for jj=1:Ny
        NN=(jj-1)*Nx+ii;
        Elements(NN,1)=(jj-1)*(Nx+1) + ii; %Node1 number
        Elements(NN,2)=(jj-1)*(Nx+1) + ii+1; %Node2 number
        Elements(NN,3)=jj*(Nx+1) + ii+1; %Node3 number
        Elements(NN,4)=jj*(Nx+1) + ii; %Node4 number
    endfor
endfor
    
    
    
endfunction     