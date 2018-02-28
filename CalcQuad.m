function KB=CalcQuad(T1,Xv,Yv)

GaussConstants=zeros(2,8);
GaussConstants(2, 1) = -0.9602898565;
GaussConstants(2, 2) = -0.7966664774;
GaussConstants(2, 3) = -0.5255324099;
GaussConstants(2, 4) = -0.1834346424;
GaussConstants(2, 5) = 0.1834346424;
GaussConstants(2, 6) = 0.5255324099;
GaussConstants(2, 7) = 0.7966664774;
GaussConstants(2, 8) = 0.9602898564;

GaussConstants(1, 1) = 0.1012285362;
GaussConstants(1, 2) = 0.2223810344;
GaussConstants(1, 3) = 0.3137066458;
GaussConstants(1, 4) = 0.3626837833;
GaussConstants(1, 5) = 0.3626837833;
GaussConstants(1, 6) = 0.3137066458;
GaussConstants(1, 7) = 0.2223810344;
GaussConstants(1, 8) = 0.1012285362;  
            KB=zeros(4,4);
Txy=inv([1,-1,-1,1;1,1,-1,-1;1,1,1,1;1,-1,1,-1]);

%Start the numerical integrration procedure
   for Xi=1:8
     xx=GaussConstants(2, Xi);
     for Yi=1:8
         yy=GaussConstants(2, Yi);
         HT=CalcH(xx,yy)*Txy;
         HTx=CalcHx(xx,yy)*Txy;
         HTy=CalcHy(xx,yy)*Txy;
         X = HT*Xv;
         Y = HT*Yv;
         Xx= HTx*Xv;
         Xy= HTy*Xv;
         Yx= HTx*Yv;
         Yy= HTy*Yv;
            %**************************
            cx=CalcHx(X,Y);
            cy=CalcHy(X,Y);
            %**************************
            Kb= ((cx'*cx)+(cy'*cy));
            %performing the weighted summation
            KB=KB+GaussConstants(1,Xi)*GaussConstants(1,Yi)*Kb*(Xx*Yy+Xy*Yx);
            %End of Calculation loop body
      end
   end
   KB= T1'*KB*T1;
