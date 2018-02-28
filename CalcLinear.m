function KB=CalcLinear(TbInv,LengthX,LengthY)

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

%Start the numerical integrration procedure
   for Xi=1:8
     X = LengthX * (GaussConstants(2, Xi) + 1) / 2;
     for Yi=1:8
         Y = LengthY * (GaussConstants(2, Yi) + 1) / 2;
            %**************************
            cx=CalcHx(X,Y);
            cy=CalcHy(X,Y);
            %**************************
            Kb= (cx'*cx)+(cy'*cy);
            %performing the weighted summation
            KB=KB+GaussConstants(1,Xi)*GaussConstants(1,Yi)*Kb;
            %End of Calculation loop body
      end
   end
   KB= TbInv'*KB*TbInv;
   %Multiplying the resulting matreces by Jacobian
   KB = KB * LengthX * LengthY/ 4;
