function TbInv=CalcTRect (Lx,Ly)
    
    T1=CalcH(0,0);
    T2=CalcH(Lx,0);
    T3=CalcH(Lx,Ly);
    T4=CalcH(0,Ly);
    
    TT=[T1;T2;T3;T4];
    TbInv=inv(TT);
endfunction