function TbInv=CalcTxy (Xv,Yv)
    
    T1=CalcH(Xv(1),Yv(1));
    T2=CalcH(Xv(2),Yv(2));
    T3=CalcH(Xv(3),Yv(3));
    T4=CalcH(Xv(4),Yv(4));
    
    TT=[T1;T2;T3;T4];
    TbInv=inv(TT);
endfunction