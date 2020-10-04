function U = haarUnitarySpinOne

U=sparse(27,27);
 
 for i=1:8
     U(i,i)=randU(1);
 end
 U(9:10,9:10)=randU(2);
 U(11,11)=randU(1);
 
 U(18:19,18:19)=randU(2);
 for i=20:27
     U(i,i)=randU(1);
 end
 
 U(14,14)=randU(1);
 U(17,17)=randU(1);
 
 U(12:13,12:13)=randU(2);
 U(15:16,15:16)=randU(2);



 