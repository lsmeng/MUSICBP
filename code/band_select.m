function [fl fh ff]=band_select(i)
    switch i
       case 1
        fl=0.1;     fh=0.25;        ff=23;
       case 2
         fl=0.25;  fh=1.0;      ff=20;         
         case 3
         fl=0.5;   fh=1;  ff=18;
         case 4
         fl=0.5;    fh=2;   ff=16;
         case 5
         fl=1;        fh=4;  ff=14; 
        otherwise 
              fprintf('error\n out of boundary');
    end
    fprintf('   fl=%f    fh=%f  ff=%f\n',fl,fh,ff);
end