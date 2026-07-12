function [fl,fh,win,range]=alignband(i)
    % [fl,fh]:  frequency band
    % win:      window length
    % range:    maximum shift allowed
    switch i
    	case 1
            fl=0.1;     fh=0.25;  	win=30 ;    range=5;
     	case 2
            fl=0.25;   	fh=0.5; 	win=15 ;    range=0.6;
      	case 3
            fl=0.5;  	fh=1.0;   	win=8 ;    range=0.1;
     	case 4
            fl=0.5;    	fh=1;    	win=8 ;     range=0.1;
       	case 5
            fl=1;    	fh=4;     	win=6;     range=0.05;
        otherwise 
        	fprintf('error\n out of boundary');
    end
    fprintf('   fl=%f    fh=%f  win=%f  range=%f\n',fl,fh,win,range);
end