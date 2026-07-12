 function [Sxy] = cmtm(X,NW)

%  midway = floor(Nx/2)+1;
%  Xf=zeros(m,Nx);
 %CMTM              Coherence function estimate using the multitaper method.
 %   Cxy = CMTM(X,Y) estimates the coherence between two equal length
 %   vectors vectors X and Y using Thomson's multitaper method.  The
 %   coherence is a complex function of frequency, with magnitude between 0
 %   and 1, that estimates
 %            E[X*(f) Y(f)] / sqrt(E[X*(f) X(f)] E[Y*(f) Y(f)])
 %   where * denotes complex conjugation.  **Note that this is different
 %   from the matlab function COHERE, which returns the magnitude squared
 %   of this value.**
 %
 %   Cxy = CMTM(X,Y,NW) specifies the time-bandwidth product for the
 %   discrete prolate spheroidal sequences (DPSS) is specified by NW; this
 %   value also determines the number of tapers as (2*NW-1).  If not
 %   given, the default NW is 4.
 %
 %   Cxy = CMTM(X,Y,NW,Fs) specifies a sampling frequency Fs.
 %
 %   [Cxy,F] = CMTM(...) also returns the vector of frequencies at which
 %   the coherence is estimated.  If Fs is specified, this vector ranges
 %   from [0,Fs/2] and the units are the same as Fs.  If Fs is not
 %   specified, F ranges from [0,1/2] and the units are relative to the
 %   sampling frequency.
 %
 %   CMTM(...) without output arguments plots the magnitude-squared
 %   and phase of the coherence (in two subplots) in the current figure.
 
 %%%%% Debugging test code:
 % Fs = 100;
 % t = 0:1/Fs:(32-1/Fs);
 % X = sin(2*pi*5*t);  Y = sin(2*pi*5*t) + 0.2*randn(size(X));
 % cmtm(X,Y,4,Fs);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (nargin < 2), error('Three arguments are required.');  end
 P = 2*NW;
 [m Nx]=size(X);
 Xf=zeros(P,Nx,m);
 Pxx=zeros(P,Nx,m);
 Sxy=zeros(Nx,m,m);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make Tapers %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [E,V] = dpss(Nx,NW,P);
%  E=E-min(E);
% if P==1
% 	tm=round(Nx*0.1);
% 	E=ones(Nx,1);
% 	E(1:tm)=cos(linspace(pi,2*pi,tm));
% 	E(end-tm+1:end)=cos(linspace(0,pi,tm));
%     E=(E+1)/2;
% 
% end

if P==1
  	E=tukeywin(Nx,0.5);
end
% E=ones(Nx,P);
%  figure(11);
%  plot(E);
 E = E';
 
 %%%%%%%%%%%%%%%%%%%%%%%%%% Calculate Coherence %%%%%%%%%%%%%%%%%%%%%%%

 for i=1:m
  X(i,:)=X(i,:)-mean(X(i,:));
  
  Xf(:,:,i) = fft(E .* repmat(X(i,:),P,1),[],2);  
%   Xf = Xf(:,1:midway);   
  Pxx(:,:,i) = conj(Xf(:,:,i)).*Xf(:,:,i); 
 end
%  disp('x')
%  X(1,1:5)
%  disp('xf')
%  Xf(1,1:5,1)
%  disp('dataout')
%  size(E(1,:))
%  size(X(1,:))
%  E(1,:).*X(1,:)
%  disp('fft')
%  fft(  E(1,:).*X(1,:))
 for i=1:m
     for j=1:m
%       Sxy(:,i,j) = conj(mean((conj(Xf(:,:,i)) .* Xf(:,:,j)) ./ sqrt((Pxx(:,:,i).*Pxx(:,:,j))),1));
        Sxy(:,i,j) = mean((Xf(:,:,i) .* conj(Xf(:,:,j))),1);

     end
 end
%  Cxy = mean(Cxy,1);
%  
 
 end
