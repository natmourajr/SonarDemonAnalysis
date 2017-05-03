function B = eeA(x,maxtau,epsilon,thrsk,condM2)
%  
%  $Date: 2001/05/07 04:34:57 $
%  $Revision: 1.4 $
%  Extreme Event Analysis on real field 
%   based on non-holonomic nested newton's methods for ica 
%     
%
%
% B:       'demixing' matrix
%
% x:       n*t matrix whose columns are observation vectors.
%
% epsilon: threshold to stop
%
% maxtau:  maximum iteration time
%
% condM2:  maximum admissible condition number
% 
% thrsk:   Small number added to stabilize the algorithm.
%          We can not detect sources 
%          with kurtoses less than (maximum kurtosis)*thrsk.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%
%
% [simplest usage]                  
%                                 
%   B=eeA(x)                      
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Please cite the article,
%
%  @InProceedings{akuzawa10,
%  author = 	 {T.Akuzawa},
%  title = 	 {Extended quasi-Newton method for the ICA},
%  booktitle = 	 {Proceedings of the International Workshop on
%  Independent Component Analysis and  Blind Signal Separation},
%  year =	 2000,
%  pages =	 {521-525}
%   }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  
%   (c) AKUZAWA, toshinao
%    akuzawa@brain.riken.go.jp
%    http://www.mns.brain.riken.go.jp/~akuzawa/  
%  

  
if nargin==0, error('Not enough input arguments.'); end
if ndims(x)>2, error('Inputs should be two dimensional!'); end

%%%%%%%%%%%%%%
% data size
[n,t] = size(x);


if nargin>6, error('Too many input arguments.'); end
if nargin<5, condM2=1e+3;  end  
if nargin<4, thrsk=1e-4;  end  
if nargin<3, epsilon=5e-6; end  
if nargin<2, maxtau=1000; end  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subtract mean from data

mzx1=x-mean(x.').'*ones(1,t);
%size(mzx1)
%%%%%%%%%%%%%%%%%
% Perhaps we have to introduce pre-scaling to reduce machine error.


co20z=cov(mzx1'); %NxN
k40z=sum(mzx1.*mzx1.*mzx1.*mzx1,2)/t-3*diag(co20z).*diag(co20z); %Nx1
zB=diag(abs(k40z).^(-1/4)); %NxN diagonal
mzx=zB*mzx1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalization by 4-th order cumulant of observed data

co20=cov(mzx');
k40=sum(mzx.*mzx.*mzx.*mzx,2)/t-3*diag(co20).*diag(co20);


%%%%%%%%%%%%%%%%%%%
% constants

fph=2; % kind of fixed-point homotopy deformation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializations

[VV,DD]=eig(cov(mzx'));
%[VV,DD]=eig(cov(mzx'),'nobalance');
%[VV, DD] = PCACOV(mzx');
C=VV'; %NxN autovalor


Delta=zeros(n);
c31=zeros(n);
c22=zeros(n);

repeat=1;
tau=1;
Cau2=100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main routine

while repeat & (tau < maxtau),

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % updating data

  cx=C*mzx; %NxT

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate fourth order cumulants
  
  co2=cx*cx.'/t; %NxN elementos fora da diag ~~ 0
  c31=(cx.*cx.*cx)*cx.'/t-3*diag(diag(co2))*co2; %NxN
  c22=(cx.*cx)*(cx.*cx).'/t-diag(co2)*diag(co2).'-2*co2.*co2; %NxN
  
  maxk=max(abs(diag(c31))); % um numero descendente
  k4b=diag(c31);
  k4=diag(c31)+sign(diag(c31))*maxk*thrsk;

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % determine each step
  %

  for ind=1:n
    fprintf('\n');
    ki=k4(ind);
    kbi=k4b(ind);
    for j=1:ind-1

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% fourth order terms

      c31ij=c31(ind,j);
      c31ji=c31(j,ind);
      c22ij=c22(ind,j);
      kj=k4(j);
      kbj=k4b(j);

      
      Sc4=1./([sqrt(abs(ki*ki*ki*kj)); sqrt(abs(ki*kj*kj*kj)); (abs(ki*kj))]); %3x1
      Sc4b=1./([sqrt(abs(kbi*kbi*kbi*kbj)); sqrt(abs(kbi*kbj*kbj*kbj));(abs(kbi*kbj))]);

      
      sFF4=Sc4.*[c31ij;c31ji;c22ij]; %3x1

      V4=[fph*c22ij, ki; ...
	  kj, fph*c22ij;  ...
	  2*c31ji, 2*c31ij]; %3x2
      
      Hosei=sFF4(1)*[6*c31ji, 5*c31ij; 5*c31ij,0]...
	    +sFF4(2)*[0, 5*c31ji; 5*c31ji,6*c31ij]...
	    +sFF4(3)*[2*kj,6*c22ij;6*c22ij,2*ki]; %2x2

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % determine whether we will use planar approximation

      bnb=V4.'*diag(Sc4b)*V4; %2x2
      bnb2=bnb+Hosei;
            
      if (trace(bnb2)>0) & (det(bnb2)>0), 
	bnb=bnb2;	  
	fprintf('s');
      else
	fprintf('z');
      end
    
      %%%%%%%%%%%%%%%%%%%%%%%%%
      % keep condition number < condM2
      
      trb=trace(bnb);
      dtb=det(bnb);
      meig=(trb-sqrt(trb^2-4*dtb))/2;
      Meig=(trb+sqrt(trb^2-4*dtb))/2;
      if meig<0 | abs(Meig/(meig+eps))>condM2,
	fprintf('X')
	bnb=bnb+(Meig-meig)/(condM2-1)*eye(2);
      end
      
      %%%%%%%%%%%%%
      % determine Delta_{ij} and Delta_{ji}
      
      ds=-inv(bnb)*(V4.'*sFF4); %2x1

            
      chone=abs(c31(ind,ind)/c31(j,j)); %um numero
      ds=ds/ceil(abs(2*log10(chone)));      

      Delta(ind,j)=ds(1);
      Delta(j,ind)=ds(2);
      
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % determine whether we can exit this function

  Dmet=diag(sqrt(abs(diag(c31))));
  Cau2=sqrt(trace(Delta*Dmet*Delta')*n/trace(Dmet))/(n*n-n);
%  Cau2=norm(Delta,'fro'); %he a sqrt(sum(diag(Delta'*Delta)))
  if  Cau2< epsilon, repeat=0;   end

  %%%%%%%%%%%%%%%%%%%%%%%
  % print to screen 
  fprintf(' ( step %d  (%e) ) ',tau, Cau2);
  %figure(9); hl=plot(tau,Cau2,'b*'); axis([1 100 0 1]);  AXIS xy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % updating 
  % I believe multiplicative updating is better  than additive one
  % at least for problems smaller than 200 channels 
  %Delta eh uma matriz NxN com delta(i,i)=0
  %figure(10); hl=mesh(expm(Delta)); cara1(hl,0); drawnow
  C=expm(Delta)*C;
  C=diag(1./sqrt(diag(C*C')))*C;
  tau=tau+1;

end

B=C*zB;

%%%%%%%%%%%%%%%% end of file %%%%%%%%%%%%%%%%%%%%%
