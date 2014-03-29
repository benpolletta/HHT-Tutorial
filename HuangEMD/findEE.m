%function [extre_IMF,zeocr_IMF,EE_IMF,ru,ru1]=findEE(EEMDIMF,chkplotEE)
%
% INPUT:
%       EEMDIMF   : Inputted data;2-d data-IMF matrix 
%       chkplotEE : when chkplotEE=1
%                   plot each IMF figure after max,min,zero-crossing are counted
% OUTPUT:
%       extre_IMF: A matrix,the number of local extrema (max+min) points for each IMF
%       zeocr_IMF: A matrix,the number of zero-crossing  points for each IMF
%       EE_IMF   : A matrix,E.E (excessive extrema) value for each IMF
%       ru       : A matrix,orthogonal index for 2 consecutive IMF
%       ru1      : A value,orthogonal index for all IMF
%
% NOTE: 
%       this function is made to find out the "goodness" of IMF sets
%       especially,for EEMD results.
%       to find out O.I(orthogonal index),E.E(Excesssive Extrema) 
%       ㄧ计琌т舱IMFぇ贺玒计,疭癸eemdぇ挡狦Τノ
%       тOrthogonal Index ---ru(ㄢㄢタユ玒计),ru1(砰タユ玒计)
%       тExcessive Extrema-E.E
%
% References:   
%  N. E Huang (2008),NCU hht class lecture 
%  2c b Ensemble EMD II.ppt 
%
% code writer: Sheng-Chung Su  
% footnote:S.C.Su 2009/04/15
%
% There are two loops coupled together in this code .
%  1.read data, find out 2 orthogonal index
%  2.find E.E
%    count local max,local min and zero-crossing numbers 
%  3.Check each IMF at a time ------------------------------------------------------loop B start
%    4.Check each data point for its value-------------------------------------loop A start
%      5.algorithm for find local-max,local-min,zero-crossing values
%        5.1-local-max
%        5.2-local-min
%        5.3-find zero-crossing
%    4.Check each data point for its value-------------------------------------loop A end  
%    6. calculate E.E  
%    7. warning for some(3) situations
%  3.Check each IMF at a time ------------------------------------------------------loop B end
%     8.Screen echo for user  
%     9.plot figure and save
%     
% Association: no
% this function ususally used for checking ensemble EMD result  
% because the sum of IMF is not a IMF,but for optimum EMD decomposition
% we check orthogonal index and E.E value 
%
% Concerned function: no
%                     

function [extre_IMF,zeocr_IMF,EE_IMF,ru,ru1]=findEE(EEMDIMF,chkplotEE)

%1.read data, find out 2 orthogonal index
%ノIMFuтbu---IMF羆计,тcu---戈羆翴计
%find out the total number of IMF & total data points
au=size(EEMDIMF);
bu=au(2);%bu should be the number of IMF
cu=au(1);%cu should be the number of total data  points
%for check if the matrix dimension need transpose
  if (cu < bu)
    EEMDIMF=EEMDIMF';
    au=size(EEMDIMF);
    bu=au(2);
    cu=au(1); 
  end

  %Exam for orthogonality for each IMF
   ru=ratioa(EEMDIMF);% neighbor IMF orthogonality 
   ru1=ratio1(EEMDIMF);% overall orthogonality


%2.find E.E
  %give initial  
  for j2=1:bu-1   %loop B start here
%3.Check each IMF at a time
   % loop B for check every IMF for EE ,one by one    
   
    extre_max=0;
    extre_min=0; 
    zero_cros=0;
    
   for j1=2:cu-1     %loop A start here
%4.Check each data point for its value
   % loop A for check every data point for extreme and zero-crossing

%5.algorithm for find local-max,local-min,zero-crossing values
     
  %talk about the values at positions  of  j1-1,j1,j1+1 
  %check which one is bigger
  %  IMF(j1-1,j2)
  %  IMF(j1,j2)
  %  IMF(j1+1,j2)
  % Find extrema
  %5.1-local-max
  % the equal sign is important and should exit only in one side --this is for flat point  
  %when (j1>=j1-1 j1>j1+1)----j1 is the local max
     if (EEMDIMF(j1-1,j2)<=EEMDIMF(j1,j2) && EEMDIMF(j1+1,j2)<EEMDIMF(j1,j2))
         extre_max=extre_max+1;
     end
  %5.2-local-min
  %when (j1<=j1-1 j1<j1+1)----j1 is the local min
     if (EEMDIMF(j1-1,j2)>=EEMDIMF(j1,j2) && EEMDIMF(j1+1,j2)>EEMDIMF(j1,j2))
       extre_min=extre_min+1;
     end
  %5.3-find zero-crossing 
  %when (j1-1*j1<0)----j1, means there is an 0
     if (EEMDIMF(j1-1,j2)*EEMDIMF(j1,j2) <= 0 )
        %add chk for bug 
        zero_cros=zero_cros+1;
     end

   end     %loop A stop here
       
%6. calculate E.E     
      %Count points for each EEMD-IMF compoments 
       extre_IMF(j2)=extre_max+extre_min;
       zeocr_IMF(j2)=zero_cros;
      %calculate EE for each EEMD-IMF compoments
       %add chk for bug line 119
       EE_IMF(j2)=1-(zeocr_IMF(j2)/extre_IMF(j2));
       
%7. warning for 3 situations
      %After EE is calculated, three conditions should be warned 
      %add a condition for extrema and zero-crossing difference=1 or -1
      
       differValue=zeocr_IMF(j2)-extre_IMF(j2);
       if abs(differValue)==1
           disp('-----When EE=0 ,this EEMD-compoment is an IMF--------,for the order of');
           disp(j2);
          EE_IMF(j2)=0;   
       end
       %add a condition for extrema=0, devided by zero
       if (extre_max==0 )
           disp('-----When EE=-999 ,this EEMD-compoment have 0 extrema point--------,for the order of' );
           disp(j2);
           EE_IMF(j2)=-999;   
       end 
       %add a condition for zero-crossing=0, not a IMF
       if (zero_cros==0 )
           disp('-----When EE=-555 ,this EEMD-compoment have 0 zero-crossing point,for the order of' );
           disp(j2);
           EE_IMF(j2)=-555;   
       end 
  end    %loop B stop here

%8. Screen echo for user  
  %do some demo on the command screen
  disp('----Calculation Report----')
  disp('Exreme numbers of IMFs are listed as follows:')
  disp(extre_IMF)
  disp('Zero-crossing numbers of IMFs are listed as follows:')
  disp(zeocr_IMF)
  disp('E.E numbers of IMFs are listed as follows:')
  disp(EE_IMF)


%9.plot figure and save
if chkplotEE==1
    destpath=input('Please give output folder paht \n','s');
    %B1=EEMDIMF';
%     %count the number of IMF
%      numberIMFA=size(EEMDIMF);
%      numberIMF=min(numberIMFA(1),numberIMFA(2));
%      numberPOI=max(numberIMFA(1),numberIMFA(2));
     %plot all IMF individually
     for g=1:bu-1
         h2=figure(2);
         g1=num2str(g);
         sE1=num2str(extre_IMF(g));
         sZ1=num2str(zeocr_IMF(g));
         s=1:cu;
         plot(s,EEMDIMF(:,g));axis tight;
         title(['EEMD-',g1,'IMF, Extrem=',sE1,' ,Zero-crossing=',sZ1]);
         filename=[destpath,'\EEMD-',g1,'.jpg'];
         saveas(h2,filename);
         close (h2);
     end 
end   

          