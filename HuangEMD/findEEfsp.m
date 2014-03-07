%function [extre_IMF,zeocr_IMF,EE_IMF,ru,ru1,F_FAP_IMF,V_FAP_IMF,F_FEP_IMF,V_FEP_IMF]=findEE(EEMDIMF,chkplotEE,samplerate,chkplotfsp,plotFSPfreq)
%
% INPUT:
%       EEMDIMF   : Inputted data;2-d data-IMF matrix 
%       chkplotEE : when chkplotEE=1
%                   plot each IMF figure after max,min,zero-crossing are counted
%       samplerate: original data sample rate,for Fourier-energy and Fourier-amplitude spectrum 
%       chkplotfsp: when chkplotfsp=1
%                   plot Fourier-energy spectrum for each IMF 
%      plotFSPfreq: max frequency value when plot Fourier spectrum 
%                   this value must lower than Nyquist frequency 
%
% OUTPUT:
%       extre_IMF: A matrix,the number of local extrema (max+min) points for each IMF
%       zeocr_IMF: A matrix,the number of zero-crossing  points for each IMF
%       EE_IMF   : A matrix,E.E (excessive extrema) value for each IMF
%       ru       : A matrix,orthogonal index for 2 consecutive IMF
%       ru1      : A value,orthogonal index for all IMF
%       F_FAP_IMF: A 1-D matrix,the frequency axis value for Fourier-Amplitude spectrum
%       V_FAP_IMF: A 2-D matrix,the Fourier-Amplitude value for Fourier-Amplitude spectrum
%                  Put value columnwise for each IMF 
%       F_FEP_IMF: A 1-D matrix,the frequency axis value for Fourier-Energy spectrum
%       V_FEP_IMF: A 2-D matrix,the Fourier-Energy value for Fourier-Energy spectrum
%                  Put value columnwise for each IMF 
%
%
% NOTE: 
%       this function is made to find out the "goodness" of IMF sets
%       especially,for EEMD results.
%       to find out O.I(orthogonal index),E.E(Excesssive Extrema) 
%       此函數是找出一組IMF之各種係數,特別對eemd後之結果有效用
%       去找出Orthogonal Index ---ru(兩兩正交係數),ru1(全體正交係數)
%       去找出Excessive Extrema-E.E
%       去畫出各個IMF之Fourier-Amplitude Spectrum 與 Fourier-Energy spectrum
%       並畫出所有IMF的能譜於一張圖上
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
%  2.find E.E and Fourier Spectrum
%    count local max,local min and zero-crossing numbers 
%    calculate Fourier Spectrum
%  3.Check each IMF at a time ------------------------------------------------------loop B start
%     4.calculate Fourier Spectrum    
%     5.Check each data point for its value-------------------------------------loop A start
%      6.algorithm for find local-max,local-min,zero-crossing values
%        6.1-local-max
%        6.2-local-min
%        6.3-find zero-crossing
%     5.Check each data point for its value-------------------------------------loop A end
%      7. calculate E.E  
%      8. warning for some(3) situations
%  3.Check each IMF at a time ------------------------------------------------------loop B end
%     9.Screen echo for user  
%     10.plot IMF figure and save 
%     11.plot Fourier Spectrum figure and save 
%     
%
% Association: no
% 1. this function ususally used for checking ensemble EMD result  
%    because the sum of IMF is not a IMF,but for optimum EMD decomposition
%    we check orthogonal index and E.E value 
% 2. this function help us know the energy distribution of ensemble EMD result
%    we can see nonlinear filter result on the spectrum
% 
% Concerned function: no
%            

function [extre_IMF,zeocr_IMF,EE_IMF,ru,ru1,F_FAP_IMF,V_FAP_IMF,F_FEP_IMF,V_FEP_IMF]=findEE(EEMDIMF,chkplotEE,samplerate,chkplotfsp,plotFSPfreq)

%1.read data, find out 2 orthogonal index
%利用IMFu去找出bu---IMF總數,找出cu---資料總點數
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


%2.find E.E and Fourier Spectrum
  %give initial  
       NmaxfN=2^nextpow2(cu)/2;
       F_FAP_IMF(1:NmaxfN)=0;
       V_FAP_IMF(1:NmaxfN,1:bu-1)=0;
       F_FEP_IMF(1:NmaxfN+1)=0;
       V_FEP_IMF(1:NmaxfN+1,1:bu-1)=0;
  
  
  
  for j2=1:bu-1   %loop B start here
%3.Check each IMF at a time
   % loop B for check every IMF for EE ,one by one    
   
    extre_max=0;
    extre_min=0; 
    zero_cros=0;

%4.calculate Fourier Spectrum    
    %find Eourier Energy spectrum for every IMF---start
      %Sample Rate
       npt=cu;
       data=EEMDIMF(:,j2);
       dt=1/samplerate;
       times=0:dt:dt*(npt-1);
       totaltu=npt*dt;
       %計算上下限 數學上FSP最大與最小頻率
       high_limit_fsp=0.5/dt;
       low_limit_fspu=1/totaltu;
       %進行兩個FFT計算
       %a----基線未修正 
       %FFT計算 
       nfft_u=2^nextpow2(npt);
       fftxu=0.5/dt*linspace(0,1,nfft_u/2);
       %基線修正
         %1.取全部長度之1/5來做base line correlation
                basemin=npt/5;  
         %U1avg1------全部長度之1/5-前1/5平均值
                U1avg1=mean(data(1:basemin));
         %arrayu
               arrayu=data-U1avg1;
         %FFT-Fourier Amplitude Spectrum 計算        
             ffty1n=2*fft(arrayu(1:npt),nfft_u)/npt;
             ffty2n=abs(ffty1n(1:nfft_u/2));
         %FFT-Fourier Energy Spectrum 計算 
             [Pxx,f]=pwelch(arrayu,[],[],nfft_u,samplerate);    
       
      %Store the FSP into matrix 
      if j2==1
       F_FAP_IMF(:)=fftxu;
       F_FEP_IMF(:)=f;
      end
       V_FAP_IMF(:,j2)=ffty2n;
       V_FEP_IMF(:,j2)=Pxx;
           
           
    %find Eourier Energy spectrum for every IMF---end
    %extrema bug should be fixed~~~~~~~~~~~~~~~~~~~
    %add chk for bug line 86

%5.Check each data point for its value
   for j1=2:cu-1     %loop A start here
   % loop A for check every data point for extreme and zero-crossing

%6.algorithm for find local-max,local-min,zero-crossing values     
  %talk about the values at positions  of  j1-1,j1,j1+1 
  %check which one is bigger
  %  IMF(j1-1,j2)
  %  IMF(j1,j2)
  %  IMF(j1+1,j2)
  % Find extrema
  %6.1-local-max
  %when (j1>=j1-1 且j1>j1+1)----j1 is the local max
     if (EEMDIMF(j1-1,j2)<=EEMDIMF(j1,j2) && EEMDIMF(j1+1,j2)<EEMDIMF(j1,j2))
        %add chk for bug 
         extre_max=extre_max+1;
     end
  %6.2-local-min
  %when (j1<=j1-1 且j1<j1+1)----j1 is the local min
     if (EEMDIMF(j1-1,j2)>=EEMDIMF(j1,j2) && EEMDIMF(j1+1,j2)>EEMDIMF(j1,j2))
         %add chk for bug line 102
         extre_min=extre_min+1;
     end
  %6.3-find zero-crossing 
  %when (j1-1*j1<0)----j1, means there is an 0
     if (EEMDIMF(j1-1,j2)*EEMDIMF(j1,j2) <= 0 )
        %add chk for bug 
        zero_cros=zero_cros+1;
     end

   end     %loop A stop here

%7. calculate E.E          
      %Count points for each EEMD-IMF compoments 
       extre_IMF(j2)=extre_max+extre_min;
       zeocr_IMF(j2)=zero_cros;
      %calculate EE for each EEMD-IMF compoments
       %add chk for bug line 119
       EE_IMF(j2)=1-(zeocr_IMF(j2)/extre_IMF(j2));

%8. warning for 3 situations       
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
  
%9. Screen echo for user
  %do some demo on the command screen
  disp('----Calculation Report----')
  disp('Exreme numbers of IMFs are listed as follows:')
  disp(extre_IMF)
  disp('Zero-crossing numbers of IMFs are listed as follows:')
  disp(zeocr_IMF)
  disp('E.E numbers of IMFs are listed as follows:')
  disp(EE_IMF)

%10.plot IMF figure and save 
if chkplotEE==1
    destpathEE=input('Please give output folder path for E.E figures \n','s');
     %plot all IMF individually
     for g=1:bu-1
         h2=figure(2);
         set(h2,'Visible','on');
         set(h2,'PaperUnits','inches');
         set(h2,'PaperPosition',[0 0 22 17]);
         set(h2,'PaperSize',[22 17]);
         set(h2,'PaperType','A2');
         g1=num2str(g);
         sE1=num2str(extre_IMF(g));
         sZ1=num2str(zeocr_IMF(g));
         plot(times,EEMDIMF(:,g));axis tight;
         title(['EEMD-',g1,' IMF, Extrem=',sE1,' ,Zero-crossing=',sZ1]);
         filename=[destpathEE,'\EEMD-',g1,'.jpg'];
         saveas(h2,filename);
         close (h2);
     end 
end    

%11.plot Fourier Spectrum figure and save     
if chkplotfsp==1
    destpathfsp=input('Please give output folder path for FSP figures\n','s');
     %plot all IMF individually
     for gg=1:bu
         
       if gg ~=bu 
         h3=figure(3);
         set(h3,'Visible','on');
         set(h3,'PaperUnits','inches');
         set(h3,'PaperPosition',[0 0 22 17]);
         set(h3,'PaperSize',[22 17]);
         set(h3,'PaperType','A2'); 
         gg1=num2str(gg);
         sE1=num2str(extre_IMF(gg));
         sZ1=num2str(zeocr_IMF(gg));
         subplot(2,2,1);plot(F_FAP_IMF,V_FAP_IMF(:,gg));axis tight;title(['EEMD-',gg1,' IMF, Fourier Amplitude Spectrum']);xlim([0 plotFSPfreq]);
         subplot(2,2,2);loglog(F_FAP_IMF,V_FAP_IMF(:,gg));axis tight;title(['EEMD-',gg1,' IMF, Fourier Amplitude Spectrum-LOG']);
         subplot(2,2,3);plot(F_FEP_IMF,V_FEP_IMF(:,gg));axis tight;title(['EEMD-',gg1,' IMF, Fourier   Energy  Spectrum']);xlim([0 plotFSPfreq]);
         subplot(2,2,4);loglog(F_FEP_IMF,V_FEP_IMF(:,gg));axis tight;title(['EEMD-',gg1,' IMF, Fourier   Energy  Spectrum-LOG']);
         filenameg=[destpathfsp,'\EEMDfsp-',gg1,'.jpg'];
         saveas(h3,filenameg);
         close (h3);
       else
         h4=figure(4); 
         set(h4,'Visible','on');
         set(h4,'PaperUnits','inches');
         set(h4,'PaperPosition',[0 0 22 17]);
         set(h4,'PaperSize',[22 17]);
         set(h4,'PaperType','A2');
         for st=1:bu-1  
         color1=rand;color2=rand;color3=rand;
         subplot(2,2,1);plot(F_FAP_IMF,V_FAP_IMF(:,st),'Color',[color1 color2 color3]);axis tight;title('EEMD-all IMF, Fourier Amplitude Spectrum');xlim([0 plotFSPfreq]);hold on;
         subplot(2,2,2);loglog(F_FAP_IMF,V_FAP_IMF(:,st),'Color',[color1 color2 color3]);axis tight;title('EEMD-all IMF, Fourier Amplitude Spectrum-LOG');hold on;
         subplot(2,2,3);plot(F_FEP_IMF,V_FEP_IMF(:,st),'Color',[color1 color2 color3]);axis tight;title('EEMD-all IMF, Fourier   Energy  Spectrum');xlim([0 plotFSPfreq]);hold on;
         subplot(2,2,4);loglog(F_FEP_IMF,V_FEP_IMF(:,st),'Color',[color1 color2 color3]);axis tight;title('EEMD-all IMF, Fourier   Energy  Spectrum-LOG');hold on;  
         end
         filenameg=[destpathfsp,'\EEMDfsp-all.jpg'];
         filenameg1=[destpathfsp,'\EEMDfsp-all.fig'];
         saveas(h4,filenameg);
         saveas(h4,filenameg1);
         close (h4);
         
          
       end
         
     end      
end   

          