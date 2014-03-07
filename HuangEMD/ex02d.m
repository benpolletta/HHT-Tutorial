%ex02
%topic-HOW  TO 
%"RUN EEMD or EMD  and use the result to 
%PLOT Hilbert-Energy or Hilbert-Energy Spectrum"
clear;close all
clc
%load data file  
  datafile=csvread('lod78.csv');
  
%==============Run EEMD=================================  
%give noise level
 noiselevel=0;
%give Number of ensemble 
 Nensemble=1;
%do eemd  and show result
 EEMDIMF=eemd(datafile,noiselevel,Nensemble);
%plot the EEMD result 
     figure(1) 
     anoiselevel=num2str(noiselevel);
     aNensemble=num2str(Nensemble);
     strips(EEMDIMF);title(['EEMD result     ,noise level=',anoiselevel,' Number of ensemble=',aNensemble]);
     clear noiselevel Nensemble
 disp('===== EEMD calculation complete! =====')    
 
% %==============Run EMD================================= 
% %give noise level
%  noiselevel=0;
% %give Number of ensemble 
%  Nensemble=1;
% %do eemd  and show result
%  EMDIMF=eemd(datafile,noiselevel,Nensemble);
% %plot the EMD result 
%      figure(4) 
%      anoiselevel=num2str(noiselevel);
%      aNensemble=num2str(Nensemble);
%      strips(EMDIMF);title(['EMD result     ,noise level=',anoiselevel,' Number of ensemble=',aNensemble]);
%  disp('===== EMD calculation complete! =====')    
     

%==============Run Hilbert-Energy Spectrum for EEMD =======================   
%give samplerate 
  samplerate=365;
%give frequency-axis resolution for hilbert-spectrum
  freqsol=400;
%give time-axis resolution for hilbert-spectrum
  timesol=800;
%give frequency-axis maximun value
  hsp_fre1=40;

% dealing with EEMD result 
     au=size(EEMDIMF);nIMF=au(2)-2;nPT=au(1)-1;
     totalt=(nPT+1)/samplerate;
     Xlow=1/totalt;
     Xhig=samplerate/2;
     
%Calculate Hilbert-Energy-Spectrum(NNSPE.m)
     [nte,tae,fae]=nnspe(EEMDIMF(1:nPT,2:nIMF+1), 0, totalt, freqsol, timesol, 0, hsp_fre1,0,totalt);
     %nnspe(EEMDIMF(1:nPT,2:nIMF-1), 0, totalt, freqsol, timesol, 0.00001,%hsp_fre1,0,totalt);
    
 %Smooth the pic ,use an Gaussian Filter "q",let nt been filtered twice,HESP spectrum will be good for reading
     q=fspecial('gaussian', 7, 0.6);
     nse=filter2(q, nte);nsue=filter2(q, nse);
     
 %Plot the HESP and MS  
     figure(2);
     imagesc(tae,fae,nsue.^.5);axis xy;set(gca,'FontSize',6);title('Hilbert Energy Spectrum- for EEMD','FontSize',9);xlabel('Time(sec)','FontSize',7,'VerticalAlignment','middle');ylabel('Freq(hz)','FontSize',7);ylim([Xlow hsp_fre1]);
 
%  %==============Run Hilbert-Amplitude Spectrum for EEMD =======================   
%  %Calculate Hilbert-Amplitude-Spectrum(NNSPE.m)
%      [nta,taa,faa]=nnspa(EEMDIMF(1:nPT,2:nIMF+1), 0, totalt, freqsol, timesol, 0, hsp_fre1,0,totalt);
%      
%  %Smooth the pic ,use an Gaussian Filter "q",let nt been filtered twice,HESP spectrum will be good for reading
%      nsa=filter2(q, nta);nsua=filter2(q, nsa);
%      
%  %Plot the HESP and MS  
%      figure(3);
%      imagesc(taa,faa,nsua.^.5);axis xy;set(gca,'FontSize',6);title('Hilbert Amplitude Spectrum- for EEMD','FontSize',9);xlabel('Time(sec)','FontSize',7,'VerticalAlignment','middle');ylabel('Freq(hz)','FontSize',7);ylim([Xlow hsp_fre1]);
%    
%  
%      
%      
%      
%     clear  samplerate freqsol timesol hsp_fre1 au totalt Xlow Xhig nta taa faa nte tae fae
%  %==============Run Hilbert-Energy Spectrum for EMD =======================   
% %give samplerate 
%   samplerate=365;
% %give frequency-axis resolution for hilbert-spectrum
%   freqsol=400;
% %give time-axis resolution for hilbert-spectrum
%   timesol=800;
% %give frequency-axis maximun value
%   hsp_fre1=40;
% 
% % dealing with EEMD result 
%      au=size(EMDIMF);nIMF=au(2)-2;nPT=au(1)-1;
%      totalt=(nPT+1)/samplerate;
%      Xlow=1/totalt;
%      Xhig=samplerate/2;   
%  %Calculate Hilbert-Amplitude-Spectrum(NNSPE.m)
%      [nte,tae,fae]=nnspe(EMDIMF(1:nPT,2:nIMF+1), 0, totalt, freqsol, timesol, 0, hsp_fre1,0,totalt);
%      
%  %Smooth the pic ,use an Gaussian Filter "q",let nt been filtered twice,HESP spectrum will be good for reading
%      nse=filter2(q, nte);nsue=filter2(q, nse);
%      
%  %Plot the HESP and MS  
%      figure(5);
%      imagesc(tae,fae,nsue.^.5);axis xy;set(gca,'FontSize',6);title('Hilbert Energy Spectrum- for EMD','FontSize',9);xlabel('Time(sec)','FontSize',7,'VerticalAlignment','middle');ylabel('Freq(hz)','FontSize',7);ylim([Xlow hsp_fre1]);
%      
%  %==============Run Hilbert-Amplitude Spectrum for EEMD =======================   
%  %Calculate Hilbert-Amplitude-Spectrum(NNSPE.m)
%      [nta,taa,faa]=nnspa(EMDIMF(1:nPT,2:nIMF+1), 0, totalt, freqsol, timesol, 0, hsp_fre1,0,totalt);
%      
%  %Smooth the pic ,use an Gaussian Filter "q",let nt been filtered twice,HESP spectrum will be good for reading
%      nsa=filter2(q, nta);nsua=filter2(q, nsa);
%      
%  %Plot the HESP and MS  
%      figure(6);
%      imagesc(taa,faa,nsua.^.5);axis xy;set(gca,'FontSize',6);title('Hilbert Amplitude Spectrum- for EMD','FontSize',9);xlabel('Time(sec)','FontSize',7,'VerticalAlignment','middle');ylabel('Freq(hz)','FontSize',7);ylim([Xlow hsp_fre1]);
%          
     
     