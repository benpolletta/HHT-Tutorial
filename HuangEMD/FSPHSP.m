%function [HESPw,HESPv,FESPw,FESPv]=FSPHSP(file,samplerate,freqsol,timesol,hsp_fre1,noiselevel,Nensemble)
%
%
%INPUT
%     file----------full path file name,for instance,'c:\example\lod78.csv'
%     samplerate----sample rate of data 
%     freqsol-------frequency axis bin number when plotting an Hilbert Energy Spectrum
%     timesol-------time  axis bin number when plotting an Hilbert Energy Spectrum 
%     hsp_fre1------maximun frequency when plotting an Hilbert Energy Spectrum
%     noiselevel----EEMD parameters--add noise level 
%     Nensemble-----EEMD parameters--number of ensemble
%
%OUTPUT
%     HESPw---------frequency values of Marginal Spectrum (from Hilbert Energy Spectrum)
%     HESPv---------Energy values of Marginal Spectrum (from Hilbert Energy Spectrum)
%     FESPw---------frequency values of Fourier Energy Spectrum
%     FESPw---------Energy values of Fourier Energy Spectrum
%
%NOTE:
%     give an 1D data set
%     1.this program performs EEMD ,do Hilbert-Spectrum, do Marginal-Spectrum
%     2.this  program performs Pwelch --Fourier Energy Spectrum
%     3.Calculate  comparison coefficients ,and multiply to those spectrum
%     4.plot 1. and 2. together 
%     
%Note:
%     the translation coefficients are formulaed by N.E.H(Norden E. Huang) in HHT class 2008,NCU
%     the code is written by Student --S.C.Su(Sheng-Chung Su)
%     footnote:S.C.Su

function [HESPw,HESPv,FESPw,FESPv]=FSPHSP(file,samplerate,freqsol,timesol,hsp_fre1,noiselevel,Nensemble)

%0.load data
  datafile=load(file);
                                
%1.this program performs EEMD ,do Hilbert-Spectrum, do Marginal-Spectrum                                
    %do eemd  and show result
     EMDIMF=eemd(datafile,noiselevel,Nensemble);
     figure(1)
     anoiselevel=num2str(noiselevel);
     aNensemble=num2str(Nensemble);
     strips(EMDIMF);title(['EEMD result     ,noise level=',anoiselevel,' Number of ensemble=',aNensemble]);
     %plot Hilbert-Energy-SPectrum
    
    % emd result counting
     au=size(EMDIMF);nIMF=au(2)-2;nPT=au(1)-1;
     totalt=(nPT+1)/samplerate;
     Xlow=1/totalt;
     Xhig=samplerate/2;
     
    %plot NNSPE
     [nt,ta,fa]=nnspe(EMDIMF(1:nPT,2:nIMF-1), 0, totalt, freqsol, timesol, 0.00001, hsp_fre1,0,totalt,'hilbtm','spline',5);
    %calculate the marginal spectrum
     ms=sum(nt,2);
    %Smooth the pic ,use an Gaussian Filter "q",let nt been filtered twice,HESP spectrum will be good for reading
     q=fspecial('gaussian', 7, 0.6);
     nsu=filter2(q, nt);nsu=filter2(q, nsu);
    %Plot the HESP and MS  
     figure(2);
     imagesc(ta,fa,nsu.^.5);axis xy;set(gca,'FontSize',6);title('Time-frequency analysis-Hilbert Spectrum','FontSize',9);xlabel('Time(sec)','FontSize',7,'VerticalAlignment','middle');ylabel('Freq(hz)','FontSize',7);ylim([Xlow hsp_fre1]);
     
% 2.this  program performs Pwelch --Fourier Energy Spectrum  
     
     %plot Fourier-Energy-SPectrum  
       stroridata='Pwelch';
      [Pxx,w] =pwelch(datafile,[],[],[],samplerate,'onesided');
      
%3.Calculate  comparison coefficients ,and multiply to those spectrum     
     
     %Calculate the specrum to bin=400,max freq=100hz
       fco=2*samplerate*(nPT+1)/2/freqsol;
       hco=samplerate/2/hsp_fre1;
%4.plot two different kind of spectrum together        
       figure(5)
       loglog(w(1:end),fco*Pxx(1:end),'g',fa(1:end),hco*ms(1:end),'r');set(gca,'FontSize',6);title('Spectrum Comparison-Energy -and -frequency','FontSize',9);xlabel('Hz(1/sec)','FontSize',7,'VerticalAlignment','middle');ylabel('Energy','FontSize',7);xlim([Xlow hsp_fre1]);legend('FESP','MESP');
     
     %return those values  
       HESPw=fa;
       HESPv=hco*ms;
       FESPw=w;
       FESPv=fco*Pxx;     