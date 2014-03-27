%function [begx1,begy1,endx1,endy1,begx2,begy2,endx2,endy2] = endprocess1(IMF, chkplot,methodtype)
%
%
%input
%   IMF-1D time series-be an IMF or not be an IMF both ok!
%   methodtype='A'-for 5 point extend
%              'B'-for first and last extrema copy once
%              'C'-for neighbor peak slope extend out
%   chkplot-when the value=100,plot end-process results
%
%output-usage 1 ,use like copy endpoints.m
% additional extend points
%   begx1-beginning points for spline control-x values for max
%   begy1-beginning points for spline control-y values for max
%   endx1-endding points for spline control-x values for max
%   endy1-endding points for spline control-y values for max
%   begx2-beginning points for spline control-x values for min
%   begy2-beginning points for spline control-y values for min
%   endx2-endding points for spline control-x values for min
%   endy2-endding points for spline control-y values for min
%
%Note:
%     using cubic spline as the envelope is widely use in EMD
%     end-effect is an unsolved problem in the beginning and ending part
%     and the spline envelope is greatly influenced by end process
%     here collectin 3 different ways to process the end 
%
%Reference:  
%
%
%  code writer:S.C.Su
% footnote:S.C.Su 2009/05/14
%
%1.find out extrema for the 1-D IMF
%2.use 3 kind of end-process to add additional data-points
%3.form the envelope with cubic spline
%4.plot the result and you can check  step by step 
%
%
%
% Association:  those procedure of HHT need this code
%1.EMD 
%2.EEMD
%3.Normalize procedure of HHT
%
% Concerned function: emax.m emin.m
%                     above mentioned m file must be put together
%                     or be put in the matlab working path
%

function [begx1,begy1,endx1,endy1,begx2,begy2,endx2,endy2] = endprocess1(IMF, chkplot,methodtype)

%

%There are three different ways of endprocess in this code  
%1.find out extrema for the 1-D IMF
   
   ggg=1:length(IMF);
   [Ymax, Xmax]=emax(IMF);%emax find local max .Return- Ymax-max_points and  Xmax-their_coordinates
    MAXnumber=size(Xmax,1); 
   [Ymin, Xmin]=emin(IMF);%emax find local min .Return- Ymin-min_points and Xmin-their_coordinates
    MINnumber=size(Xmin,1);   

%2.use 3 kind of end-process to add additional data-points  
   
   
   switch methodtype
   case{'A'}      
   %case A-from copyendpoints.m
      %local max 
       %method A-copy endpoints-beginning points-for local max
         % for x-dir
          if(MAXnumber > 1)
             delta_x1 = Xmax(2) - Xmax(1);
          else
             delta_x1 = Xmax(1);
          end
          %----- Make sure distances are greater than the distance to data end
          if((Xmax(1) - delta_x1) > 0)
             delta_x1 = Xmax(1);
          end
         % for y-dir
           y_can1 = Ymax(1);
          %----- Make sure extrema is greater than original data
            if((y_can1 < IMF(1)))
                 y_can1 = IMF(1);
            end
        %extend value-extend local max values-begx,begy-extend first time
            begx1 = Xmax(1)-delta_x1;
            begy1 = y_can1;
            
       %method A-copy endpoints-endding points-for local max      
         % for x-dir
          if(MAXnumber > 1)
             delta_x2 = Xmax(end) - Xmax(end-1);
          else
             delta_x2 = length(IMF)-Xmax(end);
          end
          %----- Make sure distances are greater than the distance to data end
          if((delta_x2 + Xmax(end)) < length(IMF))
             delta_x2 = length(IMF) - Xmax(end);
          end
         % for y-dir
           y_can2 = Ymax(end);
           %----- Make sure extrema is greater than original data
             if(y_can2 < IMF(end))
                 y_can2 = IMF(end);
             end
        %extend value-extend local max values-begx,begy-extend first time
            endx1 = Xmax(end)+delta_x2;
            endy1 = y_can2;
      
      %extend out for 4 more times
           for i=1:4
               begx1 = [begx1(1)-delta_x1; begx1 ];
               begy1 = [y_can1; begy1];
               endx1 = [endx1; endx1(end)+delta_x2];
               endy1 = [endy1; y_can2];
           end
      
      
      %local min
       %method A-copy endpoints-beginning points-for local min
         % for x-dir
          if(MINnumber > 1)
             delta_x3 = Xmin(2) - Xmin(1);
          else
             delta_x3 = Xmin(1);
          end
          %----- Make sure distances are greater than the distance to data end
          if((Xmin(1) - delta_x3) > 0)
             delta_x3 = Xmin(1);
          end
         % for y-dir
           y_can3 = Ymin(1);
          %----- Make sure extrema is greater than original data
            if((y_can3 > IMF(1)))
                 y_can3 = IMF(1);
            end
        %extend value-extend local max values-begx,begy-extend first time
            begx2 = Xmin(1)-delta_x3;
            begy2 = y_can3;
            
       %method A-copy endpoints-endding points-for local min      
         % for x-dir
          if(MINnumber > 1)
             delta_x4 = Xmin(end) - Xmin(end-1);
          else
             delta_x4 = length(IMF)-Xmin(end);
          end
          %----- Make sure distances are greater than the distance to data end
          if((delta_x4 + Xmin(end)) < length(IMF))
             delta_x4 = length(IMF) - Xmin(end);
          end
         % for y-dir
           y_can4 = Ymin(end);
           %----- Make sure extrema is greater than original data
             if(y_can4 > IMF(end))
                 y_can4 = IMF(end);
             end
        %extend value-extend local max values-begx,begy-extend first time
            endx2 = Xmin(end)+delta_x4;
            endy2 = y_can4;
      
      %extend out for 4 more times 
           for i=1:4
               begx2 = [begx2(1)-delta_x3; begx2 ];
               begy2 = [y_can3; begy2];
               endx2 = [endx2; endx2(end)+delta_x4];
               endy2 = [endy2; y_can4];
           end
        
                     
      
    case{'B'}
 
        
%----- Treat the head
    n_mx=Xmax(1);
    n_mn=Xmin(1);
    e_mx=Xmax(end);
    e_mn=Xmin(end);
    
    if (n_mx==-1) | (n_mn==-1) % no max or min
        disp('At least one component does not have Max or Min!');
        %break;
    elseif n_mn<n_mx,
        dlt=n_mx-n_mn;
    else
        dlt=n_mn-n_mx;
    end
%----- extend end point out for the head   
    begx1=Xmax(1)-2*dlt;
    begx1=[begx1(1)-2*dlt; begx1];
    begy1=Ymax(1);
    begy1=[Ymax(1);begy1];
    begx2=Xmin(1)-2*dlt;
    begx2=[begx2(1)-2*dlt; begx2];
    begy2=Ymin(1);
    begy2=[Ymin(1);begy2];

%----- Treat the tail
    if e_mn<e_mx,
        dlt=e_mx-e_mn;
    else
        dlt=e_mn-e_mx;
    end

%----- extend end point out for the end  
    endx1=Xmax(end)+2*dlt;
    endx1=[endx1; endx1(end)+2*dlt];
    endy1=Ymax(end);
    endy1=[endy1; Ymax(end)];
    endx2=Xmin(end)+2*dlt;
    endx2=[endx2; endx2(end)+2*dlt];
    endy2=Ymin(end);
    endy2=[endy2; Ymin(end)];
        
        
        
        
        
    case{'C'}

        
%End point process-please see reference about spline end effect
%extend the slpoe of neighbor 2 max value ---as extend value
%original value of end point -----as original value
%compare extend and original value 

if MAXnumber>=4
    %connect two local max value with a line
    slope1=(Ymax(1)-Ymax(2))/(Xmax(1)-Xmax(2));
    %extend out to the beginning end
    tmp1=slope1*(1-Xmax(1))+Ymax(1);
    %compare the value with original data end point
    % the bigger value is favorite
    if tmp1>IMF(1)
        begx1=1;
        begy1=tmp1;
    else
        begx1=1;
        begy1=IMF(1);
    end

    %connect two local max value with a line 
    slope2=(Ymax(end)-Ymax(end-1))/(Xmax(end)-Xmax(end-1));
    %extend out to the endding end
    tmp2=slope2*(length(IMF)-Xmax(end))+Ymax(end);
    %compare the value with original data end point
    % the bigger value is favorite
    if tmp2>IMF(end)
        endx1=length(IMF);
        endy1=tmp2;
    else
        endx1=length(IMF);
        endy1=IMF(end);
    end
else
    disp('not an IMF-for no maxinum')
end


if MINnumber>=4
    %connect two local min value with a line
    slope3=(Ymin(1)-Ymin(2))/(Ymin(1)-Ymin(2));
    %extend out to the beginning end
    tmp3=slope3*(1-Xmin(1))+Ymin(1);
    %compare the value with original data end point
    % the smaller value is favorite
    if tmp1<IMF(1)
        begx2=1;
        begy2=tmp3;
    else    
        begx2=1;
        begy2=IMF(1);
    end

    %connect two local min value with a line
    slope4=(Ymin(end)-Ymin(end-1))/(Xmin(end)-Xmin(end-1));
    %extend out to the endding end
    tmp4=slope4*(length(IMF)-Xmin(end))+Ymin(end);
    %compare the value with original data end point
    % the smaller value is favorite
    if tmp4<IMF(end)
        endx2=length(IMF);
        endy2=tmp4;
    else
        endx2=length(IMF);
        endy2=IMF(end);
    end
else
     disp('not an IMF-for no minumn')
end
        
        
    
        
        
        
   end

     
%3.form the envelope with cubic spline         
%   those end control points are used
          %envelpe control points
           Xptu=[begx1;Xmax;endx1];
           Yptu=[begy1;Ymax;endy1];
           Xptd=[begx2;Xmin;endx2];
           Yptd=[begy2;Ymin;endy2];
           gg1u=Xptu(1):1:Xptu(end);
           gg1d=Xptd(1):1:Xptd(end);
    if chkplot==100       
           %cubic spline envelpe
           uper1=spline(Xptu,Yptu,gg1u);
           down1=spline(Xptd,Yptd,gg1d);
           uper2=spline(Xptu,Yptu,ggg);
           down2=spline(Xptd,Yptd,ggg);

%4.plot the result and you can check  step by step          
           
           %plot the result
           figure(2);
           plot(begx1,begy1,'bo');hold on
           plot(endx1,endy1,'bo');
           plot(ggg,IMF,'r');
           plot(begx2,begy2,'bo');
           plot(endx2,endy2,'bo');
           plot(gg1u,uper1,'g',gg1d,down1,'g',ggg,uper2,'m',ggg,down2,'m');
           
     end      

