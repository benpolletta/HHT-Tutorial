Y='LOD78';

function allmode=eemd(Y,Nstd,NE)
xsize=length(Y);
dd=1:1:xsize;
Ystd=std(Y);
Y=Y/Ystd;

TNM=fix(log2(xsize))-1;
TNM2=TNM+2;
for kk=1:1:TNM2, 
    for ii=1:1:xsize,
        allmode(ii,kk)=0.0;
    end
end

for iii=1:1:NE,
    for i=1:xsize,
        temp=randn(1,1)*Nstd;
        X1(i)=Y(i)+temp;
    end

    for jj=1:1:xsize,
        mode(jj,1) = Y(jj);
    end
    
    xorigin = X1;
    xend = xorigin;
    
    nmode = 1;
    while nmode <= TNM,
        xstart = xend;
        iter = 1;
   
        while iter<=10,
            [spmax, spmin, flag]=extrema(xstart);
            upper= spline(spmax(:,1),spmax(:,2),dd);
            lower= spline(spmin(:,1),spmin(:,2),dd);
            mean_ul = (upper + lower)/2;
            xstart = xstart - mean_ul;
            iter = iter +1;
        end
        xend = xend - xstart;
   
   	    nmode=nmode+1;
        
        for jj=1:1:xsize,
            mode(jj,nmode) = xstart(jj);
        end
    end
   
    for jj=1:1:xsize,
        mode(jj,nmode+1)=xend(jj);
    end
   
    allmode=allmode+mode;
    
end

allmode=allmode/NE;
allmode=allmode*Ystd;