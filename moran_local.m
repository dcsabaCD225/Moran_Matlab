%%------Local Moran's I calculations-------
% for theoretical details see Dávid et al.......

%you can change parameters till the double dashed line
tic;
%Import the tiff or jpg file and rename it for kep
[M,N]=size(kep);

%different orders can be tested, must be listed here
%ords=[1;2;3;4;5;6;7;8;9;10;16;32;64];
ords=[8];
ord=numel(ords);

%Results to csv?
csvw=1;
%if it will be written in csv, then it'll be its name
%if a file with the same name exist, it will be overwritten
ne=['my_morans_csv_'];

%values in convolution matrix 
%either following Gauss(1) distribution or beeing constant (0)
Ga=1;

%number of repeats for calculating psudosignificance
ism=99;

%do not change the parameters below
%- - - - - - - - - - - -
%- - - - - - - - - - - -

%screen size
scrS=get(0,'ScreenSize');
mer=M*N;
%conversion in one column
ko=reshape(kep,mer,1);
ko=double(ko);

%mean, SD 
km=mean(ko);
ksd=std(ko,0,1);

%normalize
Zi=(ko-km)./ksd;

%M*N format again
kep_Zi=reshape(Zi,M,N);

ido=zeros(ord,1);
Iertekek=zeros(ord,1);

%loop for the orders
for oo=1:ord
    o=ords(oo);
    %size of weight mtx
    Gm=1+2*o;
    %same weight
    if Ga==0
        W=ones(Gm,Gm);
        W(o+1,o+1)=0;
    end % end of Ga==0
    %Gaussion weight distribution
    if Ga==1
        W = zeros(Gm,Gm);
        Gx0=floor(Gm./2)+1;
        Gy0=Gx0;
        %smaller number in denominator results in a wider curve 
        %larger results in a thinner curve
        Gsig=Gx0./1.7;
    
        for Gc = 1 : Gm
          for Gr = 1 : Gm
             W(Gr, Gc) = exp(-1*(((Gc-Gx0)^2)./(2*Gsig^2)+((Gr-Gy0)^2)./(2*Gsig^2)));
          end
        end
        %central weight is always 0
        W(o+1,o+1)=0;
    end %end of Ga==1
          
    %convolution
    WZj=conv2(kep_Zi,W,'same');
    
    %number of neighbors for each pixel of image
    W_2d=ones(M,N);
    nS=conv2(W_2d,W,'same');
    
    %lagged Zi
    sumWZj_n=WZj./nS;
    
    
    %Moran's I
    kep_I=kep_Zi(:,:,1).*sumWZj_n;
    
    %original image
    figure;
    imshow(kep);
    title('Original image');
    movegui([25.*(oo-1)+600, scrS(1,4)-500]);
   
    %Moran's global I scatterplot
    figure;
    xpl=reshape(kep_Zi(:,:,1),mer,1);
    ypl=reshape(sumWZj_n,mer,1);
    scatter(xpl,ypl);
    %calculating the slope fitted on all points, i.e. the global Moran's I
    Pf = polyfit(xpl,ypl,1);
    yfit = Pf(1)*xpl+Pf(2);
    hold on;
    plot(xpl,yfit,'r-.');
    axis equal
    title({'Global Moran''s I plot';['ord=' num2str(o) ', I='  num2str(Pf(1))]});
    plot([0 0], [min(ypl)-0.1*(max(ypl)-min(ypl)) max(ypl)+0.1*(max(ypl)-min(ypl))], '-- k');
    plot([min(xpl)-0.1*(max(xpl)-min(xpl)) max(xpl)+0.1*(max(xpl)-min(xpl))], [0 0], '-- k');
    xlim([min(xpl)-0.1*(max(xpl)-min(xpl)) max(xpl)+0.1*(max(xpl)-min(xpl))]);
    ylim([min(ypl)-0.1*(max(ypl)-min(ypl)) max(ypl)+0.1*(max(ypl)-min(ypl))]);
    
    movegui([25.*(oo-1), scrS(1,4)-500]);
    Iertekek(oo)=Pf(1);
    hold off
    
    %calculating pseudosignificance
    %permutations
    p_p=zeros(M,N);
    %the process above will be repeated with randomly distributed pixels
    %ism times
    for p=1:ism
        kep_I_rnd=zeros(M,N);
        d_rnd=randperm(M*N);
        Zi_rnd=zeros(M*N,1);
        Zi_rnd(:,:)=Zi(d_rnd(:,:));
        Zi_rnd=reshape(Zi_rnd,M,N);
        WZj_rnd=conv2(Zi_rnd,W,'same');
        sumWZj_rnd=WZj_rnd./nS;
        %az I értékeket rögtön a perm mtx-ba írom
       kep_I_rnd(:,:)=kep_Zi(:,:).*sumWZj_rnd;
        poz_nagy=zeros(M,N);
        neg_kics=zeros(M,N);
        poz_nagy(kep_I(:,:)>0 & kep_I_rnd(:,:)>=kep_I(:,:))=1;
        neg_kics(kep_I(:,:)<0 & kep_I_rnd(:,:)<=kep_I(:,:))=1;
        p_p=p_p+poz_nagy+neg_kics;
    end
    
   %probability of a pixel value is higher than the original calculation
    p_perm=(p_p)./(ism+1);
        
    %significance levels
    szign=[0.05;0.01;0.001;0.0001];
    
    %labeling of clasters according to their values and their neighborhood
    %0=non significant
    %1=high-high
    %2=low-low
    %3=low-high
    %4=high-low
    
    kl=zeros(mer,1);
    sig=zeros(mer,1);
    
    P_perm=reshape(p_perm,M*N,1);
    
    szi=2;
    %all pixels where p_perm<0.05 will be clasificated as cluster
    kl(P_perm<=szign(szi) & xpl>0 & ypl>0)=1;
    kl(P_perm<=szign(szi) & xpl<0 & ypl<0)=2;
    kl(P_perm<=szign(szi) & xpl<0 & ypl>0)=3;
    kl(P_perm<=szign(szi) & xpl>0 & ypl<0)=4;
    
    %to achive the same color coding in any cases, all the values must be
    %included in the upper left corner of the image
   
    kl(1)=1;
    kl(2)=2;
    kl(3)=3;
    kl(4)=4;
    
    %map of significance
    sig(P_perm<=szign(1))=1;
    sig(P_perm<=szign(2))=2;
    sig(P_perm<=szign(3))=3;
    sig(P_perm<=szign(4))=4;
    
    kep_kl=reshape(kl,M,N);
    %writing to csv
    if csvw>0
        %the name contains the order too
        nev=strcat(ne,num2str(o),'.csv');
        csvwrite(nev,kep_kl);
    end
    kep_sig=reshape(sig,M,N);
    
    figure
    imagesc(kep_sig);
    %significance map
    title({'Map of significance';['ord=' num2str(o) ', I='  num2str(Pf(1))]});
    axis equal
    ylim([0,M])
     map=[0.9 0.9 0.9; 0.7 1 0.7; 0.5 0.9 0.5; 0.3 0.9 0.3 ;0.1 0.9 0.1];
    colormap(map);
    movegui([25.*(oo-1)+600, scrS(1,4)-1000]);
       
    
    figure
    imagesc(kep_kl);
    %Cluster map
    title({'Clusters';['ord=' num2str(o) ', I='  num2str(Pf(1))]});
    axis equal
    ylim([0,M])
    map=[0.9 0.9 0.9; 1 0.2 0.2; 0.2 0.2 1; 0.7 0.7 1;1 0.7 0.7];
    colormap(map);
    movegui([25.*(oo-1), scrS(1,4)-1000]);
    
    ido(oo)=toc;

end %end oo loop

toc;
