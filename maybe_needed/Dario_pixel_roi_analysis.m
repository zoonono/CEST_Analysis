
clear all;

iptsetpref('ImshowInitialMagnification','fit')
directory_pre_all=pwd;
%dir_work=uigetdir;
%cd(dir_work);

%apri il file .mat con le analisi salvate
[file_v,path_v]=uigetfile('*.mat','File con le variabili salvate');
cd(path_v);

load([path_v, file_v]);
close all

%mostra l'immagine con la mappa ST calcolata reale
h_st=figure, imshow(image_seg, [min(image_seg(:)) max(image_seg(:))]);
colormap gray
freezeColors
hold on
if flag_st_f_t==1
   imshow((1)*(medfilt2(double(image_ST),[3 3])), [0 0.5])
end
if flag_st_f_t==-1
   imshow((-1)*(medfilt2(double(image_ST),[3 3])), [0 0.5])
end
%colormap hot
colormap(gca, mycmap)
%set(h_st,'Colormap',mycmap) 
alpha 0.7
colorbar
if size_zx-size_zy>0 | size_zx-size_zy<0
   axis square
end
title(['ST superimposed on morphological image at ', num2str(st_ppm), ' ppm' ])
%impixelinfo(h_st)


%scegli se analizzare un pixel oppure una roi
button = questdlg('Analisi sul pixel o su una roi?','','Pixel','ROI','ROI-auto','Pixel');

reg_weight=str2num(char(inputdlg('Which Regularization Factor ?' ,'',1 )));

reg_model = questdlg('Choose b-spline model: ','','sym','asym','sym');

Rsqrfilter= str2num(char(inputdlg('Which Rsqr filter ?' ,'',1, {'0.9'} )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PIXEL               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(button,'Pixel')
   
  
   title({'Seleziona un pixel - poi premi invio'; '(per uscire clicca su un pixel nero)'});
   set(h_st,'Units', 'Normalized')
   set(h_st,'Position', [0.01 0.09 0.5 0.7])
   
   flag_pixel=1;
   flag_fig=0;
   num_pixel=1;
   
   [px,py]= getpts(h_st);
   px=round(px)
   py=round(py)
   %val_sel=image_cor(py,px,1);
   text(px, py, 'o', 'FontSize', 10, 'Color','blue')
   text(px+1, py, num2str(num_pixel), 'FontSize', 16, 'Color','blue')
   
   
   while(flag_pixel==1)
      
   %calcola lo spettro zeta del pixel
   
   %fitta i valori medi della roi selezionata con le B-spline
   %le x sono le frequenze in ppm a cui sono state aquisite le immagini
   %le y sono i valori medi normalizzati delle immagini nella roi selezionata
   int_pixel=(D(:,py,px));
   B=[x;int_pixel];
   %figure, plot(x,squeeze(D(:,py,px)),'ro');
   B=B';
   B=sortrows(B,1);
   %normalizza i valori di intensitï¿½ 
   B(:,2)=B(:,2)/max(B(:,2));

   %Figura 2 SPETTRO Z
   
   if flag_fig==1
       clf(2);
   end
   figure(2); 
   subplot (2,1,1); plot(B(:,1),B(:,2),'o') %punti sperimentali originali
   
   %se il fattore di regolarizzazione ï¿½ automatico, reg_weight non esiste
   %%if exist('reg_weight','var') == 0
   %   reg_weight=reg_spectrum;
   %%end
   
   %[spline_i,reg_spectrum]=csaps(B(:,1),B(:,2));
   
   if strcmp(reg_model,'asym')
      cs = csaps(B(:,1), B(:,2), reg_weight, [], [repmat(1,1,27), repmat(50,1,14), repmat(1,1,2) ] );
   else
       [cs,p] = csaps(double(B(:,1)), double(B(:,2)), reg_weight);
   end
       
   %cs = csaps(B(:,1),B(:,2), reg_weight,[] , []);
   [minval_mean,min_absc_mean]=fnmin(cs)
   hold on, subplot (2,1,1); fnplt(cs)  %B spline
   %hold on, plot(B(:,1)-min_absc_mean,B(:,2),'x', 'Color','g', 'MarkerSize',12)
   xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
   ylabel('Normalized Intensity Values','Fontname','axial','Fontsize',14);
   
   %saveas(h_glob_z,[num2str(init-2) '_global_roi_zeta' '.jpg']);

   %utilizza il valore interpolato dalla B spline per ricavare il valore di y
   %a passi di 0.1 ppm centrando lo zero sul valore minimo assoluto della B spline
   delta=0.1;
   for i=1:(var1*10)
       y_abscissa_pos(i)=(fnval(cs,i*delta+min_absc_mean));
       abscissa_pos(i)=0.1*i;
   end
   for i=1:(var1*10)
       y_abscissa_neg(i)=(fnval(cs,-i*delta+min_absc_mean));
       abscissa_neg(i)=0.1*(-i);
   end
   
   %aggiungi alla figura i punti calcolati dal fitting della B spline
   hold on, subplot (2,1,1); plot(abscissa_neg,y_abscissa_neg,'Color','g','LineWidth',2);
   hold on, subplot (2,1,1); plot(abscissa_pos,y_abscissa_pos,'Color','g','LineWidth',2)

   %calcola Rsqr: bontï¿½ della curva interpolata di passare per i
   %punti sperimentali
           
   ss_reg=(B(:,2)-fnval(cs,B(:,1))).^2;
   ss_tot=(B(:,2)-mean(B(:,2))).^2;
   SS_reg=sum(ss_reg(:));
   SS_tot=sum(ss_tot(:));
   Rsqr_pixel=1-(SS_reg/SS_tot);

   title(['Z-spectra pixel ', num2str(num_pixel), ' (' num2str(px), ',' num2str(py), ')' ' - Rsqr= ' num2str(Rsqr_pixel)  ])
   
   %calcola il valore ST% nel pixel selezionato
   %calcola il valore di ST% come 1-[(Is/Io)*100]^sign_st
   for i=1:(var1*10)
       diff_ST(i)=(1-((y_abscissa_pos(i))/y_abscissa_neg(i))^((sign_st)*(flag_st_f_t)));
   end
   %calcola dove cade il massimo dell' ST% globale
   [st_glob_max,ppm_glob_max]=max(diff_ST);
   st_glob_max=st_glob_max*100
   ppm_glob_max=ppm_glob_max*delta

   % Figura 3 ST%
   %h_glob_st=figure, 
   subplot (2,1,2); plot(abscissa_pos,flag_st_f_t*diff_ST*100,'LineWidth',2)
   axis([0 var1 -Inf Inf]);
   xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
   ylabel('ST%','Fontname','axial','Fontsize',16);
   title(['Effetto ST% pixel ', num2str(num_pixel), ' (' num2str(px), ',' num2str(py), ')' ]);
   
   set(2,'Units', 'Normalized')
   set(2,'Position', [0.53 0.09 0.45 0.8])
   
   
   %flag_pixel=0
   
   [px,py]= getpts(h_st);
   px=round(px)
   py=round(py)
   
   
   if image_seg(py,px)==0
      flag_pixel=0;
   end
   
   flag_fig=1;
   num_pixel=num_pixel+1;
   text(px, py, 'o', 'FontSize', 10, 'Color','blue')
   text(px+1, py, num2str(num_pixel), 'FontSize', 16, 'Color','blue')
   
   end  %fine ciclo while
  
end

%forse rifittando i dati non ottengo gli stessi valori di effetto ST
%oppure ï¿½ dovuto al fatto che l'immagine con la mappa ST ï¿½ filtrata da un
%filtro di tipo mediano, posso provare a toglierlo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               ROI                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(button,'ROI')

   title({'Seleziona una ROI '; '(per uscire seleziona una roi esterna)'});
   
   set(h_st, 'Units', 'Normalized');
   set(h_st, 'Position', [0.01 0.09 0.5 0.7]);
   %set(h_st, 'KeyPressFcn',@printfig);
   
   
   flag_roi=1;
   flag_fig=0;
   num_roi=1;
   
   LoadIt = questdlg('Vuoi caricare la maschera di una roi ?','','Si','No','Si');
   if strcmp(LoadIt,'Si')
        
       [file_m,path_m]=uigetfile('*.mat','File con la maschera della ROI');
        load(file_m);
   else 
        %disegna roi e assegna a BW2 la maschera della roi
       [BW2,x_p,y_p]=roipoly;   
   end
   
   %salva la maschera
    SaveIt = questdlg('Vuoi salvare la maschera della roi ?','','Si','No','Si');
   if strcmp(SaveIt,'Si')
        
        save maschera BW2 x_p y_p
        
   end
   
   
   while(flag_roi==1)
   
    %prendi le coordinate della roi e contornala con una linea
    for n=1:(size(x_p,1)-1)
        h_line=line([x_p(n) x_p(n+1)],[y_p(n) y_p(n+1)]);
        set(h_line,'color','g');
        set(h_line,'LineWidth',2);
    end
    cx=mean(x_p);
    cy=mean(y_p);
    text(cx, cy, num2str(num_roi), 'FontSize', 16, 'Color','blue') 
     
    cont_p=1;
    
    for i=1:size_x_res
            for j=1:size_y_res
                if BW2(i,j) > 0
    
                      vect=D(:,i,j);
                     
                      P=[x;vect'];
                      P=P';
                      P=sortrows(P,1);
                      Pnonorm=P;
                      %normalizza i valori di intensità del pixel
                      maxP=max(P(:,2));
                      P(:,2)=P(:,2)/max(P(:,2));
                      dati_pixel(cont_p,:,:)=P(:,:);
                      %fitta i valori del pixel selezionato con le B-spline
                      if strcmp(reg_model,'asym')
                          cs = csaps(P(:,1), P(:,2), reg_weight, [], [repmat(1,1,27), repmat(50,1,14), repmat(1,1,2) ] );
                      else
                          cs = csaps(P(:,1),P(:,2), reg_weight,[] , []);
                      end
                     
                      [minval_mean,min_absc_mean]=fnmin(cs);
                      
                      %calcola Rsqr: bontï¿½ della curva interpolata di passare per i punti sperimentali
           
                      ss_reg=(P(:,2)-fnval(cs,P(:,1))).^2;
                      ss_tot=(P(:,2)-mean(P(:,2))).^2;
                      SS_reg=sum(ss_reg(:));
                      SS_tot=sum(ss_tot(:));
                      Rsqr_pixel=1-(SS_reg/SS_tot);
            
            
                      if Rsqr_pixel > Rsqrfilter
                                            
                      delta=0.1;
                      v=delta:delta:10;
                      y_abscissa_pos=(fnval(cs,v+min_absc_mean));
                      abscissa_pos=v;
                      y_abscissa_neg=(fnval(cs,-v+min_absc_mean));
                      abscissa_neg=-v;
                      
                      diff_ST=(1.-((y_abscissa_pos)./y_abscissa_neg).^((sign_st).*(flag_st_f_t)));
                      
                      %for v=1:(var1*10)
                      %    y_abscissa_pos(v)=(fnval(cs,v*delta+min_absc_mean));
                      %    abscissa_pos(v)=0.1*v;
                      %end
                      %for v=1:(var1*10)
                      %     y_abscissa_neg(v)=(fnval(cs,-v*delta+min_absc_mean));
                      %     abscissa_neg(v)=0.1*(-v);
                      %end

                      %for v=1:(var1*10)
                      %    diff_ST(v)=(1-((y_abscissa_pos(v))/y_abscissa_neg(v))^((sign_st)*(flag_st_f_t)));
                      %end
                      
                      st_ppm_4e2(cont_p)=diff_ST(4.2*10)*100;
                      st_ppm_5e5(cont_p)=diff_ST(5.5*10)*100;
                      clear P;
                      clear maxP;
                      clear vect;
                      clear Rsqr_pixel;
                      
                      cont_p=cont_p+1;  %contatore del numero di pixel della roi                    
                      
                      end
                      
                end
            end
     end
        
        
        if flag_fig==1
           clf(2);
        end
        figure(2);
        dati_roi=mean(dati_pixel);
        std_roi=std(dati_pixel(:,:,2));
        maxroi=max(dati_roi(1,:,2));
        dati_roi(1,:,2)=dati_roi(1,:,2)/maxroi;
        
        plot(dati_roi(1,:,1),dati_roi(1,:,2),'o', 'Color', 'b' ) %punti sperimentali originali
        
        
        if strcmp(reg_model,'asym')
           csroi = csaps(dati_roi(1,:,1), dati_roi(1,:,2), reg_weight, [], [repmat(1,1,27), repmat(50,1,14), repmat(1,1,2) ] );
        else
            csroi = csaps(dati_roi(1,:,1), dati_roi(1,:,2), reg_weight, [], []);
        end
       
        
        
        %csroi = csaps(dati_roi(1,:,1),dati_roi(1,:,2), reg_weight,[] , []);
        [minval_roi,min_absc_roi]=fnmin(csroi);
        
        delta=0.1;
        v=delta:delta:var1;
        roi_abscissa_pos=(fnval(csroi,v+min_absc_roi));
        roi_abscissa_neg=(fnval(csroi,-v+min_absc_roi));
        roi_pos=v;
        roi_neg=-v;
       
        diff_ST=(1.-((roi_abscissa_pos)./roi_abscissa_neg).^((sign_st).*(flag_st_f_t)));
       %diff_MTR=(y_abscissa_neg - y_abscissa_pos)./abs(max(y_abscissa_neg));
        diff_MTR=(roi_abscissa_neg - roi_abscissa_pos);
        
        %hold on, fnplt(csroi)  %B spline non centrata sullo zero offset
        hold on, plot([fliplr(roi_neg),0.0,roi_pos],[fliplr(roi_abscissa_neg),fnval(csroi,min_absc_roi),roi_abscissa_pos],'Color', 'g','LineWidth',2);
       
   
        mean_ST_4e2(num_roi)=mean(st_ppm_4e2);
        mean_ST_5e5(num_roi)=mean(st_ppm_5e5);
        median_ST_4e2(num_roi)=median(st_ppm_4e2);
        median_ST_5e5(num_roi)=median(st_ppm_5e5);
        std_ST_4e2(num_roi)=std(st_ppm_4e2);        
        std_ST_5e5(num_roi)=std(st_ppm_5e5);
        
    
        disp(['ROI numero ' num2str(num_roi) ])
        disp(['MEAN ST 4.2 ppm= ' num2str((mean_ST_4e2(num_roi)))])
        %disp(['MEDIAN ST 4.2 ppm= ' num2str((median_ST_4e2(num_roi)))])
        disp(['SD ST 4.2 ppm= ' num2str(std_ST_4e2(num_roi))])
        disp(['MEAN ST 5.5 ppm= ' num2str((mean_ST_5e5(num_roi)))])
        %disp(['MEDIAN ST 5.5 ppm= ' num2str((median_ST_5e5(num_roi)))])
        disp(['SD ST 5.5 ppm= ' num2str(std_ST_5e5(num_roi))])
        
       
   
   %salva i dati x e y della roi scelta
   SaveIt = questdlg('Vuoi salvare i dati x e y della roi ?','','Si','No','Si');
   if strcmp(SaveIt,'Si')
        
        fileout=[path_v char(file_v(1:2)) '_data_roi_man.txt' ]
        fid_data=fopen(fileout,'a+');
        
        fprintf(fid_data,['ROI numero ' num2str(num_roi) ]);
        
        fprintf(fid_data, '\n');
        
        %printa i punti sperimentali medi della roi
        for k = 1:length(x)          
            fprintf(fid_data, '%5.5f\t %5.5f\t %5.5f\n' , dati_roi(1,k,1), dati_roi(1,k,2), std_roi(k) );
        end
        fprintf(fid_data, '\n\n');
        
        
    
        for v=(10/delta):-1:1          
            fprintf(fid_data, '%5.5f\t %5.5f\t \n' , roi_neg(v), roi_abscissa_neg(v) );
        end
        fprintf(fid_data, '%5.5f\t %5.5f\t \n' , 0.0, fnval(csroi,(0.0*delta+min_absc_roi)) );
        for v=1:(10/delta)
            fprintf(fid_data, '%5.5f\t %5.5f\t %5.5f\t %5.5f\t \n' , roi_pos(v), roi_abscissa_pos(v), diff_ST(v), diff_MTR(v) );
        end
        fprintf(fid_data, '\n\n');
        fclose(fid_data);
       
   end
   
 

   
   %per ridisegnare la roi, devi ridare il focus alla prima figura
   %set(1,'Selected','on')
   figure(h_st);
   [BW2,x_p,y_p]=roipoly;
   
   if image_morph(round(y_p(1)), round(x_p(1)))==0
      flag_roi=0;
   end
   
   flag_fig=1;
   num_roi=num_roi+1;
   clear vect
   clear mean_roi
   
   end  %fine ciclo while
  
  
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               ROI-automatica      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(button,'ROI-auto')

    fileout=[path_v char(num2str(init)) '_data_roi_auto.txt' ]
    fid_data=fopen(fileout,'a+');

h=1;
colori = ['g' 'r' 'b' 'k' 'y' 'm' 'c' 'g' 'r' 'b' 'k' 'y' 'm' ];

dati_sperim=zeros(num_roi+1,numel(x));
dati_sperim(1,:)=B(:,1);
id_elab = waitbar(0,'Calculating ...');
cont_roi=1;

for h=1:num_roi

    fprintf(fid_data,['ROI numero ' num2str(h) ]);
    
    waitbar(h/num_roi);
    roi_number=1;
    clear vect;
    
    cont_p=1;
    
    
          
        for i=1:size_x_res
            for j=1:size_y_res
                if mask_slice_roi(i,j,h) > 0
                    
                    if image_ST(i,j) > 0.02
                      
                   if roi_number == 1 
                      figure(h_st);
                      text(j, i, num2str(h), 'FontSize', 16, 'Color','blue');
                      roi_number=0;
                   end
              
                      vect=D(:,i,j);
                     
                      P=[x;vect'];
                      P=P';
                      P=sortrows(P,1);
                      Pnonorm=P;
                      %normalizza i valori di intensità del pixel
                      maxP=max(P(:,2));
                      P(:,2)=P(:,2)/max(P(:,2));
                      dati_pixel(cont_p,:,:)=P(:,:);
                      %fitta i valori del pixel selezionato con le B-spline
                      if strcmp(reg_model,'asym')
                          
                          if strcmp(coil_1H,  '<RF RES 125 1H 059/035 QUAD TR>')   %B0 = 3T
                        
                        freqneg= sum(P(:,1)<=0);
                       
                       %[cs,p] = csaps(double(C(:,1)), double(C(:,2)), reg_weight, [], [repmat(1,1,freqneg), repmat(50,1,freqpos) ] ); 
                       
                       [cs1,p1] = csaps(double(P(1:freqneg,1)), double(P(1:freqneg,2)), 0.6, [], []);
                       [cs2,p2] = csaps(double(P(freqneg:end,1)), double(P(freqneg:end,2)), reg_weight, [], []);
                       
                       cs=struct('form', 'pp', 'breaks', zeros(1,cs_pezzi+1), 'coefs', zeros(cs_pezzi,4), 'pieces', cs_pezzi, 'order', 4, 'dim', 1);
                       cs.breaks=[cs1.breaks cs2.breaks(2:end)];
                       cs.coefs=cat(1,cs1.coefs, cs2.coefs);
                       
                       
                          end
                          
                          
                          %cs = csaps(P(:,1), P(:,2), reg_weight, [], [repmat(1,1,27), repmat(50,1,14), repmat(1,1,2) ] );
                      else
                          cs = csaps(P(:,1),P(:,2), reg_weight,[] , []);
                      end
                     
                      [minval_mean,min_absc_mean]=fnmin(cs);
                      
                      %calcola Rsqr: bontà della curva interpolata di passare per i punti sperimentali
           
                      ss_reg=(P(:,2)-fnval(cs,P(:,1))).^2;
                      ss_tot=(P(:,2)-mean(P(:,2))).^2;
                      SS_reg=sum(ss_reg(:));
                      SS_tot=sum(ss_tot(:));
                      Rsqr_pixel=1-(SS_reg/SS_tot);
            
            
                      if Rsqr_pixel > Rsqrfilter
                                            
                      delta=0.1;
                      v=delta:delta:abs(max(x));
                      y_abscissa_pos=(fnval(cs,v+min_absc_mean));
                      abscissa_pos=v;
                      y_abscissa_neg=(fnval(cs,-v+min_absc_mean));
                      abscissa_neg=-v;
                      
                      diff_ST=(1.-((y_abscissa_pos)./y_abscissa_neg).^((sign_st).*(flag_st_f_t)));
                      %%diff_ST_Rex=((1/y_abscissa_pos)-(1/y_abscissa_neg));
                      
                     % for v=1:(var1*10)
                     %     y_abscissa_pos(v)=(fnval(cs,v*delta+min_absc_mean));
                     %     abscissa_pos(v)=0.1*v;
                     % end
                     % for v=1:(var1*10)
                     %      y_abscissa_neg(v)=(fnval(cs,-v*delta+min_absc_mean));
                     %      abscissa_neg(v)=0.1*(-v);
                     % end

                      %for v=1:(var1*10)
                      %    diff_ST(v)=(1-((y_abscissa_pos(v))/y_abscissa_neg(v))^((sign_st)*(flag_st_f_t)));
                      %end
                      
                      %st_ppm_4e2(cont_p)=diff_ST(4.2*10)*100;
                      %st_ppm_5e5(cont_p)=diff_ST(5.5*10)*100;
                      zero_shift(cont_p)=min_absc_mean;
                      clear P;
                      clear maxP;
                      clear vect;
                      clear Rsqr_pixel;
                      
                      cont_p=cont_p+1;  %contatore del numero di pixel della roi                    
                      
                      end
                      
                    else
                        %dati_pixel=[x;0;0]
                    end
                      
                end
            end
        end
        
        
        
        figure(2);
        
        if exist('dati_pixel','var')
        
        if cont_p==2
            dati_roi=dati_pixel;
            std_roi=dati_pixel;
        else
             dati_roi=mean(dati_pixel);
             std_roi=std(dati_pixel(:,:,2));
             maxroi=max(dati_roi(1,:,2));
             dati_roi(1,:,2)=dati_roi(1,:,2)/maxroi;
        end
        plot(dati_roi(1,:,1),dati_roi(1,:,2),'o', 'Color', colori(h) ) %punti sperimentali originali
        dati_sperim(h+1,:)=dati_roi(1,:,2);
        
        if strcmp(reg_model,'asym')
            
            if strcmp(coil_1H,  '<RF RES 125 1H 059/035 QUAD TR>')   %B0 = 3T
                        
                        freqneg= sum(sort(x)<=0);
                       
                       %[cs,p] = csaps(double(C(:,1)), double(C(:,2)), reg_weight, [], [repmat(1,1,freqneg), repmat(50,1,freqpos) ] ); 
                       
                       [cs1,p1] = csaps(dati_roi(1,1:freqneg,1), dati_roi(1,1:freqneg,2), 0.6, [], []);
                       [cs2,p2] = csaps(dati_roi(1,freqneg:end,1), dati_roi(1,freqneg:end,2), reg_weight, [], []);
                       
                       csroi=struct('form', 'pp', 'breaks', zeros(1,cs_pezzi+1), 'coefs', zeros(cs_pezzi,4), 'pieces', cs_pezzi, 'order', 4, 'dim', 1);
                       csroi.breaks=[cs1.breaks cs2.breaks(2:end)];
                       csroi.coefs=cat(1,cs1.coefs, cs2.coefs);
                       
            end
            
           %csroi = csaps(dati_roi(1,:,1), dati_roi(1,:,2), reg_weight, [], [repmat(1,1,27), repmat(50,1,14), repmat(1,1,2) ] );
        else
            csroi = csaps(dati_roi(1,:,1), dati_roi(1,:,2), reg_weight, [], []);
        end
       
        
        
        %csroi = csaps(dati_roi(1,:,1),dati_roi(1,:,2), reg_weight,[] , []);
        [minval_roi,min_absc_roi]=fnmin(csroi);
        
        delta=0.1;
        v=delta:delta:abs(max(x));
        roi_abscissa_pos=(fnval(csroi,v+min_absc_roi));
        roi_abscissa_neg=(fnval(csroi,-v+min_absc_roi));
        roi_pos=v;
        roi_neg=-v;
       
        diff_ST=(1.-((roi_abscissa_pos)./roi_abscissa_neg).^((sign_st).*(flag_st_f_t)));
        diff_ST_Rex=((1./roi_abscissa_pos)-(1./roi_abscissa_neg));
       %diff_MTR=(y_abscissa_neg - y_abscissa_pos)./abs(max(y_abscissa_neg));
        diff_MTR=(roi_abscissa_neg - roi_abscissa_pos);
        
        
        
        %hold on, fnplt(csroi)  %B spline non centrata sullo zero offset
        hold on, plot([fliplr(roi_neg),0.0,roi_pos],[fliplr(roi_abscissa_neg),fnval(csroi,min_absc_roi),roi_abscissa_pos],'Color',colori(h),'LineWidth',2);
       
   
        figure(3);
        hold on, plot(roi_pos,diff_ST,'Color',colori(h),'LineWidth',2);
       % hold on, plot(roi_pos,diff_ST_Rex,'Color',colori(h),'LineWidth',2, 'LineStyle', '--');
        
        %mean_ST_4e2(h)=mean(st_ppm_4e2);
        %mean_ST_5e5(h)=mean(st_ppm_5e5);
        %median_ST_4e2(h)=median(st_ppm_4e2);
        %median_ST_5e5(h)=median(st_ppm_5e5);
        %std_ST_4e2(h)=std(st_ppm_4e2);        
        %std_ST_5e5(h)=std(st_ppm_5e5);
        mean_zero_shift(h)=mean(zero_shift);
        std_zero_shift(h)=std(zero_shift);
        
      
    
        disp(['ROI numero ' num2str(h) ])
        %disp(['MEAN ST 4.2 ppm = ' num2str((mean_ST_4e2(h)))])
        %disp(['MEDIAN ST 4.2 ppm = ' num2str((median_ST_4e2(h)))])
        %disp(['SD ST 4.2 ppm = ' num2str(std_ST_4e2(h))])
        %disp(['MEAN ST 5.5 ppm = ' num2str((mean_ST_5e5(h)))])
        %disp(['MEDIAN ST 5.5 ppm = ' num2str((median_ST_5e5(h)))])
        %disp(['SD ST 5.5 ppm = ' num2str(std_ST_5e5(h))])
        
        disp(['MEAN ZERO SHIFT = ' num2str((mean_zero_shift(h)))])
        disp(['SD ZERO SHIFT = ' num2str(std_zero_shift(h))])
        
        
        fprintf(fid_data, '\n');
        
        %printa i punti sperimentali medi della roi
        for k = 1:length(x)          
            fprintf(fid_data, '%8.4f\t %8.4f\t %8.4f\n' , dati_roi(1,k,1), dati_roi(1,k,2), std_roi(k) );
        end
        fprintf(fid_data, '\n\n');
        
        
        fprintf(fid_data, '\n');
    
        for v=round((abs(max(x))/delta)):-1:1          
            fprintf(fid_data, '%8.4f\t %8.4f\t \n' , roi_neg(v), roi_abscissa_neg(v) );
        end
        fprintf(fid_data, '%8.4f\t %8.4f\t \n' , 0.0, fnval(csroi,(0.0*delta+min_absc_roi)) );
        for v=1:round((abs(max(x))/delta))  
            fprintf(fid_data, '%8.4f\t %8.4f\t %8.4f\t %8.4f\t %8.4f\n' , roi_pos(v), roi_abscissa_pos(v), diff_ST(v), diff_ST_Rex(v), diff_MTR(v) );
            
        end
        fprintf(fid_data, '\n\n');
      
        
        yinterp_neg(h,:)=roi_abscissa_neg;
        yinterp_pos(h,:)=roi_abscissa_pos;
        
        clear st_ppm_4e2;
        clear st_ppm_5e5;
        clear dati_pixel;
        clear dati_roi;
        clear csroi;
        clear roi_abscissa_pos;
        clear roi_pos;
        clear roi_abscissa_neg;
        clear roi_neg;
        clear diff_ST;
        clear diff_MTR;
        
       
        legendInfo{cont_roi} = strcat(['roi ' num2str(h)]);
        legendInfo{(cont_roi+1)} = strcat(['roi ' num2str(h)]);
        cont_roi=cont_roi+2;
        end
        
        
end

figure(2);
legend(legendInfo, 'Location','SouthWest');
%legend('roi1','','roi2','','roi3','','roi4','','roi5','','roi6','','roi7','','roi8','','roi9','','roi10','Location','SouthWest');
saveas(figure(2),[num2str(init) '_Zspectra_' num2str(reg_weight) '_' reg_model '.jpg']);
saveas(figure(2),[num2str(init) '_Zspectra_' num2str(reg_weight) '_' reg_model '.fig']);

figure(3);
legend(legendInfo, 'Location','NorthEast');
%legend('roi1','','roi2','','roi3','','roi4','','roi5','','roi6','','roi7','','roi8','','roi9','','roi10','Location','NorthEast');
saveas(figure(3),[num2str(init)  '_STcurve_' num2str(reg_weight) '_' reg_model '.jpg']);
saveas(figure(3),[num2str(init)  '_STcurve_' num2str(reg_weight) '_' reg_model '.fig']);




fclose(fid_data);
close(id_elab);

dati_var=['data_' num2str(power) 'uT'];
save(dati_var,'dati_sperim');


end



%esempi di callback su pressione di tasti "KeyPressFcn"

%arg1=init;
%arg2=st_ppm;
%h_st=figure('KeyPressFcn',{@printfig,arg1,arg2}), imshow((1)*(medfilt2(double(image_cor),[3 3])), [0 0.5])
%h_st=figure('KeyPressFcn',@savefig), imshow((1)*(medfilt2(double(image_cor),[3 3])), [0 0.5])

%function printfig(src,evnt)
%      if evnt.Character == 's'
%        saveas(h_st,[num2str(init-2) '_ST_map_' num2str(st_ppm) '_ppm' '.jpg']);
%        saveas(h_st,[num2str(init-2) '_ST_map_' num2str(st_ppm) '_ppm' '.fig']);
%      %elseif length(evnt.Modifier) == 1 & %strcmp(evnt.Modifier{:},'control') & evnt.Key == 't'
%         %print ('-dtiff','-r200',['-f' num2str(src)])
%      end




