size_xy_cont=0;    
for cont_i=1:size_zx
    for cont_j=1:size_zy   
%         waitbar(size_xy_cont/(size_zx*size_zy));
        size_xy_cont=size_xy_cont+1;
        clear C
        %la variabile BW2->image_morph 
        %la variabile BW->mask_slice_roi(:,:,1) 
        if image_morph(cont_i,cont_j)>0
        
            %leggi lo stesso pixel lungo tutte le acquisizioni e metti i
            %valori in C
            for h=1:size_D(1)
                C(h,1)=x(h);
                C(h,2)=D(h,cont_i,cont_j);
            end
            size(C);
        
            %normalizza i valori asse y del pixel
            C(:,2)=C(:,2)/max(C(:,2));
            C=sortrows(C,1);
            
            %fitta con le B spline
            if reg_flag_weight==1
                [cs,p] = csaps(double(C(:,1)), double(C(:,2)));
                reg_weight=p;
            else
                cs = csaps(double(C(:,1)), double(C(:,2)), reg_weight, [],[]); 
            end
            
            %trova il minimo shiftato dell'acqua e assegna questo valore
            %all'immagine zero
            [minval,min_absc]=fnmin(cs);
            image_zero(cont_i,cont_j)=min_absc;
            
            %correggi i valori di x (RF in ppm) per lo shift trovato
            x2=C(:,1)-min_absc;

            %calcola R2: bontà della curva interpolata di passare per i
            %punti sperimentali
            ss_reg=(C(:,2)-fnval(cs,C(:,1))).^2;
            ss_tot=(C(:,2)-mean(C(:,2))).^2;
            SS_reg=sum(ss_reg(:));
            SS_tot=sum(ss_tot(:));
            R2_pixel=1-(SS_reg/SS_tot);
           
            image_r2(cont_i,cont_j)=R2_pixel;
            
            %%%%%%%%%%%%% ST puntuale classico %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %calcola l'immagine ST corretta nello shift dello zero per quel valore di ppm
            image_ST(cont_i,cont_j)=1-((fnval(cs,st_ppm+min_absc)/fnval(cs,-st_ppm+min_absc)))^(sign_st*flag_st_f_t);
          
            if flag_ratio==1
                image_ST2(cont_i,cont_j)=1-((fnval(cs,st_ppm2+min_absc)/fnval(cs,-st_ppm2+min_absc)))^(sign_st*flag_st_f_t);
            end
           
             
             %------------------------------------%
             % calcolo ST integrale bulk          %
             %------------------------------------%
            
             area_bulk=1*((min_absc-(flag_integral-1)*abs(st_ppm_bulk))-(min_absc-flag_integral*abs(st_ppm_bulk)));
             num=fnval(fnint(cs),[min_absc-flag_integral*abs(st_ppm_bulk), min_absc-(flag_integral-1)*abs(st_ppm_bulk)]);
             den=fnval(fnint(cs),[min_absc+(flag_integral-1)*abs(st_ppm_bulk), min_absc+flag_integral*abs(st_ppm_bulk)]);
            
             %ST integrale di bulk
             image_int_bulk(cont_i,cont_j)=1-((((num(2)-num(1))/(den(2)-den(1)))^flag_st_f_t)); 
          
             %------------------------------------%
             % calcolo ST integrale peak          %
             %------------------------------------%
             
             %area teorica rettangolo centrato sul picco
             area_peak=1*width_integral;
             num_p=fnval(fnint(cs),[min_absc+((-1)^flag_integral)*(st_ppm)-(width_integral/2), min_absc+((-1)^flag_integral)*(st_ppm)+(width_integral/2)]);
             den_p=fnval(fnint(cs),[min_absc-((-1)^flag_integral)*(st_ppm)-(width_integral/2), min_absc-((-1)^flag_integral)*(st_ppm)+(width_integral/2)]);
             
             %ST integrale di peak
             image_int_peak(cont_i,cont_j)=1-((((num_p(2)-num_p(1))/(den_p(2)-den_p(1)))^flag_st_f_t));
                        
             %------------------------------------%
             % calcolo ST integrale combinati     %
             %------------------------------------%
             
             %ST integrale combinato di bulk
             image_comb_bulk(cont_i,cont_j)=1-((((fnval(cs,st_ppm+min_absc)/fnval(cs,-st_ppm+min_absc)))^(sign_st*flag_st_f_t))*((num(2)-num(1))/(den(2)-den(1)))^flag_st_f_t);

             %calcolo integrale combinato ST di peak
             image_comb_peak(cont_i,cont_j)=1-((((fnval(cs,st_ppm+min_absc)/fnval(cs,-st_ppm+min_absc)))^(sign_st*flag_st_f_t))*(((num_p(2)-num_p(1))/(den_p(2)-den_p(1)))^flag_st_f_t));
 
             %salva i coefficienti della spline per ciascun pixel
             cs_coeff(cont_i,cont_j,:,:)= cs.coefs;

        else
             %se il pixel non appartiene all'immagine segmentata, assegna
             %questi valori di default
             image_zero(cont_i,cont_j)=-300;            %valore zero shift corretto
             image_ST(cont_i,cont_j)=0;                 %valore ST% puntuale
             
             if flag_ratio==1
                 image_ST2(cont_i,cont_j)=0;                 
             end 
           
             image_int_bulk(cont_i,cont_j)=0;           %valore ST integrale di bulk
             image_int_peak(cont_i,cont_j)=0;           %valore ST integrale di peak    
             image_r2(cont_i,cont_j)=0;                 %valore immagine pixel con R2 < 0.97
             
             %salva i coefficienti della spline per ciascun pixel
%             cs_coeff(cont_i,cont_j,:,:)= no_coeff;
        end
        clear cs;      
    end     %fine ciclo pixel asse y
end     %fine ciclo pixel asse x

%scegli se visualizzare le mappe solo con la maschera della roi selezionata o con la
%maschera della segmentazione
roiomorf = questdlg('mappe con maschera roi o maschera morfologica?','','roi','morfologica','morfologica');
if strcmp(roiomorf,'roi')
    image_morph=zeros(size_zx,size_zy);
   for h=1:num_roi
       image_morph=imadd(image_morph, mask_slice_roi(:,:,h));
   end
end
