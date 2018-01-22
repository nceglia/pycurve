clear all

%%P-Y curves derivation%%



%%%%% PARAMETERS TO UPDATE %%%%%%%%%%%%%%%%%%
resolution = 20;
deflection_levels = 26; %considered deflections until 3 in included
pos_deflection_levels = deflection_levels *0.5;
neg_deflection_levels = deflection_levels *0.5;
elevation_resolution = linspace(0,-93,resolution);     
midline = ones(size(elevation_resolution)) * 0;
time = linspace (1,5770,5770);

         curvature_values = zeros(length(elevation_resolution), deflection_levels);
         rotation_values = zeros(length(elevation_resolution), deflection_levels);
         deflection_values = zeros(length(elevation_resolution), deflection_levels);
         moment_values = zeros(length(elevation_resolution), deflection_levels);
         shear_values = zeros(length(elevation_resolution), deflection_levels);
         reaction_values = zeros(length(elevation_resolution), deflection_levels);
%          delta_rotation = zeros(1, deflection_levels);
%          delta_deflection = zeros(1, deflection_levels);
%          delta_shear = zeros(1, deflection_levels);
%          delta_reaction = zeros(1, deflection_levels);
         pos_curvature_values = zeros(length(elevation_resolution), pos_deflection_levels );
         pos_rotation_values = zeros(length(elevation_resolution), pos_deflection_levels);
         pos_deflection_values = zeros(length(elevation_resolution), pos_deflection_levels);
         pos_moment_values = zeros(length(elevation_resolution), pos_deflection_levels);
         pos_shear_values = zeros(length(elevation_resolution), pos_deflection_levels);
         pos_reaction_values = zeros(length(elevation_resolution), pos_deflection_levels); 

         neg_curvature_values = zeros(length(elevation_resolution), neg_deflection_levels);
         neg_rotation_values = zeros(length(elevation_resolution), neg_deflection_levels);
         neg_deflection_values = zeros(length(elevation_resolution), neg_deflection_levels);
         neg_moment_values = zeros(length(elevation_resolution), neg_deflection_levels);
         neg_shear_values = zeros(length(elevation_resolution), neg_deflection_levels);
         neg_reaction_values = zeros(length(elevation_resolution), neg_deflection_levels); 
         
set(0,'DefaultFigureVisible','on') %set to off to stop plotting everything
plot_moment_curvature = 0;
plot_curvature = 0;
plot_rotation = 0;
plot_deflection = 0;
plot_moment = 0;
plot_shear = 0;
plot_reaction = 0;
plot_pos_values = 0;
plot_neg_values = 0;
plot_all_values = 0;
plot_py = 1;


fprintf('\n\nPlot Configuration (1=PlotOn, 0=PlotOff) %i\nMomentCurvature: %i\nCurvature: %i\nRotation: %i\nDeflection: %i\nMoment: %i\nShear: %i\nReaction: %i\nPositiveValues: %i\nNegativeValues: %i\nAllValues: %i\nPY: %i\n\n', plot_moment_curvature, plot_curvature, plot_rotation, plot_deflection, plot_moment, plot_shear, plot_reaction, plot_pos_values, plot_neg_values, plot_all_values, plot_py)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Read text File for 8in PC Rebar Pile%%
load top_loads.txt
load top_deflections.txt
load elevation.txt % took B2 (z=93.5") out in order to substitute this value with soil ground surface elevation z=93"
load Ave1.txt

%%Read text File from LPILE -phi=43- Sand Model= Reese%% 
for i= 1:10
    filename = sprintf('py_%d.txt',i);
    py_lpile = load(filename);
    y_lpile= py_lpile (:,1); %first column = y[in]
    p_lpile= py_lpile (:,2); %second column = p [lbs/in]

    for j = 1 : length(y_lpile)
        y_lpile_values(j,i) = y_lpile(j,1);
        p_lpile_values(j,i) = y_lpile(j,1);
    end
end

r1_strains0 = Ave1 (:, 2:9);%Rebar 1 -B3:I5772 
r1_strains00 = Ave1 (:, 58:60);%Rebar 1 - BF3:BH5772 
r1_strains = [r1_strains0 r1_strains00]; %Rebar 1
r2_strains = Ave1 (:, 61:71); %Rebar 2 - 'BI3:BS5772'

%%Find Max & Min Deflections
        [max_deflection, max_locs] = findpeaks(top_deflections, 'MinPeakDistance', 90, 'MinPeakHeight',0.1);
        DataInv = - top_deflections;
        [min_deflection, min_locs] = findpeaks(DataInv, 'MinPeakDistance', 90, 'MinPeakHeight',0.1);
        min_deflection = - min_deflection;
    peaks_deflection = [max_deflection; min_deflection];
    peaks_deflection( [4 7 10 11 13 14 15 17 19 21 22 24 26 27 29 30 32 33 36 38 46 48 51 52 54 55 56 58 60 62 63 65 67 68 70 71 73 74 77 79], : ) = []; %all deflections
    peaks_locs = [max_locs; min_locs] ;
    peaks_locs( [4 7 10 11 13 14 15 17 19 21 22 24 26 27 29 30 32 33 36 38 46 48 51 52 54 55 56 58 60 62 63 65 67 68 70 71 73 74 77 79], : ) = []; %all deflections
    peaks_deflection ([14 15 16 17 18 19 20 21 35 36 37 38 39 40 41 42],:) = []; %until 3" included
    peaks_locs( [14 15 16 17 18 19 20 21 35 36 37 38 39 40 41 42], : ) = [];%until 3" included
    
 %%Find Max & Min Loads  
    peaks_loads = top_loads(peaks_locs);
 %% Find Max & Min Strain in Rebars
    peaks_r1_strains = r1_strains(peaks_locs,:);
    peaks_r2_strains = r2_strains(peaks_locs,:) ; 
  
%%Read text File from Response%%
load curvature_response.txt %rad/10^6 in 'L3:L88'- sheet=19
load moment_response.txt %ft-kips M3:M88'- sheet=19

curvature_response= curvature_response * 0.000001; %in
moment_response= moment_response * 12; %kips-in
response= spap2 (38, 6, curvature_response, moment_response); %spline to approximate M-Curvature profile from Response 2000

%%Experimental Moment-Curvature
    curvature_original = (peaks_r1_strains - peaks_r2_strains)/5.2; % where 5.2 in is the distance from strain gauges
    curvature_original_top = curvature_original (:,1); %at ground surface
    moment_original_top = peaks_loads * 21.8; %at ground surface
    
    if plot_moment_curvature                                                                   
    figure
    plot(curvature_original_top, moment_original_top,'x');
    hold on
    fnplt(response,'g');
    title('Moment - Curvature Profile')
    xlabel('Curvature, 1/in') 
    ylabel('Moment, kips-in')
    legend('Experimental Curve','Response 2000')
    grid on
    hold off
    end   
   
 curvature_original= curvature_original';
 
%%New Elevation Profile with artificial points
    k0 = 0;  %initial row position
    elevation11 = [elevation(1:k0,:); 93; elevation(k0+1:end,:);0]; %elevation from 93 to 0 in
    elevation111 = (-1)* elevation11;
    elevation1 = sort (elevation111,'descend');

fprintf('Everything Loaded... Starting Computation\n')
 %l= 0; %to read load and deflection value
 
%fprintf('\nRunning Sheet %i\n\n', sheet   
    for i= 1: length(peaks_deflection) 
        curvature = curvature_original (:,i);
        %%NEW CURVATURE PROFILE - add artificial point at the pile tip:
          curvature = [curvature(1:end,:); 0]; %curvature=0 @ z=-93"
        %%Curvature vs Elevation Fitting with Spline Method                                                                          
           w4 = ones(size(curvature));
           %w4(6) = 0;
           w4(5) = 0;
           w4(9) = 0;
           sp_curvature = spap2(8, 3, elevation1, curvature, w4); %spline knots, spline order, x, y, weigths
           curvature_est = fnval(sp_curvature, elevation_resolution);
           
        if plot_curvature
            figure
            plot(elevation1,curvature,'x');
            xlim([-93 0])
            grid on
            str = sprintf('Curvature for Top Deflection= %d in', peaks_deflection(i));
            title(str)
            xlabel('Elevation, in') 
            ylabel('Curvature, 1/in') 
            hold on
            fnplt(sp_curvature,'g');
        end
 
            %%1st Integral & Rotation Boundary conditions%%
                 %slope_top = xlsread(filename,sheet,'E30:F30');
                 %first_integral_bc = slope_top(:,i);
            first_integral = fnint(sp_curvature);      
            rotation_est = fnval(first_integral, elevation_resolution); 
            %rotation_est (1) = first_integral_bc
            %rotation_est2 = fnval(first_integral, elevation_resolution)- (fnval(first_integral, 0) - first_integral_bc);
%                
%             w = ones(size(rotation_est))
%             w(1)=100
%             w(20)=100
%             rotation_test = spap2(17, 4, elevation_resolution, rotation_est,w)
                        
            %%2nd Integral & Deflection Boundary conditions%%
                   %deflection_top = xlsread(filename,sheet,'C30:D30');
                   %second_integral_bc = deflection_top(:,i);
             second_integral = fnint(first_integral); 
%              w = ones(size(rotation_est))
%              w(1)=100
%              w(20)=100
             deflection_est = fnval(second_integral,elevation_resolution); 
            % deflection_est (1) = second_integral_bc
             %deflection_test = spap2(17, 4, elevation_resolution, deflection_est,w)
             %deflection_est2 = fnval(second_integral,elevation_resolution) - (fnval(second_integral,0)- second_integral_bc);
             
        if plot_rotation
            figure
            plot(elevation_resolution, rotation_est, elevation_resolution, midline);
            xlim([-93 0])
            str = sprintf('Rotation for Top Deflection= %d in', peaks_deflection(i));
            title(str)
            xlabel('Elevation, in') 
            ylabel('Rotation, rad')
            grid on
        end
        if plot_deflection
            figure
            plot(elevation_resolution,deflection_est, elevation_resolution, midline);
            xlim([-93 0])
            str = sprintf('Deflection for Top Deflection= %d in', peaks_deflection(i));
            title(str)
            xlabel('Elevation, in') 
            ylabel('Deflection, in')
            grid on
        end
        
  %fprintf('********Finished Curvature_Rotation_Deflection - Sheet %i\n', sheet)
  %%%%%%%%%%to review   
        %%Moment vs Elevation Fitting
            moment_fitted = fnval (response,curvature);
            for j = 1 : length(elevation11)
                moment_original(j,i) = moment_fitted(j,1);
            end
            
            %%Add artificial points at the top
             k0 = 0; %initial row position
                if moment_fitted (1,1) > 0
                        b1 = moment_fitted (1,1) *0.96;
                else
                        b1 = moment_fitted (1,1) * 1.04;
                end
            moment = [moment_fitted(1:k0,:); b1; moment_fitted(k0+1:end,:)];     
            elevation_moment = [elevation1(1:k0,:); 1; elevation1(k0+1:end,:)]; %elevation from +1 in (z=94") to -93 in
               w3 = ones(size(moment));
               %properties for positive deflections
               %w3(1) = 100;
               w3(7) = 100;
               %w3(6) = 100;
               w3(13) = 100;
              %w3(9) = 0;
               sp_moment = spap2(1, 5, elevation_moment, moment, w3); 
              %properties for negative deflections
%                w3(6) = 0;
%                w3(7) = 100;
%                w3(9) = 0;
%                sp_moment = spap2(2, 5, elevation_moment, moment, w3); %knots, order, x, y, weigths -- for better Moment fitting knots= 9, good with order 4
               moment_est = fnval(sp_moment,elevation_resolution);
               
            if plot_moment
                figure
                plot(elevation_moment, moment,'x');
                xlim([-93 1])
                str = sprintf('Moment for Top Deflection= %d in', peaks_deflection(i));
                title(str)
                xlabel('Elevation, in') 
                ylabel('Moment, kips-in')
                grid on
                hold on
                fnplt(sp_moment,'r');
            end
            %elevation_resolution_diff = linspace(1,-93,resolution);
            first_derivative = fnder(sp_moment);
            shear_est = fnval(first_derivative,elevation_resolution);
             %first_derivative_est = fnval(first_derivative,elevation_resolution);
            second_derivative = fnder(sp_moment,2);
            reaction_est = fnval(second_derivative,elevation_resolution);
             %second_derivative_est = -fnval(second_derivative,elevation_resolution);
                                                                                               
            if plot_shear
                figure
                plot(elevation_resolution,shear_est, elevation_resolution, midline);
                xlim([-93 0])
                str = sprintf('Shear for Top Deflection= %d in', peaks_deflection(i));
                title(str)
                xlabel('Elevation, in') 
                ylabel('Shear, kips')
                grid on
            end
            if plot_reaction
                figure
                plot(elevation_resolution,reaction_est, elevation_resolution, midline);
                xlim([-93 0])
                str = sprintf('Soil Reaction for Top Deflection= %d in', peaks_deflection(i));
                title(str)
                xlabel('Elevation, in') 
                ylabel('Soil Reaction, kips/in')
                grid on
            end

         %%Collect all the values and calculate the percentage of difference between measured and computed top boundary conditions
                                                                                       % index = ((sheet-4) * 2) + i;
           for j = 1 : length(elevation_resolution)
                curvature_values(j,i) = curvature_est(1,j);
                rotation_values(j,i)  = rotation_est(1,j);
                deflection_values(j,i)= deflection_est(1,j);
                moment_values(j,i)    = moment_est(1,j);
                shear_values(j,i)     = shear_est(1,j);
                reaction_values(j,i)  = reaction_est(1,j);
           end
            
               % delta_rotation (1, i) = abs((first_integral_bc - rotation_values (1,index))/rotation_values (1,index)) * 100;
                delta_deflection (1, i) = abs((peaks_deflection (i,1) - deflection_values (1,i))/deflection_values (1,i)) * 100;
                delta_moment (1, i) = abs((moment_original_top (i,1) - moment_values (1,i))/moment_values (1,i)) * 100;
                delta_shear (1, i) = abs((peaks_loads (i,1) - shear_values (1,i))/shear_values (1,i)) * 100;
                delta_reaction (1, i) = abs((0 - reaction_values (1,i))/reaction_values (1,i)) * 100;
    end
    %fprintf('********Finished Moment_Shear_SoilReaction - Sheet %i\n', sheet)
%    l=l+1;
    
 %end

%%Divide Positive and Negative Deflections values%%
  
for j= 1 : length(elevation_resolution)
    
        for deflection_change = 1 : pos_deflection_levels
            
            pos_curvature_values(j,deflection_change) = curvature_values(j, deflection_change);
            pos_rotation_values(j,deflection_change) = rotation_values(j,deflection_change);
            pos_deflection_values(j,deflection_change) = deflection_values(j,deflection_change);
            pos_moment_values(j,deflection_change) = moment_values(j,deflection_change);
            pos_shear_values(j,deflection_change) = shear_values(j,deflection_change);
            pos_reaction_values(j,deflection_change) = reaction_values(j,deflection_change);
           
            neg_curvature_values(j,deflection_change) = curvature_values(j,(deflection_change + pos_deflection_levels));
            neg_rotation_values(j,deflection_change) = rotation_values(j,(deflection_change + pos_deflection_levels));
            neg_deflection_values(j,deflection_change) = deflection_values(j,(deflection_change + pos_deflection_levels)); 
            neg_moment_values(j,deflection_change) = moment_values(j,(deflection_change + pos_deflection_levels));
            neg_shear_values(j,deflection_change) = shear_values(j,(deflection_change + pos_deflection_levels));
            neg_reaction_values(j,deflection_change) = reaction_values(j,(deflection_change + pos_deflection_levels));
            
        end 
          
end
     

%%summary plots- positive deflections%%
if plot_pos_values
        figure
        plot(elevation_resolution,pos_curvature_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Positive Curvature')
        xlabel('elevation, in') 
        ylabel('positive curvature, 1/in')
        for x = 1 : pos_deflection_levels
            peaks_deflection1 = peaks_deflection(x)'
            legendInfo{x} = ['Deflection [in] = ' num2str(peaks_deflection1)];
            legend(legendInfo);
        end
        grid on
        
        figure
        plot(elevation_resolution,pos_rotation_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Positive Rotation')
        xlabel('elevation, in') 
        ylabel('Rotation, rad')
        legend(legendInfo); 
        grid on
        
        figure
        plot(elevation_resolution,pos_deflection_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Positive Deflection')
        xlabel('elevation, in') 
        ylabel('deflection, in')
        legend(legendInfo); 
        grid on
                
        figure
        plot(elevation_resolution,pos_moment_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Positive Moment')
        xlabel('elevation, in') 
        ylabel('moment, kips-in')
        legend(legendInfo); 
        grid on

        figure
        plot(elevation_resolution,pos_shear_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Positive Shear')
        xlabel('elevation, in') 
        ylabel('shear, kips')
        legend(legendInfo); 
        grid on

        figure
        plot(elevation_resolution,pos_reaction_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Positive Soil Reaction')
        xlabel('elevation, in') 
        ylabel('soil reaction, kips/in')
        legend(legendInfo); 
        grid on
end 
%%Summary Plots- Negative Deflections%% 
if plot_neg_values
        figure
        plot(elevation_resolution,neg_curvature_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Negative Curvature')
        xlabel('Elevation, in') 
        ylabel('Curvature, 1/in')
        for xx = 1 : pos_deflection_levels
            peaks_deflection2 = peaks_deflection(xx + pos_deflection_levels,1)
            legendInfo{xx} = ['Deflection [in] = ' num2str(peaks_deflection2)];
            legend(legendInfo);   
        end
        grid on
    
        figure
        plot(elevation_resolution,neg_rotation_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Negative Rotation')
        xlabel('Elevation, in') 
        ylabel('Rotation, rad')
        legend(legendInfo);  
        grid on
        
        figure
        plot(elevation_resolution,neg_deflection_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Negative Deflection')
        xlabel('Elevation, in') 
        ylabel('Deflection, in')
        legend(legendInfo); 
        grid on
        
        figure
        plot(elevation_resolution,neg_moment_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Negative Moment')
        xlabel('Elevation, in') 
        ylabel('Moment, kips-in')
        legend(legendInfo); 
        grid on

        figure
        plot(elevation_resolution,neg_shear_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Negative Shear')
        xlabel('Elevation, in') 
        ylabel('Shear, kips')
        legend(legendInfo); 
        grid on

        figure
        plot(elevation_resolution,neg_reaction_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Negative Soil Reaction')
        xlabel('Elevation, in') 
        ylabel('Soil Reaction, kips/in')
        legend(legendInfo); 
        grid on
end        

%%Summary Plots- All Deflections%%
if plot_all_values
        figure
        plot(elevation_resolution,curvature_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Curvature')
        xlabel('Elevation, in') 
        ylabel('Curvature, 1/in')
        grid on

        figure
        plot(elevation_resolution,rotation_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Rotation')
        xlabel('Elevation, in') 
        ylabel('Rotation')
        %legend('Deflection= 0.4','Deflection= -0.4','Deflection= 0.6','Deflection= -0.6','Deflection= 0.8','Deflection= -0.8','Deflection= 0.9','Deflection= -0.9','Deflection= 1.0','Deflection= -1.0')
        grid on

        figure
        plot(elevation_resolution,deflection_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Deflection')
        xlabel('Elevation, in') 
        ylabel('Deflection, in')
        %legend('Deflection= 0.4','Deflection= -0.4','Deflection= 0.6','Deflection= -0.6','Deflection= 0.8','Deflection= -0.8','Deflection= 0.9','Deflection= -0.9','Deflection= 1.0','Deflection= -1.0')
        grid on

        figure
        plot(elevation_resolution,moment_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Moment')
        xlabel('Elevation, in') 
        ylabel('Moment, kips-in')
        %legend('Deflection= 0.4','Deflection= -0.4','Deflection= 0.6','Deflection= -0.6','Deflection= 0.8','Deflection= -0.8','Deflection= 0.9','Deflection= -0.9','Deflection= 1.0','Deflection= -1.0')
        grid on

        figure
        plot(elevation_resolution,shear_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Shear')
        xlabel('Elevation, in') 
        ylabel('Shear, kips')
        %legend('Deflection= 0.4','Deflection= -0.4','Deflection= 0.6','Deflection= -0.6','Deflection= 0.8','Deflection= -0.8','Deflection= 0.9','Deflection= -0.9','Deflection= 1.0','Deflection= -1.0')
        grid on

        figure
        plot(elevation_resolution,reaction_values, elevation_resolution, midline);
        xlim([-93 0])
        title('Soil Reaction')
        xlabel('Elevation, in') 
        ylabel('Soil Reaction, kips/in')
        %legend('Deflection= 0.4','Deflection= -0.4','Deflection= 0.6','Deflection= -0.6','Deflection= 0.8','Deflection= -0.8','Deflection= 0.9','Deflection= -0.9','Deflection= 1.0','Deflection= -1.0')
        grid on
end

if plot_py
    fprintf('Plotting Positive Deflection PY\n')
    figure
    hold all 
    for z = 1:10
          pos_p = abs (pos_reaction_values(z+5,:));
          pos_y = abs (pos_deflection_values(z+5,:));
         
          k0 = 0; %initial row position
          pos_p1 = [pos_p(:,1:k0) 0 pos_p(:,k0+1:end)]; 
          pos_y1 = [pos_y(:,1:k0) 0 pos_y(:,k0+1:end)]; 
       
          depth = elevation_resolution(z+5);
       for j = 1:length(pos_p1)
            pos_y_values(z ,j) = pos_y1(1,j);
            pos_p_values(z ,j) = pos_p1(1,j);
       end
        %pos_py = spap2 (4, 3, pos_y1, pos_p1); %spline to approximate p-y curves
        plot(pos_y1,pos_p1,'-o');
        %hold on
        %fnplt(pos_py);
        title('Positive p-y curves')
        xlabel('Deflection, in') 
        ylabel('Soil Reaction, kips/in')
        legendInfo{z} = ['Depth [in] = ' num2str(depth)];
        grid on
    end
    
    legend(legendInfo);
    hold off
    title('Positive Deflection PY');
    
    fprintf('Plotting Negative Deflection PY\n')
    figure
    hold all
    for z = 1:10   
            neg_p = abs (neg_reaction_values(z+5,:));
            neg_y = abs (neg_deflection_values(z+5,:));
            
            neg_p1 = [neg_p(:,1:k0) 0 neg_p(:,k0+1:end)]; 
            neg_y1 = [neg_y(:,1:k0) 0 neg_y(:,k0+1:end)]; 
            
            depth = elevation_resolution(z+5);
       for j = 1:length(neg_p1)
            neg_y_values(z ,j) = neg_y1(1,j);
            neg_p_values(z ,j) = neg_p1(1,j);
       end
       neg_py = spap2 (4, 3, neg_y1, neg_p1); %spline to approximate p-y curves
        plot(neg_y1,neg_p1,'-o');
        %hold on
        %fnplt(neg_py);
        title('Negative p-y curves')
        xlabel('Deflection, in') 
        ylabel('Soil Reaction, kips/in')
        legendInfo{z} = ['Depth [in] = ' num2str(depth)]; 
        
    end
    legend(legendInfo)
    hold off
    title('Negative Deflection PY');
    grid on

end

fprintf('Completed!\n')
