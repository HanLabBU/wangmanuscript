function dual_sensor_vm_pp_peak_trace()

    fluorescence_time = 1; % second

    plot_window_time = [-0.5 0.5]; % second
    normalize_window_time = [-0.5 0]; % second

    vm_trace_fieldname = 'refine_trace';
    timestamp_fieldname = 'estimate_frame_time';
    
    fprintf('Select mat files....')
    [selected_files,selected_folder] = uigetfile('*.mat','MultiSelect','on');    
    if class(selected_files)=='char'
        file_list(1).name = selected_files;
    else
        file_list = cell2struct(selected_files,'name',1);
    end
    fprintf([num2str(length(file_list)),' files\n'])
    
    whole_tic = tic;
    
    for file_idx=1:length(file_list)
        
        cd(selected_folder)
    
        filename = file_list(file_idx).name;
        fprintf(['Analyzing ',filename,'\n'])
        load(filename)
        
        vm_trace_data = [];
        gcamp_trace_data = [];
        
        for trial_idx=1:numel(result.vm_trial)
            
            if file_idx==1 & trial_idx==1          
                vm_timestamp = result.vm_trial(trial_idx).(timestamp_fieldname);
                sampling_freq = 1/mean(diff(vm_timestamp));
                plot_window_idx = [floor(plot_window_time(1)*sampling_freq) ceil(plot_window_time(2)*sampling_freq)];
                x_axis = [plot_window_idx(1):plot_window_idx(2)]/sampling_freq;
                normalize_window_idx = find(x_axis>=normalize_window_time(1) & x_axis<=normalize_window_time(2));
            end
            
            if trial_idx==1
                vm_fluorescence_idx = find((result.vm_trial(trial_idx).estimate_frame_time-result.vm_trial(trial_idx).estimate_frame_time(1))<=fluorescence_time);
                file_data(file_idx).vm_fluorescence = mean(result.vm_trial(trial_idx).trace(vm_fluorescence_idx));
                gcamp_fluorescence_idx = find((result.gcamp_trial(trial_idx).raw_estimate_frame_time-result.gcamp_trial(trial_idx).raw_estimate_frame_time(1))<=fluorescence_time);
                file_data(file_idx).gcamp_fluorescence = mean(result.gcamp_trial(trial_idx).raw_trace(gcamp_fluorescence_idx));
            end
            
            vm_trace = result.vm_trial(trial_idx).(vm_trace_fieldname);
            vm_trace = (vm_trace-min(vm_trace))/(max(vm_trace)-min(vm_trace));

            all_gcamp_trace = cat(2,[result.gcamp_trial.trace]);
            gcamp_trace = result.gcamp_trial(trial_idx).trace;
            gcamp_trace = (gcamp_trace-min(all_gcamp_trace(:)))/(max(all_gcamp_trace(:))-min(all_gcamp_trace(:)));
        
            current_idx_list = result.vm_trial(trial_idx).plateau_potential.peak_idx;
            
            if ~isempty(current_idx_list)
                for idx=1:numel(current_idx_list)
                    
                    current_idx = current_idx_list(idx);
                    
                    if current_idx+plot_window_idx(1)>=1 & current_idx+plot_window_idx(2)<=numel(vm_trace)
                        vm_trace_data = cat(2,vm_trace_data,vm_trace(current_idx+plot_window_idx(1):current_idx+plot_window_idx(2)));
                        gcamp_trace_data = cat(2,gcamp_trace_data,gcamp_trace(current_idx+plot_window_idx(1):current_idx+plot_window_idx(2)));
                    end
                    
                end
            end
        end
        
        vm_trace_data = bsxfun(@minus, vm_trace_data, mean(vm_trace_data(normalize_window_idx,:),1));
        gcamp_trace_data = bsxfun(@minus, gcamp_trace_data, mean(gcamp_trace_data(normalize_window_idx,:),1));
        
        file_data(file_idx).vm = vm_trace_data;
        file_data(file_idx).gcamp = gcamp_trace_data;
        
        file_data(file_idx).avg_vm = mean(vm_trace_data,2);
        file_data(file_idx).avg_gcamp = mean(gcamp_trace_data,2);
        if ~isempty(lowpass_freq)
            file_data(file_idx).avg_f_vm = mean(f_vm_trace_data,2);
        end
        
    end  

    plot_data = cat(2,[file_data.avg_vm]);
    plot_data_avg = mean(plot_data,2);
    plot_data_std = std(plot_data,[],2);
    vm_peak_amp = max(plot_data,[],1);
    figure;
    hold on
    plot(x_axis,mean(plot_data,2),'color','r');
    x_axis_shade = [x_axis, fliplr(x_axis)];
    shade_area = [reshape(plot_data_avg-plot_data_std,1,[]),fliplr(reshape(plot_data_avg+plot_data_std,1,[]))];
    fill(x_axis_shade, shade_area, 'r', 'FaceAlpha',0.1, 'EdgeColor','none');
    xlim(plot_window_time)
    xlabel('Time (s)')
    ylabel('Normalized Vm')
    title({'PP peak triggered Vm',['N=',num2str(size(plot_data,2)),' Neurons'],['Vm peak amp: ',num2str(mean(vm_peak_amp)),'+',num2str(std(vm_peak_amp))]},'Interpreter','none')

    plot_data = cat(2,[file_data.avg_gcamp]);
    plot_data_avg = mean(plot_data,2);
    plot_data_std = std(plot_data,[],2);
    [gcamp_peak_amp,gcamp_peak_amp_idx] = max(plot_data,[],1);
    gcamp_peak_amp_time = x_axis(gcamp_peak_amp_idx);
    gcamp_peak_amp_time_avg = mean(gcamp_peak_amp_time);
    gcamp_peak_amp_time_std = std(gcamp_peak_amp_time);

    figure;
    hold on
    plot(x_axis,mean(plot_data,2),'color','g');
    x_axis_shade = [x_axis, fliplr(x_axis)];
    shade_area = [reshape(plot_data_avg-plot_data_std,1,[]),fliplr(reshape(plot_data_avg+plot_data_std,1,[]))];
    fill(x_axis_shade, shade_area, 'g', 'FaceAlpha',0.1, 'EdgeColor','none');
    xlim(plot_window_time)
    xlabel('Time (s)')
    ylabel('Normalized GCaMP')
    title({'PP peak triggered GCaMP',['N=',num2str(size(plot_data,2)),' Neurons'],['GCaMP peak amp: ',num2str(mean(gcamp_peak_amp)),'+',num2str(std(gcamp_peak_amp)),' GCaMP peak time: ',num2str(gcamp_peak_amp_time_avg),'+',num2str(gcamp_peak_amp_time_std)]},'Interpreter','none')    

    fit_result = fitlm(vm_peak_amp,gcamp_peak_amp);
    figure;
    hold on
    scatter(vm_peak_amp,gcamp_peak_amp,'filled')
    plot(vm_peak_amp,fit_result.predict(reshape(vm_peak_amp,[],1)))
    axis image
    xlim([0 1])
    ylim([0 1])
    xlabel('Vm peak amp')
    ylabel('GCaMP peak amp')
    title({'PP peak',['N=',num2str(size(plot_data,2)),' Neurons'],['R square = ',num2str(fit_result.Rsquared.Ordinary)]},'Interpreter','none') 

    vm_flourence = [file_data.vm_fluorescence];
    gcamp_flourence = [file_data.gcamp_fluorescence];

    fit_result = fitlm(vm_flourence,gcamp_flourence);
    figure;
    hold on
    scatter(vm_flourence,gcamp_flourence,'filled')
    plot(vm_flourence,fit_result.predict(reshape(vm_flourence,[],1)))
    axis image
    xlabel('Vm raw flourence')
    ylabel('GCaMP raw flourence')
    title({['Raw flourence ',num2str(fluorescence_time),' second'],['N=',num2str(size(plot_data,2)),' Neurons'],['R square = ',num2str(fit_result.Rsquared.Ordinary)]},'Interpreter','none') 

    fit_result = fitlm(vm_flourence,vm_peak_amp);
    figure;
    hold on
    scatter(vm_flourence,vm_peak_amp,'filled')
    plot(vm_flourence,fit_result.predict(reshape(vm_flourence,[],1)))
    axis image
    xlabel('Vm raw flourence')
    ylabel('Vm peak amp')
    title({['Vm peak amp ',num2str(fluorescence_time),' second'],['N=',num2str(size(plot_data,2)),' Neurons'],['R square = ',num2str(fit_result.Rsquared.Ordinary)]},'Interpreter','none') 

    fit_result = fitlm(gcamp_flourence,gcamp_peak_amp);
    figure;
    hold on
    scatter(gcamp_flourence,gcamp_peak_amp,'filled')
    plot(gcamp_flourence,fit_result.predict(reshape(gcamp_flourence,[],1)))
    axis image
    xlabel('GCaMP raw flourence')
    ylabel('GCaMP peak amp')
    title({['GCaMP peak amp ',num2str(fluorescence_time),' second'],['N=',num2str(size(plot_data,2)),' Neurons'],['R square = ',num2str(fit_result.Rsquared.Ordinary)]},'Interpreter','none') 
    
    fprintf(['Total time: ',num2str(toc(whole_tic)),' seconds.\n']);

end