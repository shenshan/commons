%{
opt.BarrelMap (imported) # get a map for barrel field, with single whisker stimulation
->common.OpticalMovie
tau            : double         # time constant of regressor
-----
barrel_amp     : longblob        # amplitude map of barrel cortex, difference between on and off time period.
stimulus       : longblob        # raw stimulus signal from the piezo
mov_time       : longblob        # time of each frame of the movie, in sec
mov_ind        : longblob        # indicator of stimulus on off for each frame
%}

classdef BarrelMap < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = common.OpticalMovie & 'purpose="whisker"'
    end
    
	methods(Access=protected)

		function makeTuples(self, key)
		
            disp 'loading movie...'

            filename = fullfile(...
            fetch1(common.OpticalSession(key), 'opt_path'),...
            [fetch1(common.OpticalMovie(key), 'filename') '.h5']);
            [X, framerate, stim, pdFs] = opt.utils.getOpticalData(getLocalPath(filename),'pd');
            
            stim_ind = conv(abs(diff(stim)),gausswin(10),'same')>0.05;
            
            stim_time = (1:length(stim_ind))/pdFs;
            
            key.mov_time = (1:size(X,1))'/framerate;
            key.mov_ind = logical(interp1(stim_time, double(stim_ind)', key.mov_time, 'nearest','extrap'));
            
            onsets = key.mov_time(diff(key.mov_ind)==1);
            offsets = key.mov_time(diff(key.mov_ind)==-1);
            
            % check pairing of onsets and offsets
            if length(onsets)~=length(offsets)
                if length(onsets) == length(offsets)+1
                    onsets = onsets(1:end-1);
                else
                    error('Mismatched on off sets.')
            
                end
            end
            
            % time constant of intrinsic signal
            taus = 2;
            
            sz = size(X);
            X = reshape(X,sz(1),[]);
            for itau = 1:length(taus)
                tau = taus(itau);
                disp 'constructing design matrix...'
                G = zeros(length(key.mov_time),1);
                for iTrial = 1:length(onsets)
                    onset = onsets(iTrial);
                    offset = offsets(iTrial);
                    idx = find(key.mov_time>=onset & key.mov_time < offset);
                    G(idx) = G(idx) + 1 - exp((onset-key.mov_time(idx))/tau);
                    
                    idx = find(key.mov_time>=offset & key.mov_time < offset+5*tau);
                    G(idx) = G(idx) ...
                        + (1-exp((onset-offset)/tau))*exp((offset-key.mov_time(idx))/tau);
                    
                end
                G = G - mean(G);
                
                % regression
                disp 'regressing...'
               
                B = ne7.stats.regress(X, G, 0);
                key.barrel_amp = reshape(single(B)', sz(2), sz(3));
                key.tau = tau;
                key.stimulus = stim;
                self.insert(key)
            end
            
            
            
%             dt = 1/framerate;
%             off_period = floor(-2/dt):1:0;
%             on_period = floor(0/dt):1:floor(2/dt);
%             barrel_map = zeros(length(onset)-2,size(X,2),size(X,3));
%             
%             for iTrial = 1:length(onset)-2
%                 off_amp = mean(X(onset(iTrial+1) + off_period, :,:));
%                 on_amp = mean(X(onset(iTrial+1) + on_period, :,:));
%                 barrel_map(iTrial,:,:) = on_amp - off_amp;
%             end
%             
%             key.barrel_amp = squeeze(mean(barrel_map));
%             key.stimulus = stim;
% 			
      
		end
	end

end