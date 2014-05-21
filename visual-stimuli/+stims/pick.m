function pick    % stims.pick allows picking one of several preconfigured visual stimuli�
reload(psy.getSchema)  % preload to speed up the initialization

parentTable = common.Animal;


movingBar = struct(...
    'prompt', 'moving bar: 0.1 Hz (300 s)', ...
    'logger', stims.core.Logger(psy.Session, psy.Condition, psy.Trial, psy.MovingBar), ...
    'constants', ...
    struct(...
    'stimulus', 'moving bar', ...  % stimulus name recorded in the session table
    'monitor_distance', 10, ... (cm)
    'monitor_size', 7, ...      (inches) diagonal
    'monitor_aspect', 1.7, ...  (physical aspect ratio W/H)
    'resolution_x', 1024, ...   (pixels)
    'resolution_y',  600 ...    (pixels)
    ), ...
    'blocks', 30, ...
    'stim', {{
    setParams(stims.MovingBar,...
    'pre_blank', 0, ...   (s) blank period preceding trials
    'luminance', 30, ...    cd/m^2 mid-value luminance"
    'contrast', 0.99, ...  Michelson contrast
    'bg_color', 30, ...   0-254
    'bar_color', 254, ... 0-254
    'direction', 0, ... (degrees) 0=north, 90=east
    'bar_length', 1, ... in units of half-diagonal
    'bar_width', 0.01, ... in units of half-diagonal
    'bar_offset', 0, ... offset to the right (when facing in direction of motion) in units of half-diagonal
    'start_pos', -1, ... the starting position of the bar moviement. 1 is the distance from the center to corner"
    'end_pos', 1, ... (-1 1) the ending position of the bar movement
    'trial_duration', 10) ... (s) movement duration
    }});
    



simpleGrating = struct(...
    'prompt', 'grating: 1s ON/1s OFF (384 s)', ...
    'logger', stims.core.Logger(psy.Session, psy.Condition, psy.Trial, psy.Grating), ...
    'constants', ...
    struct(...
    'stimulus', 'grating', ...  % stimulus name recorded in the session table
    'monitor_distance', 10, ... (cm)
    'monitor_size', 7, ...      (inches) diagonal
    'monitor_aspect', 1.7, ...  (physical aspect ratio W/H)
    'resolution_x', 1024, ...   (pixels)
    'resolution_y',  600 ...    (pixels)
    ), ...
    'blocks', 30, ...
    'stim', {{    
    setParams(stims.Grating, ...
    'direction', 0:30:359, ...
    'pre_blank', 2, ...
    'trial_duration', 0.3, ...
    'aperture_radius', 2.0, ...
    'second_photodiode', [-1 1], ...
    'init_phase', 0, ...
    'temp_freq', 4, ...
    'spatial_freq', 0.08) ...
    }});

jakeGrating = struct(...
    'prompt', 'JakeGrating: 0.5s ON/1s OFF (540s)', ...
    'logger', stims.core.Logger(psy.Session, psy.Condition, psy.Trial, psy.Grating), ...
    'constants', ...
    struct(...
    'stimulus', 'grating', ...  % stimulus name recorded in the session table
    'monitor_distance', 10, ... (cm)
    'monitor_size', 7, ...      (inches) diagonal
    'monitor_aspect', 1.7, ...  (physical aspect ratio W/H)
    'resolution_x', 1024, ...   (pixels)
    'resolution_y',  600 ...    (pixels)
    ), ...
<<<<<<< Updated upstream
    'blocks', 5, ...
=======
    'blocks', 5, ... % 100 for imaging @38000 frames
>>>>>>> Stashed changes
    'stim', {{    
    setParams(stims.Grating, ...
    'grating','sqr', ...
    'direction', 0:30:359, ... % 0:45:359 for imaging @38000 frames
    'pre_blank', 1 , ...
    'trial_duration', 0.5, ...
    'second_photodiode',[-1], ...
    'aperture_radius', 2.0, ...
    'init_phase', [0]) ...
    }});

gammaGrating = struct(...
    'prompt', 'gammaGrating for luminance measurements', ...
    'logger', stims.core.Logger(psy.Session, psy.Condition, psy.Trial, psy.Grating), ...
    'constants', ...
    struct(...
    'stimulus', 'grating', ...  % stimulus name recorded in the session table
    'monitor_distance', 10, ... (cm)
    'monitor_size', 7, ...      (inches) diagonal
    'monitor_aspect', 1.7, ...  (physical aspect ratio W/H)
    'resolution_x', 1024, ...   (pixels)
    'resolution_y',  600 ...    (pixels)
    ), ...
    'blocks', 1, ...
    'stim', {{    
    setParams(stims.Grating, ...
    'grating','sqr', ...
    'direction', 0:45:359, ...
    'pre_blank', 5 , ...
    'trial_duration', 5, ...
    'second_photodiode',[-1], ...
    'aperture_radius', 2.0, ...
    'temp_freq', 2,...
    'spatial_freq', 0.04, ...
    'init_phase', [0]) ...
    }});

groupedGrating = struct(...
    'prompt', 'Jake''s grouped stim ', ...
    'logger', stims.core.Logger(psy.Session, psy.Condition, psy.Trial, psy.Grating), ...
    'constants', struct(...
    'stimulus', 'grating', ... % stimulus name recorded in the session table
    'monitor_distance', 10, ...  (cm)
    'monitor_size', 7, ...       (inches) diagonal
    'monitor_aspect', 1.7, ...   (physical aspect ratio W/H)
    'resolution_x', 1024, ...     (pixels)
    'resolution_y',  600 ...      (pixels)
    ), ...
    'blocks', 20, ...
    'stim', {{
    
    setParams(stims.Grating, ...
    'second_photodiode', [-1 1], ...
    'direction', 0:15:359, ...
    'pre_blank', 0.1, ...
    'trial_duration', 0.2, ...
    'temp_freq', 2,...
    'spatial_freq', 0.04, ...
    'aperture_radius', 2.0)
    
    stims.Pause(10)
    
%     setParams(stims.Grating, ...
%     'second_photodiode',  [-1 1], ...
%     'direction', 0:30:359,  ...
%     'pre_blank', 0, ...
%     'trial_duration', 0.1, ...
%     'temp_freq', 2,...
%     'spatial_freq', 0.04, ...
%     'aperture_radius', 2.0)
%     
%     stims.Pause(20)
    }}...
    );


quadrantGrating = struct(...
    'prompt', 'Shan''s quadrant grating (960 s)', ...
    'logger', stims.core.Logger(psy.Session, psy.Condition, psy.Trial, psy.Grating), ...
    'constants', struct(...
    'stimulus', 'grating', ... % stimulus name recorded in the session table
    'monitor_distance', 10, ...  (cm)
    'monitor_size', 7, ...       (inches) diagonal
    'monitor_aspect', 1.7, ...   (physical aspect ratio W/H)
    'resolution_x', 1024, ...     (pixels)
    'resolution_y',  600 ...      (pixels)
    ), ...
    'blocks', 8, ...
    'stim', {{
    
    
    setParams(stims.Grating, ...    
    'pre_blank', 3, ...
    'trial_duration', 2, ...
    'direction', [90, 180], ...
    'aperture_radius', 0.15, ...
    'aperture_x', [-0.4,0.2], ...
    'aperture_y', [-0.36,0.32], ...
    'temp_freq', 4,...
    'spatial_freq',0.03)...
    
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 90, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', -0.47, ...
%     'aperture_y', 0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
    
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 90, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', 0.1, ...
%     'aperture_y', -0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
%     
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 90, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', 0.1, ...
%     'aperture_y', 0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
%     
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 180, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', -0.47, ...
%     'aperture_y', -0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
%     
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 180, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', -0.47, ...
%     'aperture_y', 0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
%     
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 180, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', 0.1, ...
%     'aperture_y', -0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
%     
%     setParams(stims.Grating, ...    
%     'pre_blank', 5, ...
%     'trial_duration', 2, ...
%     'direction', 180, ...
%     'aperture_radius', 0.2, ...
%     'aperture_x', 0.1, ...
%     'aperture_y', 0.32, ...
%     'temp_freq', 4,...
%     'spatial_freq',0.03)...
        
    }}...
    );


% menu items callback
menu = [
    simpleGrating
    groupedGrating
    quadrantGrating
    jakeGrating
    gammaGrating
    movingBar
    ];

clc, disp 'Welcome to stims.pick'

% enter primary key
while true
    try
        for keyField = parentTable.primaryKey
            key.(keyField{1}) = input(sprintf('Enter %s: ', keyField{1}));
            assert(~isempty(key.(keyField{1})), 'cannot have empty key')
        end
        disp 'Entered:'
        disp(key)
        assert(count(parentTable & key)==1, 'not found in database')
        break
    catch err
        disp(err.message)
    end
end

assert(isempty(javachk('desktop')), 'no MATLAB desktop! Restart.')
fprintf '\nAt runtime, press numbers to select stimulus, "r"=run, "q"=quit:\n'
for i = 1:length(menu)
    fprintf('%d. %s\n', i, menu(i).prompt)
end
fprintf \n\n

disp 'While the screen is blanked you can:'
disp '   press 1-9 to select or change the stimulus (memorize them now)'
disp '   press "r" to run the selected stimulus'
disp '   press ESC to stop an ongoing stimulus (only while frames are flipping)'
disp '   press "q" to quit'
disp ' '
disp 'Now press any key when you are ready to blank the screen.'

pause
stims.core.run(menu, key)
end