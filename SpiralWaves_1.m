function SpiralWaves(StimProtocol)
    % Parameters
    ncols = 128; 
    nrows = 128; 
    dur = 25000; 
    h = 2.0; 
    h2 = h^2; 
    dt = 0.15; 
    Iex = 30; 
    mu = 1.0; 
    Gx = 1; 
    Gy = Gx/mu;
    a = 0.13; 
    b = 0.013; 
    c1 = 0.26; 
    c2 = 0.1; 
    d = 1.0;
    
    % Initialize variables
    v = zeros(nrows, ncols); 
    r = v;
    
    % Set initial stim current and pattern
    iex = zeros(nrows, ncols); 
    if StimProtocol == 1 
        iex(62:67, 62:67) = Iex;
    else
        iex(:, 1) = Iex;
    end
    
    % Setup figure
    figure('position', [500 600 256 256], 'color', [1 1 1], 'menubar', 'none');
    ih = imagesc(v); 
    set(ih, 'cdatamapping', 'direct');
    colormap(hot); 
    axis image off; 
    th = title('');
    
    % Create 'Quit' pushbutton
    uicontrol('style', 'pushbutton', 'string', 'Quit', ...
              'units', 'normalized', 'position', [.45 .02 .13 .07], ...
              'callback', @(src, event) set(gcf, 'userdata', 1));
    
    % Preallocate movie structure
    mov = struct('cdata', [], 'colormap', []);
    k = 0; % Frame index
    
    % Simulation loop
    for n = 1:dur
        % Two-point stimulation
        if StimProtocol == 1 && n == 3800
            iex(62:67, 49:54) = Iex;
        end
        
        % Cross-field stimulation
        if StimProtocol == 2 && n == 5400
            iex(end, :) = Iex; 
        end
        
        % Update boundary conditions
        vv = [[0 v(2, :) 0]; [v(:, 2) v v(:, end-1)]; [0 v(end-1, :) 0]];
        vxx = (vv(2:end-1, 1:end-2) + vv(2:end-1, 3:end) - 2*v) / h2; 
        vyy = (vv(1:end-2, 2:end-1) + vv(3:end, 2:end-1) - 2*v) / h2; 
        
        % FitzHugh-Nagumo equations
        dvdt = c1*v.*(v-a).*(1-v) - c2*v.*r + iex + Gx*vxx + Gy*vyy; 
        v_new = v + dvdt*dt;
        drdt = b*(v-d*r); 
        r = r + drdt*dt; 
        v = v_new; 
        
        % Update colormap
        m = 1 + round(63*v); 
        m = max(m, 1); 
        m = min(m, 64);
        set(ih, 'cdata', m);
        
        % Update title
        set(th, 'string', sprintf('%d %0.2f %0.2f', n, v(1, 1), r(1, 1))) 
        drawnow
        
        % Capture every 500th frame for movie
        if rem(n, 500) == 0 
            k = k + 1; 
            mov(k).cdata = getframe(gcf);
            mov(k).colormap = [];
        end
        
        % Check termination conditions
        if max(v(:)) < 1.0e-4 || ~isempty(get(gcf, 'userdata'))
            break;
        end
    end
    
    % Save movie
    if k > 0
        [fn, pn] = uiputfile('SpiralWaves.avi', 'Save movie as:');
        if ischar(fn)
            video_file = VideoWriter(fullfile(pn, fn), 'Motion JPEG AVI'); 
            video_file.Quality = 75; 
            open(video_file); 
            writeVideo(video_file, mov); 
            close(video_file);
        else
            disp('User cancelled movie save');
        end
    else
        disp('No frames captured for movie');
    end
    
    % Close figure
    close(gcf)
end
