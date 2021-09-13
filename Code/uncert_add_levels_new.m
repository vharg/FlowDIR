clear all close all  clc;

tic
warning('off','all')

    %% Set up input dialogue
    
prompt = {'Volcano name:','DEM file:','Default swath length:', 'Buffer (m)', 'Elevation threshold (m):',...
    'Maximum number of steps allowed:', 'Capture uncertainty in start? (0/1)', 'Start uncertainty (m)'};
dlgtitle = 'FlowDir inputs'; dims = [1 50];
definput = {'Merapi','Merapi_2015.tif','800', '50','30', '200', '0', '5'};

inputs = inputdlg(prompt,dlgtitle,dims,definput);
volcano = inputs{1};
dem = inputs{2};
defaultSW = str2double(inputs(3));
buff = str2double(inputs(4));
Thr = str2double(inputs(5));
steps = str2double(inputs(6));
uncertainty = str2double(inputs(7));
noCells = str2double(inputs(8));
swath_nb = 360;

DEM = GRIDobj(dem); 
waitfor(msgbox('Plotting DEM...draw a polygon to clip, double click inside when finished'))
MASK = createmask(DEM, 'usehillshade'); 
DEMc = crop(DEM,MASK,NaN);

% Set the number of degrees to interpolate azimuths to
interp_to = 22.5; 
swL = defaultSW;

    %% Plot DEM and select start point
    
figure(1)
title('Click +/- to zoom in/out')
imageschs(DEMc); colormap(parula); xlabel('East'); ylabel('North'); hold on
FD = FLOWobj(DEMc);
A  = flowacc(FD);
S = STREAMobj(FD, flowacc(FD)>100); % Use low threshold to define streams
plot(S,'b')
waitfor(msgbox('Click again for start point'))
[craterX_temp, craterY_temp] = ginput(1);
% craterX_temp = [438896.5004];
% craterY_temp = [9166379.859];
hold on, scatter(craterX_temp, craterY_temp, 'rx')
xlabel('East'); ylabel('North');


    %% Get coords of start points used for start point uncertainty
    
% Generate a polygon
polyout = polybuffer([craterX_temp craterY_temp],'lines',noCells,'JointType','square'); 

xcells = linspace(min(polyout.Vertices(:,1)),max(polyout.Vertices(:,1)),...
    round((max(polyout.Vertices(:,1))-min(polyout.Vertices(:,1)))/DEM.cellsize));

ycells = linspace(min(polyout.Vertices(:,2)),max(polyout.Vertices(:,2)),...
    round((max(polyout.Vertices(:,1))-min(polyout.Vertices(:,1)))/DEM.cellsize));

[XCELLS, YCELLS] = meshgrid(xcells,ycells); 
xcells = XCELLS(:); ycells = YCELLS(:);

if uncertainty == 1
scatter(xcells, ycells, 'rx')
end

xcells = [craterX_temp; xcells]; ycells = [craterY_temp; ycells];

% space to store tables output for each point in stage 1.
T_all = cell(length(xcells),1); 

figure(2)
set(gcf, 'visible','off')
sub3 = subplot(2,3,4:6);

fprintf('%s\n', 'Running FlowDir, please wait...' );

    %% Calculate for each of the start points
if uncertainty == 0
    xcells = xcells(1); ycells = ycells(1);
else
    disp('Running with uncertainty on start point')
end
    
for coord = 1:length(xcells) % for all start points
    
    fprintf('%s[%d%s%d]%s\n', 'Start point # ', coord, '/', length(xcells), '...'); 
 
    
    craterX = xcells(coord); craterY = ycells(coord);

    % Setup swaths, each swath will have 2 pixels width.
    swath_width = DEM.cellsize.*2;
    % Generate azimuths
    angles = linspace(0, 359, swath_nb); 

    % Get coords points on a circle (radius = swL) for each angle
    x = swL * cosd(angles) + craterX; 
    y = swL * sind(angles) + craterY;

    % Correct radians to geographic azimuth
    tmp = zeros(swath_nb,1); 
    
        for iS = 1:swath_nb
        tmp(iS) = 90-(iS-1)*(swath_nb/360);
        end
        
    % If less than zero add 360 to make positive 
    tmp(tmp<0) = tmp(tmp<0)+360;
    [~, angleIdx] = sort(tmp);
    
    % x and y are now points on a clockwise circle that has 1 = N
    x = x(angleIdx);
    y = y(angleIdx);

    swathSize = zeros(swath_nb, 1);
    
           %% Make swaths

        for iS = 1:swath_nb
        
     % Create swath from centre to each point
        SW = SWATHobj(DEM, [craterX, x(iS)], [craterY, y(iS)], 'width', swath_width);

        Z = SW.Z;
        X = SW.X(round(size(SW.X,1)/2),:);
        Y = SW.Y(round(size(SW.Y,1)/2),:);
    
        % DEM quality check: find values where elevation is <=0
        checkElev = find(Z<=0);
            if isempty(checkElev)== 0
                waitfor(msgbox('Swaths exist with elevation values <= 0 m, check the DEM', 'Error', 'error'));
                error('Swaths exist with elevation values <= 0 m, check the DEM');
            end
    
        % Take the mean of the pixels across the width, to smooth irregularities in the DEM
        swath = mean(SW.Z,1); 
    
            % Load the first swath to set the size of the storage matrix
            if iS == 1
                swathAll = zeros(length(swath), swath_nb);
            end
             
        % Handle swaths of different sizes
            if length(swath) == size(swathAll,1)
                swathAll(:, iS) = swath';
                elseif length(swath) < size(swathAll,1)
                swathAll(1:length(swath), iS) = swath';
            else
                swathAll(:, iS) = swath(1:size(swathAll,1));
            end
                
        % Do the same for X and Y coordinates of swath
            if iS == 1
            X_all = zeros(length(swath), swath_nb);
            Y_all = zeros(length(swath), swath_nb);
            end
        if length(swath) == size(X_all,1)
            X_all(:, iS) = X';
            Y_all(:, iS) = Y';
            elseif length(swath) < size(X_all,1)
            X_all(1:length(swath), iS) = X';
            Y_all(1:length(swath), iS) = Y';
        else
            X_all(:, iS) = X(1:size(swathAll,1));
            Y_all(:, iS) = Y(1:size(swathAll,1));

        end

        end

        %% Identify crater edges

    % For each azimuth find the maximum elevation along the swath, and clip the
    % swath using this
    swathAll_clip = zeros(size(swathAll));
    X_all_clip = zeros(size(swathAll));
    Y_all_clip = zeros(size(swathAll));
    
    % convert user input buffer to number of cells
    b = round(buff/DEM.cellsize);

        for i = 1:swath_nb
            sw = swathAll(:,i);
            maximum = max(sw);
            idMax = find(sw == maximum);
            swathAll_clip(1:idMax+b,i) = swathAll(1:idMax+b,i);
            X_all_clip(1:idMax+b,i) = X_all(1:idMax+b,i);
            Y_all_clip(1:idMax+b,i) = Y_all(1:idMax+b,i);
        end
    
    X_all_clip(X_all_clip == 0)= NaN;
    Y_all_clip(Y_all_clip == 0)= NaN;
    x_all = []; y_all = [];
    elev_all = [];
    id_clip_all = [];

        for i = 1:swath_nb
            id_clip = max(find(isnan(X_all_clip(:,i))==0));
            hold on
            id_clip_all = [id_clip_all, id_clip];
            x_all = [x_all, X_all_clip(id_clip,i)];
            y_all = [y_all, Y_all_clip(id_clip,i)];
            elev_all = [elev_all, swathAll_clip(id_clip,i)];
        end
    
        if coord == 1 % Plot the first profile (using the selected start point)
            plot(sub3,(1:360),elev_all, 'k')
            xlim([0, swath_nb])
            xlabel('Azimuth')
            ylabel('Crater elevation (m)')
            elevStart = swathAll(1,1);
            yline(elevStart, 'r--')
            text(20, double(elevStart+2),'Elevation of start point', 'color', 'r')
            set(gca, 'Color', 'none')
            X = X_all_clip(1:10,:);
            Y = Y_all_clip(1:10,:);
            set(gca, 'Position', [0.13,0.2173,0.8255,0.2613])
            set(gcf, 'Position', get(0, 'Screensize'));
            lim = ylim;

        end
        
    hold on,
    figure(1)
    plot(x_all, y_all, 'k')

        %% 1) Linear elevation gradient

    swathAll_clip(swathAll_clip==0)= NaN;
    swathAll_clip = swathAll_clip(1:id_clip,:);
    gradZn = [];  
        % Calculate the gradient between adjacent cells along swath
        % (positive grad = uphill)
        for iS = 1:swath_nb
            gradZn(:,iS) = gradient(swathAll_clip(:,iS),DEM.cellsize); 
        end

    gradZn(gradZn == 0)= NaN;
    gradZn( ~any(gradZn,2), : ) = [];

    sumAllG = cumsum(gradZn,1); % Cumulative sum of gradient
    sumAllG(isnan(sumAllG))= 0; 
    sumAllG( ~any(sumAllG,2), : ) = [];
    sumCountG = sum(sumAllG, 1); 

        
    %% Bin the summed gradients into secondary intercardinal directions
        
    rank_val_G = zeros(1,(length(sumCountG)/interp_to));
    rank_val_G = rank_val_G'; 

    LL = ((interp_to/2:interp_to:length(sumCountG)))';
    UL = ((LL+interp_to)-1);
    UL(end) = LL(1)-1;

    mid = (LL+UL)/2; mid(end) = 360;

        for i = 1:length(LL)
            rank_val_G(i) = mean(sumCountG(LL(i):UL(i)));
            
                if LL(i) > UL(i) % Wrap around the bin edges
                    rank_val_G(i) = mean([sumCountG(LL(i):end), sumCountG(1:UL(i))]);
                end    
        end
    
    T = table(LL, UL, mid, rank_val_G);
    
    % Invert so that high values are the most likely travel direction
    T.inv_rank_val_G = max(T.rank_val_G)-T.rank_val_G;
    % Sum all of the inverse rank values
    tot_SUM_G = sum(T.inv_rank_val_G);
    % Find the percentage of the total sum
    T.prob_G = (T.inv_rank_val_G/tot_SUM_G)*100;
    save('Table_strt_1', 'T');
    
    t0 = T; 
    T_all{coord} = t0; % One table is generated per start point

    
        %% 2) Path of steepest descent
        
    swathAll_clip( ~any(swathAll_clip,2), : ) = [];

    [s1, sub2] = size(swathAll_clip);
    
    % create empty matrix a bit bigger than swathAll
    mat = zeros(s1+1, sub2+2); 

    % Pad out matrix so that the calculation cant go backwards
    mat(1,:) = 10000; mat(:,1) = 10000; mat(:,sub2+2) = 10000;
    mat(2:end, 2:end-1) = swathAll_clip;
    mat_all = [];

 
    i = 2; % Start row

parfor j = 2:360  % start col
        
        i_all = []; j_all = [];
        i_all = [i_all, i]; j_all = [j_all, j];
        
        % Get cells surrounding active cell
        surround_cells = [mat(i, j-1); mat(i, j+1); mat(i+1, j-1); mat(i+1, j); mat(i+1, j+1)]; 
        
        idI_surround = [i, i, i+1, i+1, i+1]'; % store positions 
        idJ_surround = [j-1, j+1, j-1, j, j+1]';

        % Wrap the matrix edges
        % wrap azimuth 1 to 360 (moving to left through mat)
            if numel(j_all) >=2
                if j_all(end) < j_all(end-1)
                    i_all = [i_all, i];
                    j_all = [j_all, j];
                    wrap_1 = find(idJ_surround ==1);
                    idJ_surround(wrap_1) = 361;
                    surround_cells(wrap_1) = mat(idI_surround(wrap_1), 361);
                end
            end
            
            % If we are at azimuth 1
            if numel(j_all) == 1
                wrap_1 = find(idJ_surround == 1);
                idJ_surround(wrap_1) = 361;
                surround_cells(wrap_1) = mat(idI_surround(wrap_1), 361);
            end
   
            % wrap azimuth 360 to 1 (moving to right through mat)
            if numel(j_all) >=2
                if j_all(end) > j_all(end-1)
                % wrap azimuth 360 to 1 (moving to right through mat)
                wrap_360 = find(idJ_surround ==362);
                idJ_surround(wrap_360) = 2;
                surround_cells(wrap_360) = mat(idI_surround(wrap_360), 2);
                end
            end

surround = table(idI_surround, idJ_surround, surround_cells);

count = 1;

    while i_all(end) < steps %(i.e specify the max number of cells expected to be traversed)    
    
    count = count +1; 
    
        if count > steps % break if max number of steps is reached
            break
        end 
        
    surround = table(idI_surround, idJ_surround, surround_cells);
    S_surround = sortrows(surround, 'surround_cells'); %sort the table by lowest elev value

    % Change 10000 to NaN
    S_surround.surround_cells(S_surround.surround_cells==10000)=NaN;
    % delete rows with Nan
    S_surround = S_surround(~any(ismissing(S_surround),2),:);


    % Check if any cells in table have already been traversed at any point along the path, if 
    % yes delete from table for each row of the table, check the i and j values against all
    % traversed
    SidI = S_surround.idI_surround;
    SidJ = S_surround.idJ_surround;

tId = [];

        for c = 1:height(S_surround)
        id_t = find((SidI(c) == i_all(:)) & (SidJ(c) == j_all(:)));
            if isempty(id_t)==0
            tId = [tId, c];
            end
        end
        
       % IF all cells have already been traversed...
        if numel(tId) == height(S_surround)
               tId = tId(2:end-1);
        end
   
       S_surround(tId,:) = []; %delete rows

        if numel(S_surround)== 0
            break
        end 
        %% Check for cells with the same elevation
        
    %if there are no repeat values in the lowest vals in the table use the first line of the table

    if height(S_surround) > 1

        if S_surround.surround_cells(1) < S_surround.surround_cells(2)    
            iL = S_surround.idI_surround(1);
            jL = S_surround.idJ_surround(1);  
        end

    % if the two lowest values are the same...
    if S_surround.surround_cells(1) == S_surround.surround_cells(2)

        % check to see if there are more than 2    
        r = find(S_surround.surround_cells==S_surround.surround_cells(1)); 
        
    % Use the one with the highest row id (further from the crater)   
    if max(S_surround.idI_surround(r)) > min(S_surround.idI_surround(r))
     [maxI, idmaxI] = max(S_surround.idI_surround(r));
     iL = S_surround.idI_surround(idmaxI);
     jL = S_surround.idJ_surround(idmaxI);
         
        % if the row ids are the same, turn to the column     
        elseif max(S_surround.idI_surround(r)) == min(S_surround.idI_surround(r))
            [maxJ, idmaxJ] = max(S_surround.idJ_surround(r));
            jL = S_surround.idJ_surround(idmaxJ);
            iL = S_surround.idI_surround(idmaxJ);
    end
    end   
    else
        iL = S_surround.idI_surround(1);
        jL = S_surround.idJ_surround(1);  
    end

    if iL > s1
        break
    end
    
    % Add the lowest cell to the vector
    i_all = [i_all, iL]; j_all = [j_all, jL];
  
    surround_cells = [mat(iL, jL-1); mat(iL, jL+1); mat(iL+1, jL-1); mat(iL+1, jL); mat(iL+1, jL+1)];

    % Also store the positions of the surrounding cells in  mat
    idI_surround = [iL, iL, iL+1, iL+1, iL+1]';
    idJ_surround = [jL-1, jL+1, jL-1, jL, jL+1]';
    end   
    new_mat = zeros(size(mat));

    for m = 1:length(i_all)
        new_mat(i_all(m), j_all(m)) = 1;
    end  

 
    mat_all(:,:,j) = new_mat; % add the new matrix
end
    
    % Sum the matrix to see the 'most travelled paths'
    output = sum(mat_all,3);
    
    % realign matrices
    output = output(2:end, 2:end-1);
    is_nan = find(isnan(swathAll_clip));
    output(is_nan)= NaN;
    output(output ==0 )= NaN;
    % distance * 360 * number of start points sized matrix
    output_all(1:size(output,1),1:size(output,2),coord) = output;
end


    %% Process results for plotting  
    
    %% Calculate the average and standard deviation for bins (stage 1)

    % extract the probability column of each table
    all_probG = []; 
       for i = 1:length(xcells)
            act = T_all{i}.prob_G;
            all_probG = [all_probG, act];       
       end
     
    sumT = sum(all_probG,2)/length(xcells);
    
    % calculate SD
        for i = 1:size(all_probG,1)
            stdev(i) = std(all_probG(i,:)); 
        end
        % convert std into standard error
    SE = stdev/sqrt(length(xcells)); 

    
  %% Calculate the average and standard deviation for bins (stage 2)
    output_all(isnan(output_all))= 0;
    
        % sum down cols to get a score per azimuth for each start point
        % sum in 3rd dimension to get average hits per cell
        O = sum(output_all,3)/length(xcells); 
        
    sum_output_all = [];
        for j = 1:length(xcells)
        sum_output_all(:,:,j) = sum(output_all(:,:,j),1);
        end 

    avHits = sum_output_all/length(xcells); % this is the average

        % Bin avHits for sub-cardinal directions
    for i = 1:length(LL)
        for j = 1:length(xcells)
        avHits_temp = avHits(:,:,j);
        avHits_int(i,:,j) = mean(avHits_temp(LL(i):UL(i)));
        % Wrap around bin edges
        if LL(i) > UL(i)
            avHits_int(i,:,j) = mean([avHits_temp(LL(i):end), avHits_temp(1:UL(i))]);
        end    
        end
    end

  % find average per bin
av2 = mean(avHits_int,3);

  % find STD per bin
sq_avHits = squeeze(avHits_int);
std2 = std(sq_avHits,0,2);

% convert std into standard error
SE2 = std2/sqrt(length(xcells)); 


    
    %% 3) Total elevation change over swath length 
    
    % For start point 1
load('Table_strt_1.mat')
[r, c] = size(swathAll_clip);
Dif = zeros(r-1, c);

    % Calculates the difference in elevation between adjacent cells along
    % an azimuth
    for iS = 1:swath_nb
    Dif(:,iS) = diff(swathAll_clip(:,iS));
    end
    
Dif(isnan(Dif)) = 0;
sumDif = sum(Dif, 1);
Mall = [];
swathSizeIntAll = [];

for m = (1:interp_to:swath_nb)
M = mean(sumDif(m:(m+(interp_to-1))));
swathSizeInt = mean(swathSize(m:round(m+(interp_to-1))));
swathSizeIntAll = [swathSizeIntAll,swathSizeInt]; % interpolate swath size from 360 to 36
Mall = [Mall, M];
end

T.elevChange = Mall';

    % Identify elevation change above threshold 
T.isBelow = T.elevChange < Thr;
T.isAbove = T.elevChange >= Thr;
T.color(T.isBelow) = 'g';
T.color(T.isAbove) = 'r';


    %% ########################### Plotting ###########################
    
figure(2)
sub1 = subplot(2,3,1);
imageschs(DEMc), colormap(sub1, parula); 
hold on; plot(S, 'b'); scatter(xcells, ycells, 'rx')
plot(x_all, y_all, 'k')
xlabel('East')
ylabel('North')
col2 = colorbar; col2.Limits(2) = max(max(DEMc.Z));
ylabel(col2, 'Elevation (m)')
t1 = title(sprintf('%s', volcano), 'fontsize', 18);
t1.Position(2) = t1.Position(2)+1e2;
set(gca, 'Color', 'None'), set(gca, 'Fontsize', 14)

subplot(2,3,2)
sub2 = subplot(2,3,2);
for i = 1:length(T.elevChange)
    color = T.color(i);
polarscatter(T.mid(i)*pi/180,  max(T.prob_G)+1, 'd', color , 'filled')
hold on
end
edges = [(T.LL*pi/180); ((T.LL(end)+interp_to)*pi/180)]';
p = polarhistogram('BinEdges',edges,'BinCounts',sumT, 'FaceColor','blue','FaceAlpha',.7);
set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top');
thetaticklabels({'{\bf N}', 'NNE', 'NE', 'ENE', '{\bf E}', 'ESE', 'SE', 'SSE', '{\bf S}',...
    'SSW', 'SW', 'WSW', '{\bf W}', 'WNW', 'NW', 'NNW'})
thetaticks(0:interp_to:swath_nb)
title('Linear elevation gradient')
hold on
ang = T.mid.*pi/180;
polarwitherrorbar(ang,sumT,SE') % plot error bars on the data
set(gca, 'Fontsize', 14)

all = [];
rtl = rticklabels;
for r = 1:length(rtl)
    tmp = rtl{r};
    new_rtl = sprintf('%s%s', tmp, '%');
all = char(all, new_rtl);
end
hold on
rticklabels(all)
set(gca, 'Color', 'None')

subplot(2,3,3)
sub3 = subplot(2,3,3);
% plot matrix
O(O==0) = NaN;
pcolor(O)
shading('flat')
colormap(sub3, jet)
col3 = colorbar;
ylabel(col3, 'Number of hits')
xlabel('Azimuth ({\circ})')
ylabel('Distance (m)')
ticks=get(gca,'YTickLabels');%retrieve current ticks
ticks=round(str2double(ticks)*DEM.cellsize); ticks = round(ticks-DEM.cellsize);
set(gca,'YTickLabels',ticks)
title('Path of steepest descent')
set(gca, 'Color', 'None'), set(gca, 'Fontsize', 14);

subplot(2,3,4:6)
plot((1:360),elev_all, 'k')
xlim([0, swath_nb])
xlabel('Azimuth ({\circ})')
ylabel('Crater elevation (m)')
elevStart = swathAll(1,1);
yline(elevStart, 'r--')
text(20, double(elevStart+2),'Elevation of start point', 'color', 'r')
set(gca, 'Color', 'none')
X = X_all_clip(1:10,:);
Y = Y_all_clip(1:10,:);
set(gca, 'Position', [0.13,0.217,...
    0.825,0.261])
set(gcf, 'Position', get(0, 'Screensize'));
set(gca, 'Color', 'None'), set(gca, 'Fontsize', 14);
lim = ylim;
set(figure(2), 'PaperType', 'A1', 'PaperOrientation', 'landscape')

figure(3)
title('Path of steepest descent')
O(:,swath_nb+1) = NaN; O(O==0) = NaN;
R = (1:(size(O,1))); % number of rows in output
Az = (0:swath_nb);
N = 360;
temp = round(DEM.cellsize*(linspace(1,size(O, 1), size(O, 1))));
temp = round(temp - DEM.cellsize);
Rticks = string(temp);
pl = polarPcolor(R,Az,O, 'colormap', 'jet','Nspokes',17, 'Ncircles',size(O, 1),...
    'autoOrigin','off', 'RtickLabel', Rticks);
col4 = colorbar;
col4.Location = 'eastOutside';
ylabel(col4, 'Average number of hits')
t_tmp = get(gca, 'children');
t = findobj(t_tmp, 'type', 'text');
lab = {'NNW', 'NW','WNW','{\bf W}','WSW','SW','SSW','{\bf S}', 'SSE','SE','ESE','{\bf E}',...
    'ENE','NE','NNE', '{\bf N}'};
for i = 1:numel(lab)
    set(t(i), 'String', lab(i))
end

fprintf('%s\n', 'Saving figures...' ); 
    
% Save figures and workspace
if not(isfolder(sprintf('%s/%s','Out', volcano)))
    mkdir(sprintf('%s/%s','Out', volcano))
end

if not(isfolder(sprintf('%s/%s/%d', 'Out', volcano, 0)))
    mkdir(sprintf('%s/%s/%d', 'Out', volcano, 0))
end
    
runs = dir(sprintf('%s/%s', 'Out', volcano));
r = struct2cell(runs);
r_temp = r(1,:);
r_db = str2double(r_temp);
maxr = max(r_db);
    
    mkdir(sprintf('%s/%s/%d', 'Out', volcano, maxr+1))
    
    save(sprintf('%s/%s/%d/%s', 'Out', volcano, maxr+1, 'workspace'))
    savefig(figure(1), sprintf('%s/%s/%d/%s%s', 'Out',volcano, maxr+1, volcano, '_roi'))
    print(figure(1), '-dpdf', fullfile(sprintf('%s/%s/%d/%s%s', 'Out', volcano, maxr+1, volcano, '_roi')))
    savefig(figure(2), sprintf('%s/%s/%d/%s','Out', volcano, maxr+1, volcano))
    print(figure(2), '-dpdf', fullfile(sprintf('%s/%s/%d/%s','Out', volcano, maxr+1, volcano)))
    savefig(figure(3), sprintf('%s/%s/%d/%s%s', 'Out', volcano, maxr+1, volcano, '_polar'))
    print(figure(3), '-dpdf', fullfile(sprintf('%s/%s/%d/%s%s', 'Out', volcano, maxr+1, volcano, '_polar')))
    
delete Table_strt_1.mat;

str = sprintf('%s', 'Finished!' ); 
    disp(str)



