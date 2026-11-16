


%% Mass over time tracks

%%
% matches requires Matlab 2019b or later
drugloc   = matches(tracksColHeaders,'condition_drugID');
% drugloc   = 11;
concloc   = matches(tracksColHeaders,'condition_concentration (um)');
% concloc   = 12;
timeloc   = matches(tracksColHeaders,'Time (h)');
% timeloc   = 4;
frameloc  = matches(tracksColHeaders,'Frame ID');
% frameloc  = 9;
cellsloc  = matches(tracksColHeaders,'id');
% cellsloc  = 5;

%% Founder Cells Filter 

StartTime = 0;
EndTime = 72; 
minp = 10;

% find the frame id that exist between starttime and endtime
target      = tracks(:,timeloc)>= StartTime & tracks(:,timeloc)<= EndTime;
% find unique frames that are encompassed by StartTime and EndTime
frames      = unique(tracks(target,frameloc));
% get first and last frame
StartFrame  = min(frames);
EndFrame    = max(frames);

clearvars target frames
%%
% Find cells that start at frame 1
startFrame1Cells = unique(tracks(tracks(:, frameloc) == 1, cellsloc));

% Create a logical array to identify cells that start at frame 1
startFrame1Logical = ismember(tracks(:, cellsloc), startFrame1Cells);

% Check if cells that start at frame 1 have all frames after frame 1
target = startFrame1Logical & all(ismember(1:EndFrame, tracks(:, frameloc)), 2);

%[cellsID, ia, ~] = unique(tracks(target, cellsloc));
cellsID       = unique(tracks(target,cellsloc));
Ncount = histc(tracks(target, cellsloc), cellsID);

% throws out cells whose number of measurements inside of StartFrame and
% EndFrame is less than min path length (minp). This cuts off cells that
% would otherwise meet min path length, but due to a shortened StartTime or
% EndTime are excluded.
targetcells   = cellsID(Ncount >= minp);
newtracks  = tracks(target,:);

%%

for z = 1:length(targetcells) % parfor
    
    target_tracks  = newtracks(newtracks(:,5) == targetcells(z),:);
%     cells{z}      =[target_tracks(:,timeloc) target_tracks(:,3)];
    cells{z}      =[target_tracks(:,timeloc) medfilt1(target_tracks(:,3),5,'truncate') target_tracks(:,frameloc)];
    cells_loc(z)  = target_tracks(1,8);
    cells_drug(z) = target_tracks(1,drugloc);
    cells_conc(z) = target_tracks(1,concloc);    
    cells_ID(z)   = target_tracks(1,5);

    %     

end

%%
for i = 1:numel(cells)
    % column 1 of cells{i} is time
    [t_sorted, order] = sort(cells{i}(:,1));
    cells{i} = cells{i}(order, :);
end

    %% EXPONENTIAL FITTING

    leen = cellfun(@length, cells);

for i = 1:length(cells)
    % Set each cell initial time to zero
    cells_mod{i} = [cells{i}(:,1) - cells{i}(1,1), cells{i}(:,2)];
    
    % Log-transform the y-values for exponential fit
    log_y = log(cells_mod{i}(:,2));
    
    % Calculate parameters used for regression (linear fit on log-transformed data)
    xav(i) = mean(cells_mod{i}(:,1));
    x = cells_mod{i}(:,1) - xav(i);
    x2(i) = sum(x.^2);
    yav_log(i) = mean(log_y);
    y_log = log_y - yav_log(i);
    xy_log(i) = sum(x .* y_log);
    b_log(i) = xy_log(i) / x2(i);
    % Exponential growth rate of each cell
    a_log(i) = yav_log(i) - b_log(i) * xav(i);
    % Convert back to exponential parameters
    a_exp(i) = exp(a_log(i));  % Initial value (intercept)
    b_exp(i) = b_log(i);       % Growth rate
    
    % Calculate r-squared for the exponential fit
    fitted_log_y = a_log(i) + b_log(i) * cells_mod{i}(:,1);
    SS_tot = sum((log_y - mean(log_y)).^2);
    SS_res = sum((log_y - fitted_log_y).^2);
    r2_log(i) = 1 - (SS_res / SS_tot);
    
    % Additional statistics (optional)
    TSS_log(i) = sum(y_log.^2);
    LRSS_log(i) = (xy_log(i)^2 / x2(i));
    syx_log(i) = sqrt((TSS_log(i) - LRSS_log(i)) / (leen(i) - 2));
    syxy_log(i) = syx_log(i) / yav_log(i);
    sb_log(i) = sqrt(syx_log(i)^2 / x2(i));
end

%% plot tracks and exp fits 


c_idx = 1; 
% Create a single figure
f = figure;

% Maximize the figure window
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

hold on % Keep plotting on the same figure
for jj = 1:length(c_idx)
    j = c_idx(jj); 
    %x = cells_drug == j & b_exp > 0.02 & r2_log >= 0.5;
    x = cells_drug == j & r2_log >= 0.8;
    ccc = find(x);
    for i = 1:length(ccc)
        shifted_x = cells{ccc(i)}(:,1); % time
        shifted_y = cells{ccc(i)}(:,2);     % Mass

        t_shift = shifted_x - shifted_x(1);         % same as in fitting


        fitted = a_exp(ccc(i)) * exp(b_exp(ccc(i)) * t_shift);
        cell_idx = ccc(i);

    % Plot the track
        plot(shifted_x, shifted_y, 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
        plot(shifted_x, fitted);

    end
end
hold off
title(drugLabels{c_idx})
xlabel('time (h)')
ylabel('mass (pg)')


%% EXP + C FITTING

model = @(p, t) p(1)*exp(p(2)*t) + p(3);

nCells  = numel(cells);
a_ec   = nan(1,nCells);
b_ec   = nan(1,nCells);
c_ec   = nan(1,nCells);
r2_ec   = nan(1,nCells);   % R² for exp+C

for i = 1:nCells
    t = cells{i}(:,1);
    y = cells{i}(:,2);

    t = t - t(1);   

    p0 = [max(y)-min(y), 0.05, min(y)];
    lb = [0, -Inf, 0];
    ub = [Inf, Inf, Inf];

    p = lsqcurvefit(model, p0, t, y, lb, ub);

    a_ec(i) = p(1);
    b_ec(i) = p(2);
    c_ec(i) = p(3);

    y_fit = model(p, t);
    SS_tot = sum((y - mean(y)).^2);
    SS_res = sum((y - y_fit).^2);
    r2_ec(i) = 1 - SS_res/SS_tot;
end



%% plot tracks and fits for exp+C

c_idx = 1; 
% Create a single figure
f = figure;

% Maximize the figure window
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

hold on % Keep plotting on the same figure
    for jj = 1:length(c_idx)
        j = c_idx(jj); 
    
        x = cells_loc == j; %& r2_ec >= 0.75 & cells_loc == 57 ;
    
        ccc = find(x);
        for i = 1:length(ccc)
        cell_idx = ccc(i);   % index into cells{} and cells_ID
    
        shifted_x = cells{cell_idx}(:,1); % time
        shifted_y = cells{cell_idx}(:,2); % mass
        t_shift   = shifted_x - shifted_x(1); 
        %fitted    = a_ec(cell_idx) * exp(b_ec(cell_idx) * t_shift) + c_ec(cell_idx);
    
        % Plot the track
        plot(shifted_x, shifted_y, 'LineWidth', 1);% 'Color', [0.6 0.6 0.6]);
        %plot(t_shift, fitted);
    
        % Label with the correct cell ID
        text(t_shift(end), shifted_y(end), ...
            sprintf('%d', cells_ID(cell_idx)), ...   % <-- use cell_idx here
            'FontSize', 8, ...
            'Color', 'k', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','middle');
        end 
    end

hold off
title(drugLabels{c_idx})
%xlim([0 50])
xlabel('time (h)')
ylabel('mass (pg)')


%% Compute AIC for both models 
nCells = numel(cells);
RSS_exp  = nan(1,nCells);
AIC_exp  = nan(1,nCells);

for i = 1:nCells
    t = cells{i}(:,1);
    y = cells{i}(:,2);

    t_shift = t - t(1);
    y_fit = a_exp(i).*exp(b_exp(i).*t_shift);

    RSS_exp(i) = nansum((y - y_fit).^2);
end


RSS_ec = nan(1,nCells);
AIC_ec = nan(1,nCells);

for i = 1:nCells
    t = cells{i}(:,1);
    y = cells{i}(:,2);

    t_shift = t - t(1);
    y_fit = a_ec(i).*exp(b_ec(i).*t_shift) + c_ec(i);

    RSS_ec(i) = nansum((y - y_fit).^2);
end


for i = 1:nCells
    n = numel(cells{i}(:,2));

    AIC_exp(i) = n * log(RSS_exp(i)/n) + 2*2;   % 2 parameters
    AIC_ec(i)  = n * log(RSS_ec(i)/n)  + 2*3;   % 3 parameters
end


            


%% Save params

% basic "good fit" mask
good_fit = isfinite(b_ec);

% MAD-based outlier filter on b across all cells
b_med = median(b_ec(good_fit));
b_mad = mad(b_ec(good_fit), 1);   % median absolute deviation

k = 5;   % 3 = strict, 5 = lenient
mad_keep = abs(b_ec - b_med) <= k*b_mad;

good_keep = mad_keep & r2_ec >= 0.75;

% indices of cells that pass both r2 and MAD
gc = find(good_keep);
n  = numel(gc);


valid_cells = struct( ...
    'id',           cell(1, n), ...
    'conditionNum', cell(1, n), ...
    'condition',    cell(1, n), ...
    'time_of_div',  cell(1, n), ...
    'mass_at_div',  cell(1, n), ...
    'init_mass',    cell(1, n), ...
    'a',            cell(1, n), ...
    'b',            cell(1, n), ...
    'c',            cell(1, n), ...
    'y_intercept',  cell(1, n), ...
    'r2',           cell(1, n),...
    'RSS',          cell(1, n));   

for ii = 1:n
    cell_idx = gc(ii);


    cells_mod{cell_idx} = [cells{cell_idx}(:,1) - cells{cell_idx}(1,1), cells{cell_idx}(:,2)];

    t = cells_mod{cell_idx}(:,1);   % times
    m = cells_mod{cell_idx}(:,2);   % mass

    % core info
    valid_cells(ii).id           = targetcells(cell_idx);
    valid_cells(ii).conditionNum = cells_drug(cell_idx);
    valid_cells(ii).condition    = drugLabels{cells_drug(cell_idx)};
    valid_cells(ii).time_of_div  = max(t);
    [~, k] = max(t);
    valid_cells(ii).mass_at_div  = m(k); %m(end);
    valid_cells(ii).init_mass    = m(1);

    % exp+constant parameters (MAD-filtered cells only)
    valid_cells(ii).a  = a_ec(cell_idx);
    valid_cells(ii).b  = b_ec(cell_idx);
    valid_cells(ii).c  = c_ec(cell_idx);
    valid_cells(ii).y_intercept = a_ec(cell_idx) + c_ec(cell_idx);
    valid_cells(ii).r2 = r2_ec(cell_idx);
    valid_cells(ii).RSS = RSS_ec(cell_idx);
end

%%

% number of valid cells
n = numel(valid_cells);

% extract condition labels from struct
conds = categorical({valid_cells.condition});

% impose custom order
conds = categorical(conds, ...
    {'MZ CpG','MZ LPS','FO CpG','FO LPS'}, ...
    {'MZ CpG','MZ LPS','FO CpG','FO LPS'});

init_mass = [valid_cells.init_mass];

figure;
violinplot(init_mass, conds);
ylabel('Initial mass (pg)');
xlabel('Condition');
title('Initial mass from imaging');


%% Initial mass measured vs. fitted

y0_fit =  [valid_cells.y_intercept];
y0_meas = [valid_cells.init_mass];

figure;
scatter(y0_meas, y0_fit, 10, 'filled');
xlabel('Measured initial mass');
ylabel('Fitted intercept (a + c)');
grid on;
axis equal;
hold on;

lims = axis;                    
minv = min(lims([1 3]));
maxv = max(lims([2 4]));
plot([minv maxv], [minv maxv], 'k--', 'LineWidth', 1.5);

R = corr(y0_meas(:), y0_fit(:), 'rows','complete');
p = polyfit(y0_meas, y0_fit, 1);        % linear fit: y = p(1)*x + p(2)
xx = linspace(minv, maxv, 100);
yy = polyval(p, xx);
plot(xx, yy, 'r-', 'LineWidth', 1.5);
txt = sprintf('R = %.2f', R);
text(0.05, 0.9, txt, 'Units','normalized', ...
     'FontSize', 12, 'FontWeight','bold');
text(0.05, 0.82, sprintf('y = %.3f x + %.3f', p(1), p(2)), ...
     'Units','normalized','FontSize',12,'FontWeight','bold');









