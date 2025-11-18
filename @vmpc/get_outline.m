function [region_cells, peak_indices, peak_values] = get_outline(obj, varargin)
    % 1. Find all peaks above mean + 1.5*std
    % 2. Filter out peaks below threshold (percentile or % of max)
    % 3. Process remaining peaks from largest to smallest
    % 4. Extract regions stopping at 0.5*(peak-mean) distance
    
    Args = struct('threshold', 0.75, 'Smooth', 1, 'CellIndex', 1, ...
                  'std_multiplier', 1.0, 'use_percentile', 1);
    [Args, ~] = getOptArgs(varargin, Args);

    if Args.Smooth
        map_1d = obj.data.maps_adsm(Args.CellIndex, :);
    else
        map_1d = obj.data.maps_raw(Args.CellIndex, :);
    end

    grid_size = sqrt(length(map_1d));
    map_2d = reshape(map_1d, grid_size, grid_size);

    % Remove NaN values for statistics
    valid_values = map_2d(~isnan(map_2d));
    
    if isempty(valid_values)
        region_cells = {};
        peak_indices = [];
        peak_values = [];
        return;
    end
    
    % Calculate statistics
    map_mean = mean(valid_values);
    map_std = std(valid_values);
    map_max = max(valid_values);
    
    % Statistical threshold: mean + 1.5*std
    statistical_threshold = map_mean + Args.std_multiplier * map_std;
    
    % Detection threshold: percentile OR percentage of max
    if Args.use_percentile
        % Use PERCENTILE (e.g., 90th percentile)
        detection_threshold = prctile(valid_values, Args.threshold * 100);
        threshold_description = sprintf('%.0fth percentile', Args.threshold * 100);
    else
        % Use PERCENTAGE of max (e.g., 90% of max)
        detection_threshold = Args.threshold * map_max;
        threshold_description = sprintf('%.0f%% of max', Args.threshold * 100);
    end
    
    fprintf('Map stats: mean=%.2f, std=%.2f, max=%.2f\n', map_mean, map_std, map_max);
    fprintf('Statistical threshold (mean+%.1f*std): %.2f Hz\n', Args.std_multiplier, statistical_threshold);
    fprintf('Detection threshold (%s): %.2f Hz\n', threshold_description, detection_threshold);
    fprintf('Region extraction: stops at 0.5 distance from peak to mean\n\n');
    
    % Initialize outputs
    region_cells = {};
    peak_indices = [];
    peak_values = [];
    
    % Create working copy
    working_map = map_2d;
    already_assigned = false(grid_size, grid_size);
    
    % Find all peaks iteratively
    peak_count = 0;
    while true
        % Find the maximum in remaining map
        masked_map = working_map;
        masked_map(already_assigned) = -inf;
        masked_map(isnan(masked_map)) = -inf;
        
        [current_peak_val, peak_idx] = max(masked_map(:));
        
        % Check 1: Must be above statistical threshold
        if current_peak_val < statistical_threshold || isinf(current_peak_val)
            fprintf('No more peaks above statistical threshold (%.2f Hz)\n', statistical_threshold);
            break;
        end
        
        [peak_row, peak_col] = ind2sub([grid_size, grid_size], peak_idx);
        
        % Check 2: Must be above detection threshold
        if current_peak_val < detection_threshold
            fprintf('Excluding peak at (%d,%d): %.2f Hz (below %s = %.2f Hz)\n', ...
                peak_row, peak_col, current_peak_val, threshold_description, detection_threshold);
            already_assigned(peak_row, peak_col) = true;
            continue;
        end
        
        % Extract region
        peak_count = peak_count + 1;
        
        fprintf('Peak %d at (%d,%d): %.2f Hz\n', peak_count, peak_row, peak_col, current_peak_val);
        
        % Region threshold: 0.5 distance from peak to mean
        distance_from_mean = current_peak_val - map_mean;
        region_threshold = current_peak_val - 0.3 * distance_from_mean;
        
        fprintf('  Distance from mean: %.2f Hz\n', distance_from_mean);
        fprintf('  Region stops at: %.2f Hz\n', region_threshold);
        
        % Find connected region
        above_threshold = (map_2d >= region_threshold) & ~isnan(map_2d);
        region_mask = false(grid_size, grid_size);
        region_mask(peak_row, peak_col) = true;
        
        % Flood fill
        queue = [peak_row, peak_col];
        head = 1;
        
        while head <= size(queue, 1)
            r = queue(head, 1);
            c = queue(head, 2);
            head = head + 1;
            
            neighbors = [r-1, c; r+1, c; r, c-1; r, c+1];
            
            for i = 1:4
                nr = neighbors(i, 1);
                nc = neighbors(i, 2);
                
                if nr >= 1 && nr <= grid_size && nc >= 1 && nc <= grid_size
                    if ~region_mask(nr, nc) && above_threshold(nr, nc) && ~already_assigned(nr, nc)
                        region_mask(nr, nc) = true;
                        queue = [queue; nr, nc];
                    end
                end
            end
        end
        
        % Store region
        region_idx = find(region_mask);
        if ~isempty(region_idx)
            region_values = map_2d(region_idx);
            
            region_cells{end+1} = region_idx;
            peak_indices(end+1) = peak_idx;
            peak_values(end+1) = current_peak_val;
            
            fprintf('  Region size: %d pixels (%.2f-%.2f Hz)\n\n', ...
                length(region_idx), min(region_values), max(region_values));
            
            already_assigned = already_assigned | region_mask;
        end
    end
    
    fprintf('\n=== Summary ===\n');
    fprintf('Total regions extracted: %d\n', peak_count);
    
    % Convert to column vectors
    peak_indices = peak_indices(:);
    peak_values = peak_values(:);
end