function imageFiles = imageOrderStruct(imagePath)
% AIM: Reads .tif files and sorts them naturally based on numeric sequences in filenames


%% 1. Retrieve file list
imageFiles = dir(fullfile(imagePath, '*.tif')); 
num_files = length(imageFiles);
file_indices = zeros(num_files, 1);

%% 2. Generate sorting indices
for i = 1:num_files
    % Extract all consecutive digits from the filename
    digits = regexp(imageFiles(i).name, '\d+', 'match');
    
    if ~isempty(digits)
        % Concatenate digits to form the sorting key (e.g., "Img_1_Ph_2.tif" -> 12)
        file_indices(i) = str2double(strjoin(digits, ''));
    end
end

%% 3. Execute sorting
[~, sort_order] = sort(file_indices);
imageFiles = imageFiles(sort_order);
end