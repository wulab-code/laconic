function filteredim = filtercellareas(cellareavector,filter,im)

idx = cellareavector > filter;
% n_cells = sum(idx);

for j = 1:length(idx)
    if idx(j) == 0
        im(im == j) = 0;        
    end
end

% get cells 
im(im > 0) = 1;
filteredim = im;