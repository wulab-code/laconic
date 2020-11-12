function interior = getcellinterior(imlabel)

b_dC = imlabel == 2; % on dC, get all borders
o_dC = imlabel >=2; % get all borders+ cell+nuclei
interior = o_dC-b_dC; % subtract off borders to create cells