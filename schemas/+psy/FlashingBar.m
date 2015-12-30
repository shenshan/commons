%{
psy.FlashingBar (manual) # flashing bar
-> psy.Condition
---
luminance              : float                         # (cd/m^2) mid-value luminance
contrast               : float                         # (0-1) Michelson contrast of values 0..255
orientation            : decimal(4,1)                  # (degrees) 0=horizontal,  90=vertical
offset                 : float                         # normalized by half-diagonal
width                  : float                         # normalized by half-diagonal
flash_on               : float                         # (s) 
flash_frequency        : float                         # (Hz) will be rounded to the nearest fraction of fps 
%}

classdef FlashingBar < dj.Relvar
end