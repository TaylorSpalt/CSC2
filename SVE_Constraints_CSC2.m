function [cos_d_ang_SVE] = ...
            SVE_Constraints_CSC2...
                (cos_d_ang_SVE,SVE_sin_theta)


%------------------------------------------------------------------------  
% Refer to Section 3.2 of "Constrained Spectral Conditioning for spatial
%   sound level estimation", Spalt et al., 2016 for details. 
%------------------------------------------------------------------------

%___Save the sign of the cross-spectra cosine projections onto the 
%   orthogonal phase vector for later use____
s_cos_d_ang = sign(cos_d_ang_SVE);

%___Modify the cosine projections by the sine projection of the
%   uncertainty in the steering vector phase (Fig. 4 & Eq. (20))_____
cos_d_ang_SVE = abs(cos_d_ang_SVE) - SVE_sin_theta;   

%___Negative magnitudes disqualify those ref sensors (Eq. (21))___
cos_d_ang_SVE(cos_d_ang_SVE<0) = NaN;            
   
%___Reintroduce sign________________________
cos_d_ang_SVE = s_cos_d_ang.*cos_d_ang_SVE;	


end