function [G,G_diag,Y_true] = ...
            CSM_and_True_Source...
                (nFFT,c,S,sx,sy,sz,mx,my,mz,z_M_M,K,two_invK,M,v_G_diag,...
                    z_S,SNR,spert,sz0,v_S,v_M,v_K,ftrue,e,cTe,wY,wG) 
                    
    
%% Create frequency domain data that would be generated from user-defined,
%   computer-generated point sources "measured" at user-defined array
%   sensor locations using the free-field Green's function_________________
[G,G_diag,Gs,nsource,xsource,ysource,zsource,v_nsource] = ...
    Generate_FFT_Array_Data...
        (K,SNR,spert,sz0,c,two_invK,nFFT,ftrue,mx,my,mz,M,v_G_diag,...
            v_M,v_K);
           
%% Calculate true source distribution via frequency-domain beamforming
%   using the individual source Gs at the grid points defined by the 
%   scanning grid (as opposed to their true locations)_____________________
[Y_true] = ...
    True_Source_Distribution...
        (S,nsource,cTe,e,Gs,wY,xsource,ysource,zsource,sx,sy,sz,wG,...
            z_M_M,z_S,v_S,v_nsource);

        
end