%% This function generates surrogate data on which to apply total
% activation for the thresholding process (selecting relevant innovations).
% The phase of the surrogate data has been scrambled
%
% Inputs:
% - TC2 is a 2D (n_ret_voxels x n_time_points) matrix with the data to
% scramble
% - param is the structure containing TA-relevant parameters; here, we need
% the 'Dimension' (size of X, Y, Z and T), the 'NbrVoxels' (number of
% retained voxels), and the 'date_TA' (date and time when TA was launched
% for a particular trial) fields
%
% Outputs:
% - Surrogate is a n_ret_voxels x n_time_points 2D matrix with the
% surrogate data
function [Surrogate] = GenerateSurrogate(TC,Path,param,fid)

    % output matrix
    Surrogate = zeros(size(TC));

    for iter_tc=1:param.NbrVoxels,

        % phase_signal is a time x 1 vector filled with random phase
        % information (in rad, from -pi to pi)
        rand_signal=fft(rand(param.Dimension(4),1),param.Dimension(4));
        phase_signal=angle(rand_signal);

        % We multiply the magnitude of the original data with random phase
        % information to generate surrogate data
        Surrogate(iter_tc,:) = real(ifft(exp(1i*phase_signal).*abs(fft(TC(iter_tc,:)',param.Dimension(4))),param.Dimension(4)));
    end
    
    WriteInformation(fid,['Surrogate data generated and saved at: ',...
        fullfile(Path,'TA_results',param.title,'Surrogate'),'...']);
end