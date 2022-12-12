%% function to compute the quality indices for the comparison of different k
%
% 15.6.2017 - Daniela Z?ller: created based on function by Thomas
%               adapted for large sparse matrices structures

function [CDF,AUC] = ComputeClusteringQuality(Consensus,K_range)

    
    % Number of K values to check
    nK = length(K_range);
    n_items = size(Consensus,1);
    
    % Creates the CDF range
    c = 0:0.01:1;
    
    for k = 1:nK
        
        % Sorted consensus entries
        Cons_val = jUpperTriMatToVec(squeeze(Consensus(:,:,k)));
        Cons_val = sort(Cons_val,'ascend');
        
        % Computation of CDF
        for i = 1:length(c)
            CDF(k,i) = nnz(Cons_val <= c(i));
        end
        CDF(k,:)=CDF(k,:)./length(Cons_val);
        
        
        % vectorized computation of AUC
        AUC(k)=diff(c)*CDF(k,2:end)';
        
%         % Computation of the AUC
%         AUC(k) = 0;
%         for i = 2:length(Cons_val)%(n_items*(n_items-1)/2)
%             AUC(k) = AUC(k) + (Cons_val(i)-Cons_val(i-1))* interp1q(c',CDF(k,:)',Cons_val(i));
%         end

        clear Cons_val

    end
    
%     for k = 2:nK
%         % Computation of Delta
%         tmp_max_AUC = max(AUC(1:k-1));
%         if K_range(k) == 2
%             Delta(k) = 0;
%         else
%             Delta(k) = (AUC(k) - tmp_max_AUC)/tmp_max_AUC;
%         end
%     end
%     
%     % Computation of Delta AUC
%     max_AUC = AUC(1);
%     Delta(1) = AUC(1);
%     for k = 2:nK
%         Delta_percent(k) = (AUC(k) - max_AUC)/max_AUC;
%         Delta(k) = AUC(k) - max_AUC;
%         if AUC(k)>max_AUC
%             max_AUC=AUC(k);
%         end
%     end
    
end