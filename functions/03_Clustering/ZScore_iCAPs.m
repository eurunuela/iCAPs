function [iCAPs_z] = ZScore_iCAPs(iCAPs,I_sig,IDX)

    for i = 1:size(iCAPs,1)
        % 1. Temporal z-scoring with all frames belonging to the cluster
%         iCAPs(i,:) = iCAPs(i,:)./std(I_sig(IDX==i,:),[],1);
        
        % 2. Spatial z-scoring
        [a,b] = hist(iCAPs(i,:), 100);
        aind = find(a == max(a));
        med  = b(aind(1));
%         iCAPs(i,:) = (iCAPs(i,:)-med)/std(squeeze(iCAPs(i,:)));
        iCAPs_z(i,:) = (iCAPs(i,:)-med)/sqrt((sum((iCAPs(i,:)-med).^2))/length(iCAPs(i,:)));    % normalization copied from Isik
    end
end