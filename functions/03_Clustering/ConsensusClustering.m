%% This function performs consensus clustering over a range of K values
% The goal is to provide a measure of how good each value of K is
% 
% Inputs:
% - X is the data matrix (n_DP x n_DIM)
% - K_range is the range of K values to examine
% - Subsample_type defines how subsampling is done: across items (data
% points) if 'items', and across dimensions if 'dims', across subjects (all
% items of one subject) if 'subjects'
% - Subsample_fraction is the fraction of the original data points, or
% dimensions, to keep for a given fold
% - n_folds is the number of folds over which to run
% - DistType: 
% - subject_labels is the assignment of every item in X to a subject, only
% needed for Subsample_type='subjects'
%
% Outputs:
%   saves the results for every K in the consensus subfolder:
%   Consensus_ordered - Consensus matrix for every fold (n_items x n_items)
%
%
% 14.6.2017 - Daniela:
%       created based on Thomas' version, added subject subsampling,
%       changed outputs to structs
% 29.05.2018 - Daniela:
%       changed data saving, updated Build_Connectivity_Matrix, removed
%       dimensions sampling



function ConsensusClustering(X,subject_labels,param)

    K_range=param.K_vect;
    n_folds=param.cons_n_folds;
    n_rep=param.n_folds;
    DistType=param.DistType;
    Subsample_fraction=param.Subsample_fraction;
    Subsample_type=param.Subsample_type;
    if isfield(param,'MaxIter')
        MaxIter=param.MaxIter;
    else
        MaxIter=100;
    end
    outDir_cons=param.outDir_cons;
    
    % Number of data points
    n_items = size(X,1);
    disp(['n_items ',num2str(n_items),'...']);
    
    % Number of dimensions
    n_dims = size(X,2);
    disp(['n_dims ',num2str(n_dims),'...']);
    
    % Loop over all K values to assess
    for k = 1:length(K_range)
        disp(['Running consensus clustering for K = ',num2str(K_range(k)),'...']);
        
        if exist(fullfile(outDir_cons,['consensusResults_' num2str(K_range(k)) '.mat']), 'file') && ...
            ~param.force_ConsensusClustering
            disp('consensus clustering already done, skipping ...')
            continue;
        end
        
        if exist(fullfile(outDir_cons,['Consensus_' num2str(K_range(k)) '.mat']),'file')
            load(fullfile(outDir_cons,['Consensus_' num2str(K_range(k)) '.mat']));
        else
        
        % Connectivity matrix that will contain 0s or 1s depending on whether
        % elements are clustered together or not
        M_sum=zeros(n_items,n_items);
        I_sum=zeros(n_items,n_items);
        if strcmp(Subsample_type,'dims')
            M_sum=zeros(n_dims,n_dims);
            I_sum=zeros(n_dims,n_dims);
        end
        
        % Loops over all the folds to perform clustering for
        for h = 1:n_folds
            disp(['Fold ' num2str(h) ':'])
            switch Subsample_type
                case 'items'
                    
                    % Number of items to subsample
                    n_items_ss = floor(Subsample_fraction*n_items);
                    
                    % Does the subsampling
                    [~,tmp_ss] = datasample(X,n_items_ss,1,'Replace',false);
                    
                    % Vector
                    I_vec = zeros(n_items,1);
                    I_vec(tmp_ss) = 1;
                    
                    %Constructs the indicator matrix
                    I_sum=I_sum+I_vec*I_vec';
                    
                    X_ss=X(I_vec>0,:);
                    tmp_ss=[];
%                 case 'dims'
%                     error('Subsampling according to dimensions is not implemented yet')
                    % this part of the routine is not testet, modifications
                    % may be needed before using it
%                     % Number of dimensions to subsample
%                     n_dims_ss = floor(Subsample_fraction*n_dims);
%                     
%                     % Does the subsampling
%                     [X_ss,tmp_ss] = datasample(X,n_dims_ss,2,'Replace',false);
%                     
%                     % Vector
%                     I_vec = zeros(n_dims,1);
%                     I_vec(tmp_ss) = 1;
%                     
%                     %Constructs the indicator matrix
%                     I_sum=I_sum+I_vec_s*I_vec_s';
                    
                case 'subjects'
                    
                    if ~exist('subject_labels','var') || isempty(subject_labels)
                        error('subject labels needed for subjects subsampling!')
                    end
                    
                    subject_list=unique(subject_labels);
                    n_subjects=length(subject_list);
                    
                    % Number of items to subsample
                    n_subjects_ss = floor(Subsample_fraction*n_subjects);
                    
                    % Does the subject subsampling
                    disp('Subject subsampling...'); tic
                    [subjects_ss,~] = datasample(subject_list,n_subjects_ss,1,'Replace',false);
                    
                    % Vector
                    I_vec = zeros(n_items,1);
                    for iS=1:n_subjects_ss
                        I_vec(subject_labels==subjects_ss(iS)) = 1;
                    end
                    
                    %Constructs the indicator matrix
                    I_sum=I_sum+I_vec*I_vec';
                    
                    % subsampled data
                    disp('Data subsampling...');tic
                    X_ss=X(I_vec>0,:);
                otherwise
                    errordlg('PROBLEM IN TYPE OF SUBSAMPLING');
            end
            
            % Does the clustering (for now, only with k-means), so that IDX
            % contains the indices for each datapoint
            disp('Clustering ...');tic
            IDX = kmeans(X_ss,K_range(k),'Distance',DistType,'Replicates',n_rep,'MaxIter',MaxIter);
            
            clear X_ss
            
            % Builds the connectivity matrix
            disp('Buiding connectivity matrix M ...');tic
            M_sum=M_sum+Build_Connectivity_Matrix(IDX,find(I_vec>0),Subsample_type,n_items);
            
            clear I_vec
            clear X_ss
            clear tmp_ss
            clear IDX
        end
        
        
        disp('Computing Consensus Matrix ...');tic
        
        % Constructs the consensus matrix for the considered K
        Consensus = M_sum./I_sum;
        if any(I_sum(:)==0)
            warning([num2str(nnz(isnan(Consensus(:)))) ' (' ...
                num2str(nnz(isnan(Consensus(:)))/length(Consensus(:))*100) ...
                '%) items have not been selected during subsampling, you should increase the number of folds!']);
            Consensus(isnan(Consensus))=0;
        end
        
        clear M_sum I_sum M I
        
        disp('Saving consensus results (not ordered) ...');
        save(fullfile(outDir_cons,['Consensus_' num2str(K_range(k))]),'Consensus','-v7.3');
        
        end
        
        disp('Ordering consensus matrix ...'); tic
        tree = linkage(1-Consensus,'average');
        toc
        % Leaf ordering to create a nicely looking matrix
        leafOrder = optimalleaforder(tree,1-Consensus);
        toc
        % Ordered consensus matrix
        Consensus_ordered = Consensus(leafOrder,leafOrder);
        toc
        
        clear Consensus leafOrder
        
        % computing CDF and AUC of consensus matrix
        [CDF(k,:),AUC(k,:)] = ComputeClusteringQuality(Consensus_ordered,K_range(k));
        toc
        disp('Saving consensus results ...');
        save(fullfile(outDir_cons,['Consensus_ordered_' num2str(K_range(k))]),'Consensus_ordered','-v7.3');
        
        figure;imagesc(Consensus_ordered,[0 1]);colorbar;title(['k= ' num2str(K_range(k))]);
%         savefig(fullfile(outDir_cons,['Consensus_ordered_' num2str(K_range(k))]));
        print(fullfile(outDir_cons,['Consensus_ordered_' num2str(K_range(k))]), '-dpng','-painters');
        close gcf
        
        clearvars Consensus_ordered
    end
    save(fullfile(outDir_cons,'CDF'),'CDF');
    save(fullfile(outDir_cons,'AUC'),'AUC');
    
    figure;hold on;
    for k = 1:length(K_range)
        pl(k)=plot(0:0.01:1,CDF(k,:),'linewidth',2);
    end
    legend(pl,[num2str(K_range')],'location','southeast');
    title('CDF');
    print(fullfile(outDir_cons,'CDF'),'-depsc2','-painters');
    
    figure;hold on;
    plot(K_range,AUC,'-o','linewidth',2);
    xlabel('K');
    title('AUC');
    print(fullfile(outDir_cons,'AUC'),'-depsc2','-painters');
    
end