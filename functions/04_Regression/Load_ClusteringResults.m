function [out] = Load_ClusteringResults(param)

load(fullfile(param.outDir_main,'AI'));
load(fullfile(param.outDir_main,'AI_subject_labels'));
load(fullfile(param.outDir_main,'time_labels'));
load(fullfile(param.outDir_main,'subject_labels'));

load(fullfile(param.outDir_iCAPs,'iCAPs'));
load(fullfile(param.outDir_iCAPs,'IDX'));
load(fullfile(param.outDir_iCAPs,'dist_to_centroid'));


out.iCAPs=iCAPs;
out.AI=AI;
out.AI_subject_labels=AI_subject_labels;
out.IDX=IDX;
out.dist_to_centroid=dist_to_centroid;
out.subject_labels=subject_labels;
out.time_labels=time_labels;

end