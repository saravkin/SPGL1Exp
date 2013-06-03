j = batch('My_script_name','configuration','torque4x1','matlabpool',3,'CaptureDiary',true)
% I use 3 workers, but use your own
% replace 'torque4x1' with the configuration you want
% also in 'matlabpool' specify the max number of workers minus 1
% (1 worker will be used as head node)

% IMPORTANT: remember your job id!
jl.state                         % shows what state the job is in
jl.display                       % shows job details
jl.id                            % this give the Job ID. Rememeber the Job ID!

% Quit and wait, and when coming back
jA = getClusterJob('torque2x1for24', 2)  % in this example my job ID is 5 (getClusterJob is a script I wrote)
j.state                         % check that it says 'finished'
data = j.getAllOutputArguments  % retrieve all the remote job variables into 'data'
jl.diary                         % displays all the command line window output of the remote job

% finally after you are all done:
j.destroy                       % cleans up the job


%%

j1  = batch('SPGL1_lowrank_HSS_l1','configuration','torque1x2for24','matlabpool',1,'CaptureDiary',true)
j2  = batch('SPGL1_lowrank_HSS_l2','configuration','torque1x2for24','matlabpool',1,'CaptureDiary',true)


j2  = batch('interpolation4d_LowRank50_BG_highfreq','configuration','torque1x2for24','matlabpool',1,'CaptureDiary',true)
j3  = batch('interpolation4d_LowRank30_BG_highfreq','configuration','torque1x2for24','matlabpool',1,'CaptureDiary',true)
j4  = batch('interpolation4d_LowRank10_BG_highfreq','configuration','torque1x2for24','matlabpool',1,'CaptureDiary',true)

