function [goodCA_trace env_CSUS] = movingtimetraining(trace, US_timestoconvert, CA_timestamps, pos, put_1_for_a_2_for_b)
%outputs two matrices -- one vector of calcium traces that occur during movement OR cs/us
%and a second matrix that has two rows-- one row indicating environment (1 for A, 2 for B)
%and the second row indicating CS/US/none (CS=10, US=20, none =0)




pos = convertpostoframe(pos, CA_timestamps);

velthreshold = 9;
vel = ca_velocity(pos);
goodvel = find(vel(1,:)>=velthreshold); %want these for their velocity

%outputs train for CS/US
CSUStrain = converttoframe(US_timestoconvert, CA_timestamps);
goodCSUS = find(CSUStrain>0); %want these for US/US
allwant = [goodvel, goodCSUS];
allwant = sort(allwant);


%outputs wanted trace
goodCA_trace = trace(:,allwant);
%outputs wanted CS/US train
CSUStrain = CSUStrain(allwant);
%makes env train
env_train = ones(1,length(CSUStrain)).*put_1_for_a_2_for_b;
%concat

env_CSUS = [env_train; CSUStrain];
