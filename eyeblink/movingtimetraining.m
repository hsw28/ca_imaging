function [wanted goodCA_trace env_CSUS vel pos] = movingtimetraining(trace, CS_timestoconvert, US_timestoconvert, CA_timestamps, pos, put_1_for_a_2_for_b, velthreshold)
%outputs two matrices -- one vector of calcium traces that occur during movement OR cs/us
%and a second matrix that has two rows-- one row indicating environment (1 for A, 2 for B)
%and the second row indicating CS/US/none (CS=10, US=20, none =0)




pos = convertpostoframe(pos, CA_timestamps);

%velthreshold = 9;
vel = ca_velocity(pos);
goodvel = find(vel(1,:)>=velthreshold); %want these for their velocity

%outputs train for CS/US
CSUStrain = converttoframe(CS_timestoconvert, US_timestoconvert, CA_timestamps);
goodCSUS = find(CSUStrain>0); %want these for US/US

%if only want fast moving then:
first = min(goodCSUS);
last = max(goodCSUS);
[C,IA,IB] = intersect(goodvel,goodCSUS);
allwant = goodvel(setdiff(1:end,IA));

%%%%%%

%if want moving and CS/US
%first = min(goodCSUS);
%allwant = [goodvel, goodCSUS];
%allwant = sort(allwant);
%allwant = unique(allwant);

%outputs wanted trace
trace(:,1:first) = NaN;

goodCA_trace = trace(:,allwant);
nowant1 = max(find(isnan(goodCA_trace(1,:))==1));
nowant2 = min(find((goodCA_trace(1,:))==10000000000));

  if length(nowant1)>0 & length(nowant2)>0
    goodCA_trace = goodCA_trace(:,nowant1+1:nowant2-1);
  elseif length(nowant1)>0 & length(nowant2)==0
    goodCA_trace = goodCA_trace(:,nowant1+1:end);
  elseif length(nowant1)==0 & length(nowant2)>0
    goodCA_trace = goodCA_trace(:,1:nowant2-1);
  elseif length(nowant1)==0 & length(nowant2)==0
    goodCA_trace = goodCA_trace;
  end

%outputs wanted CS/US train
CSUStrain = CSUStrain(allwant);
%makes env train
env_train = ones(1,length(CSUStrain)).*put_1_for_a_2_for_b;
%concat


wanted1 = find(allwant>first);
wanted2 = find(allwant<last);
wanted3 = intersect(wanted1, wanted2);
wanted = allwant(wanted3);
vel = vel(:,wanted);
pos = pos(wanted, :);
env_CSUS = [env_train; CSUStrain];
