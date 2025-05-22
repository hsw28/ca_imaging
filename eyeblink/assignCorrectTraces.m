function [trimmedCSUS_ID, trimmed_traces] = assignCorrectTraces(animalName)
%takes CSUS_id and changes numbering to either be 1 for correct trial or 0 for incorrect trial

animal = evalin('base', animalName);
CSUSid = animal.CSUS_id;
CRs = animal.CRs;
traces = animal.Ca_traces;
dateList = autoDateList(animal);
nDays = (numel(dateList));

trimmedCSUS_ID = struct();
trimmed_traces = struct();

for d = 1:nDays
    dateStr = dateList{d}
    CRs_Var_name = ['EMGts_' dateStr];
    CSUSid_Var_name = ['CSUS_id_' dateStr];
    traces_var_name = ['CA_traces_' dateStr];

    CRs_Var = animal.CRs.(CRs_Var_name);
    CSUSid_Var = animal.CSUS_id.(CSUSid_Var_name);
    traces_var = animal.Ca_traces.(traces_var_name);


    if length(CRs_Var) < 2
      fprintf('no trials on')
      trimmedCSUS_ID.(CSUSid_Var_name) = NaN;
      trimmed_traces.(traces_var_name) = NaN;
      continue;
    end


    %trials = find(CSUSid_Var(1,:)>0);
    trials = find(CSUSid_Var(1,:)>0 & CSUSid_Var(1,:)<7);
    if max(trials)>length(traces_var)
      warning('your trial data is longer than your cell traces for day')
      dateStr
      trimmedCSUS_ID.(CSUSid_Var_name) = NaN;
      trimmed_traces.(traces_var_name) = NaN;
      continue;
    end


    CSUSid_Var = CSUSid_Var(:,trials);
    traces_var = traces_var(:,trials);
    start=(find(CSUSid_Var(1,:)>1 & CSUSid_Var(1,:)<3));



    for k = 1:length(start)
      curtrial_start = start(k);
      if curtrial_start>1
        curtrial_start=curtrial_start-1;
          if CRs_Var(k) == 1
            CSUSid_Var(1,curtrial_start:curtrial_start+5) = 1;
          elseif CRs_Var(k) == 0
            CSUSid_Var(1,curtrial_start:curtrial_start+5) = 0;
          elseif isnan(CRs_Var(k)) == 1
            CSUSid_Var(1,curtrial_start:curtrial_start+5) = NaN;
          end
      else
          if CRs_Var(k) == 1
            CSUSid_Var(1,curtrial_start:curtrial_start+4) = 1;
          elseif CRs_Var(k) == 0
            CSUSid_Var(1,curtrial_start:curtrial_start+4) = 0;
          elseif isnan(CRs_Var(k)) == 1
            CSUSid_Var(1,curtrial_start:curtrial_start+4) = NaN;
          end
      end

    end

    legittrial = find(isnan(CSUSid_Var(1,:))==0);
    CSUSid_Var = CSUSid_Var(:,legittrial);
    traces_var = traces_var(:, legittrial);

    trimmedCSUS_ID.(CSUSid_Var_name) =CSUSid_Var;
    trimmed_traces.(traces_var_name) = traces_var;
end
